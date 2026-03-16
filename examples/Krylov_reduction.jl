using FlexibleMultibody
using Gridap, GridapGmsh
using Arpack
using LinearAlgebra
using SparseArrays
using Printf
using Plots

# ==============================================================================
#  МЕТОД КРЫЛОВСКИХ ПОДПРОСТРАНСТВ
# ==============================================================================

# ==============================================================================
#  Шаг 1: Строим модель
# ==============================================================================
mat = LinearElasticMaterial(2.1e11, 0.3, 7.85e3)

rbe_info = [
    RBEInfo("left_hole",  VectorValue(0.0, 0.0, 0.0)),
    RBEInfo("right_hole", VectorValue(0.7, 0.0, 0.0)),
]

println("Сборка модели...")
fl_rod  = FlexibleComponent("link_mesh.msh", mat; rbe_info)
K       = fl_rod.mK
M       = fl_rod.mM
n_total = size(K, 1)
println("Всего DOF: $n_total")

# ==============================================================================
#  Шаг 2: DOF возбуждения
# ==============================================================================
DOF_EXCITE  = 1
DOF_RESPOND = 1

b = zeros(n_total)
b[DOF_EXCITE] = 1.0

# ==============================================================================
#  Шаг 3: Построение базиса Крылова через процесс Арнольди
# ==============================================================================
KRYLOV_PER_DOF = 3    # векторов на каждый граничный DOF
                       # итого: 12 граничных × 3 = 36 векторов

println("\nПостроение блочного базиса Крылова...")
println("  Граничных DOF: $n_b  ×  $KRYLOV_PER_DOF векторов = $(n_b * KRYLOV_PER_DOF) векторов")

K_dense = Matrix(K)
M_dense = Matrix(M)

println("  Факторизация K...")
K_fact = factorize(K_dense)

# Собираем все базисные векторы сюда
cols = Vector{Vector{Float64}}()

# Для каждого граничного DOF строим цепочку Крылова
for dof in boundary_dofs
    # Стартовый вектор — единичная нагрузка в этом DOF
    b_dof = zeros(n_total)
    b_dof[dof] = 1.0

    v = K_fact \ b_dof   # статический отклик K⁻¹·b

    for k in 1:KRYLOV_PER_DOF
        # Ортогонализация Грама-Шмидта относительно всех уже добавленных векторов
        for u in cols
            v = v - dot(v, u) * u
        end

        nv = norm(v)
        nv < 1e-12 && break   # вектор линейно зависим — пропускаем

        v = v / nv
        push!(cols, copy(v))

        # Следующий вектор в цепочке: K⁻¹·M·v
        k < KRYLOV_PER_DOF && (v = K_fact \ (M_dense * v))
    end
end

V   = hcat(cols...)
n_V = size(V, 2)
println("  Итоговый размер базиса: $n_V векторов")

ortho_err = norm(V' * V - I)
println("  Ошибка ортогональности ‖VᵀV - I‖ = $(round(ortho_err, digits=8))")

# ==============================================================================
#  Шаг 4: Проецируем матрицы на подпространство Крылова
# ==============================================================================
println("\nПроекция матриц на подпространство Крылова...")
K_kr = Symmetric(V' * K_dense * V)
M_kr = Symmetric(V' * M_dense * V)

n_red = size(K_kr, 1)
println("Размер редуцированной системы: $n_red × $n_red")
println("Степень редукции: $(round(100*(1 - n_red/n_total), digits=1))%")

# ==============================================================================
#  Шаг 5: Собственные значения полной модели (для сравнения)
# ==============================================================================
NUM_MODES = 20
println("\nВычисление мод полной модели...")

λ_fb, Φ_fb = eigs(M, K; nev=NUM_MODES+6, which=:LM)
order_fb   = sortperm(real.(λ_fb), rev=true)
λ_fb       = real.(λ_fb)[order_fb]
Φ_fb       = real.(Φ_fb)[:, order_fb]

ω2_fb   = 1.0 ./ abs.(λ_fb)
mask_fb = ω2_fb .> (2π)^2
ω2_fb   = ω2_fb[mask_fb]
Φ_fb    = Φ_fb[:, mask_fb]

# Нормировка
for i in axes(Φ_fb, 2)
    mn = sqrt(abs(Φ_fb[:,i]' * M_dense * Φ_fb[:,i]))
    mn > 1e-12 && (Φ_fb[:,i] ./= mn)
end

freq_full = sqrt.(ω2_fb) ./ (2π)

# ==============================================================================
#  Шаг 6: Собственные значения редуцированной модели Крылова
# ==============================================================================
println("Вычисление мод редуцированной модели Крылова...")
λ_kr, Φ_kr = eigen(K_kr, M_kr)
order_kr   = sortperm(real.(λ_kr))
λ_kr       = real.(λ_kr)[order_kr]
Φ_kr       = real.(Φ_kr)[:, order_kr]

ω2_kr   = λ_kr
mask_kr = ω2_kr .> (2π)^2
ω2_kr   = ω2_kr[mask_kr]
Φ_kr    = Φ_kr[:, mask_kr]

freq_kr = sqrt.(ω2_kr) ./ (2π)

# Переводим моды в физическое пространство: Φ_phys = V · Φ_kr
Φ_kr_phys = V * Φ_kr

# Нормировка
for i in axes(Φ_kr_phys, 2)
    mn = sqrt(abs(Φ_kr_phys[:,i]' * M_dense * Φ_kr_phys[:,i]))
    mn > 1e-12 && (Φ_kr_phys[:,i] ./= mn)
end

# ==============================================================================
#  Шаг 7: Таблица сравнения частот
# ==============================================================================
n_cmp = min(length(freq_full), length(freq_kr), NUM_MODES)
println("\n" * "="^62)
println("    СРАВНЕНИЕ ЧАСТОТ (Гц) — Крыловское подпространство")
println("="^62)
@printf("%-6s  %-16s  %-16s  %-10s\n", "Мода", "Полная модель", "Krylov", "Ошибка %")
println("-"^62)
for i in 1:n_cmp
    f1  = freq_full[i]
    f2  = freq_kr[i]
    err = abs(f1 - f2) / f1 * 100
    @printf("%-6d  %-16.2f  %-16.2f  %-10.4f\n", i, f1, f2, err)
end
println("="^62)
println("Редуцированная модель: $n_red DOF вместо $n_total DOF")

# ==============================================================================
#  Шаг 8: FRF — полная модель vs Крылов
# ==============================================================================
println("\nПостроение FRF...")

ζ     = 0.02
f_min = 50.0
f_max = freq_full[end] * 1.3
freqs = range(f_min, f_max, length=3000)
ωs    = 2π .* collect(freqs)

# Полная модель
H_full = zeros(ComplexF64, length(ωs))
for (k, ω) in enumerate(ωs)
    for i in eachindex(ω2_fb)
        ωi = sqrt(ω2_fb[i])
        H_full[k] += Φ_fb[DOF_EXCITE, i] * Φ_fb[DOF_RESPOND, i] /
                     (ω2_fb[i] - ω^2 + 2im*ζ*ωi*ω)
    end
end

# Крыловская модель
H_kr = zeros(ComplexF64, length(ωs))
for (k, ω) in enumerate(ωs)
    for i in eachindex(ω2_kr)
        ωi = sqrt(ω2_kr[i])
        H_kr[k] += Φ_kr_phys[DOF_EXCITE, i] * Φ_kr_phys[DOF_RESPOND, i] /
                   (ω2_kr[i] - ω^2 + 2im*ζ*ωi*ω)
    end
end

# ==============================================================================
#  График
# ==============================================================================
H_full_dB = 20 .* log10.(abs.(H_full) .+ 1e-40)
H_kr_dB   = 20 .* log10.(abs.(H_kr)   .+ 1e-40)

p = plot(
    freqs, H_full_dB,
    label     = "Полная модель",
    color     = :blue,
    linewidth = 2,
    xlabel    = "Частота (Гц)",
    ylabel    = "Амплитуда |H(ω)| (дБ)",
    title     = "FRF: Krylov ($n_red векторов) vs Полная модель  |  DOF $DOF_EXCITE → $DOF_RESPOND",
    legend    = :topright,
    grid      = true,
    size      = (1100, 550)
)

plot!(p,
    freqs, H_kr_dB,
    label     = "Krylov ($n_red векторов)",
    color     = :green,
    linewidth = 1.5,
    linestyle = :dash
)

for f in freq_full
    vline!(p, [f], color=:gray, alpha=0.4, linewidth=1, label="")
end

savefig(p, "bode_krylov.png")
println("График сохранён: bode_krylov.png")
display(p)
