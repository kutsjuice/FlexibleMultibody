using FlexibleMultibody
using Gridap, GridapGmsh
using Arpack
using LinearAlgebra
using SparseArrays
using Printf
using Plots

# Редукция методом Крылова: H(ω) раскладывается в ряд по ω² около σ=0
# Базис строится из статических откликов на нагрузки в граничных DOF
# + инерционные поправки K⁻¹Mv для каждого вектора

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

# DOF входа и выхода для FRF
DOF_EXCITE  = 1
DOF_RESPOND = 1

b = zeros(n_total)
b[DOF_EXCITE] = 1.0

# Граничные DOF — два интерфейса RBE2 по 6 DOF каждый
boundary_dofs = Int[]
for (_, iface) in fl_rod.interfaces
    append!(boundary_dofs, collect(iface.dofs))
end
sort!(boundary_dofs)
interior_dofs = setdiff(1:n_total, boundary_dofs)
n_b = length(boundary_dofs)
n_i = length(interior_dofs)

# Число векторов Крылова на каждый граничный DOF:
# k=1 → только статика (= моды ограничений CB)
# k=2 → + инерция первого порядка
# k=3 → + инерция второго порядка
KRYLOV_PER_DOF = 3

println("\nПостроение базиса Крылова ($n_b граничных DOF × $KRYLOV_PER_DOF векторов = $(n_b*KRYLOV_PER_DOF) итого)...")

K_dense = Matrix(K)
M_dense = Matrix(M)

# Факторизуем один раз — используем для всех итераций
K_fact = factorize(K_dense)

cols = Vector{Vector{Float64}}()

for dof in boundary_dofs
    b_dof = zeros(n_total)
    b_dof[dof] = 1.0

    v = K_fact \ b_dof  # v₁ = K⁻¹b, статический отклик

    for k in 1:KRYLOV_PER_DOF
        # Грам-Шмидт: убираем компоненты вдоль уже добавленных векторов
        for u in cols
            v = v - dot(v, u) * u
        end

        nv = norm(v)
        nv < 1e-12 && break  # вектор линейно зависим — подпространство исчерпано

        v = v / nv
        push!(cols, copy(v))

        k < KRYLOV_PER_DOF && (v = K_fact \ (M_dense * v))  # следующий: K⁻¹Mv
    end
end

V   = hcat(cols...)
n_V = size(V, 2)
println("  Размер базиса: $n_V векторов")
println("  Ошибка ортогональности ‖VᵀV - I‖ = $(round(norm(V'*V - I), digits=8))")

# Проекция: K_red = VᵀKV, M_red = VᵀMV
println("\nРедукция матриц...")
K_kr = Symmetric(V' * K_dense * V)
M_kr = Symmetric(V' * M_dense * V)

n_red = size(K_kr, 1)
println("Размер редуцированной системы: $n_red DOF  (редукция $(round(100*(1-n_red/n_total), digits=1))%)")

# Собственные моды полной модели для сравнения
NUM_MODES = 10
println("\nМоды полной модели...")

λ_fb, Φ_fb = eigs(M, K; nev=NUM_MODES+6, which=:LM)
order_fb   = sortperm(real.(λ_fb), rev=true)
λ_fb       = real.(λ_fb)[order_fb]
Φ_fb       = real.(Φ_fb)[:, order_fb]

ω2_fb   = 1.0 ./ abs.(λ_fb)
mask_fb = ω2_fb .> (2π)^2  # убираем жёсткотельные моды
ω2_fb   = ω2_fb[mask_fb]
Φ_fb    = Φ_fb[:, mask_fb]

# Нормировка по массе: φᵀMφ = 1
for i in axes(Φ_fb, 2)
    mn = sqrt(abs(Φ_fb[:,i]' * M_dense * Φ_fb[:,i]))
    mn > 1e-12 && (Φ_fb[:,i] ./= mn)
end

freq_full = sqrt.(ω2_fb) ./ (2π)

# Собственные моды редуцированной модели
println("Моды редуцированной модели...")
λ_kr, Φ_kr = eigen(K_kr, M_kr)
order_kr   = sortperm(real.(λ_kr))
λ_kr       = real.(λ_kr)[order_kr]
Φ_kr       = real.(Φ_kr)[:, order_kr]

ω2_kr   = λ_kr
mask_kr = ω2_kr .> (2π)^2
ω2_kr   = ω2_kr[mask_kr]
Φ_kr    = Φ_kr[:, mask_kr]

freq_kr = sqrt.(ω2_kr) ./ (2π)

# Возврат мод в физическое пространство: φ_phys = V·φ_kr
Φ_kr_phys = V * Φ_kr

for i in axes(Φ_kr_phys, 2)
    mn = sqrt(abs(Φ_kr_phys[:,i]' * M_dense * Φ_kr_phys[:,i]))
    mn > 1e-12 && (Φ_kr_phys[:,i] ./= mn)
end

# Сравнение частот
n_cmp = min(length(freq_full), length(freq_kr), NUM_MODES)
println("\n" * "="^62)
println("        СРАВНЕНИЕ ЧАСТОТ (Гц)")
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

# FRF методом модальной суперпозиции
println("\nПостроение FRF...")

ζ     = 0.02
f_min = 50.0
f_max = freq_full[end] * 1.3
freqs = range(f_min, f_max, length=3000)
ωs    = 2π .* collect(freqs)

H_full = zeros(ComplexF64, length(ωs))
for (k, ω) in enumerate(ωs)
    for i in eachindex(ω2_fb)
        ωi = sqrt(ω2_fb[i])
        H_full[k] += Φ_fb[DOF_EXCITE, i] * Φ_fb[DOF_RESPOND, i] /
                     (ω2_fb[i] - ω^2 + 2im*ζ*ωi*ω)
    end
end

H_kr = zeros(ComplexF64, length(ωs))
for (k, ω) in enumerate(ωs)
    for i in eachindex(ω2_kr)
        ωi = sqrt(ω2_kr[i])
        H_kr[k] += Φ_kr_phys[DOF_EXCITE, i] * Φ_kr_phys[DOF_RESPOND, i] /
                   (ω2_kr[i] - ω^2 + 2im*ζ*ωi*ω)
    end
end

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
