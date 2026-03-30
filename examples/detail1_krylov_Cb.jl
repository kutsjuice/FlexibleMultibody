using Pkg;
Pkg.activate(".");

using FlexibleMultibody
using Gridap, GridapGmsh
using Arpack
using LinearAlgebra
using SparseArrays
using Printf
using Plots

# ==============================================================================
#  СРАВНЕНИЕ МЕТОДОВ РЕДУКЦИИ: Craig-Bampton vs Krylov
#  Деталь 1 четырёхзвенника (U-образная скоба)
# ==============================================================================

DOF_EXCITE  = 1          # DOF возбуждения (1-й DOF первого RBE2)
DOF_RESPOND = 1          # DOF отклика
NUM_MODES   = 10         # число мод для сравнения
ζ           = 0.02       # коэффициент демпфирования

# ==============================================================================
#  Модель
# ==============================================================================
# Единицы: мм, тонны, секунды → МПа
#   E = 210 000 МПа = 2.1e5
#   ν = 0.3
#   ρ = 7850 кг/м³ = 7.85e-9 тонн/мм³
# ==============================================================================
mat = LinearElasticMaterial(2.1e5, 0.3, 7.85e-9)

# Центры RBE2 — середина между парами отверстий на каждом конце
# Конец A (z ≈ -6.6 мм): отверстия в x = ±8.517, y = 6.5
# Конец B (z ≈ -63.4 мм): отверстия в x = ±8.517, y = 6.5
rbe_info = [
    RBEInfo("hole_A", VectorValue(0.0, 6.5, -6.608)),
    RBEInfo("hole_B", VectorValue(0.0, 6.5, -63.392)),
]

println("Сборка модели из detail1_mesh.msh...")
fl_part  = FlexibleComponent("detail1_mesh.msh", mat; rbe_info)
K        = fl_part.mK
M        = fl_part.mM
n_total  = size(K, 1)
println("  Полная система: $n_total DOF")

# Разделение на граничные и внутренние DOF
boundary_dofs = Int[]
for (_, iface) in fl_part.interfaces
    append!(boundary_dofs, collect(iface.dofs))
end
sort!(boundary_dofs)
interior_dofs = setdiff(1:n_total, boundary_dofs)
n_b = length(boundary_dofs)
n_i = length(interior_dofs)
println("  Граничные DOF: $n_b,  Внутренние DOF: $n_i")

K_dense = Matrix(K)
M_dense = Matrix(M)

# ==============================================================================
#  CRAIG-BAMPTON
# ==============================================================================
println("\n── Craig-Bampton ──────────────────────────────────────────")
NUM_CB_MODES = 10

K_ii = K[interior_dofs, interior_dofs]
K_ib = K[interior_dofs, boundary_dofs]
M_ii = M[interior_dofs, interior_dofs]

println("  Моды ограничений (constraint modes)...")
Ψ = Matrix(-K_ii \ Matrix(K_ib))

println("  Фиксированные моды ($NUM_CB_MODES)...")
λ_fix, Φ_fix = eigs(M_ii, K_ii; nev=NUM_CB_MODES, which=:LM)
order_fix    = sortperm(real.(λ_fix), rev=true)
Φ_fix        = real.(Φ_fix)[:, order_fix]

T_cb = [Φ_fix                        Ψ                  ;
        zeros(n_b, NUM_CB_MODES)      Matrix(I, n_b, n_b)]

K_cb = Symmetric(T_cb' * K_dense * T_cb)
M_cb = Symmetric(T_cb' * M_dense * T_cb)
n_cb = size(K_cb, 1)
println("  Размер CB системы: $n_cb DOF  (редукция $(round(100*(1-n_cb/n_total),digits=1))%)")

λ_cb, Φ_cb = eigen(K_cb, M_cb)
order_cb   = sortperm(real.(λ_cb))
λ_cb       = real.(λ_cb)[order_cb];  Φ_cb = real.(Φ_cb)[:, order_cb]
ω2_cb      = filter(x -> x > (2π)^2, λ_cb)
Φ_cb       = Φ_cb[:, λ_cb .> (2π)^2]
Φ_cb_phys  = T_cb * Φ_cb
for i in axes(Φ_cb_phys, 2)
    mn = sqrt(abs(Φ_cb_phys[:,i]' * M_dense * Φ_cb_phys[:,i]))
    mn > 1e-12 && (Φ_cb_phys[:,i] ./= mn)
end
freq_cb = sqrt.(ω2_cb) ./ (2π)

# ==============================================================================
#  KRYLOV
# ==============================================================================
println("\n── Krylov ─────────────────────────────────────────────────")
KRYLOV_PER_DOF = 3

println("  Факторизация K...")
K_fact = factorize(K_dense)

cols = Vector{Vector{Float64}}()
for dof in boundary_dofs
    b_dof = zeros(n_total);  b_dof[dof] = 1.0
    v = K_fact \ b_dof
    for k in 1:KRYLOV_PER_DOF
        for u in cols;  v = v - dot(v, u) * u;  end
        nv = norm(v);  nv < 1e-12 && break
        v = v / nv;  push!(cols, copy(v))
        k < KRYLOV_PER_DOF && (v = K_fact \ (M_dense * v))
    end
end
V_kr  = hcat(cols...)
n_kr  = size(V_kr, 2)
println("  Размер Krylov системы: $n_kr DOF  (редукция $(round(100*(1-n_kr/n_total),digits=1))%)")

K_kr = Symmetric(V_kr' * K_dense * V_kr)
M_kr = Symmetric(V_kr' * M_dense * V_kr)

λ_kr, Φ_kr = eigen(K_kr, M_kr)
order_kr   = sortperm(real.(λ_kr))
λ_kr       = real.(λ_kr)[order_kr];  Φ_kr = real.(Φ_kr)[:, order_kr]
ω2_kr      = filter(x -> x > (2π)^2, λ_kr)
Φ_kr       = Φ_kr[:, λ_kr .> (2π)^2]
Φ_kr_phys  = V_kr * Φ_kr
for i in axes(Φ_kr_phys, 2)
    mn = sqrt(abs(Φ_kr_phys[:,i]' * M_dense * Φ_kr_phys[:,i]))
    mn > 1e-12 && (Φ_kr_phys[:,i] ./= mn)
end
freq_kr = sqrt.(ω2_kr) ./ (2π)

# ==============================================================================
#  ПОЛНАЯ МОДЕЛЬ (моды для FRF и таблицы)
# ==============================================================================
println("\n── Полная модель ──────────────────────────────────────────")
λ_fb, Φ_fb = eigs(M, K; nev=NUM_MODES+6, which=:LM)
order_fb   = sortperm(real.(λ_fb), rev=true)
λ_fb       = real.(λ_fb)[order_fb];  Φ_fb = real.(Φ_fb)[:, order_fb]
ω2_fb      = 1.0 ./ abs.(λ_fb)
mask_fb    = ω2_fb .> (2π)^2
ω2_fb      = ω2_fb[mask_fb];  Φ_fb = Φ_fb[:, mask_fb]
for i in axes(Φ_fb, 2)
    mn = sqrt(abs(Φ_fb[:,i]' * M_dense * Φ_fb[:,i]))
    mn > 1e-12 && (Φ_fb[:,i] ./= mn)
end
freq_full = sqrt.(ω2_fb) ./ (2π)

# ==============================================================================
#  ТАБЛИЦА СРАВНЕНИЯ ЧАСТОТ
# ==============================================================================
n_cmp = min(length(freq_full), length(freq_cb), length(freq_kr), NUM_MODES)
println("\n" * "="^80)
println("               СРАВНЕНИЕ ЧАСТОТ (Гц) — Деталь 1 четырёхзвенника")
println("="^80)
@printf("%-6s  %-14s  %-14s  %-10s  %-14s  %-10s\n",
        "Мода", "Полная", "CB", "Ошибка CB%", "Krylov", "Ошибка Kr%")
println("-"^80)
for i in 1:n_cmp
    f0  = freq_full[i]
    fcb = freq_cb[i]
    fkr = freq_kr[i]
    @printf("%-6d  %-14.2f  %-14.2f  %-10.4f  %-14.2f  %-10.4f\n",
            i, f0, fcb, abs(f0-fcb)/f0*100, fkr, abs(f0-fkr)/f0*100)
end
println("="^80)
@printf("%-30s  %d DOF\n", "Craig-Bampton:", n_cb)
@printf("%-30s  %d DOF\n", "Krylov:", n_kr)
@printf("%-30s  %d DOF\n", "Полная модель:", n_total)

# ==============================================================================
#  FRF
# ==============================================================================
println("\nПостроение FRF...")
f_min = 50.0
f_max = freq_full[end] * 1.3
freqs = range(f_min, f_max, length=3000)
ωs    = 2π .* collect(freqs)

function compute_frf(ω2_modes, Φ_phys, dof_in, dof_out, ζ, ωs)
    H = zeros(ComplexF64, length(ωs))
    for (k, ω) in enumerate(ωs)
        for i in eachindex(ω2_modes)
            ωi = sqrt(ω2_modes[i])
            H[k] += Φ_phys[dof_in, i] * Φ_phys[dof_out, i] /
                    (ω2_modes[i] - ω^2 + 2im*ζ*ωi*ω)
        end
    end
    return H
end

H_full = compute_frf(ω2_fb,  Φ_fb,       DOF_EXCITE, DOF_RESPOND, ζ, ωs)
H_cb   = compute_frf(ω2_cb,  Φ_cb_phys,  DOF_EXCITE, DOF_RESPOND, ζ, ωs)
H_kr   = compute_frf(ω2_kr,  Φ_kr_phys,  DOF_EXCITE, DOF_RESPOND, ζ, ωs)

# ==============================================================================
#  ГРАФИК
# ==============================================================================
to_dB(H) = 20 .* log10.(abs.(H) .+ 1e-40)

p = plot(
    freqs, to_dB(H_full),
    label     = "Полная модель ($n_total DOF)",
    color     = :blue,
    linewidth = 2.5,
    xlabel    = "Частота (Гц)",
    ylabel    = "Амплитуда |H(ω)| (дБ)",
    title     = "FRF: Деталь 1 — Craig-Bampton vs Krylov  |  DOF $DOF_EXCITE → $DOF_RESPOND",
    legend    = :topright,
    grid      = true,
    size      = (1200, 600),
)

plot!(p, freqs, to_dB(H_cb),
    label     = "Craig-Bampton ($n_cb DOF)",
    color     = :red,
    linewidth = 1.5,
    linestyle = :dash
)

plot!(p, freqs, to_dB(H_kr),
    label     = "Krylov ($n_kr DOF)",
    color     = :green,
    linewidth = 2.5,
    linestyle = :dot
)

for f in freq_full
    vline!(p, [f], color=:gray, alpha=0.4, linewidth=1, label="")
end

savefig(p, "detail1_compare_methods.png")
println("График сохранён: detail1_compare_methods.png")
display(p)

# ==============================================================================
#  ЭКСПОРТ МОД В VTK (опционально)
# ==============================================================================
println("\nЭкспорт собственных форм в VTK...")
using FlexibleMultibody: write_solution_to_vtk

cellfields = Vector{Pair{String, Vector{Float64}}}(undef, length(freq_full))
for (i, f) in enumerate(freq_full)
    freq_str = round(f, digits=2)
    cellfields[i] = "freq=$(freq_str)Hz" => real.(Φ_fb[:, i])
end
write_solution_to_vtk(fl_part, "detail1_eigen_modes", cellfields)
println("VTK сохранён: detail1_eigen_modes.vtu")
