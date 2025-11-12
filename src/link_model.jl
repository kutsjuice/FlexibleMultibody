using Gridap
using GridapGmsh
using Gridap.Geometry
using SparseArrays
using LinearAlgebra
using StaticArrays


function applydirichelet!(
    mK::AbstractMatrix{<:Number},
    vF::AbstractVector{<:Number},
    indices::AbstractVector{<:Integer},
    dirsmask::BitVector,
    func::Function,
    nodes_coords::Vector{VectorValue{3,Float64}}
)

    dims = size(dirsmask, 1)
    dirs = findall(dirsmask)
    for ind in indices

        p = nodes_coords[ind]
        dofsvalue = func(p)

        for dir in dirs
            dof = dims * (ind - 1) + dir
            mK[dof, dof] = 1
            vF[dof] = dofsvalue[dir]

            for k in eachindex(vF)
                if dof != k
                    if mK[dof, k] != 0
                        vF[k] -= dofsvalue[dir] * mK[dof, k]
                        mK[dof, k] = mK[k, dof] = 0
                    end
                end
            end
        end
    end

end

function applydirichelet!(
    mK::AbstractMatrix{<:Number},
    vF::AbstractVector{<:Number},
    dofs::AbstractVector{<:Integer},
    dofsvalue::AbstractVector{<:Number}
)

    for (id, dof) in enumerate(dofs)
        mK[dof, dof] = 1
        vF[dof] = dofsvalue[id]
        for k in eachindex(vF)
            if dof != k
                if mK[dof, k] != 0
                    vF[k] -= dofsvalue[id] * mK[dof, k]
                    mK[dof, k] = mK[k, dof] = 0
                end
            end
        end
    end

end

function mises(s)
    return 0.5 * sqrt((s[1, 1] - s[2, 2])^2 + (s[2, 2] - s[3, 3])^2 + (s[3, 3] - s[1, 1])^2 + 6 * (s[2, 3]^2 + s[3, 1]^2 + s[1, 2]^2))
end

model = GmshDiscreteModel("link_mesh.msh")
writevtk(model, "link")

# Материал
E = 2.1e5
ν = 0.3
λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
μ = E / (2 * (1 + ν))

σ(ε) = λ * tr(ε) * one(ε) + 2 * μ * ε

# Пространства
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, 1)
V = TestFESpace(model, reffe, conformity=:H1)
U = TrialFESpace(V) # без граничных условий на узлах

# Интегрирование
EL_ORDER = 1
degree = EL_ORDER * 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Билинейная и линейная формы
a(u, v) = ∫(ε(v) ⊙ (σ ∘ ε(u))) * dΩ
b(v) = 0

op = AffineFEOperator(a, b, U, V)
mK = get_matrix(op)
vF = get_vector(op)

# Получаем узлы
left_hole_nodes = get_nodes_by_tag(model, "left_hole")
right_hole_nodes = get_nodes_by_tag(model, "right_hole")

# Определяем координаты центров RBE2
left_rbe_coords = VectorValue(0.0, 0.0, 0.0)  # Пример координат для левой дыры
right_rbe_coords = VectorValue(0.7, 0.0, 0.0) # Пример координат для правой дыры

# Создаем вектор RBE2 элементов
rbe2_elements = [
    RBE2Element(left_hole_nodes, left_rbe_coords),
    RBE2Element(right_hole_nodes, right_rbe_coords)
]

# Вычисляем матрицу конденсации
RBE2_mat = compute_rbe2_condensation_mat(
    rbe2_elements,
    V.metadata.node_and_comp_to_dof,
    model.grid.node_coordinates,
    V.nfree
)

# Преобразуем матрицу и вектор
mK = RBE2_mat' * mK * RBE2_mat
vF = RBE2_mat' * vF


# Граничные условия теперь задаются на DOF master узлов
# Предположим, что в RBE2_mat первые 6 DOF соответствуют первому RBE2 (левому), следующие 6 — второму (правому)
left_rbe_dofs = (size(RBE2_mat, 2)-11):(size(RBE2_mat, 2)-6) # первые 6 DOF в RBE2_mat (левый RBE2)
right_rbe_dofs = (size(RBE2_mat, 2)-5):size(RBE2_mat, 2)      # последние 6 DOF в RBE2_mat (правый RBE2)
n̂ = size(RBE2_mat, 2)
n_rbe = 2*6
n_free = size(RBE2_mat, 2) - n_rbe
rbe_dofs = [(n_free+1+(i-1)*6:n_free+i*6) for i in 1:2]
# g1(x) = VectorValue(0.0, 0.0, 0.0) -> фиксируем все 6 DOF левого RBE2
# g2(x) = VectorValue(0.0, 0.0, 0.1) -> задаем смещение 0.1 по Z на правом RBE2

# Применяем условия: левый RBE2 полностью фиксирован (0)
applydirichelet!(mK, vF, left_rbe_dofs[1:3], zeros(6)[1:3])

# Применяем условия: правый RBE2 - смещение 0.1 по Z (индекс 3)
# Устанавливаем смещение: 0 по X, 0 по Y, 0.1 по Z, 0 по вращениям
disp_vals = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0]

applydirichelet!(mK, vF, right_rbe_dofs, disp_vals)

# Решение
x0 = mK \ vF
x0 = RBE2_mat * x0
uh_lin = FEFunction(U, x0)

# Запись результата
res_file = "results_new"
writevtk(Triangulation(model), res_file, cellfields=["uh" => uh_lin, "sigma" => σ ∘ ε(uh_lin)])