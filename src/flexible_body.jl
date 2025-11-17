
struct FlexibleInterface
    coords::VectorValue{3,Float64}
    rotation::RotMatrix{3,Float64}
    dofs::SVector{6,Int64}
end

struct LinearElasticMaterial
    young_modulus::Float64
    poisson_ratio::Float64
    density::Float64
end


function generate_stress_law(material::LinearElasticMaterial)
    # Материал
    E = material.young_modulus
    ν = material.poisson_ratio

    λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
    μ = E / (2 * (1 + ν))
    return (ε) -> λ * tr(ε) * one(ε) + 2 * μ * ε
end
mutable struct FlexibleComponent
    model::Gridap.Geometry.UnstructuredDiscreteModel
    trial::FESpace
    mK::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}
    mM::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}
    mT::SparseMatrixCSC{Float64,Int64}
    interfaces::Dict{String,FlexibleInterface}
end

function FlexibleComponent(mesh::String, material::LinearElasticMaterial; rbe_info::Vector{RBEInfo}=RBEInfo[])
    model = GmshDiscreteModel(mesh)
    writevtk(model, "link")

    σ = generate_stress_law(material)

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
    blf_stiff(u, v) = ∫(ε(v) ⊙ (σ ∘ ε(u))) * dΩ
    blf_mass(u, v) = material.density * ∫(v ⋅ u)dΩ
    # b(v) = 0

    # op = AffineFEOperator(blf_stiff, b, U, V)
    mK = assemble_matrix(blf_stiff, V, U)
    mM = assemble_matrix(blf_mass, V, U)

    rbe_nodes_groups = [get_nodes_by_tag(model, info.tag) for info in rbe_info]

    rbe_coords = [info.coords for info in rbe_info]

    rbe2_elements = [RBE2Element(rbe_nodes_groups[i], rbe_coords[i]) for i in eachindex(rbe_info)]

    RBE2_mat = compute_rbe2_condensation_mat(
        rbe2_elements,
        V.metadata.node_and_comp_to_dof,
        model.grid.node_coordinates,
        V.nfree
    )

    # Преобразуем матрицу и вектор
    mK = RBE2_mat' * mK * RBE2_mat

    n̂ = size(RBE2_mat, 2)
    n_rbe = length(rbe_info) * 6
    n_free = n̂ - n_rbe
    rbe_dofs = [(n_free+1+(i-1)*6:n_free+i*6) for i in eachindex(rbe_info)]

    interfaces = Dict{String,FlexibleInterface}()
    for (i, info) in enumerate(rbe_info)
        interfaces[info.tag] = FlexibleInterface(
            info.coords,
            RotMatrix{3,Float64}(I),
            rbe_dofs[i]
        )
    end
    return FlexibleComponent(model,U, mK, mM, RBE2_mat, interfaces)
end




