using Rotations, StaticArrays, SparseArrays, Gridap, GridapGmsh
struct FlexibleInterface
    coords::SVector{3,Float64}
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
    model::UnstructuredDiscreteModel
    mK::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}
    mM::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}
    mT::SparseMatrixCSC{Float64,Int64}
    interfaces::Dict{String,FlexibleInterface}
end

function FlexibleComponent(mesh::String, material::LinearElasticMaterial; rbe_labels = [])
    model = GmshDiscreteModel("link_mesh.msh")
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
    blf_mass(u, v) = material.density * ∫(v⋅u)dΩ
    # b(v) = 0

    # op = AffineFEOperator(blf_stiff, b, U, V)
    mK = assemble_matrix(blf_stiff, V, U)
    mM = assemble_matrix(blf_mass, V, U)
    mT = spzeros(size(mK)) + I(size(mK, 1))


    #TODO: add automatic rbe condensation



    return FlexibleComponent(model, mK, mM, mT, Dict{String, FlexibleInterface}())
end

mat = LinearElasticMaterial(2.1e11, 0.3, 7.85e3)
FlexibleComponent("link_mesh.msh", mat)