using Gridap
using GridapGmsh

using SparseArrays
using LinearAlgebra
using StaticArrays
function skew_symmetric(r)
    return SA[
        0.0 -r[3] r[2];
        r[3] 0.0 -r[1];
        -r[2] r[1] 0.0
    ]
end
function compute_rbe2_condensation_mat(slave_nodes, node_and_comp_to_dof, node_coords::Vector{VectorValue{3,Float64}}, rbe_coords::VectorValue{3,Float64}, num_dofs_glob)
    num_rbe_dofs = 6
    num_slave_dofs = 0
    for node in slave_nodes
        num_slave_dofs += length(node_and_comp_to_dof[node])
    end
    dofs_glob = 1:num_dofs_glob
    dofs_slave = Vector{Int64}(undef, num_slave_dofs)
    ind = 1
    for node in slave_nodes
        dofs_slave[ind] = node_and_comp_to_dof[node][1]
        dofs_slave[ind+1] = node_and_comp_to_dof[node][2]
        dofs_slave[ind+2] = node_and_comp_to_dof[node][3]
        ind += 3
    end
    dofs_indep = setdiff(dofs_glob, dofs_slave)
    RBE2_mat = spzeros(num_dofs_glob, num_dofs_glob - num_slave_dofs + num_rbe_dofs)
    RBE2_mat[dofs_indep, 1:length(dofs_indep)] = I(length(dofs_indep))

    for slave_node in slave_nodes
        if length(node_and_comp_to_dof[slave_node]) != 3
            error("Each slave node must have 3 DOFs.")
        end
        r = node_coords[slave_node] - rbe_coords
        RBE2_mat[Vector(node_and_comp_to_dof[slave_node]), end-5:end-3] = I(3)
        RBE2_mat[Vector(node_and_comp_to_dof[slave_node]), end-2:end] = skew_symmetric(r)
    end
    return RBE2_mat
end



model = GmshDiscreteModel("link_mesh.msh")
writevtk(model, "link")

const E = 2.1e5
const ν = 0.3
const λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
const μ = E / (2 * (1 + ν))

σ(ε) = λ * tr(ε) * one(ε) + 2 * μ * ε


reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, 1)
V = TestFESpace(model, reffe, conformity=:H1)

g1(x) = VectorValue(0.0, 0.0, 0.0) # boundary condition on the inner face of the left hole – displacement 0.0.
g2(x) = VectorValue(0.0, 0.0, 0.1) # boundary condition on the inner face of the right hole – displacement -0.1.


U = TrialFESpace(V) # create trial space with boundary conditions
EL_ORDER = 1
degree = EL_ORDER * 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# neumanntags = []#["loaded"]
# Γ = BoundaryTriangulation(model,tags=neumanntags)
# dΓ = Measure(Γ,degree)


# f(x) = VectorValue(0.0, 0.0, 0.0);

a(u, v) = ∫(ε(v) ⊙ (σ ∘ ε(u))) * dΩ
b(v) = 0; #∫( dot(v, f))*dΓ 


op = AffineFEOperator(a, b, U, V)

mK = get_matrix(op)
vF = get_vector(op)

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

    for (id,dof) in enumerate(dofs)
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

using Gridap.Geometry

function get_nodes_by_tag(model::DiscreteModel, tag::String)::Vector{Int64}


    labels = get_face_labeling(model)
    tag_id = get_tag_from_name(labels, tag)

    dim = 0 # we will look for nodes which has dimension of 0

    dface_to_entity = get_face_entity(labels, dim)

    dface_to_isontag = BitVector(undef, num_faces(labels, dim))

    tag_entities = get_tag_entities(labels, tag)

    for i in eachindex(dface_to_entity)
        buf = false
        for entity in tag_entities
            buf += dface_to_entity[i] == entity
        end
        dface_to_isontag[i] = buf
    end

    return findall(dface_to_isontag)

end

left_hole_nodes = get_nodes_by_tag(model, "left_hole")
right_hole_nodes = get_nodes_by_tag(model, "right_hole")

applydirichelet!(mK, vF, left_hole_nodes, trues(3), g1, model.grid.node_coordinates)

right_rbe_coords = VectorValue(0.7, 0.0, 0.0)
RBE2_mat = compute_rbe2_condensation_mat(right_hole_nodes, V.metadata.node_and_comp_to_dof, model.grid.node_coordinates, right_rbe_coords, V.nfree)
rbe2_dofs = (size(RBE2_mat, 2) - 5):size(RBE2_mat, 2)
mK = RBE2_mat' * mK * RBE2_mat
vF = RBE2_mat' * vF
applydirichelet!(mK, vF, rbe2_dofs[[1,2,3]], [0,0, -0.1])

x0 = mK \ vF
x0 = RBE2_mat * x0
uh_lin = FEFunction(U, x0)

res_file = "results_new"

function mises(s)
    return 0.5 * sqrt((s[1, 1] - s[2, 2])^2 + (s[2, 2] - s[3, 3])^2 + (s[3, 3] - s[1, 1])^2 + 6 * (s[2, 3]^2 + s[3, 1]^2 + s[1, 2]^2))
end

writevtk(Ω, res_file, cellfields=["uh" => uh_lin, "sigma" => σ ∘ ε(uh_lin)])