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



r = VectorValue(0.35, 0.0, 0.0)
th = VectorValue(0.0, 0.0, 1.0)
skew_R = [
    0.0 -r[3] r[2];
    r[3] 0.0 -r[1];
    -r[2] r[1] 0.0
]
cross(r, th)

[I(3) skew_symmetric(r)]