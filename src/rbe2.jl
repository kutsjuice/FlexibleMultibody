struct RBE2Element
    slave_nodes::Vector{Int}
    master_coords::VectorValue{3,Float64}
end

struct RBEInfo
    tag::String
    coords::VectorValue{3,Float64}
end


function compute_rbe2_condensation_mat(
    rbe2_elements::Vector{RBE2Element},
    node_and_comp_to_dof::Vector{VectorValue{3,Int32}}, # например, [[1,2,3], [4,5,6], ...]
    node_coords::Vector{VectorValue{3,Float64}},
    num_dofs_glob::Int
)
    num_rbe_dofs_total = 6 * length(rbe2_elements)
    num_slave_dofs_total = 0

    for rbe in rbe2_elements
        for node in rbe.slave_nodes
            num_slave_dofs_total += length(node_and_comp_to_dof[node])
        end
    end

    dofs_glob = 1:num_dofs_glob
    dofs_slave_all = Vector{Int64}()

    for rbe in rbe2_elements
        for node in rbe.slave_nodes
            append!(dofs_slave_all, node_and_comp_to_dof[node])
        end
    end

    dofs_indep = setdiff(dofs_glob, dofs_slave_all)

    RBE2_mat = spzeros(num_dofs_glob, num_dofs_glob - num_slave_dofs_total + num_rbe_dofs_total)

    # Просто копируем независимые DOF
    RBE2_mat[dofs_indep, 1:length(dofs_indep)] = I(length(dofs_indep))

    # Заполняем зависимости для каждого RBE2
    current_master_dof_start = length(dofs_indep) + 1

    for rbe in rbe2_elements
        master_coords = rbe.master_coords
        for slave_node in rbe.slave_nodes
            if length(node_and_comp_to_dof[slave_node]) != 3
                error("Each slave node must have 3 DOFs.")
            end

            r = node_coords[slave_node] - master_coords
            slave_dofs = node_and_comp_to_dof[slave_node]

            # 3 поступательных DOF master узла
            RBE2_mat[Vector(slave_dofs), current_master_dof_start:current_master_dof_start+2] = I(3)

            # 3 вращательных DOF master узла
            RBE2_mat[Vector(slave_dofs), current_master_dof_start+3:current_master_dof_start+5] = -skew_symmetric(r)
        end
        current_master_dof_start += 6
    end

    return RBE2_mat
end