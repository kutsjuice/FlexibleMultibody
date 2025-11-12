

function skew_symmetric(r::VectorValue{3,Float64})
    return SA[
        0.0 -r[3] r[2];
        r[3] 0.0 -r[1];
        -r[2] r[1] 0.0
    ]
end

function get_nodes_by_tag(model::DiscreteModel, tag::String)::Vector{Int64}


    labels = get_face_labeling(model)

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