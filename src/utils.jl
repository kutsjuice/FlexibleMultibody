

function skew_symmetric(r::VectorValue{3, Float64})
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

function applydirichelet!(
	mK::AbstractMatrix{<:Number},
	vF::AbstractVector{<:Number},
	dofs::AbstractVector{<:Integer},
	dofsvalue::AbstractVector{<:Number},
)
	for (id, dof) in enumerate(dofs)
		mK[dof, dof] = 1
		mM[dof, dof] = 1
		vF[dof] = dofsvalue[id]
		for k in eachindex(vF)
			if dof != k
				if mK[dof, k] != 0
					vF[k] -= dofsvalue[id] * mK[dof, k]
					mK[dof, k] = mK[k, dof] = 0
					mM[dof, k] = mM[k, dof] = 0
				end
			end
		end
	end
end

function applydirichelet!(
	comp::FlexibleComponent,
	dofs::AbstractVector{<:Integer},
)
	for (id, dof) in enumerate(dofs)
		comp.mK[dof, dof] = 1
		comp.mM[dof, dof] = 1
		for k in 1:size(comp.mK, 1)
			if dof != k
				comp.mK[dof, k] = comp.mK[k, dof] = 0
				comp.mM[dof, k] = comp.mM[k, dof] = 0
			end
		end
	end
end
