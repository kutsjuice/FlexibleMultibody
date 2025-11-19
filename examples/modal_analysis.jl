using FlexibleMultibody
using Gridap, GridapGmsh
using FlexibleMultibody: applydirichelet!, write_solution_to_vtk
using Arpack
mat = LinearElasticMaterial(2.1e11, 0.3, 7.85e3)


rbe_info = [
    RBEInfo("left_hole", VectorValue(0.0, 0.0, 0.0)),
    RBEInfo("right_hole", VectorValue(0.7, 0.0, 0.0)),
]

fl_rod = FlexibleComponent("link_mesh.msh", mat; rbe_info);
f = zeros(size(fl_rod.mK, 1))

left_interface = fl_rod.interfaces["left_hole"]
applydirichelet!(fl_rod, left_interface.dofs)

# right_interface = fl_rod.interfaces["right_hole"]
# disp_vals = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0]
# applydirichelet!(fl_rod.mK, f, right_interface.dofs, disp_vals)
num_modes = 10
λ, Φ = eigs( fl_rod.mM, fl_rod.mK; nev=num_modes, which=:LM)

1 ./ λ
cellfields = Vector{Pair{String, Vector{Float64}}}(undef, length(λ))
for (i, λᵢ) in enumerate(λ)
    freq = round( sqrt(1 / abs(λᵢ)) / 2 /pi , digits=2)
    cellfields[i] = "freq=$(freq)" => Φ[:,i]
end

write_solution_to_vtk(fl_rod, "eigen_modes", cellfields)
