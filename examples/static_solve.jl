using FlexibleMultibody
using Gridap, GridapGmsh
using FlexibleMultibody: applydirichelet!, write_solution_to_vtk
mat = LinearElasticMaterial(2.1e11, 0.3, 7.85e3)

rbe_info = [
    RBEInfo("left_hole", VectorValue(0.0, 0.0, 0.0)),
    RBEInfo("right_hole", VectorValue(0.7, 0.0, 0.0)),
]

fl_rod = FlexibleComponent("link_mesh.msh", mat; rbe_info);
f = zeros(size(fl_rod.mK, 1))

left_interface = fl_rod.interfaces["left_hole"]
applydirichelet!(fl_rod.mK, f, left_interface.dofs[1:3], zeros(3))

right_interface = fl_rod.interfaces["right_hole"]
disp_vals = [0.0, 0.0, 0.1, 0.0, 0.0, 0.0]
applydirichelet!(fl_rod.mK, f, right_interface.dofs, disp_vals)

u_condensed = fl_rod.mK \ f

cellfields = Dict(
    "displacement" => u_condensed,
)

write_solution_to_vtk(fl_rod, "static_solution", cellfields)
