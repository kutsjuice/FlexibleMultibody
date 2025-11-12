using FlexibleMultibody
using Gridap, GridapGmsh
mat = LinearElasticMaterial(2.1e11, 0.3, 7.85e3)

rbe_info = [
    RBEInfo("left_hole", VectorValue(0.0, 0.0, 0.0)),
    RBEInfo("right_hole", VectorValue(0.7, 0.0, 0.0)),
]

fl_rod = FlexibleComponent("link_mesh.msh", mat; rbe_info)