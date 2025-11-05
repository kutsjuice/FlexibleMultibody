using Gridap
using GridapGmsh

model = GmshDiscreteModel("link_mesh.msh")
writevtk(model, "link")