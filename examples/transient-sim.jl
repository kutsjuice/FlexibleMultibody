using Gridap
using Gmsh
using GridapGmsh
gmsh.initialize()
gmsh.model.add("transient")
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.occ.importShapes("examples/model.stp")

gmsh.model.occ.synchronize()


gmsh.option.setNumber("Mesh.MeshSizeMax", 2.5)
EPS = 1e-5

entities = gmsh.model.occ.get_entities(3)
bounding_box = gmsh.model.occ.get_bounding_box(entities[1]...)

xmin, ymin, zmin, xmax, ymax, zmax = bounding_box

xhin = [9.75847, -9.75847]
zhin = [0.0, -70.0]
refinement = []
for i in 1:2
    for j in 1:2

        push!(refinement, gmsh.model.mesh.field.add("Box"))
        gmsh.model.mesh.field.setNumber(refinement[end], "VIn", 0.1)
        gmsh.model.mesh.field.setNumber(refinement[end], "VOut", 2.5)
        gmsh.model.mesh.field.setNumber(refinement[end], "Thickness", 15.0)
        gmsh.model.mesh.field.setNumber(refinement[end], "XMin", xhin[i] - 0.25)
        gmsh.model.mesh.field.setNumber(refinement[end], "XMax", xhin[i] + 0.25)
        gmsh.model.mesh.field.setNumber(refinement[end], "YMin", ymin - EPS)
        gmsh.model.mesh.field.setNumber(refinement[end], "YMax", ymax + EPS)
        gmsh.model.mesh.field.setNumber(refinement[end], "ZMin", zhin[j] - 0.5)
        gmsh.model.mesh.field.setNumber(refinement[end], "ZMax", zhin[j] + 0.5)

    end
end

min_field = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", refinement)
gmsh.model.mesh.field.setAsBackgroundMesh(min_field)


bot_bnd = gmsh.model.getEntitiesInBoundingBox(
    bounding_box[1] - EPS, -Inf, -Inf,
    bounding_box[1] + EPS, +Inf, +Inf
)

top_bnd = gmsh.model.getEntitiesInBoundingBox(
    bounding_box[4] - EPS, -Inf, -Inf,
    bounding_box[4] + EPS, +Inf, +Inf
)

for dim in (0, 1, 2)
    gmsh.model.addPhysicalGroup(dim, [dimtag[2] for dimtag in bot_bnd if dimtag[1] == dim], -1, "bottom")
    gmsh.model.addPhysicalGroup(dim, [dimtag[2] for dimtag in top_bnd if dimtag[1] == dim], -1, "top")
end

gmsh.model.occ.synchronize()
gmsh.model.addPhysicalGroup(3, [entities[1][2]], -1, "domain")
gmsh.model.mesh.generate(3)

gmsh.fltk.run()

msh_file = "examples/transient.msh"



gmsh.write(msh_file)

gmsh.finalize()
##
model = GmshDiscreteModel(msh_file)
writevtk(model, "examples/model")

order = 1