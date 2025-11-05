using Gmsh


w = 0.1
h = 0.08
l = 0.7

hole_diam = 0.05
try
    gmsh.finalize()
catch
end
gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("link_mesh")
lc = .03


#              
#       p4-------------p3
#     p5                 p6
#       p1-------------p2
#

p1 = gmsh.model.occ.add_point(0.0, -w/2, 0.0, lc)
p2 = gmsh.model.occ.add_point(l, -w/2, 0.0, lc)
p3 = gmsh.model.occ.add_point(l, w/2, 0.0, lc)
p4 = gmsh.model.occ.add_point(0.0, w/2, 0.0, lc)

p5 = gmsh.model.occ.add_point(0.0, 0.0, 0.0, lc)
p6 = gmsh.model.occ.add_point(l, 0.0, 0.0, lc)

p7 = gmsh.model.occ.add_point(hole_diam/2, 0.0, 0.0, lc)
p8 = gmsh.model.occ.add_point(-hole_diam/2, 0.0, 0.0, lc)

p9 = gmsh.model.occ.add_point(l - hole_diam/2, 0.0, 0.0, lc)
p10 = gmsh.model.occ.add_point(l + hole_diam/2, 0.0, 0.0, lc)


c1 = gmsh.model.occ.add_circle_arc(p1, p5, p4,-1,true)
c2 = gmsh.model.occ.add_circle_arc(p3, p6, p2,-1,true)

c3 = gmsh.model.occ.add_circle_arc(p8,p5,p7,-1,true)
c4 = gmsh.model.occ.add_circle_arc(p7,p5,p8,-1,true)
c5 = gmsh.model.occ.add_circle_arc(p9,p6,p10,-1,true)
c6 = gmsh.model.occ.add_circle_arc(p10,p6,p9,-1,true)

l1 = gmsh.model.occ.add_line(p1, p2)
l2 = gmsh.model.occ.add_line(p3, p4)

curve_loop = gmsh.model.occ.add_curve_loop([-c1, l1, -c2, l2])
hole1 = gmsh.model.occ.add_curve_loop([c3, c4])
hole2 = gmsh.model.occ.add_curve_loop([c5, c6])

surface = gmsh.model.occ.add_plane_surface([curve_loop, hole1, hole2])

volume = gmsh.model.occ.extrude([(2, surface)], 0, 0, h)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)

# gmsh.fltk.run()
gmsh.write("link_mesh.msh")
gmsh.finalize()