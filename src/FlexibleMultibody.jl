module FlexibleMultibody

using Rotations, StaticArrays, SparseArrays, Gridap, GridapGmsh, LinearAlgebra, Gmsh
using Gridap.Geometry

# Экспортируем основные типы и функции
export FlexibleComponent, FlexibleInterface, RBEInfo, LinearElasticMaterial
# Включаем модули
include("utils.jl")
include("rbe2.jl")
include("flexible_body.jl")

end