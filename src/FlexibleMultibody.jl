module FlexibleMultibody

using Rotations, StaticArrays, SparseArrays, Gridap, GridapGmsh, LinearAlgebra, Gmsh
using Arpack
using Gridap.Geometry

# Экспортируем основные типы и функции
export FlexibleComponent, FlexibleInterface, RBEInfo, LinearElasticMaterial
# Включаем модули
include("rbe2.jl")
include("flexible_body.jl")
include("modal_analysis.jl")
include("utils.jl")
include("export.jl")
end