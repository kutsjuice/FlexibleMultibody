function write_solution_to_vtk(comp::FlexibleComponent, filename::String, solutions::Dict{String, Vector{Float64}})
    #TODO: generalize for additional FEFunction fields (σ(mode_shape), etc.)
    trian = Triangulation(comp.model)
    cellfields = Dict{String, Any}()
    for (name, sol) in solutions
        mode_shape = FEFunction(comp.trial, comp.mT*sol)
        cellfields[name] = mode_shape
    end
    writevtk(trian, filename; cellfields)
end