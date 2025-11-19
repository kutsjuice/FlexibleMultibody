function write_solution_to_vtk(comp::FlexibleComponent, filename::String, solutions::Vector{Pair{String, Vector{Float64}}})
    #TODO: generalize for additional FEFunction fields (σ(mode_shape), etc.)
    trian = Triangulation(comp.model)
    cellfields = Dict{String, Any}()
    for (i,sol) in enumerate(solutions)
        name = "sol #$i - $(sol[1])"
        mode_shape = FEFunction(comp.trial, comp.mT*sol[2])
        cellfields[name] = mode_shape
    end
    writevtk(trian, filename; cellfields)
end