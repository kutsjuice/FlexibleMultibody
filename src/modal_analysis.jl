function compute_modes(comp::FlexibleComponent, num_modes::Int; fixed_interfaces::Vector{String}=String[])
    # Создаем копии матриц для модального анализа
    # copying matrices for modal analysis
    K_modal = copy(comp.mK)
    M_modal = copy(comp.mM)
    
    # Применяем граничные условия к закрепленным интерфейсам
    for interface_name in fixed_interfaces
        if haskey(comp.interfaces, interface_name)
            interface = comp.interfaces[interface_name]
            # Закрепляем все 6 степеней свободы интерфейса
            applydirichelet!(K_modal, zeros(size(K_modal, 1)), interface.dofs, zeros(6))
            applydirichelet!(M_modal, zeros(size(M_modal, 1)), interface.dofs, zeros(6))
        else
            @warn "Interface $interface_name not found"
        end
    end
    
    # Решаем проблему собственных значений
    # Для небольших задач используем встроенный решатель, для больших - ARPACK
    λ, Φ_condensed = eigs(K_modal, M_modal; nev=num_modes, which=:SM)
        
    # Берем только запрошенное количество мод
    λ = λ[1:num_modes]
    Φ_condensed = Φ_condensed[:, 1:num_modes]
    
    # Преобразуем в полное пространство
    Φ_full = comp.mT * Φ_condensed
    
    return λ, Φ_full
end

function write_modes_to_vtk(comp::FlexibleComponent, Φ_full::Matrix, λ::Vector, prefix::String="mode")
    for i in 1:size(Φ_full, 2)
        # Создаем FE функцию для i-й формы колебаний
        mode_shape = FEFunction(comp.U, Φ_full[:, i])
        
        # Вычисляем частоту в Гц
        freq_hz = sqrt(λ[i]) / (2π)
        
        # Записываем в VTK
        writevtk(comp.model, "$(prefix)_$(i)_$(round(freq_hz, digits=2))Hz", 
                cellfields=["mode_shape" => mode_shape])
    end
end