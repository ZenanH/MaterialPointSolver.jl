const MPMplugins = ["Mechanics"] #, "Saturated", "GlenNye", "FiniteStrain"] 
const lib_path = joinpath(@__DIR__, "../lib")

@inbounds for module_name in MPMplugins
    include(joinpath(lib_path, module_name, "src/main.jl"))
end