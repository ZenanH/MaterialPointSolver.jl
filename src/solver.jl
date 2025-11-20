export mpmsolver!

include(joinpath(@__DIR__, "solver/datatransfer.jl"))
include(joinpath(@__DIR__, "solver/shapefuncs.jl"))
include(joinpath(@__DIR__, "solver/materials.jl"))
include(joinpath(@__DIR__, "solver/calculator.jl"))
include(joinpath(@__DIR__, "solver/randomfields.jl"))
include(joinpath(@__DIR__, "solver/timestep.jl"))

mpmsolver!(solver::Function, args...) = solver(args...)