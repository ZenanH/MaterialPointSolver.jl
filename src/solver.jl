#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Main solver                                                                |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export mpmsolver!

include(joinpath(@__DIR__, "solver/datatransfer.jl"))
include(joinpath(@__DIR__, "solver/interpolation.jl"))
include(joinpath(@__DIR__, "solver/calculator.jl"))

# constitutive models
include(joinpath(@__DIR__, "solver/materials/linearelastic.jl"))
include(joinpath(@__DIR__, "solver/materials/druckerprager.jl"))

mpmsolver!(solver::Function, args...) = solver(args...)