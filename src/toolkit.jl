#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : MaterialPointSolver.jl                                                     |
|  Description: Some helper functions for MPMSolver.jl                                     |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

include(joinpath(@__DIR__, "toolkits/devicehelpfunc.jl"))
include(joinpath(@__DIR__, "toolkits/mpbasisfunc.jl"   ))
include(joinpath(@__DIR__, "toolkits/terminaltxt.jl"   ))
include(joinpath(@__DIR__, "toolkits/modelinfo.jl"     ))
include(joinpath(@__DIR__, "toolkits/warmup.jl"        ))
include(joinpath(@__DIR__, "toolkits/p2nindex.jl"      ))
include(joinpath(@__DIR__, "toolkits/hardwareinfo.jl"  ))
include(joinpath(@__DIR__, "toolkits/randomfield.jl"   ))