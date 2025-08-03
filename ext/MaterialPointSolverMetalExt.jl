#==========================================================================================+
|                OpenMPM.jl: High-performance MPM Solver for Geomechanics                  |
+------------------------------------------------------------------------------------------+
|  Description: Metal extension for OpenMPM.jl                                             |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

module MaterialPointSolverMetalExt

using BenchmarkTools, Metal, KernelAbstractions, Printf, MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device

include(joinpath(@__DIR__, "MetalExt/datatransfer.jl"))

end