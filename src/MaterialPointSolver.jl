#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Module file                                                                |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

module MaterialPointSolver

using Adapt, BenchmarkTools, Dates, DelimitedFiles, HDF5, KernelAbstractions, Logging, 
      PrecompileTools, Pkg, PrettyTables, Printf, SysInfo

import Adapt.adapt as KAupload
import KernelAbstractions.synchronize as KAsync
import KernelAbstractions.Extras: @unroll as @KAunroll
macro Σ(expr)
    esc(quote
        KernelAbstractions.@atomic :monotonic $expr
    end)
end

include(joinpath(@__DIR__, "type.jl"))
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "solver.jl"))

export @Σ, @KAunroll, KAsync, KAupload

quiet(f) = redirect_stdout(devnull) do
    redirect_stderr(devnull) do
        with_logger(NullLogger()) do
            f()
        end
    end
end

include(joinpath(@__DIR__, "precompile.jl"))

end