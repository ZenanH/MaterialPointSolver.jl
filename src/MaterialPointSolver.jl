module MaterialPointSolver

using Adapt, BenchmarkTools, Dates, DelimitedFiles, HDF5, KernelAbstractions, Printf, 
      ProgressMeter

import Adapt.adapt as KAupload
import Adapt.@adapt_structure as @KAadapt
import KernelAbstractions.synchronize as KAsync
import KernelAbstractions.Extras: @unroll as @KAunroll
macro Σ(expr)
    esc(quote
        KernelAbstractions.@atomic :monotonic $expr
    end)
end

export KAupload, KAsync
export @Σ, @KAadapt, @KAunroll

include(joinpath(@__DIR__, "type/conf.jl"))
include(joinpath(@__DIR__, "type/grid.jl"))
include(joinpath(@__DIR__, "type/particle.jl"))
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "solver.jl"))

end