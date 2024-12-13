module MaterialPointSolverMetalExt

using BenchmarkTools
using Metal
using KernelAbstractions
using Printf
using MaterialPointSolver

# rewrite with CUDA
import MaterialPointSolver: host2device, device2host!, clean_device!, Tpeak, getBackend,
       warmup, grf_gc!, grf_ec!, getArray

#=-------------------------------------------------------------------------------------=#
# here we need something like CUDA.allowscalar(false) # disable scalar operation in GPU #
#--------------------------------------------------------------------------------------=#
include(joinpath(@__DIR__, "MetalExt/devicehelpfunc_metal.jl"))
include(joinpath(@__DIR__, "MetalExt/warmup_metal.jl"        ))
include(joinpath(@__DIR__, "MetalExt/randomfield_metal.jl"   ))

end