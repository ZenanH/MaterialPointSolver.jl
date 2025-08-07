module MaterialPointSolverMetalExt

using BenchmarkTools
using Metal
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device

dev_backend(::Val{:metal}) = MetalBackend()

host2device(::MetalBackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(MtlArray, grid), KAupload(MtlArray, mpts)
host2device(::MetalBackend, data) = KAupload(MtlArray, data)
end