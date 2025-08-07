module MaterialPointSolveroneAPIExt

using BenchmarkTools
using oneAPI
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device

dev_backend(::Val{:oneapi}) = oneAPIBackend()

host2device(::oneAPIBackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(oneArray, grid), KAupload(oneArray, mpts)
host2device(::oneAPIBackend, data) = KAupload(oneArray, data)

end