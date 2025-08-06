module MaterialPointSolverAMDGPUExt

using BenchmarkTools
using AMDGPU
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device

dev_backend(::Val{:rocm}) = ROCBackend()

host2device(::ROCBackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(ROCArray, grid), KAupload(ROCArray, mpts)

end