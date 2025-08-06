module MaterialPointSolverCUDAExt

using BenchmarkTools
using CUDA
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device
CUDA.allowscalar(false) # disable scalar operation in GPU

dev_backend(::Val{:cuda}) = CUDABackend()

host2device(::CUDABackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(CuArray, grid), KAupload(CuArray, mpts)

end