module MaterialPointSolverCUDAExt

using BenchmarkTools
using CUDA
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device
CUDA.allowscalar(false) # disable scalar operation in GPU

dev_backend(::Val{:cuda}) = CUDABackend()

function host2device(::CUDABackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} 
    dev_grid = KAupload(CuArray, grid)
    dev_mpts = KAupload(CuArray, mpts)
    memsize = memorysize(dev_grid) + memorysize(dev_mpts)
    content = "uploaded $(@sprintf("%.2f", memsize)) GiB data to CUDA"
    println("\e[1;32m[▲ I/O:\e[0m \e[0;32m$(content)\e[0m")
    return dev_grid, dev_mpts
end

host2device(::CUDABackend, data) = KAupload(CuArray, data)

end