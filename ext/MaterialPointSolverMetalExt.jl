module MaterialPointSolverMetalExt

using BenchmarkTools
using Metal
using KernelAbstractions
using Printf
using MaterialPointSolver

import MaterialPointSolver: dev_backend, host2device

dev_backend(::Val{:metal}) = MetalBackend()

function host2device(::MetalBackend, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2}
    dev_grid = KAupload(MtlArray, grid)
    dev_mpts = KAupload(MtlArray, mpts)
    memsize = memorysize(dev_grid) + memorysize(dev_mpts)
    content = "uploaded $(@sprintf("%.2f", memsize)) GiB data to Metal"
    println("\e[1;32m[▲ I/O:\e[0m \e[0;32m$(content)\e[0m")
    return dev_grid, dev_mpts
end

host2device(::MetalBackend, data) = KAupload(MtlArray, data)
end