# Backend-agnostic solver

MaterialPointSolver.jl is a high-performance solver that is backend-agnostic, meaning we only need to write a unified codebase that can run on different accelerator backends. The main computational process of the solver relies on kernel functions implemented by [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), so users do not need to worry about the data migration process from the CPU to other hardware backends. 

Users should use MaterialPointSolver.jl, and then selects the backend based on your hardware, for example:

```julia-repl
julia> using MaterialPointSolver
julia> using CUDA

# using AMDGPU # for AMD   gpu
# using Metal  # for Apple gpu
# using oneAPI # for Intel gpu
```

Taking the CUDA backend as an example, we define automatic conversion methods in `/ext/CUDAExt/devicehelpfunc_cuda.jl`: `host2device` for transferring data to the device, and `device2host!` for downloading data back to the CPU from the backend.

!!! note

    The reason this solution can succeed is that the structs' (including custom structs) fields only contain `scalar` and `array` types for the kernel function, and currently do not support other forms.

Although Julia has a garbage collection mechanism to automatically reclaim memory, we hope to clean up GPU memory usage promptly after each simulation. Currently, in `/ext/CUDAExt/devicehelpfunc_cuda.jl`, we use the function `clean_device!` to do this.

!!! note
    
    For all functions related to the accelerated backend that are not part of the MPM process, a dedicated implementation may be required in the `/ext` folder.

!!! todo

    Currently,
    1) Memory clean can only be achieved on the Nvidia platform, and other platforms are waiting for it to provide the relevant API.
    2) Metal on the Apple platform does not support FP64 computation precision.  
    3) Intel's GPU devices have not been tested, but the code implementation should be relatively easy. Ideally, it should only require filling in some parts of the `/ext` folder.
