# Backend-agnostic solver

***MaterialPointSolver.jl*** is a high-performance solver that is backend-agnostic, meaning we only need to write a unified codebase that can run on different accelerator backends. The main computational process of the solver relies on kernel functions implemented by [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), so users do not need to worry about the data migration process from the CPU to other hardware backends. 

Users should use ***MaterialPointSolver.jl***, and then selects the backend based on your hardware, for example:

```julia-repl
julia> using MaterialPointSolver
julia> using CUDA

# using AMDGPU # for AMD   gpu
# using Metal  # for Apple gpu
# using oneAPI # for Intel gpu
```

!!! note

    In addition to selecting the appropriate backend package, users must also specify the backend when instantiating `Config`, which will be explained in later sections.