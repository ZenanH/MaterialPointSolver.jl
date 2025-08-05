# Kernel functions

The main computational processes of the MPM solver, such as `P2G`, `G2P`, stress updating, etc., are accomplished through kernel functions written in `KernelAbstractions.jl`, which are independent of the backend. Writing these kernel functions requires some additional knowledge, and an introduction to them is beyond the scope of the documentation. Here, only some useful, step-by-step materials are provided for reference.

- [Book: The Art of HPC](https://github.com/VictorEijkhout/TheArtofHPC_pdfs)
- [Multi-Threading in Julia](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [CUDA Parallel Programming Model](https://docs.nvidia.com/cuda/cuda-c-programming-guide/)
- [CUDA.jl documentation](https://github.com/JuliaGPU/CUDA.jl)
- [`KernelAbstractions.jl` documentation](https://github.com/JuliaGPU/KernelAbstractions.jl)
- [ask ChatGPT to provide a code example 🤖](https://openai.com)

Here we mainly introduce some key points in writing kernel functions for `MaterialPointSolver.jl`, with most kernel functions located under the `/src/solvers/`. Taking the one-phase single-point MPM as an example, which is in the file `utils_OS.jl`, it includes most of the computation process within a time step (stress updates/constitutive equations are located in the `/src/materials`).

## Writing a kernel function

!!! note
    
    1) Kernel functions generally use `@inbounds=true` to prevent boundary checks for performance improvement, which requires careful examination of the code indices to ensure there are no out-of-bounds accesses.
    2) Currently, we have mainly developed kernel functions for one-phase single-point MPM and two-phase single-point MPM, but users can develop new algorithms on their own, including two/three-point methods.

A basic kernel function's writing and invocation are as follows:

```math
c = a + b
```

```julia
using KernelAbstractions

@kernel inbounds=true function kernel1!(a, b, c)
    ix = @index(Global)
    if ix ≤ length(a)
        c[ix] = a[ix] + b[ix]
    end
end

a = rand(10)
b = rand(10)
c = zeros(10)

kernel1!(CPU())(ndrange=10, a, b, c)
```

## Kernel performance

The code of the kernel function almost determines the performance of the solver. If performance investigation is needed, we recommend testing these kernel functions one by one. A useful example can be found in `/test/c2_perf/t1_wallclock_time`.

!!! note

    In Julia, especially when measuring computation time on non-CPU hardware accelerators, it is recommended to use [`BenchmarkTools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl) and to consider synchronization time. On NVIDIA platforms, Nsight System/Compute can be used for profiling, especially in multi-GPU situations.

## Atomic operations

In the `P2G` process, we need to map the physical information of particles to the related nodes of the background grid. If we parallelize based on particles, it can lead to data races because different particle threads will try to write to the same node memory address. To avoid this issue, we use atomic locks to ensure thread safety. In our paper, we performed a performance analysis of the solver and found that the performance overhead of atomic locks has gradually decreased on recent GPU architectures, along with design variations that facilitate hardware. To find a balance between code maintainability and performance, we chose this particle-loop mode based on atomic locks.

## Basis function optimization

In our paper, we found that the solver spends most of its time calculating basis functions. To address this, we used shared memory for basis function calculations in uGIMP to enhance performance. Currently, this feature is only available on NVIDIA platform GPUs. Additionally, we offer an option for users to define a global computation precision of `FP64`, but during shape function calculations, we convert this computation to `FP32` and finally promote it to FP64 to significantly improve speed.

!!! danger

    This optional setting of partial `FP64-FP32` will cause the calculation results to differ from those of full `FP64`.
   
!!! tip

    This option has a very noticeable effect on consumer-grade graphics cards, as the `FP32` throughput on consumer platforms is significantly higher than `FP64`.