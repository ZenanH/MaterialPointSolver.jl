## Introduction

This package provides a high-performance, backend-agnostic implementation of the Material Point Method (MPM) using the [Julia Language](https://julialang.org). It is lightweight and user-friendly, allowing efficient execution on various hardware accelerators with a single codebase.

Supported platform:

[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)

!!! warning
    
    This is the latest version of `MaterialPointSover.jl`, if you want to see the examples in the paper, please move to [https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper](https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper).

In the following, we will introduce the code structure and design logic of the solver, which we call as the `Interface`. After understanding the code logic, we provide several practical examples in the `Tutorial` to demonstrate how users should use the solver. Next, we discuss the advanced usage of the solver, mainly involving customizing the MPM process/algorithm. This is a very powerful feature, allowing users to implement their algorithms at any level and seamlessly integrate them into any accelerator backend. Finally, we list some very useful functions (`Useful Tools`) in this package. 

In Julia, if you have any questions about a function, you can type `?` in the REPL followed by the function name, for example:

```julia-repl
julia>?
help?> materialpointsolver!
search: materialpointsolver! MaterialPointSolver

  materialpointsolver!(args::DeviceArgs{T1, T2}, grid::DeviceGrid{T1, T2}, 
      mp::DeviceParticle{T1, T2}, attr::DeviceProperty{T1, T2}, 
      bc::DeviceVBoundary{T1, T2}; workflow::Function=procedure!)

  Description:
  ============

  This function is the main function of the MPM solver, user has to pre-define the data of args, grid, mp, attr and bc,
  they are the model configuration, background grid, material points, particle property and boundary conditions (2/3D).
```

If you have any questions about this solver, please do not hesitate to file an issue on [GitHub Issue](https://github.com/LandslideSIM/MaterialPointSolver.jl/issues).

## Installation

Just type `]` in Julia's `REPL` to enter the built-in Pkg manager:

```julia-repl
julia> ]
(@1.11) Pkg> add MaterialPointSolver
```

## Citation
If you use `MaterialPointSolver.jl` in your research, please consider to cite this paper:

```bib
@article{index,
  title={Here is the title},
  author={authors},
  journal={journal},
  year={year}
}
```

## Acknowledgement

This project is sponsored by [Risk Group | Université de Lausanne](https://wp.unil.ch/risk/) and [China Scholarship Council [中国国家留学基金管理委员会]](https://www.csc.edu.cn/).
