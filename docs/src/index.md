## Introduction

This package provides a high-performance, backend-agnostic implementation of the Material Point Method (MPM) using the [Julia Language](https://julialang.org). It is lightweight and user-friendly, allowing efficient execution on various hardware accelerators with a single codebase.

Supported platform:

[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)

!!! warning
    
    This is the latest version of `MaterialPointSover.jl`, if you want to see the examples in the paper, please move to [https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper](https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper).

If you have any questions about this solver, please do not hesitate to file an issue on [GitHub Issue](https://github.com/LandslideSIM/MaterialPointSolver.jl/issues).

## Installation

Just type `]` in Julia's `REPL`:

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
