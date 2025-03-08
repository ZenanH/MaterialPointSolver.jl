# *MaterialPointSolver* <img src="docs/src/assets/logo.png" align="right" height="126" />

[![CI](https://github.com/LandslideSIM/MaterialPointSolver.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/LandslideSIM/MaterialPointSolver.jl/actions/workflows/ci.yml) 
[![](https://img.shields.io/badge/docs-stable-blue.svg?logo=quicklook)](https://landslidesim.github.io/MaterialPointSolver.jl/stable/)
[![](https://img.shields.io/badge/version-v0.4.1-pink)]()

[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)

<p>
This package provides a high-performance, backend-agnostic implementation of the Material Point Method (MPM) using the <a href="https://julialang.org" target="_blank"><img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em"> Julia Language</a>. It is lightweight and user-friendly, allowing efficient execution on various hardware accelerators with a single codebase. Please check here for the <a href="https://landslidesim.github.io/MaterialPointSolver.jl/stable/" target="_blank">documentation</a>.
</p>

<img src="docs/assets/readme.gif" width="100%" align="center">

## Installation ⚙️

Just type <kbd>]</kbd> in Julia's `REPL`:

```julia
julia> ]
(@1.11) Pkg> add MaterialPointSolver
```

## Features 💪

*These features can be combined in any way, but MLS-MPM can only use quadratic b-spline for the speed*

- Basis function:

  - ✅ standard MPM
  - ✅ uGIMP (uniformed Generalized interpolation MPM)
  - ✅ quadratic B-spline
  - ✅ cubic B-spline (boundary modified)

- Stress update scheme:

  - ✅ USL (update stress last)
  - ✅ USF (update stress first)
  - ✅ MUSL (modified USL)

- MPM formulation:

  - ✅ one-phase single-point
  - 🚧 two-phase single-point (saturated/unsaturated)

- Constitutive model:

  - ✅ linear elastic
  - ✅ hyper elastic (Neo-Hookean)
  - ✅ Drucker-Prager (with softening/harding)
  - 🚧 Mohr-Coulomb
  - ✅ Bingham
    
    …

- Others:

  - ✅ Affine/MLS-MPM
  - ✅ $\bar{F}$-based volumetric locking elimination
  - ✅ Gaussian random field
  - ✅ one-click switch between `FP64` and `FP32`
  - ✅ user-defined algorithms/extensions at any level

## Citation 🔥

If you find `MaterialPointSolver.jl` useful or have used it in your research, please cite it as follows:

```bib
@article{index,
  title={Here is the title},
  author={authors},
  journal={journal},
  year={year}
}
```
> [!CAUTION]
> This is the latest version of `MaterialPointSover.jl`, if you want to see the examples in the paper, please move to [https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper](https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper).

## Acknowledgement 👍

This project is sponsored by [Risk Group | Université de Lausanne](https://wp.unil.ch/risk/) and [China Scholarship Council [中国国家留学基金管理委员会]](https://www.csc.edu.cn/).

## MPM ➕ Julia

* [[package]: elastoPlasm.jl](https://github.com/ewyser/elastoPlasm.jl) is fully witten in Julia, it solves explicit elasto-plastic problems within a finite deformation framework.

* [[package]: Tesserae.jl](https://github.com/KeitaNakamura/Tesserae.jl) is an MPM-related Julia package, it provides some useful functions that can be used for MPM, such as convenient macros for transferring data between grids and particles.

* [[code]: MPM-Julia](https://github.com/vinhphunguyen/MPM-Julia) is the code for the paper: Sinai, V.P. Nguyen, C.T. Nguyen and S. Bordas. Programming the Material Point Method in Julia. Advances in Engineering Software,105: 17--29, 2017.

* [[code]: jump](https://github.com/vinhphunguyen/jump) is for the theory of the MPM described in the book 'The Material Point Method: Theory, Implementations and Applications (Scientific Computation) 1st ed. 2023 Edition'. [https://link.springer.com/book/10.1007/978-3-031-24070-6](https://link.springer.com/book/10.1007/978-3-031-24070-6)
