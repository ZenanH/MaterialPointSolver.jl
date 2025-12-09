# Get Started

## Installation

Just type `]` in Julia's `REPL` to enter the built-in Pkg manager:

```julia-repl
julia> ]
(@1.12) Pkg> add MaterialPointSolver
```

## Supported platforms

::: code-group

```julia [CPU x86/arm]
using MaterialPointSolver
```

```julia [CUDA | NVIDIA]
using MaterialPointSolver
using CUDA
```

```julia [ROCm | AMD]
using MaterialPointSolver
using AMDGPU
```

```julia [Metal | Apple]
using MaterialPointSolver
using Metal
```

```julia [oneAPI | Intel]
using MaterialPointSolver
using oneAPI
```

:::

## Citation

If you find ***MaterialPointSolver.jl*** useful or have used it in your research, please cite it as follows:

```latex
@article{HUO2025107189,
title = {A high-performance backend-agnostic Material Point Method solver in Julia},
journal = {Computers and Geotechnics},
volume = {183},
pages = {107189},
year = {2025},
issn = {0266-352X},
doi = {https://doi.org/10.1016/j.compgeo.2025.107189},
url = {https://www.sciencedirect.com/science/article/pii/S0266352X25001387},
author = {Zenan Huo and Yury Alkhimenkov and Michel Jaboyedoff and Yury Podladchikov and Ludovic Räss and Emmanuel Wyser and Gang Mei},
keywords = {MPM, Julia language, Heterogeneous computing, Effective memory throughput}
```

::: warning
    
This is the latest version of ***MaterialPointSover.jl***, if you want to see the examples in the paper, please move to [https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper](https://github.com/LandslideSIM/Archive_MaterialPointSolver.jl_paper).

:::

## Acknowledgement

This project is sponsored by [Risk Group | Université de Lausanne](https://wp.unil.ch/risk/) and [China Scholarship Council [中国国家留学基金管理委员会]](https://www.csc.edu.cn/).