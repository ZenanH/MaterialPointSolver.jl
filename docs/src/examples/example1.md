# 2D Soil Collapse Test

!!! note

    1) We use a classic example: 2D & 3D granular collapse[^1], to demonstrate how to establish computational models and solve them. 
    2) To successfully run the code, we will install the dependencies at the beginning of the code provided below. If you are already familiar with the Julia Pkg ENV and have installed the necessary packages, you can ignore this part.
    3) We default to using Nvidia's GPU or x86/ARM CPUs. If you want to use other acceleration backends, please make modifications in the appropriate places in the code.
    4) We use Unicode to enhance readability (when comparing formulas), but sometimes it may be confused with regular letters, such as ``\nu`` and `v`. If you are using VSCode, you can enable the following in the settings:
    
       ```json
       "editor.unicodeHighlight.ambiguousCharacters": true,
       ```

    5) You can copy and save it in the `your_file.jl` and run the file directly using Julia in `REPL`:

       ```julia-repl
       julia> include("path/to/your_file.jl")
       ```

        Or run this file in the terminal:

       ```bash
       bash> julia path/to/your_file.jl
       ```

[^1]: Bui, H.H., Fukagawa, R., Sako, K., Ohno, S., 2008. Lagrangian meshfree particles method (SPH) for large deformation and failure flows of geomaterial using elastic–plastic soil constitutive model. Int. J. Numer. Anal. Methods Geomech. 32, 1537–1570. https://doi.org/10.1002/nag.688

::: details Here is the complete 2D code

```julia
using Pkg
Pkg.add(["MaterialPointSolver", "MaterialPointGenerator", "CairoMakie", "CUDA"])
using MaterialPointSolver
using MaterialPointGenerator
using CairoMakie
using CUDA
```

:::

## 2D model description

```@raw html
<br>
<img src="./figures/exp_model.png" width=50%>
<br>
```

In this example, aluminum bars are used to model the non-cohesive soil collapse. We use uGIMP to simulate the failure process of soil collapse. The geometry of the numerical model is depicted in the figure, with a length ``l`` of ``0.2\ m`` and a height ``h`` of ``0.1\ m``. The parameters of the numerical model are provided in the Table.

| Parameter   | Value            | Unit       | Description            |
|:-----------:|:----------------:|:----------:|:----------------------:|
| ``\rho``    | ``2650``         | ``kg/m^3`` | density                |
| ``\mu``     | ``0.3``          | -          | Poisson's ratio        |
| ``K``       | ``7\times 10^5`` | ``Pa``     | bulk modulus           |
| ``t_{sim}`` | ``1.0``          | second     | simulation time        |
| ``\phi``    | ``19.8``         | degree     | friction angle         |
| ``\alpha``  | ``1.0``          |  -         | FLIP-PIC mixing factor |

## 2D code explaination