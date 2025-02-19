# MPM processes

!!! note

    This package provides an explicit MPM solver.

We classify MPM processes based on stress update schemes. Generally, the stress update schemes include: 

- `USL`
- `USF`
- `MUSL`
- `AFFINE` 

Here, `AFFINE`[^1] is an exception; it does not represent a stress update scheme, but since we have only implemented one stress update scheme for `AFFINE`, similar to `USF`, it is also an option.

[^1]: He, K. Y., Liang, W., Yin, Z. Y., & Jin, Y. F. (2023). An efficient material point method framework based on the affine matrix. Computers and Geotechnics, 163, 105712. https://doi.org/10.1016/j.compgeo.2023.105712

This implementation is located under the `/src/solvers/` path, you can see some files like `**_MUSL.jl`. We just need to rearrange/modifiy the kernel functions, super easy! 🎸