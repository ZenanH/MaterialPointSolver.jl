# Ecosystem

Completing a high-quality numerical simulation task is not easy. 

Generally, we divide it into three parts: pre-processing, solving, and post-processing. Taking MPM as an example, pre-processing involves discretizing the problem domain into particles (structured particles) and addressing issues such as classifying particle attributes and applying boundary conditions [***MaterialPointGenerator.jl***](https://github.com/LandslideSIM/MaterialPointGenerator.jl). Solving refers to the part handled by [***MaterialPointSolver.jl***](https://github.com/LandslideSIM/MaterialPointSolver.jl). Post-processing is responsible for exporting the solution results and visualizing them [***MaterialPointVisualizer.jl***](https://github.com/LandslideSIM/MaterialPointVisualizer.jl).

The MPM algorithm is not as mature as FEM, so the related software for these three parts may not be very user-friendly. Our **LandslideSIM** project provides many relevant tools; for specific information, you can check at [LandslideSIM organization](https://github.com/LandslideSIM).