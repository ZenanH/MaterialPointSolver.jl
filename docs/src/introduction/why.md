# Why?

This project began as a research initiative where we aimed to study the run-out process of landslides after instability through numerical simulations. However, this nonlinear large deformation problem is difficult to address using mesh-based methods because the mesh can become distorted. At this point, the Material Point Method (MPM), which combines a Lagrangian perspective with an Eulerian perspective, will be quite advantageous for this type of problem.

![MPM process](./figures/mpmprocess.png)

However, MPM does not have as many available software/codes as FEM, and there is also less work focused on computational performance. Notable examples include Anura3D, which is written in Fortran, and CB-Geo MPM, which is written in C++. Installing, configuring, and using these libraries is not easy. We hope to leverage the advantages of the Julia language to allow users to install with one click and achieve better performance.

These software and code, due to being written a long time ago, are no longer compatible with the rapidly developing computing resources of today. Currently, we mainly utilize the powerful performance of GPUs, and it is important to note that GPU vendors are not limited to NVIDIA. In this case, whether as researchers or related professionals, it would be a fantastic experience if we could write just one set of code that can run on different backends.

For these reasons, we have implemented ***MaterialPointSolver.jl***. 

> Since this is a personal project and maintainability is a priority, we only provide the core infrastructure and utilities for MPM. The governing equations and algorithms are implemented through a plugin system. If you want to develop your own methods, the plugin system allows full customization of computational details without handling mesh topology, data transfer, post-processing, etc. If you prefer not to modify anything, you can use the existing features in the lib folder. Your work can be open-sourced under the MIT license (PRs welcome) or kept as private code.