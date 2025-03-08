# Tips

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