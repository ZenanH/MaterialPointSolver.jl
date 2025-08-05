# Solver's workflow

The solver's workflow is defined in the `/src/solver.jl` file. It will call the corresponding MPM processes based on the configuration and automatically update the time steps. If the user specifies that the computation results be exported as the `HDF5` file, then at the corresponding time steps, the data will be downloaded from the hardware accelerator to the CPU and then written to disk.

In addition to using a fixed time step, we also support automatically updating the time step (adaptive time step) through simulating state variables (CFL conditions).

In order to facilitate users in viewing the calculation progress, we have also introduced a progress bar through `ProgressMeter.jl` and developed a suitable algorithm to check the ETA, even with adaptive time steps.

!!! tip

    The configurations we mentioned can be set when instantiating `Args2D`/`Args3D`.

The current `submit_work!` function is very powerful, as its input parameters support a custom MPM processes function, greatly enhancing the scalability of this solver. We will demonstrate how to use it in the **Advanced Topics** section.