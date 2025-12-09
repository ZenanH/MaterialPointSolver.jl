# Structs design

The solver provides three main struct types: 1) **Configurations**, 2) **Background grid**, 3) **Material particles**.

**Configurations** use the `Config` type and store settings such as simulation time, adaptive time stepping, whether to save results as HDF5, etc.

The **Background grid** uses the `DeviceGrid` type and stores information such as grid size, spacing, and number of nodes. Only structured grids are supported.

**Material particles** use the `DeviceParticle` type and represent the discrete particles of the model, storing quantities such as particle count, velocity, stress, and more.

Both `DeviceGrid` and `DeviceParticle` support the plugin system. If the default fields are not sufficient (which is true in most cases), users may add custom fields.

!!! note

    Custom fields can only be scalar or array types of Real, e.g., `(custom1 = 1.2, custom2 = rand(10, 2))`.

!!! note

    You are responsible for specifying the numerical precision of custom fields.