#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 2d_slumpblock.jl                                                           |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : Self-weight consolidation test                                             |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using HDF5
using CairoMakie
using Printf
using CUDA
device_type = :CPU

include(joinpath(@__DIR__, "2d_funcs.jl"))
MaterialPointSolver.warmup(Val(device_type))

# model configuration
init_grid_space_x = 0.125
init_grid_space_y = 0.125
init_grid_range_x = [0-init_grid_space_x, 3+init_grid_space_x]
init_grid_range_y = [0-init_grid_space_y, 2+init_grid_space_y]
init_mp_in_space  = 2
init_ρs           = 2650
init_ρw           = 1000
init_n            = 0.4
init_k            = 1e-4
init_ν            = 0.3
init_Es           = 1e5
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_Ks           = init_Es / (3 * (1 - 2 * init_ν))
init_Kw           = 2.2e9
init_T            = 0.5
init_Te           = 0.1
init_ΔT           = 1e-6
init_step         = floor(init_T / init_ΔT / 150)
init_basis        = :bspline2
init_NIC          = 9
init_phase        = 2
init_ϵ            = "FP64"

# parameters setup
args = UserArgs2D(
    Ttol         = init_T,
    Te           = init_Te,
    ΔT           = init_ΔT,
    time_step    = :fixed,
    constitutive = :linearelastic,
    basis        = init_basis,
    hdf5         = true,
    hdf5_step    = init_step,
    MVL          = true,
    device       = device_type,
    coupling     = :TS,
    scheme       = :MUSL,
    gravity      = -9.8,
    ζs           = 0.,
    ζw           = 0.,
    project_name = "2d_slumpblock",
    project_path = @__DIR__,
    ϵ            = init_ϵ
)

# background grid setup
grid = UserGrid2D(
    ϵ     = init_ϵ,
    phase = init_phase,
    x1    = init_grid_range_x[1],
    x2    = init_grid_range_x[2],
    y1    = init_grid_range_y[1],
    y2    = init_grid_range_y[2],
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    NIC   = init_NIC
)

# material points setup
dx = grid.dx / init_mp_in_space
dy = grid.dy / init_mp_in_space
ξ0 = meshbuilder(0 + dx / 2 : dx : 2 - dx / 2, 0 + dy / 2 : dy : 2.0 - dy / 2)
mp = UserParticle2D(
    ϵ     = init_ϵ,
    phase = init_phase,
    NIC   = init_NIC,
    dx    = dx,
    dy    = dy,
    ξ     = ξ0,
    ρs    = ones(size(ξ0, 1)) .* init_ρs,
    ρw    = ones(size(ξ0, 1)) .* init_ρw,
    n     = ones(size(ξ0, 1)) .* init_n,
)

# particle property setup
nid  = ones(mp.np)
attr = UserProperty(
    ϵ   = init_ϵ,
    nid = nid,
    ν   = [init_ν],
    Es  = [init_Es],
    Gs  = [init_Gs],
    Ks  = [init_Ks],
    Kw  = [init_Kw],
    k   = [init_k]
)

# boundary condition nodes index
vx_idx  = zeros(grid.ni)
vy_idx  = zeros(grid.ni)
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0, 1:grid.ni)
tmp_idy = findall(i -> grid.ξ[i, 2] ≤ 0, 1:grid.ni)
vx_idx[tmp_idx] .= 1
vy_idx[tmp_idy] .= 1
bc = UserVBoundary2D(
    vx_s_idx = copy(vx_idx),
    vx_s_val = zeros(grid.ni),
    vy_s_idx = copy(vy_idx),
    vy_s_val = zeros(grid.ni),
    vx_w_idx = copy(vx_idx),
    vx_w_val = zeros(grid.ni),
    vy_w_idx = copy(vy_idx),
    vy_w_val = zeros(grid.ni)
)

# MPM solver
materialpointsolver!(args, grid, mp, attr, bc, workflow=tprocedure!)

# post-processing
animation(args)

let 
    # theme setup
    set_theme!(theme_latexfonts())
    # figure setup
    fig = Figure(size=(500, 300))
    # axis setup
    ax = Axis(fig[1, 1], aspect=DataAspect(), xticks=(0:0.5:3), yticks=(0:1:2))
    ax.xtrimspine = true
    ax.ytrimspine = true
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.rightspinevisible = false
    ax.topspinevisible = false
    limits!(ax, -0.1, 3.1, -0.1, 2.1)
    # plot setup
    plt1 = scatter!(ax, mp.ξ, color=-mp.σw/1e3, marker=:rect, markersize=10)
    plt1.colorrange = (0, 4e4)
    # colorbar setup
    Colorbar(fig[1, 2], plt1, spinewidth=0, width=8)
    # visualize
    display(fig)
end