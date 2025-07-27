#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 2d_liakopoulos.jl                                                          |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : Liakopoulos unsaturated soil test                                          |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using KernelAbstractions
using HDF5
using CairoMakie
using Printf

include(joinpath(@__DIR__, "2d_funcs.jl"))

MaterialPointSolver.warmup(Val(:CPU))

# model configuration
init_grid_space_x = 0.04
init_grid_space_y = 0.04
init_grid_range_x = [-init_grid_space_x*2, 0.04+init_grid_space_x*2]
init_grid_range_y = [-init_grid_space_y*2, 1.00+init_grid_space_y*2]
init_mp_in_space  = 2
init_ρs           = 2e3
init_ρw           = 1e3
init_n            = 0.2975
init_k            = 1
init_ν            = 0.4
init_Es           = 1.3e6
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_Ks           = init_Es / (3 * (1 - 2 * init_ν))
init_Kw           = 2e9
init_ΔT           = 1e-5
init_T            = 2.0
init_Te           = 0.0
init_step         = floor(init_T / init_ΔT / 200)
init_basis        = :linear
init_NIC          = 4
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
    FLIP         = 0.98,
    PIC          = 0.02,
    hdf5_step    = init_step,
    MVL          = false,
    device       = :CPU,
    coupling     = :TS,
    scheme       = :MUSL,
    gravity      = -9.8,
    ζs           = 0.,
    ζw           = 0.,
    project_name = "2d_terzaghi",
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
ξ0 = meshbuilder(0 + dx / 2 : dx : 0.1 - dx / 2, 0 + dy / 2 : dy : 1.0 - dy / 2)
mp = UserParticle2D(
    ϵ     = init_ϵ,
    phase = init_phase,
    NIC   = init_NIC,
    dx    = dx,
    dy    = dy,
    ξ     = ξ0,
    S     = ones(size(ξ0, 1)),
    ρs    = ones(size(ξ0, 1)) .* init_ρs,
    ρw    = ones(size(ξ0, 1)) .* init_ρw,
    n     = ones(size(ξ0, 1)) .* init_n,
    k     = ones(size(ξ0, 1)) .* init_k
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
    Kw  = [init_Kw]
)

# boundary condition nodes index
vx_idx  = zeros(grid.ni)
vy_idx  = zeros(grid.ni)
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 0.04, 1:grid.ni)
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
materialpointsolver!(args, grid, mp, attr, bc, workflow=Tprocedure!)

# post-processing
animation(args)