#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 2d_Liakopoulos.jl                                                          |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : 2D infiltration test                                                       |
|  Reference  : Ceccato, F., Yerro, A., Girardi, V., Simonini, P., 2021. Two-phase dynamic |
|               MPM formulation for unsaturated soil. Computers and Geotechnics 129,       |
|               103876, https://doi.org/10.1016/j.compgeo.2020.103876                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using CairoMakie
using Printf

include(joinpath(@__DIR__, "2d_utils.jl"))
include(joinpath(@__DIR__, "2d_funcs.jl"))

MaterialPointSolver.warmup(Val(:CPU))

# model configuration
init_grid_space_x = 0.04
init_grid_space_y = 0.04
init_grid_range_x = [0-init_grid_space_x*8, 0.04+init_grid_space_x*8]
init_grid_range_y = [0-init_grid_space_y*8, 1.00+init_grid_space_y*8]
init_mp_in_space  = 2
init_ρs           = 2e3
init_ρw           = 1e3
init_n            = 0.2975
init_k            = 4.41e-6
init_ν            = 0.4
init_S            = 1.0
init_S_max        = 1.0
init_S_min        = 0.125
init_P_ref        = 3e3
init_λ            = 0.7
init_α            = 2e-4
init_β            = 2.214
init_Es           = 1.3e6
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_Ks           = init_Es / (3 * (1 - 2 * init_ν))
init_Kw           = 2e9
init_T            = 600#7200
init_Te           = 0.0
init_ΔT           = 1e-5
init_step         = floor(init_T / init_ΔT / 200)
init_basis        = :linear
init_NIC          = 4
init_phase        = 2
init_Ŵ            = -1e-4
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
    project_name = "2d_liakopoulos",
    project_path = @__DIR__,
    ϵ            = init_ϵ
)

# background grid setup
tmp = UserGrid2D(
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
ext = NEWGrid(zeros(tmp.ni))
user_adapt(Array, ext)
grid = UserGrid2D(
    ϵ     = init_ϵ,
    phase = init_phase,
    x1    = init_grid_range_x[1],
    x2    = init_grid_range_x[2],
    y1    = init_grid_range_y[1],
    y2    = init_grid_range_y[2],
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    NIC   = init_NIC,
    ext   = ext
)

# material points setup
dx = grid.dx / init_mp_in_space
dy = grid.dy / init_mp_in_space
ξ0 = meshbuilder(0 + dx / 2 : dx : 0.04 - dx / 2, 0 + dy / 2 : dy : 1.0 - dy / 2)
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
    S     = ones(size(ξ0, 1)) .* init_S,
    k     = ones(size(ξ0, 1)) .* init_k
)

# particle property setup
nid  = ones(mp.np)
tmp = UserProperty(
    ϵ   = init_ϵ,
    nid = nid,
    ν   = [init_ν],
    Es  = [init_Es],
    Gs  = [init_Gs],
    Ks  = [init_Ks],
    Kw  = [init_Kw]
)
ext = NEWProperty(init_S_min, init_S_max, init_P_ref, init_λ, init_α, init_β, init_Ŵ)
attr = UserProperty(
    ϵ   = init_ϵ,
    nid = nid,
    ν   = [init_ν],
    Es  = [init_Es],
    Gs  = [init_Gs],
    Ks  = [init_Ks],
    Kw  = [init_Kw],
    ext = ext
)

# boundary condition nodes index
vx_idx  = zeros(grid.ni)
vy_idx  = zeros(grid.ni)
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 0.04, 1:grid.ni)
tmp_idy = findall(i -> grid.ξ[i, 2] ≤ 0, 1:grid.ni)
vx_idx[tmp_idx] .= 1
vy_idx[tmp_idy] .= 1
idx = findall(i -> 0 ≤ grid.ξ[i, 1] ≤ 0.11 && grid.ξ[i, 2] ≈ 1.0, 1:grid.ni)
ext = TractionBoundary(idx)
user_adapt(Array, ext)
bc = UserVBoundary2D(
    vx_s_idx = copy(vx_idx),
    vx_s_val = zeros(grid.ni),
    vy_s_idx = copy(vy_idx),
    vy_s_val = zeros(grid.ni),
    vx_w_idx = copy(vx_idx),
    vx_w_val = zeros(grid.ni),
    vy_w_idx = zeros(grid.ni),
    vy_w_val = zeros(grid.ni),
    ext      = ext
)

# MPM solver
materialpointsolver!(args, grid, mp, attr, bc, workflow=Tprocedure!)

# post-processing
Tanimation(args)