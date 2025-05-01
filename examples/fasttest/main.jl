#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_druckerprager.jl                                                        |
|  Description: Case used to vaildate the functions                                        |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using WGLMakie
using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using KernelAbstractions
using HDF5, ProgressMeter, WriteVTK, NearestNeighbors
using CUDA
backend_name = :CUDA

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "funcs.jl"))

MaterialPointSolver.warmup(Val(backend_name))

raw_pts = readxyz(joinpath(@__DIR__, "output/slide.xyz"))
raw_low = readxyz(joinpath(@__DIR__, "output/low20.xyz"))
vmin = minimum(raw_low, dims=1)
raw_pts[:, 1] .-= vmin[1]
raw_pts[:, 2] .-= vmin[2]
raw_pts[:, 3] .-= vmin[3]
raw_low[:, 1] .-= vmin[1]
raw_low[:, 2] .-= vmin[2]
raw_low[:, 3] .-= vmin[3]

savexyz(joinpath(@__DIR__, "output/low20surface.xyz"), raw_low)

init_grid_space_x = 20
init_grid_space_y = 20
init_grid_space_z = 20
init_grid_range_x = [0-init_grid_space_x*6, 3330+init_grid_space_x*6]
init_grid_range_y = [0-init_grid_space_y*6, 5610+init_grid_space_y*6]
init_grid_range_z = [0-init_grid_space_z*6, 1790+init_grid_space_z*6]
init_mp_in_space  = 2
init_T            = 40
init_Te           = 10
init_ρs           = 2650
init_ν            = 0.3
init_Es           = 2e8
init_Ks           = init_Es / (3 * (1 - 2 * init_ν))
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_ΔT           = 0.0005#0.5 * init_grid_space_x / sqrt(init_Es / init_ρs)
init_step         = floor(init_T / init_ΔT / 100)
init_ϕ            = deg2rad(19)
init_FP           = "FP64"
init_basis        = :uGIMP
init_NIC          = 27

# args setup
args = UserArgs3D(
    Ttol         = init_T,
    Te           = 0,
    ΔT           = init_ΔT,
    time_step    = :fixed,
    FLIP         = 1,
    PIC          = 0,
    constitutive = :druckerprager,
    basis        = init_basis,
    hdf5         = true,
    hdf5_step    = init_step,
    MVL          = false,
    device       = backend_name,
    coupling     = :OS,
    scheme       = :MUSL,
    gravity      = -9.8,
    ζs           = 0,
    project_name = "3d_druckerprager",
    project_path = @__DIR__,
    ϵ            = init_FP
)

# grid setup
tmp_grid = UserGrid3D(
    ϵ     = init_FP,
    phase = 1,
    x1    = init_grid_range_x[1],
    x2    = init_grid_range_x[2],
    y1    = init_grid_range_y[1],
    y2    = init_grid_range_y[2],
    z1    = init_grid_range_z[1],
    z2    = init_grid_range_z[2],
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    dz    = init_grid_space_z,
    NIC   = init_NIC
)
offset    = getoffset(tmp_grid, raw_pts)
raw_pts .-= offset'
raw_low .-= offset'
tmp_ext   = get_ext(raw_low, tmp_grid, 6*init_grid_space_x, offset)
grid = UserGrid3D(
    ϵ     = init_FP,
    phase = 1,
    x1    = init_grid_range_x[1],
    x2    = init_grid_range_x[2],
    y1    = init_grid_range_y[1],
    y2    = init_grid_range_y[2],
    z1    = init_grid_range_z[1],
    z2    = init_grid_range_z[2],
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    dz    = init_grid_space_z,
    NIC   = init_NIC,
    ext   = tmp_ext
)

# material point setup
dx = grid.dx / init_mp_in_space
dy = grid.dy / init_mp_in_space
dz = grid.dz / init_mp_in_space
pts = raw_pts
mpρs = ones(size(pts, 1)) * init_ρs
mp = UserParticle3D(
    ϵ      = init_FP,
    phase  = 1,
    NIC    = init_NIC,
    dx     = dx,
    dy     = dy,
    dz     = dz,
    ξ      = pts,
    ρs     = mpρs
)

# property setup
nid = ones(mp.np)
attr = UserProperty(
    ϵ   = init_FP,
    nid = nid,
    ν   = [init_ν],
    Es  = [init_Es],
    Gs  = [init_Gs],
    Ks  = [init_Ks],
    σt  = [0],
    ϕ   = [init_ϕ],
    ϕr  = [0],
    ψ   = [0],
    c   = [0],
    cr  = [0],
    Hp  = [0]
)

# boundary setup
vx_idx  = zeros(grid.ni)
vy_idx  = zeros(grid.ni)
vz_idx  = zeros(grid.ni)
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 3330 || 
                       grid.ξ[i, 2] ≤ 0 || grid.ξ[i, 2] ≥ 5610 || 
                       grid.ξ[i, 3] ≤ 0 || grid.ξ[i, 3] ≥ 1790, 1:grid.ni)
tmp_idy = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 3330 || 
                       grid.ξ[i, 2] ≤ 0 || grid.ξ[i, 2] ≥ 5610 || 
                       grid.ξ[i, 3] ≤ 0 || grid.ξ[i, 3] ≥ 1790, 1:grid.ni)
tmp_idz = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 3330 || 
                       grid.ξ[i, 2] ≤ 0 || grid.ξ[i, 2] ≥ 5610 || 
                       grid.ξ[i, 3] ≤ 0 || grid.ξ[i, 3] ≥ 1790, 1:grid.ni)
vx_idx[tmp_idx] .= 1
vy_idx[tmp_idy] .= 1
vz_idx[tmp_idz] .= 1
bc = UserVBoundary3D(
    ϵ        = init_FP,
    vx_s_idx = vx_idx,
    vx_s_val = zeros(grid.ni),
    vy_s_idx = vy_idx,
    vy_s_val = zeros(grid.ni),
    vz_s_idx = vz_idx,
    vz_s_val = zeros(grid.ni)
)

# solver setup
materialpointsolver!(args, grid, mp, attr, bc, workflow=Fprocedure!)
Tanimation(args, offset)
# let 
#     fig=Figure()
#     ax = LScene(fig[1, 1])
#     scatter!(ax, mp.ξ, color=:red, markersize=1)
#     vid = findall(grid.ext.id .≠ 0)
#     scatter!(ax, grid.ξ[vid, :], color=:blue, markersize=1)
#     display(fig)
# end


# let
#     figfont = MaterialPointSolver.tnr
#     fig = Figure(size=(1200, 700), fonts=(; regular=figfont, bold=figfont), fontsize=30)
#     ax = Axis3(fig[1, 1], xlabel=L"x\ (m)", ylabel=L"y\ (m)", zlabel=L"z\ (m)", 
#         aspect=:data, azimuth=0.2*π, elevation=0.1*π, xlabeloffset=60, zlabeloffset=80,
#         protrusions=100, xticks=(0:0.04:0.04), height=450, width=950)
#     pl1 = scatter!(ax, mp.ξ, color=log10.(mp.ϵq.+1), colormap=:jet, markersize=3,
#         colorrange=(0, 1))
#     Colorbar(fig[1, 1], limits=(0, 1), colormap=:jet, size=16, ticks=0:0.5:1, spinewidth=0,
#         label=L"log_{10}(\epsilon_{II}+1)", vertical=false, tellwidth=false, width=200,
#         halign=:right, valign=:top, flipaxis=false)
#     display(fig)
# end
# rm(joinpath(abspath(args.project_path), args.project_name), recursive=true, force=true)