#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_case.jl                                                                 |
|  Description: Case used to vaildate the functions                                        |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointSolver
using MaterialPointGenerator
using KernelAbstractions
using CUDA
using BenchmarkTools
using UnicodePlots
using DelimitedFiles
using Dates

MaterialPointSolver.warmup(Val(:CUDA))

# 0.01000:   8000 pts
# 0.00660:  27000 pts
# 0.00500:  64000 pts
# 0.00400: 125000 pts | grid: x[-0.008, 0.808] y[-0.008, 0.056] z[-0.08, 0.108]
# 0.00333: 216000 pts
# 0.00250: 512000 pts
init_grid_space_x = 0.0025
init_grid_space_y = 0.0025
init_grid_space_z = 0.0025
init_grid_range_x = [-0.02, 0.07]
init_grid_range_y = [-0.02, 0.75]
init_grid_range_z = [-0.02, 0.12]
init_mp_in_space  = 2
init_T            = 1
init_ρs           = 2650
init_ν            = 0.3
init_Ks           = 7e5
init_Es           = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_ΔT           = 0.5 * init_grid_space_x / sqrt(init_Es / init_ρs)
init_step         = floor(init_T / init_ΔT / 50)
init_ϕ            = deg2rad(19.8)
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
    hdf5         = false,
    hdf5_step    = init_step,
    MVL          = false,
    device       = :CUDA,
    coupling     = :OS,
    scheme       = :MUSL,
    progressbar  = true,
    gravity      = -9.8,
    ζs           = 0,
    project_name = "3d_profile",
    project_path = @__DIR__,
    ϵ            = init_FP
)

# grid setup
grid = UserGrid3D(
    phase =  1,
    x1    = -0.02,
    x2    =  0.07,
    y1    = -0.02,
    y2    =  0.75,
    z1    = -0.02,
    z2    =  0.12,
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    dz    = init_grid_space_z,
    NIC   = init_NIC,
    ϵ     = init_FP
)

# material point setup
dx = grid.dx / init_mp_in_space
dy = grid.dy / init_mp_in_space
dz = grid.dz / init_mp_in_space
ξ0 = meshbuilder(0 + dx / 2 : dx : 0.05 - dx / 2,
                 0 + dy / 2 : dy : 0.20 - dy / 2,
                 0 + dz / 2 : dz : 0.10 - dz / 2)
mpρs = ones(size(ξ0, 1)) * init_ρs
mp = UserParticle3D(
    ϵ     = init_FP,
    phase = 1,
    NIC   = init_NIC,
    dx    = dx,
    dy    = dy,
    dz    = dz,
    ξ     = ξ0,
    ρs    = mpρs
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
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 0.05 ||
                       grid.ξ[i, 3] ≤ 0 || grid.ξ[i, 2] ≤ 0, 1:grid.ni)
tmp_idy = findall(i -> grid.ξ[i, 2] ≤ 0 || grid.ξ[i, 3] ≤ 0, 1:grid.ni)
tmp_idz = findall(i -> grid.ξ[i, 3] ≤ 0, 1:grid.ni)
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

initmpstatus!(CPU())(ndrange=mp.np, grid, mp, Val(args.basis))
# variables setup for the simulation 
T1 = Int64
T2 = Float64
Ti = T2(0.0)
ΔT = args.ΔT
dev_grid, dev_mp, dev_attr, dev_bc = host2device(grid, mp, attr, bc, Val(args.device))
G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
dev = getBackend(Val(args.device))
results = ["resetgridstatus_OS!"   0.0;
           "resetmpstatus_OS_CPU!" 0.0;
           "P2G_OS!"               0.0;
           "solvegrid_a_OS!"       0.0;
           "doublemapping1_a_OS!"  0.0;
           "doublemapping2_OS!"    0.0;
           "doublemapping3_OS!"    0.0;
           "G2P_OS!"               0.0;
           "liE!"                  0.0;
           "dpP!"                  0.0;]

results[1, 2] = @belapsed begin 
    resetgridstatus_OS!($dev)(ndrange=$dev_grid.ni, $dev_grid)
    KAsync($dev)
end
if args.basis == :uGIMP && args.device == :CUDA
    results[2, 2] = @belapsed begin 
        resetmpstatus_OS_CUDA!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp, $Val(args.basis)) 
        KAsync($dev)
    end
else
    results[2, 2] = @belapsed begin 
        resetmpstatus_OS!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp, $Val(args.basis))
        KAsync($dev)
    end
end
results[3, 2] = @belapsed begin 
    P2G_OS!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp, $G)
    KAsync($dev)
end
results[4, 2] = @belapsed begin 
    solvegrid_OS!($dev)(ndrange=$dev_grid.ni, $dev_grid, $dev_bc, $ΔT, $args.ζs)
    KAsync($dev)
end
results[5, 2] = @belapsed begin 
    doublemapping1_OS!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp, $dev_attr, $ΔT, $args.FLIP, $args.PIC)
    KAsync($dev)
end
results[6, 2] = @belapsed begin 
    doublemapping2_OS!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp)
    KAsync($dev)
end
results[7, 2] = @belapsed begin 
    doublemapping3_OS!($dev)(ndrange=$dev_grid.ni, $dev_grid, $dev_bc, $ΔT)
    KAsync($dev)
end
results[8, 2] = @belapsed begin 
    G2P_OS!($dev)(ndrange=$dev_mp.np, $dev_grid, $dev_mp, $ΔT)
    KAsync($dev)
end
results[9, 2] = @belapsed begin 
    liE!($dev)(ndrange=$dev_mp.np, $dev_mp, $dev_attr)
    KAsync($dev)
end
results[10, 2] = @belapsed begin 
    dpP!($dev)(ndrange=$dev_mp.np, $dev_mp, $dev_attr)
    KAsync($dev)
end

tol_time = sum(results[:, 2])
func_value = round.(results[:, 2] ./ tol_time * 100, digits=1)
func_names = results[:, 1]
tol_time2 = round(tol_time * 1000, digits=1)
fig = barplot(func_names, func_value, title="Execution Time Proportion (%)")
display(fig)

file_name = string(Dates.format(now(), "yyyy-mm-dd-HH-MM-SS"), ".csv")
open(joinpath(args.project_path, args.project_name, file_name), "w") do io
    results[:, 2] .*= 1000
    header = ["function_name" "execution_time(ms)"]
    rst = vcat(header, results)
    writedlm(io, rst, ',')
end