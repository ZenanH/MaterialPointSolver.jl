begin
using MaterialPointGenerator  # for pre-processing
using MaterialPointSolver
using MaterialPointVisualizer # for post-processing
using CairoMakie              # for visualization
using SysInfo                 # for querying the system information
# make sure you already installed the required packages
active_dev = :CPU             # :CPU, :CUDA, :ROCm, :oneAPI, :Metal
assetsdir  = MaterialPointSolver.assets_dir
sysinfo()

init_grid_space_x = 0.04
init_grid_space_y = 0.04
init_grid_range_x = [0-init_grid_space_x*2, 0.04+init_grid_space_x*2]
init_grid_range_y = [0-init_grid_space_y*2, 1.00+init_grid_space_y*2]
init_mp_in_space  = 2
init_T            = 1
init_ρs           = 2e3
init_n            = 0#0.2975
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
init_T            = 1e-5*20#20#7200
init_Te           = 0.0
init_ΔT           = 1e-5
init_step         = floor(init_T / init_ΔT / 200)
init_basis        = :linear
init_NIC          = 4
init_Ŵ            = -1e-4
init_ϵ            = "FP64"

args = UserArgs2D(
    Ttol         = init_T,
    Te           = 0,
    ΔT           = init_ΔT,
    time_step    = :fixed,
    FLIP         = 0.98,
    PIC          = 0.02,
    constitutive = :linearelastic,
    basis        = init_basis,
    hdf5         = false,
    hdf5_step    = init_step,
    MVL          = false,
    device       = active_dev,
    coupling     = :OS,
    scheme       = :MUSL,
    gravity      = -9.8,
    ζs           = 0,
    project_name = "2d_druckerprager",
    project_path = @__DIR__,
    ϵ            = init_ϵ
)

grid = UserGrid2D(
    ϵ     = init_ϵ,
    phase = 1,
    x1    = init_grid_range_x[1],
    x2    = init_grid_range_x[2],
    y1    = init_grid_range_y[1],
    y2    = init_grid_range_y[2],
    dx    = init_grid_space_x,
    dy    = init_grid_space_y,
    NIC   = init_NIC
)

dx = grid.dx / init_mp_in_space
dy = grid.dy / init_mp_in_space
ξ0 = meshbuilder(0 + dx / 2 : dx : 0.04 - dx / 2,
                 0 + dy / 2 : dy : 1.00 - dy / 2)
mp = UserParticle2D(
    ϵ     = init_ϵ,
    phase = 1,
    NIC   = init_NIC,
    dx    = dx,
    dy    = dy,
    ξ     = ξ0,
    ρs    = ones(size(ξ0, 1)) .* init_ρs
); @. mp.σij[:, 2] = 9800 * (mp.ξ0[:, 2] - 1)

nid = ones(mp.np)
attr = UserProperty(
    ϵ   = init_ϵ,
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

vx_idx = zeros(grid.ni)
vy_idx = zeros(grid.ni)
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0.0 || grid.ξ[i, 1] ≥ 0.04, 1:grid.ni)
tmp_idy = findall(i -> grid.ξ[i, 2] ≤ 0.0, 1:grid.ni)
vx_idx[tmp_idx] .= 1
vy_idx[tmp_idy] .= 1 
bc = UserVBoundary2D(
    ϵ        = init_ϵ,
    vx_s_idx = vx_idx,
    vx_s_val = zeros(grid.ni),
    vy_s_idx = vy_idx,
    vy_s_val = zeros(grid.ni)
)
end

ΔT, Ti, T1, T2 = 1e-5, 0, Int64, Float64
Ti < args.Te ? G = args.gravity / args.Te * Ti : G = args.gravity
dev = getBackend(Val(args.device))
resetgridstatus_OS!(dev)(ndrange=grid.ni, grid)
resetmpstatus_OS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
P2G_OS!(dev)(ndrange=mp.np, grid, mp, G)
solvegrid_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs)
doublemapping1_OS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
doublemapping2_OS!(dev)(ndrange=mp.np, grid, mp)
doublemapping3_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT)