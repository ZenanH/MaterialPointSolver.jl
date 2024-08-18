#==========================================================================================+
|                MPMSolver.jl: High-performance MPM Solver for Geomechanics                |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_collapse.jl                                                             |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : 3D collapse in Drucker-Prager constitutive equation                        |
|  Notes      : We need to ignore the CFL condition, grid {σm, vol} and volume locking     |
|               when we do the MTeff test.                                                 |
+==========================================================================================#

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
init_grid_range_y = [-0.02, 0.82]
init_grid_range_z = [-0.02, 0.12]
init_mp_in_space  = 2
init_project_name = "3d_collapse"
init_project_path = joinpath(rtsdir, init_project_name)
init_constitutive = :druckerprager
init_gravity      = -9.8
init_ζ            = 0
init_ρs           = 2700
init_ν            = 0
init_E            = 1e6
init_Ks           = init_E/(3*(1-2*init_ν))
init_G            = init_E/(2*(1+  init_ν))
init_T            = 1
init_Te           = 0
init_ΔT           = 0.5*init_grid_space_x/sqrt((init_Ks+4/3*init_G)/init_ρs)
init_step         = floor(init_T/init_ΔT/20) |> Int64
init_step<10 ? init_step=1 : nothing
init_σt           = 0
init_ϕ            = 19.8*π/180
init_c            = 0
init_ψ            = 0
init_basis        = :uGIMP
init_phase        = 1
init_NIC          = 64
iInt              = Int64
iFloat            = Float64

# parameters setup
args = Args3D{iInt, iFloat}(
    Ttol         = init_T,
    ΔT           = init_ΔT,
    Te           = init_Te,
    time_step    = :fixed,
    FLIP         = 1,
    PIC          = 0,
    ζ            = init_ζ,
    gravity      = init_gravity,
    project_name = init_project_name,
    project_path = init_project_path,
    constitutive = init_constitutive,
    animation    = false,
    hdf5         = false,
    hdf5_step    = init_step,
    vollock      = true,
    device       = :CUDA,
    coupling     = :OS,
    progressbar  = true,
    basis        = init_basis
)

# background grid setup
grid = Grid3D{iInt, iFloat}(
    NIC      = init_NIC,
    range_x1 = init_grid_range_x[1],
    range_x2 = init_grid_range_x[2],
    range_y1 = init_grid_range_y[1],
    range_y2 = init_grid_range_y[2],
    range_z1 = init_grid_range_z[1],
    range_z2 = init_grid_range_z[2],
    space_x  = init_grid_space_x,
    space_y  = init_grid_space_y,
    space_z  = init_grid_space_z,
    phase    = init_phase
)

# material points setup
range_x    = [0+grid.space_x/init_mp_in_space/2, 0.05-grid.space_x/init_mp_in_space/2]
range_y    = [0+grid.space_y/init_mp_in_space/2, 0.20-grid.space_y/init_mp_in_space/2]
range_z    = [0+grid.space_z/init_mp_in_space/2, 0.10-grid.space_z/init_mp_in_space/2]
space_x    = grid.space_x/init_mp_in_space
space_y    = grid.space_y/init_mp_in_space
space_z    = grid.space_z/init_mp_in_space
vx         = range_x[1]:space_x:range_x[2] |> collect
vy         = range_y[1]:space_y:range_y[2] |> collect
vz         = range_z[1]:space_z:range_z[2] |> collect
m, n, o    = length(vy), length(vx), length(vz)
vx         = reshape(vx, 1, n, 1)
vy         = reshape(vy, m, 1, 1)
vz         = reshape(vz, 1, 1, o)
om         = ones(Int, m)
on         = ones(Int, n)
oo         = ones(Int, o)
x_tmp      = vec(vx[om, :, oo])
y_tmp      = vec(vy[:, on, oo])
z_tmp      = vec(vz[om, on, :])
mp_num     = length(x_tmp)
mp_ρs      = ones(mp_num).*init_ρs
mp_layer   = ones(mp_num)
mp_ν       = [init_ν]
mp_E       = [init_E]
mp_G       = [init_G]
mp_σt      = [init_σt]
mp_ϕ       = [init_ϕ]
mp_c       = [init_c]
mp_ψ       = [init_ψ]
mp_Ks      = [init_Ks]
mp         = Particle3D{iInt, iFloat}(space_x=space_x, space_y=space_y, space_z=space_z,
    pos=[x_tmp y_tmp z_tmp], ν=mp_ν, E=mp_E, G=mp_G, σt=mp_σt, ϕ=mp_ϕ, c=mp_c, ψ=mp_ψ, 
    ρs=mp_ρs, Ks=mp_Ks, NIC=init_NIC, layer=mp_layer, phase=init_phase)

# boundary condition nodes index
vx_idx  = zeros(iInt, grid.node_num)
vy_idx  = zeros(iInt, grid.node_num)
vz_idx  = zeros(iInt, grid.node_num)
tmp_idx = findall(i->grid.pos[i, 1]≤0||grid.pos[i, 1]≥0.05||
                     grid.pos[i, 3]≤0||grid.pos[i, 2]≤0, 1:grid.node_num)
tmp_idy = findall(i->grid.pos[i, 2]≤0||grid.pos[i, 3]≤0, 1:grid.node_num)
tmp_idz = findall(i->grid.pos[i, 3]≤0, 1:grid.node_num)
vx_idx[tmp_idx] .= iInt(1)
vy_idx[tmp_idy] .= iInt(1)
vz_idx[tmp_idz] .= iInt(1)
bc = VBoundary3D{iInt, iFloat}(
    Vx_s_Idx = vx_idx,
    Vx_s_Val = zeros(grid.node_num),
    Vy_s_Idx = vy_idx,
    Vy_s_Val = zeros(grid.node_num),
    Vz_s_Idx = vz_idx,
    Vz_s_Val = zeros(grid.node_num)
)

# MPM solver
# CUDA.@profile external=true simulate!(args, grid, mp, bc)
simulate!(args, grid, mp, bc)
# savevtu(args, mp)

# effective memory throughput
nio   = mp.num*197+mp.num*35*mp.NIC+grid.node_num*33+2*grid.cell_num
Aeff  = nio*sizeof(typeof(mp).parameters[2])*(1/1024^3)
MTeff = @sprintf "%.2e" (Aeff*args.iter_num)/(args.end_time-args.start_time)
@info "MTeff: $(MTeff) GB/s"

#=
let
    fig = Figure()
    figtitle = join([args.project_name, "_", args.device])
    axis = Axis3(fig[1, 1], aspect=:data, azimuth=0.2π, elevation=0.08π)
    scatter!(axis, mp.pos, color=log10.(mp.epII.+1), colormap=:jet)
    display(fig)
    save(joinpath(rtsdir, args.project_name, "$(figtitle).png"), fig)
end
=#