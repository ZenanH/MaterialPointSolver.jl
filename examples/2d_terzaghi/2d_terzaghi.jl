#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 2d_terzaghi.jl                                                             |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : Terzaghi consolidation test                                                |
|  Reference  : Formulation of a Dynamic Material Point Method (MPM) for Geomechanical     |
|               Problems ISSAM K. J. AL-KAFAJI, 2013, Institut für Geotechnik der          |
|               Universität Stuttgart                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using HDF5
using CairoMakie
using Printf

include(joinpath(@__DIR__, "2d_funcs.jl"))

MaterialPointSolver.warmup(Val(:CPU))

# model configuration
init_grid_space_x = 0.05
init_grid_space_y = 0.05
init_grid_range_x = [0-init_grid_space_x*2, 0.1+init_grid_space_x*2]
init_grid_range_y = [0-init_grid_space_y*2, 1.0+init_grid_space_y*2]
init_mp_in_space  = 2
init_ρs           = 2650.0
init_ρw           = 1000.0
init_n            = 0.3
init_k            = 1e-3
init_ν            = 0
init_Es           = 1e7
init_Gs           = init_Es / (2 * (1 +     init_ν))
init_Ks           = init_Es / (3 * (1 - 2 * init_ν))
init_Kw           = 2.2e9
init_ΔT           = satΔt(init_Gs, init_Ks, init_Kw, init_ρs, init_ρw, init_n, init_k, init_grid_space_x)
init_T            = 2.0
init_Te           = 0.0
init_ΔT           = 1e-5
init_step         = floor(init_T / init_ΔT / 200)
init_basis        = :uGIMP
init_NIC          = 9
init_phase        = 2
init_σw           = -1e4
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
    gravity      = 0,
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
    ρs    = ones(size(ξ0, 1)) .* init_ρs,
    ρw    = ones(size(ξ0, 1)) .* init_ρw,
    n     = ones(size(ξ0, 1)) .* init_n,
)
mp.σw .= init_σw

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
tmp_idx = findall(i -> grid.ξ[i, 1] ≤ 0 || grid.ξ[i, 1] ≥ 0.1 || 
                       grid.ξ[i, 2] ≤ 0, 1:grid.ni)
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
    vy_w_idx = copy(vy_idx),
    vy_w_val = zeros(grid.ni),
    ext      = ext
)

# MPM solver
materialpointsolver!(args, grid, mp, attr, bc, workflow=Tprocedure!)

# post-processing
animation(args)

let 
    # figure setup =========================================================================
    figfont = MaterialPointSolver.tnr
    fig = Figure(size=(800, 350), fontsize=12, fonts=(; regular=figfont, bold=figfont))
    ax1 = Axis(fig[1, 1], aspect=1, xlabel="Normalized pore pressure 𝑝 [-]", 
        ylabel="Normalized depth 𝐻 [-]", xticks=(0:0.2:1), yticks=(0:0.2:1))
    text!(0.42, 0.75, text=L"T_v = 0.1")
    text!(0.55, 0.30, text=L"T_v = 0.3")
    text!(0.37, 0.10, text=L"T_v = 0.5")
    text!(0.00, 0.10, text=L"T_v = 0.7")
    ax2 = Axis(fig[1, 2], aspect=DataAspect(), ylabel=L"y\ (m)",
        xticks=(11:12),yticks=(0:0.2:1.0), ytrimspine=true, xlabel=L"T_v = 0.1",
        xgridvisible=false, ygridvisible=false, rightspinevisible=false,
        topspinevisible=false, bottomspinevisible=false)
    ax3 = Axis(fig[1, 3], aspect=DataAspect(), xlabel=L"T_v = 0.3",
        xticks=(11:12), xgridvisible=false, ygridvisible=false, rightspinevisible=false,
        topspinevisible=false, leftspinevisible=false, bottomspinevisible=false)
    ax4 = Axis(fig[1, 4], aspect=DataAspect(), xlabel=L"T_v = 0.5",
        xticks=(11:12), xgridvisible=false, ygridvisible=false, rightspinevisible=false,
        topspinevisible=false, leftspinevisible=false, bottomspinevisible=false)
    ax5 = Axis(fig[1, 5], aspect=DataAspect(), xlabel=L"T_v = 0.7",
        xticks=(11:12), xgridvisible=false, ygridvisible=false, rightspinevisible=false,
        topspinevisible=false, leftspinevisible=false, bottomspinevisible=false)
    limits!(ax1, -0.10, 1.1, -0.1, 1.1)
    limits!(ax2, -0.05, 0.1, -0.1, 1.1)
    limits!(ax3, -0.00, 0.1, -0.1, 1.1)
    limits!(ax4, -0.00, 0.1, -0.1, 1.1)
    limits!(ax5, -0.00, 0.1, -0.1, 1.1)
    colsize!(fig.layout, 1, Relative(0.4))
    colgap!(fig.layout, 2, 0)
    colgap!(fig.layout, 3, 0)
    colgap!(fig.layout, 4, 0)
    hideydecorations!(ax3, grid=false)
    hideydecorations!(ax4, grid=false)
    hideydecorations!(ax5, grid=false)

    # load data ============================================================================
    fid = h5open(joinpath(args.project_path, args.project_name, "$(args.project_name).h5"), "r")
    num = mp.np/length(unique(mp.ξ0[:, 1])) |> Int
    mp_rst = zeros(num, 2, 4)
    timeset = [11, 30, 50, 70]
    for i in eachindex(timeset)
        c_pp = fid["group$(timeset[i])/pressure_w"] |> read
        mp_rst[:, 1, i] .= reverse(c_pp[1:num])./init_σw
        mp_rst[:, 2, i] .= reverse(mp.ξ0[1:num, 2])
    end
    ξ1 = fid["group$(timeset[1])/coords"] |> read
    ξ2 = fid["group$(timeset[2])/coords"] |> read
    ξ3 = fid["group$(timeset[3])/coords"] |> read
    ξ4 = fid["group$(timeset[4])/coords"] |> read
    cvalue1 = (fid["group$(timeset[1])/pressure_w"] |> read) * 1e-3
    cvalue2 = (fid["group$(timeset[2])/pressure_w"] |> read) * 1e-3
    cvalue3 = (fid["group$(timeset[3])/pressure_w"] |> read) * 1e-3
    cvalue4 = (fid["group$(timeset[4])/pressure_w"] |> read) * 1e-3
    close(fid)

    # plot setup ===========================================================================
    p11 = lines!(ax1, terzaghi(init_σw, 0.1), color=:black, linewidth=1, 
        label="Analytical solution")
        
    p12 = lines!(ax1, terzaghi(init_σw, 0.3), color=:black, linewidth=1)
    p13 = lines!(ax1, terzaghi(init_σw, 0.5), color=:black, linewidth=1)
    p14 = lines!(ax1, terzaghi(init_σw, 0.7), color=:black, linewidth=1)
    p15 = scatterlines!(ax1, mp_rst[:, :, 1], linewidth=0.5, markersize=4, color=:red, 
        marker=:star8, strokewidth=0, label="MPM solution")
    p16 = scatterlines!(ax1, mp_rst[:, :, 2], linewidth=0.5, markersize=4, color=:red, 
        marker=:star8, strokewidth=0)
    p17 = scatterlines!(ax1, mp_rst[:, :, 3], linewidth=0.5, markersize=4, color=:red, 
        marker=:star8, strokewidth=0)
    p18 = scatterlines!(ax1, mp_rst[:, :, 4], linewidth=0.5, markersize=4, color=:red, 
        marker=:star8, strokewidth=0)
    axislegend(ax1, merge=true, padding=(10, 6, 0, 0))

    plt1 = scatter!(ax2, ξ1, color=cvalue1, markersize=10, marker=:rect, colormap=:viridis, 
        colorrange=(-10, 0))
    plt2 = scatter!(ax3, ξ2, color=cvalue2, markersize=10, marker=:rect, colormap=:viridis, 
        colorrange=(-10, 0))
    plt3 = scatter!(ax4, ξ3, color=cvalue3, markersize=10, marker=:rect, colormap=:viridis, 
        colorrange=(-10, 0))
    plt4 = scatter!(ax5, ξ4, color=cvalue4, markersize=10, marker=:rect, colormap=:viridis, 
        colorrange=(-10, 0))

    Colorbar(fig[1, 6], plt4, label=L"p\ \text{(kPa)}", size=5, spinewidth=0, vertical=true, 
        height=Relative(1/2.5))

    Label(fig[1, 1, Bottom()], "(a) Excess pore pressure isochrones", 
        halign=:center, padding=(0,0,-40,0))
    Label(fig[1, 2:6, Bottom()], "(b) Pore pressure distribution", 
        halign=:center, padding=(0,0,-40,0))
    display(fig)
    save(joinpath(args.project_path, args.project_name, "$(args.project_name).pdf"), fig)
end