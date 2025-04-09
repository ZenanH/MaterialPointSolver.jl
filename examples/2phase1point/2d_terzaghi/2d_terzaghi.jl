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
init_basis        = :linear
init_NIC          = 16
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
    hdf5_step    = init_step,
    MVL          = true,
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

#cv = init_k / (9800 * (1 / init_Es + init_n / init_Kw))
#t = 0.1 / cv

# MPM solver
materialpointsolver!(args, grid, mp, attr, bc, workflow=Tprocedure!)

# post-processing
animation(args)
begin
    # cv = init_k / (9800 * (1 / init_Es + init_n / init_Kw))
    # t = 0.1 / cv
    # helper functions =====================================================================
    function terzaghi(p0, Tv)
        num = 100
        H = 1
        Z = range(0, 1, length=num)
        data = zeros(num, 2)
        data[:, 2] .= Z
        @inbounds for i in 1:num
            p = 0.0
            for m in 1:2:1e4
                p += 4*p0/π*(1/m)*sin((m*π*data[i, 2])/(2*H))*exp((-m^2)*((π/2)^2)*Tv)
            end
            data[num+1-i, 1] = p/p0
        end
        return data
    end

    function consolidation()
        num = 1000
        dat = zeros(num, 2)    
        dat[:, 1] .= collect(range(0, 10, length=num))
        @inbounds for i in 1:num
            tmp = 0.0
            for m in 1:2:1e4
                tmp += (8/π^2)*(1/m^2)*exp(-(m*π/2)^2*dat[i, 1]) 
            end
            dat[i, 2] = 1-tmp
        end
        return dat
    end
    # figure setup =========================================================================
    figfont = MaterialPointSolver.tnr
    fig = Figure(size=(1000, 450), fontsize=15, fonts=(; regular=figfont, bold=figfont))
    titlay = fig[0, :] = GridLayout()
    Label(titlay[1, :], "Two-phase single-point MPM (𝑣-𝑤)\n2D Terzaghi Consolidation Test", 
        fontsize=20, tellwidth=false, halign=:center, justification=:center, lineheight=1.2)
    layout = fig[1, 1] = GridLayout()
    gd1 = layout[1, 1]
    gd2 = layout[1, 2]
    gd3 = layout[1, 4]
    cb1 = layout[1, 3]
    cb2 = layout[1, 5]
    colsize!(layout, 1, 348)
    # axis setup ===========================================================================
    ax1 = Axis(gd1, xlabel="Normalized pore pressure 𝑝 [-]", 
    ylabel="Normalized depth 𝐻 [-]", xticks=0:0.2:1, yticks=0:0.2:1, 
    title="Excess pore pressure isochrones", aspect=1)
    ax2 = Axis(gd2, aspect=DataAspect(), xlabel=L"x\ (m)", ylabel=L"y\ (m)", 
        xticks=0:0.1:0.1, yticks=0:0.2:1.2, title="Pore pressure distribution")
    ax3 = Axis(gd3, aspect=DataAspect(), xlabel=L"x\ (m)", ylabel=L"y\ (m)", 
        xticks=0:0.1:0.1, yticks=0:0.2:1.2, title="Isotropic stress distribution")
    limits!(ax1, -0.1, 1.1, -0.1, 1.1)
    limits!(ax2, -0.1, 0.3, -0.1, 1.3)
    limits!(ax3, -0.1, 0.3, -0.1, 1.3)
    # plot setup ===========================================================================
    p11 = lines!(ax1, terzaghi(init_σw, 0.1), color=:black, linewidth=1, 
        label="Analytical solution")
    p12 = lines!(ax1, terzaghi(init_σw, 0.3), color=:black, linewidth=1)
    p13 = lines!(ax1, terzaghi(init_σw, 0.5), color=:black, linewidth=1)
    p14 = lines!(ax1, terzaghi(init_σw, 0.7), color=:black, linewidth=1)
    prj = joinpath(args.project_path, args.project_name)
    fid = h5open(joinpath(prj, "$(args.project_name).h5"), "r")
    num = mp.np/length(unique(mp.ξ0[:, 1])) |> Int
    mp_rst = zeros(num, 2, 4)

    cv = init_k / (9800 * (1 / init_Es + init_n / init_Kw))
    tv = [0.1, 0.3, 0.5, 0.7]
    t  = tv ./ cv
    timeset = @. round(Int, t / (init_step * init_ΔT))

    timeset = [11, 30, 51, 70]

    #timeset = [8, 23, 38, 53]
    for i in eachindex(timeset)
        c_pp = fid["group$(timeset[i])/pressure_w"] |> read
        mp_rst[:, 1, i] .= reverse(c_pp[1:num])./init_σw
        mp_rst[:, 2, i] .= reverse(mp.ξ0[1:num, 2])
    end
    close(fid)
    p15 = scatterlines!(ax1, mp_rst[:, :, 1], linewidth=0.5, markersize=6, color=:red, 
        marker=:star8, strokewidth=0, label="MPM solution")
    p16 = scatterlines!(ax1, mp_rst[:, :, 2], linewidth=0.5, markersize=6, color=:red, 
        marker=:star8, strokewidth=0)
    p17 = scatterlines!(ax1, mp_rst[:, :, 3], linewidth=0.5, markersize=6, color=:red, 
        marker=:star8, strokewidth=0)
    p18 = scatterlines!(ax1, mp_rst[:, :, 4], linewidth=0.5, markersize=6, color=:red, 
        marker=:star8, strokewidth=0)
    axislegend(ax1, merge=true, labelsize=14, padding=(10, 6, 0, 0))
    #---------------------------------------------------------------------------------------
    p21 = scatter!(ax2, mp.ξ, color=mp.σw./1000 , markersize=10, marker=:rect, 
        colormap=:viridis, colorrange=(-1, 0))
    limits!(ax2, -0.05, 0.15, -0.05, 1.10)
    #---------------------------------------------------------------------------------------
    p31 = scatter!(ax3, mp.ξ, color=mp.vs[:, 2], markersize=10, marker=:rect, 
        colormap=:viridis )
    limits!(ax3, -0.05, 0.15, -0.05, 1.10)
    # colorbar setup =======================================================================
    Colorbar(cb1, p21, label=L"\sigma_w\ \text{(kPa)}", size=5, spinewidth=0, vertical=true, 
        height=Relative(1/1.5))
    Colorbar(cb2, p31, label=L"\sigma\ \text{(Pa)}", size=5, spinewidth=0, vertical=true, 
        height=Relative(1/1.5))
    # save figure ==========================================================================
    display(fig)
    save(joinpath(prj, "$(args.project_name).png"), fig)
    @info "Figure saved in project path"
end