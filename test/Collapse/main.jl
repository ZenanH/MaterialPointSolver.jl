#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: 3D granular collapse test                                                  |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using CUDA
# using CairoMakie

init_h     = 0.0025
init_ϵ     = :double
init_basis = :bspline2
init_NIC   = 27
init_FLIP  = 1.0
init_G     = -9.8
init_ρs    = 2650
init_ν     = 0.3
init_Ks    = 7e5
init_Es    = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs    = init_Es / (2 * (1 +     init_ν))
init_σt    = 0.0
init_ϕ     = deg2rad(19.8)
init_dev   = :cuda
init_T     = 0.6
init_Tcur  = 0.0
init_ΔT    = 0.5 * init_h / sqrt(init_Es / init_ρs)
init_h5    = floor(Int, init_T / init_ΔT / 200)
init_var   = (:ξ, :ϵq)

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T, h5_int=init_h5, varnames=init_var,
    prjpath=@__DIR__, prjname="Collapse", basis=init_basis)

# grid and boundary conditions setup
rangex  = -0.02:init_h:0.07
rangey  = -0.02:init_h:0.75
rangez  = -0.02:init_h:0.12
bg      = meshbuilder(rangex, rangey, rangez)
ni      = size(bg, 1)
tmp_idx = findall(i -> bg[i, 1]≤0||bg[i, 1]≥0.05||bg[i, 3]≤0||bg[i, 2]≤0, 1:ni)
tmp_idy = findall(i -> bg[i, 2]≤0||bg[i, 3]≤0.00, 1:ni)
tmp_idz = findall(i -> bg[i, 3]≤0, 1:ni)
vx_idx  = zeros(Int, ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, ni); vy_idx[tmp_idy] .= 1
vz_idx  = zeros(Int, ni); vz_idx[tmp_idz] .= 1 
grid    = init_grid(rangex, rangey, rangez, ϵ=init_ϵ, NIC=init_NIC,
    vsxi=vx_idx, vsxv=zeros(ni), vsyi=vy_idx, vsyv=zeros(ni), vszi=vz_idx, vszv=zeros(ni))

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder(mh*0.5 : mh : 0.05-mh*0.5, 
                 mh*0.5 : mh : 0.20-mh*0.5, 
                 mh*0.5 : mh : 0.10-mh*0.5) # xrange: 0:0.05, yrange: 0:0.2, zrange: 0:0.1
np = size(ξ0, 1)
mpts = init_mpts(ϵ=init_ϵ, NIC=init_NIC,
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = ones(Int, np),
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks],
    Es   = [init_Es],
    Gs   = [init_Gs],
    ϕ    = [init_ϕ]
)

# solver setup
mpmsolver!(procedure!, conf, grid, mpts)
h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
animation(h5conf)

# let
#     set_theme!(theme_latexfonts())
#     fig = Figure(size=(1200, 700), fontsize=30)
#     ax = Axis3(fig[1, 1], xlabel=L"x\ (m)", ylabel=L"y\ (m)", zlabel=L"z\ (m)", 
#         aspect=:data, azimuth=0.2*π, elevation=0.1*π, xlabeloffset=60, zlabeloffset=80,
#         protrusions=100, xticks=(0:0.04:0.04), height=450, width=950)
#     pl1 = scatter!(ax, mpts.ξ, color=log10.(mpts.ϵq.+1), colormap=:jet, markersize=3,
#         colorrange=(0, 1))
#     Colorbar(fig[1, 1], limits=(0, 1), colormap=:jet, size=16, ticks=0:0.5:1, spinewidth=0,
#         label=L"log_{10}(\epsilon_{II}+1)", vertical=false, tellwidth=false, width=200,
#         halign=:right, valign=:top, flipaxis=false)
#     display(fig)
#     save("Collapse.png", fig, px_per_unit=2.0)
# end