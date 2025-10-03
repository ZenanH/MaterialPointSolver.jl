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
using CairoMakie

include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "funcs.jl"))

init_h     = 0.8
base, pts  = loadmodel(init_h)
init_ϵ     = :double
init_basis = :bspline2
init_mtl   = :glen_nye
init_NIC   = 27
init_FLIP  = 1.0
init_G     = -9.8
init_ρs    = 9e2
init_ν     = 0.3
init_Ks    = 1.3e7
init_Es    = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs    = init_Es / (2 * (1 +     init_ν))
init_μ     = 0.3
init_dev   = :cuda
init_T     = 7
init_Te    = 1
init_ΔT    = 0.5 * init_h / sqrt(init_Ks / init_ρs)
init_h5    = floor(Int, init_T / init_ΔT / 100)
init_var   = (:ξ, :vs, (:σm,))

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T, t_eld=init_Te,h5_int=init_h5, varnames=init_var,
    prjpath=@__DIR__, prjname="Glacier", basis=init_basis, material=init_mtl)

# grid and boundary conditions setup
rangex  = -2*init_h:init_h:200+2*init_h
rangey  = -2*init_h:init_h:200+2*init_h
rangez  = -6*init_h:init_h: 50+2*init_h
bg      = gridbuilder(rangex, rangey, rangez)
n, id   = get_ext(base, bg, 5)
grid    = init_grid(rangex, rangey, rangez, NIC = init_NIC,
    ext = (n=Array{Float64}(n), id=Array{Int64}(id), μ=init_μ))

# material point setup
np = size(pts, 1)
ϵ0̇ = 1e-2
Am = ϵ0̇^2 / (8*(0.5*init_ρs*init_h^2*inv(init_ΔT))^3)
mpts = init_mpts(ϵ=init_ϵ, NIC=init_NIC,
    ξ    = pts,
    h    = grid.h * 0.5,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = ones(Int, np),
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks],
    Es   = [init_Es],
    Gs   = [init_Gs],
    ext  = (σm = zeros(np), ϵ0̇=Float64(ϵ0̇), A=Float64(1e-14))
)

# solver setup
mpmsolver!(tprocedure!, conf, grid, mpts)
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