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

using MaterialPointSolver.Mechanics # import the extension, `mpm_mechanics!` is from this extension

init_h     = 0.0025
init_FLIP  = 1.0
init_ζs    = 0.0
init_G     = -9.8
init_ρs    = 2650
init_ν     = 0.3
init_Ks    = 7e5
init_Es    = init_Ks * (3 * (1 - 2 * init_ν))
init_Gs    = init_Es / (2 * (1 +     init_ν))
init_σt    = 0.0
init_c     = 0.0
init_ϕ     = deg2rad(19.8)
init_ψ     = deg2rad(0.0)
init_dev   = :cuda
init_T     = 0.6
init_Tcur  = 0.0
init_ΔT    = 0.5 * init_h / sqrt(init_Es / init_ρs)
init_h5    = 200
init_var   = (:ξ, (:ϵq, :σij))

# args setup
conf = init_conf(
    dev      = init_dev,
    Δt       = init_ΔT,
    t_tol    = init_T,
    h5_int   = init_h5,
    varnames = init_var,
    prjpath  = joinpath(@__DIR__, ".output"),
    prjname  = "Collapse",
    adaptive = false,
)

# grid and boundary conditions setup
margin = 0.
rangex = [-0.02-margin, 0.07+margin]
rangey = [-0.02-margin, 0.75+margin]
rangez = [-0.02-margin, 0.12+margin]
bg     = gridbuilder(rangex, rangey, rangez, init_h)
grid   = init_grid(bg, ζs=init_ζs, 
    ext = (
        ps = zeros(bg.ni, 3),
        vs = zeros(bg.ni, 3),
        vsT= zeros(bg.ni, 3),
        fs = zeros(bg.ni, 3),
        ms = zeros(bg.ni)
    )
)

# material point setup
mh = grid.h * 0.5
ξ0 = meshbuilder([mh*0.5, 0.05-mh*0.5], 
                 [mh*0.5, 0.20-mh*0.5], 
                 [mh*0.5, 0.10-mh*0.5], mh) # xrange: 0:0.05, yrange: 0:0.2, zrange: 0:0.1
np = size(ξ0, 1)
mpts = init_mpts(
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = ones(Int, np),
    ext  = ( 
        σij = zeros(np, 6),
        Ks  = [init_Ks],
        Es  = [init_Es],
        Gs  = [init_Gs],
        ϕ   = [init_ϕ],
        ψ   = [init_ψ],
        σt  = [init_σt],
        c   = [init_c],
        ϵq  = zeros(np),
        ϵk  = zeros(np),
        ρs  = fill(Float64(init_ρs), np),
        ρs0 = Float64(init_ρs),
        Ω   = fill(Float64(mh^3), np),
        Ω0  = Float64(mh^3),
        vs  = zeros(np, 3),)
)

# solver setup
mpmsolver!(mpm_mechanics!, conf, grid, mpts)
h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
animation(h5conf)