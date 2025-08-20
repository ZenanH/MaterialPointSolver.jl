#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_druckerprager.jl                                                        |
|  Description: Case used to vaildate the functions                                        |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

using MaterialPointGenerator
using MaterialPointSolver
using MaterialPointVisualizer
using CUDA

include(joinpath(@__DIR__, "funcs.jl"))

init_h     = 0.5
init_ϵ     = :double
init_mtl   = :druckerprager
init_basis = :bspline2
init_NIC   = 27
init_FLIP  = 1.0
init_G     = -9.8
init_ρs    = 2.5e3
init_ρw    = 1e3
init_ν     = 0.3
init_Es    = 1e6
init_Gs    = init_Es / (2 * (1 + init_ν))
init_Ks    = init_Es / (3 * (1 - 2 * init_ν))
init_Kw    = 1e6
init_c     = 2e4
init_cr    = 1e4
init_n     = 0.3
init_k     = 1e-8 * init_ρw * 9.8 / 1e-3
init_ϕ1    = deg2rad(20)
init_ϕ2    = deg2rad(10)
init_Hp    = 6e4
init_dev   = :cuda
init_T     = 15.0
init_Te    = 8.0
init_Tcur  = 0.0
Eu         = init_Es + init_Kw / init_n
cv         = init_k/(init_ρw*9.8*(1/init_Es+init_n)/init_Kw)
c1         = sqrt(Eu/(init_ρs*(1-init_n) + init_ρw*init_n))
init_ΔT    = 0.5 * init_h / c1
# init_ΔT    = 0.5 * init_h / sqrt(init_Es / init_ρs)
init_h5    = floor(Int, init_T / init_ΔT / 50)
init_var   = (:ξ, :ϵq, :σij, :Ω, :ρs, (:σw,))

# args setup
conf = init_conf(dev=init_dev, Δt=init_ΔT, t_tol=init_T, t_eld=init_Te, h5_int=init_h5, varnames=init_var,
    prjpath=@__DIR__, prjname="Slope", basis=init_basis, material=init_mtl)

# grid and boundary conditions setup
rangex  = -2*init_h:init_h:65+2*init_h
rangey  = -2*init_h:init_h:10+2*init_h
rangez  = -2*init_h:init_h:13+2*init_h
bg      = gridbuilder(rangex, rangey, rangez)
tmp_idx = findall(i -> bg.ξ[i, 1]≤0 || bg.ξ[i, 1]≥65 || bg.ξ[i, 3]≤0, 1:bg.ni)
tmp_idy = findall(i -> bg.ξ[i, 2]≤0 || bg.ξ[i, 2]≥10 || bg.ξ[i, 3]≤0, 1:bg.ni)
tmp_idz = findall(i -> bg.ξ[i, 3]≤0, 1:bg.ni)
vx_idx  = zeros(Int, bg.ni); vx_idx[tmp_idx] .= 1
vy_idx  = zeros(Int, bg.ni); vy_idx[tmp_idy] .= 1
vz_idx  = zeros(Int, bg.ni); vz_idx[tmp_idz] .= 1 
grid    = init_grid(rangex, rangey, rangez, NIC = init_NIC,
    vsxi=vx_idx, vsxv=zeros(bg.ni), vsyi=vy_idx, vsyv=zeros(bg.ni), vszi=vz_idx, vszv=zeros(bg.ni),
    ext = (Ω   = zeros(Float64, bg.nc),
           avg = zeros(Float64, bg.nc),
           mw  = zeros(Float64, bg.ni),
           mi  = zeros(Float64, bg.ni), 
           pw  = zeros(Float64, bg.ni, 3),
           fw  = zeros(Float64, bg.ni, 3),
           vw  = zeros(Float64, bg.ni, 3),
           vwT = zeros(Float64, bg.ni, 3),
           vwxi=vx_idx, vwxv=zeros(bg.ni), 
           vwyi=vy_idx, vwyv=zeros(bg.ni), 
           vwzi=vz_idx, vwzv=zeros(bg.ni))
)

# material point setup
mh = grid.h * 0.5
tmp1 = meshbuilder(mh*0.5 : mh : 65-mh*0.5,
                   mh*0.5 : mh : 10-mh*0.5,
                   mh*0.5 : mh : 13-mh*0.5)
idx = findall(i->(22≤tmp1[i, 1]≤32 && tmp1[i, 3]≥35-tmp1[i, 1]) ||
                 (32≤tmp1[i, 1] && tmp1[i, 3]≥3), 1:size(tmp1, 1))
tmpx = deleteat!(tmp1[:, 1], idx)
tmpy = deleteat!(tmp1[:, 2], idx)
tmpz = deleteat!(tmp1[:, 3], idx)
ξ0 = hcat(tmpx, tmpy, tmpz)
np = size(ξ0, 1)
nid = ones(Int, np)
bot = findall(i -> ξ0[i, 3] ≤ 3, 1:np)
nid[bot] .= 2
mpts = init_mpts(ϵ=init_ϵ, NIC = init_NIC,
    ξ    = ξ0,
    h    = mh,
    G    = init_G,
    FLIP = init_FLIP,
    nid  = nid,
    ρs   = fill(init_ρs, np),
    Ks   = [init_Ks, init_Ks],
    Es   = [init_Es, init_Es],
    Gs   = [init_Gs, init_Gs],
    ϕ    = [init_ϕ1, init_ϕ2],
    ψ    = [0.0    , 0.0    ],
    σt   = [0.0    , 0.0    ],
    c    = [init_c , init_c ],
    cr   = [init_cr, init_cr],
    Hp   = [init_Hp, init_Hp],
    ext = (
        σw = zeros(Float64, np),
        vw = zeros(Float64, np, 3),
        k = init_k,
        n = fill(init_n, np),
        ρw = init_ρw,
        Kw = init_Kw,)
)

# solver setup
mpmsolver!(tprocedure!, conf, grid, mpts)
h5conf = (prjdst=conf.prjdst, prjname=conf.prjname)
animation(h5conf)

# let     
#     custom_dark = theme_dark()
#     custom_dark.textcolor = :white
#     custom_dark.linecolor = :white
#     custom_dark.backgroundcolor = "#181818"
#     Makie.set_theme!(WGLMakie=(resize_to=:body,), custom_dark)
#     fig = Figure()
#     ax = LScene(fig[1, 1], show_axis=true)
#     canvas = ax.scene.plots[1]
#     canvas.ticks[:textcolor] = :white
#     canvas.frame[:axiscolor] = "#818181"
#     canvas.names[:textcolor] = :white
#     vcol = mpts.ext.σw #@. (mpts.σij[:, 1] + mpts.σij[:, 2] + mpts.σij[:, 3]) * 0.333333
#     scatter!(ax, mpts.ξ, color=vcol, colormap=:jet,
#         inspector_label = (self, i, pos) -> (
#         "x=" * string(pos[1]) *
#         "\ny=" * string(pos[2]) *
#         "\nz=" * string(pos[3]) *
#         "\nval=" * string(vcol[i])
#     ))
#     #scatter!(ax, bg, markersize=4, alpha=0.5)
#     DataInspector(fig, textcolor = :black,)
#     display(fig)
# end