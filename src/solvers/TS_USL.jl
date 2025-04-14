#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : TS_USL.jl                                                                  |
|  Description: Defaule USL (update stress last) implementation for two-phase              | 
|               single-point MPM                                                           |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
+==========================================================================================#

export solvegrid_USL_TS!

@kernel inbounds=true function solvegrid_USL_TS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2,
    ζw  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 2] = bc.vy_w_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) * 
                            (grid.fs[ix, 2] - grid.fw[ix, 2]))
        # compute node acceleration
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) - grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) - grid.fd[ix, 2])
        asx = ms_denom * (-grid.mi[ix] * awx + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * awy + grid.fs[ix, 2] + 
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        # update nodal temp velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        # apply boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vsT[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vsT[ix, 2] * ΔT
        grid.Δuw[ix, 1] = grid.vwT[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vwT[ix, 2] * ΔT
    end
end

@kernel inbounds=true function solvegrid_USL_TS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2,
    ζw  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.pw[ix, 3] = bc.vz_w_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mi_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mi_denom
        grid.vw[ix, 3] = grid.pw[ix, 3] * mi_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2]  +
                             grid.fw[ix, 3] * grid.fw[ix, 3] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) *
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) +
                            (grid.fs[ix, 3] - grid.fw[ix, 3]) *
                            (grid.fs[ix, 3] - grid.fw[ix, 3]))
        # compute node acceleration
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) - grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) - grid.fd[ix, 2])
        awz = mw_denom * (grid.fw[ix, 3] + dampvw * sign(grid.vw[ix, 3]) - grid.fd[ix, 3])
        asx = ms_denom * (-grid.mi[ix] * awx + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * awy + grid.fs[ix, 2] +
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        asz = ms_denom * (-grid.mi[ix] * awz + grid.fs[ix, 3] +
            dampvw * sign(grid.vw[ix, 3]) + dampvs * sign(grid.vs[ix, 3]))
        # update nodal temp velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + asx * ΔT
        grid.vsT[ix, 2] = grid.vs[ix, 2] + asy * ΔT
        grid.vsT[ix, 3] = grid.vs[ix, 3] + asz * ΔT
        grid.vwT[ix, 1] = grid.vw[ix, 1] + awx * ΔT
        grid.vwT[ix, 2] = grid.vw[ix, 2] + awy * ΔT
        grid.vwT[ix, 3] = grid.vw[ix, 3] + awz * ΔT
        # apply boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.vwT[ix, 3] = bc.vz_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vsT[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vsT[ix, 2] * ΔT
        grid.Δus[ix, 3] = grid.vsT[ix, 3] * ΔT
        grid.Δuw[ix, 1] = grid.vwT[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vwT[ix, 2] * ΔT
        grid.Δuw[ix, 3] = grid.vwT[ix, 3] * ΔT
    end
end

function procedure!(
    args::     DeviceArgs{T1, T2}, 
    grid::     DeviceGrid{T1, T2}, 
    mp  :: DeviceParticle{T1, T2}, 
    attr:: DeviceProperty{T1, T2},
    bc  ::DeviceVBoundary{T1, T2},
    ΔT  ::T2,
    Ti  ::T2,
        ::Val{:TS},
        ::Val{:USL}
) where {T1, T2}
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    # MPM procedure
    resetgridstatus_TS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_TS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_TS!(dev)(ndrange=mp.np, grid, mp, attr, G)
    solvegrid_USL_TS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs, args.ζw)
    doublemapping1_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    G2P_TS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT)
    # update stress status
    if args.constitutive == :hyperelastic
        hyE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :linearelastic
        liE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :druckerprager
        liE!(dev)(ndrange=mp.np, mp, attr)
        if Ti ≥ args.Te
            dpP!(dev)(ndrange=mp.np, mp, attr)
        end
    elseif args.constitutive==:mohrcoulomb
        liE!(dev)(ndrange=mp.np, mp, attr)
        if Ti ≥ args.Te
            mcP!(dev)(ndrange=mp.np, mp, attr)
        end
    end
    # cell-averaged volumetric locking elimination
    if args.MVL == true
        vollock1_TS!(dev)(ndrange=mp.np, grid, mp)
        vollock2_TS!(dev)(ndrange=mp.np, grid, mp)
    end
    return nothing
end