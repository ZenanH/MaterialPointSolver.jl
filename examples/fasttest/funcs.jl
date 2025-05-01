@kernel inbounds = true function tsolvegrid_OS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        μ = T2(0.36)
        ΔT_1 = inv(ΔT)
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.ps[ix, 3] = bc.vz_s_val[ix] : nothing
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        # damping force for solid
        dampvs = -ζs * sqrt(grid.fs[ix, 1] * grid.fs[ix, 1] + 
                            grid.fs[ix, 2] * grid.fs[ix, 2] + 
                            grid.fs[ix, 3] * grid.fs[ix, 3])
        # compute nodal total force for mixture
        Fs_x = grid.fs[ix, 1] + dampvs * sign(grid.vs[ix, 1])
        Fs_y = grid.fs[ix, 2] + dampvs * sign(grid.vs[ix, 2])
        Fs_z = grid.fs[ix, 3] + dampvs * sign(grid.vs[ix, 3])
        # update nodal velocity
        grid.vsT[ix, 1] = grid.vs[ix, 1] + Fs_x * ΔT * ms_denom
        grid.vsT[ix, 2] = grid.vs[ix, 2] + Fs_y * ΔT * ms_denom
        grid.vsT[ix, 3] = grid.vs[ix, 3] + Fs_z * ΔT * ms_denom
        # add friction boundary condition
        nv_id = grid.ext.id[ix]
        if nv_id ≠ T1(0)
            n2x = -grid.ext.n[nv_id, 1]
            n2y = -grid.ext.n[nv_id, 2]
            n2z = -grid.ext.n[nv_id, 3]
            if grid.vsT[ix, 1] * n2x > T2(0.0) || grid.vsT[ix, 2] * n2y > T2(0.0) ||
               grid.vsT[ix, 3] * n2z > T2(0.0)
                fc2x = grid.ms[ix] * -grid.vsT[ix, 1] * ΔT_1
                fc2y = grid.ms[ix] * -grid.vsT[ix, 2] * ΔT_1
                fc2z = grid.ms[ix] * -grid.vsT[ix, 3] * ΔT_1
                fn2x = fc2x * n2x * n2x + fc2y * n2y * n2x + fc2z * n2z * n2x
                fn2y = fc2x * n2x * n2y + fc2y * n2y * n2y + fc2z * n2z * n2y
                fn2z = fc2x * n2x * n2z + fc2y * n2y * n2z + fc2z * n2z * n2z
                ft2x, ft2y, ft2z = fc2x - fn2x, fc2y - fn2y, fc2z - fn2z
                ft_m = sqrt(ft2x * ft2x + ft2y * ft2y + ft2z * ft2z)
                fn_m = sqrt(fn2x * fn2x + fn2y * fn2y + fn2z * fn2z)
                if ft_m ≥ μ * fn_m
                    denom = ft_m < eps(T2) ? T2(0.0) : inv(ft_m)
                    fc2x = fn2x + (μ * fn_m * ft2x) * denom
                    fc2y = fn2y + (μ * fn_m * ft2y) * denom
                    fc2z = fn2z + (μ * fn_m * ft2z) * denom
                end
                grid.vsT[ix, 1] += fc2x * ΔT * ms_denom
                grid.vsT[ix, 2] += fc2y * ΔT * ms_denom
                grid.vsT[ix, 3] += fc2z * ΔT * ms_denom
            end
        end
        # boundary condition
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vsT[ix, 3] = bc.vz_s_val[ix] : nothing
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
    end
end

function Fprocedure!(
    args::     DeviceArgs{T1, T2}, 
    grid::     DeviceGrid{T1, T2}, 
    mp  :: DeviceParticle{T1, T2}, 
    attr:: DeviceProperty{T1, T2},
    bc  ::DeviceVBoundary{T1, T2},
    ΔT  ::T2,
    Ti  ::T2,
        ::Val{:OS},
        ::Val{:MUSL}
) where {T1, T2}
    G = Ti < args.Te ? args.gravity / args.Te * Ti : args.gravity
    dev = getBackend(Val(args.device))
    # MPM procedure
    resetgridstatus_OS!(dev)(ndrange=grid.ni, grid)
    resetmpstatus_OS!(dev)(ndrange=mp.np, grid, mp, Val(args.basis))
    P2G_OS!(dev)(ndrange=mp.np, grid, mp, G)
    tsolvegrid_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT, args.ζs)
    doublemapping1_OS!(dev)(ndrange=mp.np, grid, mp, attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_OS!(dev)(ndrange=mp.np, grid, mp)
    doublemapping3_OS!(dev)(ndrange=grid.ni, grid, bc, ΔT)
    # F-bar based volumetric locking elimination approach
    if args.MVL == false
        G2P_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    else
        G2Pvl1_OS!(dev)(ndrange=mp.np, grid, mp)
        fastdiv_OS!(dev)(ndrange=grid.ni, grid)
        G2Pvl2_OS!(dev)(ndrange=mp.np, grid, mp, ΔT)
    end
    # update stress status
    if args.constitutive == :hyperelastic
        hyE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :linearelastic
        liE!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :druckerprager
        liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && dpP!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :mohrcoulomb
        liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && mcP!(dev)(ndrange=mp.np, mp, attr)
    elseif args.constitutive == :bingham
        Ti < args.Te && liE!(dev)(ndrange=mp.np, mp, attr)
        Ti ≥ args.Te && bhP!(dev)(ndrange=mp.np, mp, attr, inv(ΔT))
    end
    return nothing
end