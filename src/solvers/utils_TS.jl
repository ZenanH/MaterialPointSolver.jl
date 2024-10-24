#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : TS.jl                                                                      |
|  Description: Basic computing functions for two-phase single-point MPM                   |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. resetgridstatus_TS! [2D]                                               |
|               02. resetgridstatus_TS! [3D]                                               |
|               03. resetmpstatus_TS!   [2D, linear basis]                                 |
|               04. resetmpstatus_TS!   [3D, linear basis]                                 |
|               05. resetmpstatus_TS!   [2D,  uGIMP basis]                                 |
|               06. resetmpstatus_TS!   [3D,  uGIMP basis]                                 |
|               07. P2G_TS!             [2D]                                               |
|               08. P2G_TS!             [3D]                                               |
|               09. solvegrid_TS!       [2D]                                               |
|               10. solvegrid_TS!       [3D]                                               |
|               11. doublemapping1_TS!  [2D]                                               |
|               12. doublemapping1_TS!  [3D]                                               |
|               13. doublemapping2_TS!  [2D]                                               |
|               14. doublemapping2_TS!  [3D]                                               |
|               15. doublemapping3_TS!  [2D]                                               |
|               16. doublemapping3_TS!  [3D]                                               |
|               17. G2P_TS!             [2D]                                               |
|               18. G2P_TS!             [3D]                                               |
|               19. vollock1_TS!        [2D]                                               |
|               20. vollock1_TS!        [3D]                                               |
|               21. vollock2_TS!        [2D]                                               |
|               22. vollock2_TS!        [3D]                                               |
+==========================================================================================#

export resetgridstatus_TS! 
export resetmpstatus_TS!
export P2G_TS!
export solvegrid_TS! 
export doublemapping1_TS!
export doublemapping2_TS!
export doublemapping3_TS!
export G2P_TS! 
export vollock1_TS!
export vollock2_TS!

"""
    resetgridstatus_TS!(grid::DeviceGrid2D{T1, T2})

Description:
---
Reset some variables for the grid.
"""
@kernel inbounds=true function resetgridstatus_TS!(
    grid::DeviceGrid2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        if ix ≤ grid.nc
            grid.σm[ix] = T2(0.0)
            grid.σw[ix] = T2(0.0)
            grid.Ω[ix]  = T2(0.0)
        end
        grid.ms[ix]    = T2(0.0)
        grid.mi[ix]    = T2(0.0)
        grid.mw[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
        grid.fw[ix, 1] = T2(0.0)
        grid.fw[ix, 2] = T2(0.0)
        grid.fd[ix, 1] = T2(0.0)
        grid.fd[ix, 2] = T2(0.0)
    end
end

"""
    resetgridstatus_TS!(grid::DeviceGrid3D{T1, T2})

Description:
---
Reset some variables for the grid.
"""
@kernel inbounds=true function resetgridstatus_TS!(
    grid::DeviceGrid3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        if ix ≤ grid.nc
            grid.σm[ix] = T2(0.0)
            grid.σw[ix] = T2(0.0)
            grid.Ω[ix]  = T2(0.0)
        end
        grid.ms[ix]    = T2(0.0)
        grid.mi[ix]    = T2(0.0)
        grid.mw[ix]    = T2(0.0)
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.pw[ix, 3] = T2(0.0)
        grid.fs[ix, 1] = T2(0.0)
        grid.fs[ix, 2] = T2(0.0)
        grid.fs[ix, 3] = T2(0.0)
        grid.fw[ix, 1] = T2(0.0)
        grid.fw[ix, 2] = T2(0.0)
        grid.fw[ix, 3] = T2(0.0)
        grid.fd[ix, 1] = T2(0.0)
        grid.fd[ix, 2] = T2(0.0)
        grid.fd[ix, 3] = T2(0.0)
    end
end

@kernel inbounds=true function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        # compute particle to cell and particle to node index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            p2n = getP2N_linear(grid, mp.p2c[ix], iy)
            mp.p2n[ix, iy] = p2n
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            # compute basis function
            Nx, dNx = linearBasis(Δdx, grid.dx)
            Ny, dNy = linearBasis(Δdy, grid.dy)
            mp.Nij[ix, iy] =  Nx * Ny
            mp.∂Nx[ix, iy] = dNx * Ny # x-gradient shape function
            mp.∂Ny[ix, iy] = dNy * Nx # y-gradient shape function
        end
    end
end

@kernel inbounds=true function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:linear}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update momentum and mass
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
        # compute particle to cell and particle to node index
        mp.p2c[ix] = unsafe_trunc(T1, 
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 3] - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            p2n = getP2N_linear(grid, mp.p2c[ix], iy)
            mp.p2n[ix, iy] = p2n
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            Δdz = mp.ξ[ix, 3] - grid.ξ[p2n, 3]
            # compute basis function
            Nx, dNx = linearBasis(Δdx, grid.dx)
            Ny, dNy = linearBasis(Δdy, grid.dy)
            Nz, dNz = linearBasis(Δdz, grid.dz)
            mp.Nij[ix, iy] =  Nx * Ny * Nz
            mp.∂Nx[ix, iy] = dNx * Ny * Nz # x-gradient shape function
            mp.∂Ny[ix, iy] = dNy * Nx * Nz # y-gradient shape function
            mp.∂Nz[ix, iy] = dNz * Nx * Ny # z-gradient shape function
        end
    end
end

"""
    resetmpstatus_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2},
       ::Val{:uGIMP})

Description:
---
1. Get topology between particle and grid.
2. Compute the value of basis function (uGIMP).
3. Update particle mass and momentum.
"""
@kernel inbounds=true function resetmpstatus_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update mass and momentum
        mp.ms[ix] = mp.Ω[ix] * mp.ρs[ix]
        mp.mw[ix] = mp.Ω[ix] * mp.ρw[ix]
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mp.ξ[ix, 2] - grid.y1, grid.dy) +
            fld(mp.ξ[ix, 1] - grid.x1, grid.dx) * grid.ncy)
        # reset Nij, ∂Nx, ∂Ny
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(16)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance between particle and related nodes
            Δdx = mp.ξ[ix, 1] - grid.ξ[p2n, 1]
            Δdy = mp.ξ[ix, 2] - grid.ξ[p2n, 2]
            # compute basis function
            if abs(Δdx) < (grid.dx + T2(0.5) * mp.dx) &&
               abs(Δdy) < (grid.dy + T2(0.5) * mp.dy)
                Nx, dNx = uGIMPbasis(Δdx, grid.dx, mp.dx)
                Ny, dNy = uGIMPbasis(Δdy, grid.dy, mp.dy)
                mp.Nij[ix, viy] =  Nx * Ny
                mp.∂Nx[ix, viy] = dNx * Ny # x-gradient shape function
                mp.∂Ny[ix, viy] = dNy * Nx # y-gradient shape function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

"""
    resetmpstatus_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, 
        ::Val{:uGIMP})

Description:
---
1. Get topology between particle and grid.
2. Compute the value of basis function (uGIMP).
3. Update particle mass and momentum.
"""
@kernel inbounds=true function resetmpstatus_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
        ::Val{:uGIMP}
) where {T1, T2}
    ix = @index(Global)   
    if ix ≤ mp.np
        mpξ1 = mp.ξ[ix, 1]
        mpξ2 = mp.ξ[ix, 2]
        mpξ3 = mp.ξ[ix, 3]
        mpms = mp.Ω[ix] * mp.ρs[ix]
        mpmw = mp.Ω[ix] * mp.ρw[ix]
        # update particle mass and momentum
        mp.ms[ix] = mpms
        mp.mw[ix] = mpmw
        mp.mi[ix] = mp.Ω[ix] * ((T2(1.0) - mp.n[ix]) * mp.ρs[ix] + mp.n[ix] * mp.ρw[ix])
        mp.ps[ix, 1] = mpms * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mpms * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mpms * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mpmw * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mpmw * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mpmw * mp.vw[ix, 3] *            mp.n[ix]
        # p2c index
        mp.p2c[ix] = unsafe_trunc(T1,
            cld(mpξ2 - grid.y1, grid.dy) +
            fld(mpξ3 - grid.z1, grid.dz) * grid.ncy * grid.ncx +
            fld(mpξ1 - grid.x1, grid.dx) * grid.ncy)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            mp.Nij[ix, iy] = T2(0.0)
            mp.∂Nx[ix, iy] = T2(0.0)
            mp.∂Ny[ix, iy] = T2(0.0)
            mp.∂Nz[ix, iy] = T2(0.0)
        end
        viy = T1(1)
        @KAunroll for iy in Int32(1):Int32(64)
            # p2n index
            p2n = getP2N_uGIMP(grid, mp.p2c[ix], iy)
            # compute distance betwe en particle and related nodes
            Δdx = mpξ1 - grid.ξ[p2n, 1]
            Δdy = mpξ2 - grid.ξ[p2n, 2]
            Δdz = mpξ3 - grid.ξ[p2n, 3]
            # compute basis function
            if abs(Δdx) < (grid.dx + T2(0.5) * mp.dx) &&
               abs(Δdy) < (grid.dy + T2(0.5) * mp.dy) &&
               abs(Δdz) < (grid.dz + T2(0.5) * mp.dz)
                Nx, dNx = uGIMPbasis(Δdx, grid.dx, mp.dx)
                Ny, dNy = uGIMPbasis(Δdy, grid.dy, mp.dy)
                Nz, dNz = uGIMPbasis(Δdz, grid.dz, mp.dz)
                mp.Nij[ix, viy] =  Nx * Ny * Nz
                mp.∂Nx[ix, viy] = dNx * Ny * Nz # x-gradient basis function
                mp.∂Ny[ix, viy] = dNy * Nx * Nz # y-gradient basis function
                mp.∂Nz[ix, viy] = dNz * Nx * Ny # z-gradient basis function
                mp.p2n[ix, viy] = p2n
                viy += T1(1)
            end
            viy > mp.NIC && break
        end
    end
end

"""
    P2G_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, gravity::T2)

Description:
---
P2G procedure for scattering the mass, momentum, and forces from particles to grid.
"""
@kernel inbounds=true function P2G_TS!(
    grid   ::    DeviceGrid2D{T1, T2},
    mp     ::DeviceParticle2D{T1, T2},
    attr   ::  DeviceProperty{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni  = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                vol = mp.Ω[ix]
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                p2n = mp.p2n[ix, iy]
                tmp_drag = (mp.n[ix] * mp.mw[ix] * gravity) / attr.k[attr.nid[ix]]
                # compute nodal mass
                @KAatomic grid.ms[p2n] += Ni * mp.ms[ix] * (T2(1.0) - mp.n[ix])
                @KAatomic grid.mi[p2n] += Ni * mp.mw[ix] *            mp.n[ix]
                @KAatomic grid.mw[p2n] += Ni * mp.mw[ix]
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mp.ps[ix, 1]
                @KAatomic grid.ps[p2n, 2] += Ni * mp.ps[ix, 2]
                @KAatomic grid.pw[p2n, 1] += Ni * mp.pw[ix, 1]
                @KAatomic grid.pw[p2n, 2] += Ni * mp.pw[ix, 2]
                # compute nodal total force
                @KAatomic grid.fw[p2n, 1] += -vol * (∂Nx * mp.σw[ix])
                @KAatomic grid.fw[p2n, 2] += -vol * (∂Ny * mp.σw[ix]) +
                                               Ni * mp.mw[ix] * gravity
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (mp.σij[ix, 1]  + mp.σw[ix]) + 
                                                     ∂Ny *  mp.σij[ix, 4])
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (mp.σij[ix, 2]  + mp.σw[ix]) + 
                                                     ∂Nx *  mp.σij[ix, 4]) +
                                               Ni * mp.mi[ix] * gravity
                @KAatomic grid.fd[p2n, 1] += Ni * tmp_drag * (mp.vw[ix, 1] - mp.vs[ix, 1])
                @KAatomic grid.fd[p2n, 2] += Ni * tmp_drag * (mp.vw[ix, 2] - mp.vs[ix, 2])
            end
        end
    end
end

"""
    P2G_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, gravity::T2)

Description:
---
P2G procedure for scattering the mass, momentum, and forces from particles to grid.
"""
@kernel inbounds=true function P2G_TS!(
    grid   ::    DeviceGrid3D{T1, T2},
    mp     ::DeviceParticle3D{T1, T2},
    attr   ::  DeviceProperty{T1, T2},
    gravity::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni  = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                ∂Nz = mp.∂Nz[ix, iy]
                p2n = mp.p2n[ix, iy]
                vol = mp.Ω[ix]
                tdg = (mp.n[ix] * mp.mw[ix] * gravity) / attr.k[attr.nid[ix]]
                # compute nodal mass
                @KAatomic grid.ms[p2n] += Ni * mp.ms[ix] * (T2(1.0) - mp.n[ix])
                @KAatomic grid.mi[p2n] += Ni * mp.mw[ix] *            mp.n[ix]
                @KAatomic grid.mw[p2n] += Ni * mp.mw[ix]
                # compute nodal momentum
                @KAatomic grid.ps[p2n, 1] += Ni * mp.ps[ix, 1]
                @KAatomic grid.ps[p2n, 2] += Ni * mp.ps[ix, 2]
                @KAatomic grid.ps[p2n, 3] += Ni * mp.ps[ix, 3]
                @KAatomic grid.pw[p2n, 1] += Ni * mp.pw[ix, 1]
                @KAatomic grid.pw[p2n, 2] += Ni * mp.pw[ix, 2]
                @KAatomic grid.pw[p2n, 3] += Ni * mp.pw[ix, 3]
                # compute nodal total force
                @KAatomic grid.fw[p2n, 1] += -vol * (∂Nx * mp.σw[ix])
                @KAatomic grid.fw[p2n, 2] += -vol * (∂Ny * mp.σw[ix])
                @KAatomic grid.fw[p2n, 3] += -vol * (∂Nz * mp.σw[ix]) + 
                                              Ni  * mp.mw[ix] * gravity
                @KAatomic grid.fs[p2n, 1] += -vol * (∂Nx * (mp.σij[ix, 1]  + mp.σw[ix]) + 
                                                     ∂Ny *  mp.σij[ix, 4]  + 
                                                     ∂Nz *  mp.σij[ix, 6])
                @KAatomic grid.fs[p2n, 2] += -vol * (∂Ny * (mp.σij[ix, 2]  + mp.σw[ix]) + 
                                                     ∂Nx *  mp.σij[ix, 4]  + 
                                                     ∂Nz *  mp.σij[ix, 5])
                @KAatomic grid.fs[p2n, 3] += -vol * (∂Nz * (mp.σij[ix, 3]  + mp.σw[ix]) + 
                                                     ∂Nx *  mp.σij[ix, 6]  + 
                                                     ∂Ny *  mp.σij[ix, 5]) + 
                                              Ni  * mp.mi[ix] * gravity
                @KAatomic grid.fd[p2n, 1] += Ni * tdg * (mp.vw[ix, 1] - mp.vs[ix, 1])
                @KAatomic grid.fd[p2n, 2] += Ni * tdg * (mp.vw[ix, 2] - mp.vs[ix, 2])
                @KAatomic grid.fd[p2n, 3] += Ni * tdg * (mp.vw[ix, 3] - mp.vs[ix, 3])
            end
        end
    end
end


"""
    solvegrid_TS!(grid::DeviceGrid2D{T1, T2}, bc::DeviceVBoundary2D{T1, T2}, ΔT::T2, ζs::T2,
        ζw::T2)

Description:
---
1. Solve equations on grid.
2. Add boundary condition.
"""
@kernel inbounds=true function solvegrid_TS!(
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
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        # compute damping force
        dampvw = -ζw * sqrt( grid.fw[ix, 1] * grid.fw[ix, 1]  + 
                             grid.fw[ix, 2] * grid.fw[ix, 2] )
        dampvs = -ζs * sqrt((grid.fs[ix, 1] - grid.fw[ix, 1]) * 
                            (grid.fs[ix, 1] - grid.fw[ix, 1]) +
                            (grid.fs[ix, 2] - grid.fw[ix, 2]) * 
                            (grid.fs[ix, 2] - grid.fw[ix, 2]))
        # compute node acceleration
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) + grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) + grid.fd[ix, 2])
        asx = ms_denom * (-grid.mi[ix] * grid.aw[ix, 1] + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * grid.aw[ix, 2] + grid.fs[ix, 2] + 
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
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
    end
end

"""
    solvegrid_TS!(grid::DeviceGrid3D{T1, T2}, bc::DeviceVBoundary3D{T1, T2}, ΔT::T2, ζs::T2,
        ζw::T2)

Description:
---
1. Solve equations on grid.
2. Add boundary condition.
"""
@kernel inbounds=true function solvegrid_TS!(
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
        # compute nodal velocity
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        grid.vw[ix, 3] = grid.pw[ix, 3] * mw_denom
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
        awx = mw_denom * (grid.fw[ix, 1] + dampvw * sign(grid.vw[ix, 1]) + grid.fd[ix, 1])
        awy = mw_denom * (grid.fw[ix, 2] + dampvw * sign(grid.vw[ix, 2]) + grid.fd[ix, 2])
        awz = mw_denom * (grid.fw[ix, 3] + dampvw * sign(grid.vw[ix, 3]) + grid.fd[ix, 3])
        asx = ms_denom * (-grid.mi[ix] * grid.aw[ix, 1] + grid.fs[ix, 1] + 
            dampvw * sign(grid.vw[ix, 1]) + dampvs * sign(grid.vs[ix, 1]))
        asy = ms_denom * (-grid.mi[ix] * grid.aw[ix, 2] + grid.fs[ix, 2] +
            dampvw * sign(grid.vw[ix, 2]) + dampvs * sign(grid.vs[ix, 2]))
        asz = ms_denom * (-grid.mi[ix] * grid.aw[ix, 3] + grid.fs[ix, 3] +
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
        # reset grid momentum
        grid.ps[ix, 1] = T2(0.0)
        grid.ps[ix, 2] = T2(0.0)
        grid.ps[ix, 3] = T2(0.0)
        grid.pw[ix, 1] = T2(0.0)
        grid.pw[ix, 2] = T2(0.0)
        grid.pw[ix, 3] = T2(0.0)
    end
end

"""
    doublemapping1_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, ΔT::T2, 
        FLIP::T2, PIC::T2)

Description:
---
1. Solve equations on grid.
2. Compute CFL conditions.
"""
@kernel inbounds=true function doublemapping1_TS!(
    grid    ::    DeviceGrid2D{T1, T2},
    mp      ::DeviceParticle2D{T1, T2},
    ΔT      ::T2,
    FLIP    ::T2,
    PIC     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξsx = ξsy = ξwx = ξwy = vsx = vsy = vwx = vwy = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξsx += Ni *  grid.vsT[p2n, 1]
                ξsy += Ni *  grid.vsT[p2n, 2]
                ξwx += Ni *  grid.vwT[p2n, 1]
                ξwy += Ni *  grid.vwT[p2n, 2]
                vsx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vsy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
                vwx += Ni * (grid.vwT[p2n, 1] - grid.vw[p2n, 1])
                vwy += Ni * (grid.vwT[p2n, 2] - grid.vw[p2n, 2])
            end
        end
        # update particle ξition
        mp.ξ[ix, 1] += ΔT * ξsx
        mp.ξ[ix, 2] += ΔT * ξsy
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vsx) + PIC * ξsx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vsy) + PIC * ξsy
        mp.vw[ix, 1] = FLIP * (mp.vw[ix, 1] + vwx) + PIC * ξwx
        mp.vw[ix, 2] = FLIP * (mp.vw[ix, 2] + vwy) + PIC * ξwy
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
    end
end

"""
    doublemapping1_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, ΔT::T2, 
        FLIP::T2, PIC::T2)

Description:
---
1. Solve equations on grid.
2. Compute CFL conditions.
"""
@kernel inbounds=true function doublemapping1_TS!(
    grid    ::    DeviceGrid3D{T1, T2},
    mp      ::DeviceParticle3D{T1, T2},
    ΔT      ::T2,
    FLIP    ::T2,
    PIC     ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ξsx = ξsy = ξsz = vsx = vsy = vsz = T2(0.0)
        ξwx = ξwy = ξwz = vwx = vwy = vwz = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ξsx += Ni *  grid.vsT[p2n, 1]
                ξsy += Ni *  grid.vsT[p2n, 2]
                ξsz += Ni *  grid.vsT[p2n, 3]
                ξwx += Ni *  grid.vwT[p2n, 1]
                ξwy += Ni *  grid.vwT[p2n, 2]
                ξwz += Ni *  grid.vwT[p2n, 3]
                vsx += Ni * (grid.vsT[p2n, 1] - grid.vs[p2n, 1])
                vsy += Ni * (grid.vsT[p2n, 2] - grid.vs[p2n, 2])
                vsz += Ni * (grid.vsT[p2n, 3] - grid.vs[p2n, 3])
                vwx += Ni * (grid.vwT[p2n, 1] - grid.vw[p2n, 1])
                vwy += Ni * (grid.vwT[p2n, 2] - grid.vw[p2n, 2])
                vwz += Ni * (grid.vwT[p2n, 3] - grid.vw[p2n, 3])
            end
        end
        # update particle ξition
        mp.ξ[ix, 1] += ΔT * ξsx
        mp.ξ[ix, 2] += ΔT * ξsy
        mp.ξ[ix, 3] += ΔT * ξsz
        # update particle velocity
        mp.vs[ix, 1] = FLIP * (mp.vs[ix, 1] + vsx) + PIC * ξsx
        mp.vs[ix, 2] = FLIP * (mp.vs[ix, 2] + vsy) + PIC * ξsy
        mp.vs[ix, 3] = FLIP * (mp.vs[ix, 3] + vsz) + PIC * ξsz
        mp.vw[ix, 1] = FLIP * (mp.vw[ix, 1] + vwx) + PIC * ξwx
        mp.vw[ix, 2] = FLIP * (mp.vw[ix, 2] + vwy) + PIC * ξwy
        mp.vw[ix, 3] = FLIP * (mp.vw[ix, 3] + vwz) + PIC * ξwz
        # update particle momentum
        mp.ps[ix, 1] = mp.ms[ix] * mp.vs[ix, 1] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 2] = mp.ms[ix] * mp.vs[ix, 2] * (T2(1.0) - mp.n[ix])
        mp.ps[ix, 3] = mp.ms[ix] * mp.vs[ix, 3] * (T2(1.0) - mp.n[ix])
        mp.pw[ix, 1] = mp.mw[ix] * mp.vw[ix, 1] *            mp.n[ix]
        mp.pw[ix, 2] = mp.mw[ix] * mp.vw[ix, 2] *            mp.n[ix]
        mp.pw[ix, 3] = mp.mw[ix] * mp.vw[ix, 3] *            mp.n[ix]
    end
end

"""
    doublemapping2_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2})

Description:
---
Scatter momentum from particles to grid.
"""
@kernel inbounds=true function doublemapping2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update grid momentum
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
                @KAatomic grid.pw[p2n, 1] += mp.pw[ix, 1] * Ni
                @KAatomic grid.pw[p2n, 2] += mp.pw[ix, 2] * Ni
            end
        end
    end
end

"""
    doublemapping2_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2})

Description:
---
Scatter momentum from particles to grid.
"""
@kernel inbounds=true function doublemapping2_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        # update grid momentum
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            Ni = mp.Nij[ix, iy]
            if Ni ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                @KAatomic grid.ps[p2n, 1] += mp.ps[ix, 1] * Ni
                @KAatomic grid.ps[p2n, 2] += mp.ps[ix, 2] * Ni
                @KAatomic grid.ps[p2n, 3] += mp.ps[ix, 3] * Ni
                @KAatomic grid.pw[p2n, 1] += mp.pw[ix, 1] * Ni
                @KAatomic grid.pw[p2n, 2] += mp.pw[ix, 2] * Ni
                @KAatomic grid.pw[p2n, 3] += mp.pw[ix, 3] * Ni
            end
        end
    end
end

"""
    doublemapping3_TS!(grid::DeviceGrid2D{T1, T2}, bc::DeviceVBoundary2D{T1, T2}, ΔT::T2)

Description:
---
Solve equations on grid.
"""
@kernel inbounds=true function doublemapping3_TS!(
    grid::     DeviceGrid2D{T1, T2},
    bc  ::DeviceVBoundary2D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni 
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 2] = bc.vy_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
        grid.Δuw[ix, 1] = grid.vw[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vw[ix, 2] * ΔT
    end
end

"""
    doublemapping3_TS!(grid::DeviceGrid3D{T1, T2}, bc::DeviceVBoundary3D{T1, T2}, ΔT::T2)

Description:
---
Solve equations on grid.
"""
@kernel inbounds=true function doublemapping3_TS!(
    grid::     DeviceGrid3D{T1, T2},
    bc  ::DeviceVBoundary3D{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ grid.ni
        ms_denom = grid.ms[ix] < eps(T2) ? T2(0.0) : inv(grid.ms[ix])
        mi_denom = grid.mi[ix] < eps(T2) ? T2(0.0) : inv(grid.mi[ix])
        mw_denom = grid.mw[ix] < eps(T2) ? T2(0.0) : inv(grid.mw[ix])
        # compute nodal velocities
        grid.vs[ix, 1] = grid.ps[ix, 1] * ms_denom
        grid.vs[ix, 2] = grid.ps[ix, 2] * ms_denom
        grid.vs[ix, 3] = grid.ps[ix, 3] * ms_denom
        grid.vw[ix, 1] = grid.pw[ix, 1] * mw_denom
        grid.vw[ix, 2] = grid.pw[ix, 2] * mw_denom
        grid.vw[ix, 3] = grid.pw[ix, 3] * mw_denom
        # fixed Dirichlet nodes
        bc.vx_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 1] = bc.vx_s_val[ix] : nothing
        bc.vy_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 2] = bc.vy_s_val[ix] : nothing
        bc.vz_s_idx[ix] ≠ T1(0) ? grid.vs[ix, 3] = bc.vz_s_val[ix] : nothing
        bc.vx_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 1] = bc.vx_w_val[ix] : nothing
        bc.vy_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 2] = bc.vy_w_val[ix] : nothing
        bc.vz_w_idx[ix] ≠ T1(0) ? grid.vw[ix, 3] = bc.vz_w_val[ix] : nothing
        # compute nodal displacement
        grid.Δus[ix, 1] = grid.vs[ix, 1] * ΔT
        grid.Δus[ix, 2] = grid.vs[ix, 2] * ΔT
        grid.Δus[ix, 3] = grid.vs[ix, 3] * ΔT
        grid.Δuw[ix, 1] = grid.vw[ix, 1] * ΔT
        grid.Δuw[ix, 2] = grid.vw[ix, 2] * ΔT
        grid.Δuw[ix, 3] = grid.vw[ix, 3] * ΔT
    end
end

"""
    G2P_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, ΔT::T2)

Description:
---
Update particle information.
"""
@kernel inbounds=true function G2P_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = T2(1.0) / ΔT
        dfs1 = dfs2 = dfs3 = dfs4 = dfw1 = dfw2 = dfw3 = dfw4 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]
                ∂Ny = mp.∂Ny[ix, iy]
                # compute solid incremental deformation gradient
                dfs1 += grid.Δus[p2n, 1] * ∂Nx
                dfs2 += grid.Δus[p2n, 1] * ∂Ny
                dfs3 += grid.Δus[p2n, 2] * ∂Nx
                dfs4 += grid.Δus[p2n, 2] * ∂Ny
                dfw1 += grid.Δuw[p2n, 1] * ∂Nx
                dfw2 += grid.Δuw[p2n, 1] * ∂Ny
                dfw3 += grid.Δuw[p2n, 2] * ∂Nx
                dfw4 += grid.Δuw[p2n, 2] * ∂Ny
            end
        end
        mp.ΔFs[ix, 1] = dfs1
        mp.ΔFs[ix, 2] = dfs2
        mp.ΔFs[ix, 3] = dfs3
        mp.ΔFs[ix, 4] = dfs4
        # strain rate (Second Invariant of Strain Rate Tensor)
        ϵvxx = dfs1 * ΔT_1
        ϵvyy = dfs4 * ΔT_1
        ϵvxy = T2(0.5) * (dfs2 + dfs3) * ΔT_1
        mp.ϵv[ix] = sqrt(ϵvxx * ϵvxx + ϵvyy * ϵvyy + T2(2.0) * ϵvxy * ϵvxy)
        # compute strain increment 
        mp.Δϵijs[ix, 1] = dfs1
        mp.Δϵijs[ix, 2] = dfs4
        mp.Δϵijs[ix, 4] = dfs2 + dfs3
        mp.Δϵijw[ix, 1] = dfw1
        mp.Δϵijw[ix, 2] = dfw4
        mp.Δϵijw[ix, 4] = dfw2 + dfw3
        # update strain tensor
        mp.ϵijs[ix, 1] += dfs1
        mp.ϵijs[ix, 2] += dfs4
        mp.ϵijs[ix, 4] += dfs2 + dfs3
        mp.ϵijw[ix, 1] += dfw1
        mp.ϵijw[ix, 2] += dfw4
        mp.ϵijw[ix, 4] += dfw2 + dfw3
        # deformation gradient matrix
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]; F4 = mp.F[ix, 4]      
        mp.F[ix, 1] = (dfs1 + T2(1.0)) * F1 + dfs2 * F3
        mp.F[ix, 2] = (dfs1 + T2(1.0)) * F2 + dfs2 * F4
        mp.F[ix, 3] = (dfs4 + T2(1.0)) * F3 + dfs3 * F1
        mp.F[ix, 4] = (dfs4 + T2(1.0)) * F4 + dfs3 * F2
        # update jacobian value and particle Ωume
        J         = mp.F[ix, 1] * mp.F[ix, 4] - mp.F[ix, 2] * mp.F[ix, 3]
        ΔJ        = J * mp.Ω0[ix] / mp.Ω[ix]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs_init[ix] / J
        mp.ρw[ix] = mp.ρw_init[ix] / J
        # update pore pressure and n
        mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
            (T2(1.0) - mp.n[ix]) * (dfs1 + dfs4) + 
                       mp.n[ix]  * (dfw1 + dfw4))
        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
    end
end

"""
    G2P_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2}, 
        attr::DeviceProperty{T1, T2}, ΔT::T2)

Description:
---
Update particle information.
"""
@kernel inbounds=true function G2P_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        ΔT_1 = T2(1.0) / ΔT
        dfs1 = dfs2 = dfs3 = dfs4 = dfs5 = dfs6 = dfs7 = dfs8 = dfs9 = T2(0.0)
        dfw1 = dfw2 = dfw3 = dfw4 = dfw5 = dfw6 = dfw7 = dfw8 = dfw9 = T2(0.0)
        @KAunroll for iy in Int32(1):Int32(mp.NIC)
            if mp.Nij[ix, iy] ≠ T2(0.0)
                p2n = mp.p2n[ix, iy]
                ∂Nx = mp.∂Nx[ix, iy]; ds1 = grid.Δus[p2n, 1]; dw1 = grid.Δuw[p2n, 1]
                ∂Ny = mp.∂Ny[ix, iy]; ds2 = grid.Δus[p2n, 2]; dw2 = grid.Δuw[p2n, 2]
                ∂Nz = mp.∂Nz[ix, iy]; ds3 = grid.Δus[p2n, 3]; dw3 = grid.Δuw[p2n, 3]
                # compute incremental deformation gradient
                dfs1 += ds1 * ∂Nx; dfs2 += ds1 * ∂Ny; dfs3 += ds1 * ∂Nz
                dfs4 += ds2 * ∂Nx; dfs5 += ds2 * ∂Ny; dfs6 += ds2 * ∂Nz
                dfs7 += ds3 * ∂Nx; dfs8 += ds3 * ∂Ny; dfs9 += ds3 * ∂Nz
                dfw1 += dw1 * ∂Nx; dfw2 += dw1 * ∂Ny; dfw3 += dw1 * ∂Nz
                dfw4 += dw2 * ∂Nx; dfw5 += dw2 * ∂Ny; dfw6 += dw2 * ∂Nz
                dfw7 += dw3 * ∂Nx; dfw8 += dw3 * ∂Ny; dfw9 += dw3 * ∂Nz
            end
        end
        mp.ΔFs[ix, 1] = dfs1; mp.Δfs[ix, 2] = dfs2; mp.Δfs[ix, 3] = dfs3
        mp.ΔFs[ix, 4] = dfs4; mp.Δfs[ix, 5] = dfs5; mp.Δfs[ix, 6] = dfs6
        mp.ΔFs[ix, 7] = dfs7; mp.Δfs[ix, 8] = dfs8; mp.Δfs[ix, 9] = dfs9
        # strain rate (Second Invariant of Strain Rate Tensor)
        ϵvxx = dfs1 * ΔT_1
        ϵvyy = dfs5 * ΔT_1
        ϵvzz = dfs9 * ΔT_1
        ϵvxy = T2(0.5) * (dfs2 + dfs4) * ΔT_1
        ϵvyz = T2(0.5) * (dfs6 + dfs8) * ΔT_1
        ϵvxz = T2(0.5) * (dfs3 + dfs7) * ΔT_1
        mp.ϵv[ix] = sqrt(ϵvxx * ϵvxx + ϵvyy * ϵvyy + ϵvzz * ϵvzz + 
            T2(2.0) * (ϵvxy * ϵvxy + ϵvyz * ϵvyz + ϵvxz * ϵvxz))
        # compute strain increment
        mp.Δϵijs[ix, 1] = dfs1
        mp.Δϵijs[ix, 2] = dfs5
        mp.Δϵijs[ix, 3] = dfs9
        mp.Δϵijs[ix, 4] = dfs2 + dfs4
        mp.Δϵijs[ix, 5] = dfs6 + dfs8
        mp.Δϵijs[ix, 6] = dfs3 + dfs7
        mp.Δϵijw[ix, 1] = dfw1
        mp.Δϵijw[ix, 2] = dfw5
        mp.Δϵijw[ix, 3] = dfw9 
        mp.Δϵijw[ix, 4] = dfw2 + dfw4
        mp.Δϵijw[ix, 5] = dfw6 + dfw8
        mp.Δϵijw[ix, 6] = dfw3 + dfw7
        # update strain tensor
        mp.ϵijs[ix, 1] += dfs1
        mp.ϵijs[ix, 2] += dfs5
        mp.ϵijs[ix, 3] += dfs9
        mp.ϵijs[ix, 4] += dfs2 + dfs4
        mp.ϵijs[ix, 5] += dfs6 + dfs8
        mp.ϵijs[ix, 6] += dfs3 + dfs7
        mp.ϵijw[ix, 1] += dfs1
        mp.ϵijw[ix, 2] += dfs5
        mp.ϵijw[ix, 3] += dfs9
        mp.ϵijw[ix, 4] += dfs2 + dfs4
        mp.ϵijw[ix, 5] += dfs6 + dfs8
        mp.ϵijw[ix, 6] += dfs3 + dfs7
        # deformation gradient matrix
        F1 = mp.F[ix, 1]; F2 = mp.F[ix, 2]; F3 = mp.F[ix, 3]
        F4 = mp.F[ix, 4]; F5 = mp.F[ix, 5]; F6 = mp.F[ix, 6]
        F7 = mp.F[ix, 7]; F8 = mp.F[ix, 8]; F9 = mp.F[ix, 9]        
        mp.F[ix, 1] = (dfs1 + T2(1.0)) * F1 + dfs2 * F4 + dfs3 * F7
        mp.F[ix, 2] = (dfs1 + T2(1.0)) * F2 + dfs2 * F5 + dfs3 * F8
        mp.F[ix, 3] = (dfs1 + T2(1.0)) * F3 + dfs2 * F6 + dfs3 * F9
        mp.F[ix, 4] = (dfs5 + T2(1.0)) * F4 + dfs4 * F1 + dfs6 * F7
        mp.F[ix, 5] = (dfs5 + T2(1.0)) * F5 + dfs4 * F2 + dfs6 * F8
        mp.F[ix, 6] = (dfs5 + T2(1.0)) * F6 + dfs4 * F3 + dfs6 * F9
        mp.F[ix, 7] = (dfs9 + T2(1.0)) * F7 + dfs8 * F4 + dfs7 * F1
        mp.F[ix, 8] = (dfs9 + T2(1.0)) * F8 + dfs8 * F5 + dfs7 * F2
        mp.F[ix, 9] = (dfs9 + T2(1.0)) * F9 + dfs8 * F6 + dfs7 * F3
        # update jacobian value and particle Ωume
        J = mp.F[ix, 1] * mp.F[ix, 5] * mp.F[ix, 9] + 
            mp.F[ix, 2] * mp.F[ix, 6] * mp.F[ix, 7] +
            mp.F[ix, 3] * mp.F[ix, 4] * mp.F[ix, 8] - 
            mp.F[ix, 7] * mp.F[ix, 5] * mp.F[ix, 3] -
            mp.F[ix, 8] * mp.F[ix, 6] * mp.F[ix, 1] - 
            mp.F[ix, 9] * mp.F[ix, 4] * mp.F[ix, 2]
        ΔJ        = J * mp.Ω0[ix] / mp.Ω[ix]
        mp.Ω[ix]  = J * mp.Ω0[ix]
        mp.ρs[ix] = mp.ρs_init[ix] / J
        mp.ρw[ix] = mp.ρw_init[ix] / J
        # update pore pressure and n
        mp.σw[ix] += (attr.Kw[attr.nid[ix]] / mp.n[ix]) * (
            (T2(1.0) - mp.n[ix]) * (dfs1 + dfs5 + dfs9) + 
                       mp.n[ix]  * (dfw1 + dfw5 + dfw9))
        mp.n[ix] = clamp(T2(1.0) - (T2(1.0) - mp.n[ix]) / ΔJ, T2(0.0), T2(1.0))
    end
end

"""
    vollock1_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2})

Description:
---
Mapping mean stress and Ωume from particle to grid.
"""
@kernel inbounds=true function vollock1_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = mp.p2c[ix]
        vol = mp.Ω[ix]
        @KAatomic grid.σm[p2c] += vol * mp.σm[ix]
        @KAatomic grid.σw[p2c] += vol * mp.σw[ix]
        @KAatomic grid.Ω[p2c]  += vol
    end
end

"""
    vollock1_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2})

Description:
---
Mapping mean stress and Ωume from particle to grid.
"""
@kernel inbounds=true function vollock1_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = mp.p2c[ix]
        vol = mp.Ω[ix]
        @KAatomic grid.σm[p2c] += vol * mp.σm[ix]
        @KAatomic grid.σw[p2c] += vol * mp.σw[ix]
        @KAatomic grid.Ω[p2c]  += vol
    end
end

"""
    vollock2_TS!(grid::DeviceGrid2D{T1, T2}, mp::DeviceParticle2D{T1, T2})

Description:
---
Mapping back mean stress and Ωume from grid to particle.
"""
@kernel inbounds=true function vollock2_TS!(
    grid::    DeviceGrid2D{T1, T2},
    mp  ::DeviceParticle2D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = mp.p2c[ix]
        σm  = grid.σm[p2c] / grid.Ω[p2c]
        σw  = grid.σw[p2c] / grid.Ω[p2c]
        mp.σij[ix, 1] = mp.sij[ix, 1] + σm
        mp.σij[ix, 2] = mp.sij[ix, 2] + σm
        mp.σij[ix, 3] = mp.sij[ix, 3] + σm
        mp.σij[ix, 4] = mp.sij[ix, 4]
        # update mean stress tensor
        σm = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.σm[ix] = σm
        mp.σw[ix] = σw
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - σm
        mp.sij[ix, 2] = mp.σij[ix, 2] - σm
        mp.sij[ix, 3] = mp.σij[ix, 3] - σm
        mp.sij[ix, 4] = mp.σij[ix, 4]
    end
end

"""
    vollock2_TS!(grid::DeviceGrid3D{T1, T2}, mp::DeviceParticle3D{T1, T2})

Description:
---
Mapping back mean stress and Ωume from grid to particle.
"""
@kernel inbounds=true function vollock2_TS!(
    grid::    DeviceGrid3D{T1, T2},
    mp  ::DeviceParticle3D{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        p2c = mp.p2c[ix]
        σm  = grid.σm[p2c] / grid.Ω[p2c]
        σw  = grid.σw[p2c] / grid.Ω[p2c]
        mp.σij[ix, 1] = mp.sij[ix, 1] + σm
        mp.σij[ix, 2] = mp.sij[ix, 2] + σm
        mp.σij[ix, 3] = mp.sij[ix, 3] + σm
        mp.σij[ix, 4] = mp.sij[ix, 4]
        mp.σij[ix, 5] = mp.sij[ix, 5]
        mp.σij[ix, 6] = mp.sij[ix, 6]
        # update mean stress tensor
        σm = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.σm[ix] = σm
        mp.σw[ix] = σw
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - σm
        mp.sij[ix, 2] = mp.σij[ix, 2] - σm
        mp.sij[ix, 3] = mp.σij[ix, 3] - σm
        mp.sij[ix, 4] = mp.σij[ix, 4]
        mp.sij[ix, 5] = mp.σij[ix, 5]
        mp.sij[ix, 6] = mp.σij[ix, 6]
    end
end