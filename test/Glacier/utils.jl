using NearestNeighbors

function loadmodel(h)
    #h = 0.8
    surf = meshbuilder(0:h:200, 0:h:200)
    surf = hcat(surf, zeros(size(surf, 1)))
    surf[findall(i->surf[i, 1]<50, axes(surf, 1)), 3] .=50
    idx = findall(i -> 100>surf[i, 1]≥50, axes(surf, 1))
    surf[idx, 3] .= 100 .- surf[idx, 1]


    idx = findall(i->(surf[i, 1]-62.5)^2 / 37.5^2 + (surf[i, 2]-100)^2 / 20^2 ≤ 1, axes(surf, 1))
    poly = get_polygon(surf[idx, 1:2], ratio=0.1)

    bottom = SLBL3D(surf, poly, 40, 75)

    up = surf[pip_query(poly, surf[:, 1:2]), :]
    down = bottom[pip_query(poly, bottom[:, 1:2]), :]

    pts = dem2particle(up, h, down)
    return bottom, pts
end

@views function get_ext(pts::Matrix, grid::NamedTuple, thickness::Int)
    id = zeros(Int, grid.ni)
    @inbounds for i in axes(pts, 1)
        mξx, mξy, mξz = pts[i, 1], pts[i, 2], pts[i, 3]
        bnx = unsafe_trunc(Int, fld(mξx - grid.x1, grid.h))
        bny = unsafe_trunc(Int, fld(mξy - grid.y1, grid.h))
        bnz = unsafe_trunc(Int, fld(mξz - grid.z1, grid.h))
        bi1 = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
        (bi1 > grid.ni || bi1 < 1) && error("Invalid particle id $(bi1), exceed the grid range")
        bi2 = bi1 + 1      ; id[bi1] = 1; id[bi2] = 1
        bi3 = bi1 + grid.ny; id[bi3] = 1
        bi4 = bi3 + 1      ; id[bi4] = 1
    end
    surf_id   = findall(id .== 1)
    grid_base = grid.ξ[surf_id, :]
    surf_norm = get_normals(pts, k=8)
    kdtree    = KDTree(pts')
    grid_nvec = zeros(size(grid_base, 1), 3)
    idx, _    = nn(kdtree, grid_base')
    grid_nvec = surf_norm[idx, :]
    @inbounds for (i, sid) in enumerate(surf_id)
        for j in 0:1:thickness
            gid = sid - j * grid.ny * grid.nx
            id[gid] = i
        end
    end
    return grid_nvec, id
end

@inline function _get_grid_coord(grid::DeviceGrid{T1, T2}, mpx::T2, mpy::T2, mpz::T2) where {T1, T2}
    bnx = unsafe_trunc(T1, floor((mpx - grid.x1) * grid.invh))
    bny = unsafe_trunc(T1, floor((mpy - grid.y1) * grid.invh))
    bnz = unsafe_trunc(T1, floor((mpz - grid.z1) * grid.invh))
    bid = grid.nx * grid.ny * bnz + grid.ny * bnx + bny + 1
    gξx = grid.x1 + bnx * grid.h
    gξy = grid.y1 + bny * grid.h
    gξz = grid.z1 + bnz * grid.h
    return bid, gξx, gξy, gξz
end

@inline function _get_grid_coord(grid::DeviceGrid{T1, T2}, bid) where {T1, T2}
    P   = grid.nx * grid.ny
    k   = bid - 1
    zn  = k ÷ P
    r   = k % P
    xn  = r ÷ grid.ny
    yn  = r % grid.ny
    x0  = grid.x1 + xn * grid.h
    y0  = grid.y1 + yn * grid.h
    z0  = grid.z1 + zn * grid.h
    return x0, y0, z0
end