export get_grid_ξ

"""
    get_grid_ξ(grid::DeviceGrid{T1, T2}, bid) where {T1, T2}

Get the coordinates of the grid node with index `bid` in the grid `grid`.
The grid is defined by its dimensions and the spacing between nodes.
The function returns the coordinates `(x0, y0, z0)` of the grid node.
"""
@inline function get_grid_ξ(grid::DeviceGrid{T1, T2}, bid) where {T1, T2}
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