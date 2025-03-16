function gettheme()
    custom_dark = theme_dark()
    custom_dark.textcolor = :white
    custom_dark.linecolor = :white
    custom_dark.backgroundcolor = "#181818"
    return custom_dark
end

function getbox(grid::DeviceGrid2D{T1, T2}) where {T1, T2}
    vmin = minimum(grid.ξ, dims=1)
    vmax = maximum(grid.ξ, dims=1)
    
    # 直接创建最终的8个Point2f元组数组（4条边 × 2个点）
    return Point2f[
        (vmin[1], vmin[2]), (vmax[1], vmin[2]),  # 底边
        (vmax[1], vmin[2]), (vmax[1], vmax[2]),  # 右边
        (vmax[1], vmax[2]), (vmin[1], vmax[2]),  # 顶边
        (vmin[1], vmax[2]), (vmin[1], vmin[2])   # 左边
    ]
end

function getbox(grid::DeviceGrid3D{T1, T2}) where {T1, T2}
    vmin = minimum(grid.ξ, dims=1)
    vmax = maximum(grid.ξ, dims=1)
    
    # 直接创建24个Point3f元组数组（12条边 × 2个点）
    return Point3f[
        # 底面 (4条边)
        (vmin[1], vmin[2], vmin[3]), (vmax[1], vmin[2], vmin[3]),
        (vmax[1], vmin[2], vmin[3]), (vmax[1], vmax[2], vmin[3]),
        (vmax[1], vmax[2], vmin[3]), (vmin[1], vmax[2], vmin[3]),
        (vmin[1], vmax[2], vmin[3]), (vmin[1], vmin[2], vmin[3]),
        
        # 顶面 (4条边)
        (vmin[1], vmin[2], vmax[3]), (vmax[1], vmin[2], vmax[3]),
        (vmax[1], vmin[2], vmax[3]), (vmax[1], vmax[2], vmax[3]),
        (vmax[1], vmax[2], vmax[3]), (vmin[1], vmax[2], vmax[3]),
        (vmin[1], vmax[2], vmax[3]), (vmin[1], vmin[2], vmax[3]),
        
        # 侧边 (4条边)
        (vmin[1], vmin[2], vmin[3]), (vmin[1], vmin[2], vmax[3]),
        (vmax[1], vmin[2], vmin[3]), (vmax[1], vmin[2], vmax[3]),
        (vmax[1], vmax[2], vmin[3]), (vmax[1], vmax[2], vmax[3]),
        (vmin[1], vmax[2], vmin[3]), (vmin[1], vmax[2], vmax[3])
    ]
end

