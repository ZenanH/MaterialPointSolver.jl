struct NV3D{
    T1 <: AbstractArray, # Array{T2, 2}
    T2 <: AbstractArray, # Array{T1, 1}
} <: UserGridExtra
    n ::T1
    id::T2
end

@user_struct NV3D

function get_ext(
    base_data::Matrix, 
    grid     ::DeviceGrid3D{T1, T2}, 
    thickness, 
    offset
) where {T1, T2}
    surf   = base_data .- offset'
    nvec   = getnormals(surf)
    id     = zeros(T1, grid.ni)
    θ      = T2(thickness)
    kdtree = KDTree(surf')

    @inbounds for i in axes(grid.ξ, 1)
        x, y, z = grid.ξ[i, 1], grid.ξ[i, 2], grid.ξ[i, 3]
        point = @view grid.ξ[i, :]
        index, _ = knn(kdtree, point, 1)
        surf_id = index[1]
        nx, ny, nz = nvec[surf_id, 1], nvec[surf_id, 2], nvec[surf_id, 3]
        sx, sy, sz = surf[surf_id, 1], surf[surf_id, 2], surf[surf_id, 3]
        vx, vy, vz = x - sx, y - sy, z - sz
        dot_val = vx * nx + vy * ny + vz * nz
        if dot_val < 0
            sz - θ ≤ z ≤ sz ? id[i] = surf_id : nothing
        end
    end

    tmp = NV3D{AbstractArray, AbstractArray}(nvec, id)
    ext = user_adapt(Array, tmp)
    return ext
end

@views function getoffset(grid::DeviceGrid3D, pts::Matrix)
    pts_min = minimum(pts, dims=1)
    offset  = zeros(3)

    vid = findall(pts_min[1] .- grid.ξ[:, 1] .≥ 0)
    _, id1 = findmin(pts_min[1] .- grid.ξ[vid, 1])
    vid = findall(pts_min[2] .- grid.ξ[:, 2] .≥ 0)
    _, id2 = findmin(pts_min[2] .- grid.ξ[vid, 2])
    vid = findall(pts_min[3] .- grid.ξ[:, 3] .≥ 0)
    _, id3 = findmin(pts_min[3] .- grid.ξ[vid, 3])

    offset[1] = pts_min[1] - grid.ξ[id1, 1] + grid.dx * 0.25
    offset[2] = pts_min[2] - grid.ξ[id2, 2] + grid.dy * 0.25
    offset[3] = pts_min[3] - grid.ξ[id3, 3] + grid.dz * 0.25
    return offset
end

@views function Tanimation(args::DeviceArgs3D{T1, T2}, offset) where {T1, T2}
    proj_path = joinpath(args.project_path, args.project_name)
    h5filedir = joinpath(proj_path, "$(args.project_name).h5")
    anim_path = mkpath(joinpath(proj_path, "animation"))
    nds_path  = joinpath(proj_path, "grid")
    mps_path  = joinpath(proj_path, args.project_name)
    fid       = h5open(h5filedir, "r")
    itr       = (read(fid["FILE_NUM"])-1) |> Int64
    nid       = fid["nid"        ] |> read
    ξ0        = fid["mp_coords0" ] |> read
    ξ0      .+= offset'
    grid_ξ    = fid["grid_coords"] |> read
    ni        = size(grid_ξ, 1)
    p         = Progress(length(1:1:itr) - 1; 
        desc      = "\e[1;36m[ Info:\e[0m $(lpad("ani_vtp", 7))",
        color     = :white,
        barlen    = 12,
        barglyphs = BarGlyphs(" ■■  ")
    )
    vtkfs = Vector{WriteVTK.DatasetFile}(undef, itr)
    times = Vector{Float64}(undef, itr)
    # generate files for particles
    @inbounds Threads.@threads for i in 1:itr
        # read data from HDF5 file
        time = fid["group$(i)/time"      ] |> read
        σij  = fid["group$(i)/stress"    ] |> read
        ϵijs = fid["group$(i)/strain_s"  ] |> read
        ϵq   = fid["group$(i)/eqstrain"  ] |> read
        ϵk   = fid["group$(i)/ekstrain"  ] |> read
        ϵv   = fid["group$(i)/eqrate"    ] |> read
        vs   = fid["group$(i)/velocity_s"] |> read
        ms   = fid["group$(i)/mass_s"    ] |> read
        Ω    = fid["group$(i)/volume"    ] |> read
        ξ    = fid["group$(i)/coords"    ] |> read
        ξ  .+= offset'
        args.coupling==:TS ? (
            σw   = fid["group$(i)/pressure_w"] |> read;
            ϵijw = fid["group$(i)/strain_w"  ] |> read;
            vw   = fid["group$(i)/velocity_w"] |> read;
            n    = fid["group$(i)/porosity"  ] |> read;
        ) : nothing            
        # write data
        vtp_cls = [MeshCell(PolyData.Verts(), [j]) for j in 1:size(ξ0, 1)]
        let vtk = vtk_grid(joinpath(anim_path, "iter_$(i)"), ξ', vtp_cls, ascii=false)
            vtk["stress"      ] = σij'
            vtk["strain_s"    ] = ϵijs'
            vtk["eqstrain"    ] = ϵq
            vtk["ekstrain"    ] = ϵk
            vtk["eqrate"      ] = ϵv
            vtk["mass_s"      ] = ms
            vtk["volume"      ] = Ω
            vtk["velocity_s"  ] = vs'
            vtk["nid"         ] = nid
            vtk["displacement"] = abs.(ξ .- ξ0)'
            args.coupling==:TS ? (
                vtk["strain_w"  ] = ϵijw';
                vtk["pressure_w"] = σw   ;
                vtk["velocity_w"] = vw'  ;
                vtk["porosity"  ] = n    ;
            ) : nothing
            times[i] = time
            vtkfs[i] = vtk
        end
        next!(p)
    end
    idx = sortperm(times)
    times_s = times[idx]
    vtkfs_s = vtkfs[idx]
    paraview_collection(mps_path) do pvd
        @inbounds for i in 1:itr
            pvd[times_s[i]] = vtkfs_s[i]
        end
    end
    # generate vtp files for nodes
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:ni]
    vtk_grid(nds_path, grid_ξ', vtp_cls, ascii=false) do vtk end
    close(fid)
end