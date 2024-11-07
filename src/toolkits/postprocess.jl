#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : postprocess.jl                                                             |
|  Description: Post-process functions                                                     |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 1. fastvtp()                                                               |
|               2. savevtp()   [2D & 3D]                                                   |
|               3. animation() [2D & 3D]                                                   |
+==========================================================================================#

export fastvtp, savevtp, animation

"""
    fastvtp(coords; vtppath="output", data::T=NamedTuple())

Description:
---
Generates a `.vtp` file by passing custom fields.
"""
function fastvtp(coords; vtppath="output", data::T=NamedTuple()) where T <: NamedTuple
    pts_num = size(coords, 1)
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:pts_num]
    vtk_grid(vtppath, coords', vtp_cls, ascii=false) do vtk
        keys(data) ≠ () && for vtp_key in keys(data)
            vtk[string(vtp_key)] = getfield(data, vtp_key)
        end
    end
    return nothing
end

"""
    savevtp(args::DeviceArgs2D{T1, T2}, grid::DeviceGrid2D{T1, T2}, 
        mp::DeviceParticle2D{T1, T2}, attr::DeviceProperty{T1, T2})

Description:
---
Generates the final geometry and properties in `.vtp` format (2D).
"""
@views function savevtp(
    args::    DeviceArgs2D{T1, T2}, 
    grid::    DeviceGrid2D{T1, T2}, 
    mp  ::DeviceParticle2D{T1, T2}, 
    attr::  DeviceProperty{T1, T2}
) where {T1, T2}
    prj_path = joinpath(args.project_path, args.project_name)
    mps_path = joinpath(prj_path, args.project_name)
    nds_path = joinpath(prj_path, "grid")
    # generate vtp files for particles
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:mp.np]
    vtk_grid(mps_path, mp.ξ', vtp_cls, ascii=false) do vtk
        vtk["stress"      ] = mp.σij[:, [1, 2, 4]]'
        vtk["strain_s"    ] = mp.ϵijs[:, [1, 2, 4]]'
        vtk["eqstrain"    ] = mp.ϵq
        vtk["eqrate"      ] = mp.ϵv
        vtk["ekstrain"    ] = mp.ϵk
        vtk["mass_s"      ] = mp.ms
        vtk["nid"         ] = attr.nid
        vtk["volume"      ] = mp.Ω
        vtk["velocity_s"  ] = mp.vs'
        vtk["displacement"] = abs.(mp.ξ .- mp.ξ0)'
        args.coupling==:TS ? (
            vtk["strain_w"  ] = mp.ϵijw[:, [1, 2, 4]]';
            vtk["pressure_w"] = mp.σw                 ;
            vtk["velocity_w"] = mp.vw'                ;
            vtk["porosity"  ] = mp.n                  ;
        ) : nothing
    end
    # generate vtp files for nodes
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:grid.ni]
    vtk_grid(nds_path, grid.ξ', vtp_cls, ascii=false) do vtk end
    @info "final vtp file is saved in project path"
    return nothing
end

"""
    savevtp(args::DeviceArgs3D{T1, T2}, grid::DeviceGrid3D{T1, T2}, 
        mp::DeviceParticle3D{T1, T2}, attr::DeviceProperty{T1, T2})

Description:
---
Generates the final geometry and properties in `.vtp` format (3D).
"""
@views function savevtp(
    args::    DeviceArgs3D{T1, T2}, 
    grid::    DeviceGrid3D{T1, T2}, 
    mp  ::DeviceParticle3D{T1, T2}, 
    attr::  DeviceProperty{T1, T2}
) where {T1, T2}
    prj_path = joinpath(args.project_path, args.project_name)
    mps_path = joinpath(prj_path, args.project_name)
    nds_path = joinpath(prj_path, "grid")
    # generate vtp files for particles
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:mp.np]
    vtk_grid(mps_path, mp.ξ', vtp_cls, ascii=false) do vtk
        vtk["stress"      ] = mp.σij'
        vtk["strain_s"    ] = mp.ϵijs'
        vtk["eqstrain"    ] = mp.ϵq
        vtk["eqrate"      ] = mp.ϵv
        vtk["ekstrain"    ] = mp.ϵk
        vtk["mass_s"      ] = mp.ms
        vtk["nid"         ] = attr.nid
        vtk["volume"      ] = mp.Ω
        vtk["velocity_s"  ] = mp.vs'
        vtk["displacement"] = abs.(mp.ξ .- mp.ξ0)'
        args.coupling==:TS ? (
            vtk["strain_w"  ] = mp.ϵijw';
            vtk["pressure_w"] = mp.σw   ;
            vtk["velocity_w"] = mp.vw'  ;
            vtk["porosity"  ] = mp.n    ;
        ) : nothing
    end
    # generate vtp files for nodes
    vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:grid.ni]
    vtk_grid(nds_path, grid.ξ', vtp_cls, ascii=false) do vtk end
    @info "final vtp file is saved in project path"
    return nothing
end

"""
    animation(args::DeviceArgs2D{T1, T2})

Description:
---
Generates animation by using the data from HDF5 file (2D).
"""
@views function animation(args::DeviceArgs2D{T1, T2}) where {T1, T2}
    proj_path = joinpath(args.project_path, args.project_name)
    h5filedir = joinpath(proj_path, "$(args.project_name).h5")
    anim_path = mkpath(joinpath(proj_path, "animation"))
    nds_path  = joinpath(proj_path, "grid")
    mps_path  = joinpath(proj_path, args.project_name)
    fid       = h5open(h5filedir, "r")
    itr       = (read(fid["FILE_NUM"])-1) |> Int64
    nid       = fid["nid"        ] |> read
    ξ0        = fid["mp_coords0" ] |> read
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
        args.coupling==:TS ? (
            σw   = fid["group$(i)/pressure_w"] |> read;
            ϵijw = fid["group$(i)/strain_w"  ] |> read;
            vw   = fid["group$(i)/velocity_w"] |> read;
            n    = fid["group$(i)/porosity"  ] |> read;
        ) : nothing
        # write data
        vtp_cls = [MeshCell(PolyData.Verts(), [j]) for j in 1:size(ξ0, 1)]
        let vtk = vtk_grid(joinpath(anim_path, "iter_$(i)"), ξ', vtp_cls, ascii=false)
            vtk["stress"      ] = σij[:, [1, 2, 4]]'
            vtk["strain_s"    ] = ϵijs[:, [1, 2, 4]]'
            vtk["eqstrain"    ] = ϵq
            vtk["ekstrain"    ] = ϵk
            vtk["eqrate"      ] = ϵv
            vtk["mass_s"      ] = ms
            vtk["volume"      ] = Ω
            vtk["velocity_s"  ] = vs'
            vtk["nid"         ] = nid
            vtk["displacement"] = abs.(ξ .- ξ0)'
            args.coupling==:TS ? (
                vtk["strain_w"  ] = ϵijw[:, [1, 2, 4]]';
                vtk["pressure_w"] = σw                 ;
                vtk["velocity_w"] = vw'                ;
                vtk["porosity"  ] = n                  ;
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
# @views function animation(args::DeviceArgs2D{T1, T2}) where {T1, T2}
#     proj_path = joinpath(args.project_path, args.project_name)
#     h5filedir = joinpath(proj_path, "$(args.project_name).h5")
#     anim_path = mkpath(joinpath(proj_path, "animation"))
#     nds_path  = joinpath(proj_path, "grid")
#     mps_path  = joinpath(proj_path, args.project_name)
#     fid       = h5open(h5filedir, "r")
#     itr       = (read(fid["FILE_NUM"])-1) |> Int64
#     nid       = fid["nid"        ] |> read
#     ξ0        = fid["mp_coords0" ] |> read
#     grid_ξ    = fid["grid_coords"] |> read
#     ni        = size(grid_ξ, 1)
#     p         = Progress(length(1:1:itr) - 1; 
#         desc      = "\e[1;36m[ Info:\e[0m $(lpad("ani_vtp", 7))",
#         color     = :white,
#         barlen    = 12,
#         barglyphs = BarGlyphs(" ■■  ")
#     )
#     # generate files for particles
#     paraview_collection(mps_path) do pvd
#         @inbounds for i in 1:itr
#             # read data from HDF5 file
#             time = fid["group$(i)/time"      ] |> read
#             σij  = fid["group$(i)/stress"    ] |> read
#             ϵijs = fid["group$(i)/strain_s"  ] |> read
#             ϵq   = fid["group$(i)/eqstrain"  ] |> read
#             ϵk   = fid["group$(i)/ekstrain"  ] |> read
#             ϵv   = fid["group$(i)/eqrate"    ] |> read
#             vs   = fid["group$(i)/velocity_s"] |> read
#             ms   = fid["group$(i)/mass_s"    ] |> read
#             Ω    = fid["group$(i)/volume"    ] |> read
#             ξ    = fid["group$(i)/coords"    ] |> read
#             args.coupling==:TS ? (
#                 σw   = fid["group$(i)/pressure_w"] |> read;
#                 ϵijw = fid["group$(i)/strain_w"  ] |> read;
#                 vw   = fid["group$(i)/velocity_w"] |> read;
#                 n    = fid["group$(i)/porosity"  ] |> read;
#             ) : nothing
#             # write data
#             vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:size(ξ0, 1)]
#             let vtk = vtk_grid(joinpath(anim_path, "iter_$(i)"), ξ', vtp_cls, ascii=false)
#                 vtk["stress"      ] = σij[:, [1, 2, 4]]'
#                 vtk["strain_s"    ] = ϵijs[:, [1, 2, 4]]'
#                 vtk["eqstrain"    ] = ϵq
#                 vtk["ekstrain"    ] = ϵk
#                 vtk["eqrate"      ] = ϵv
#                 vtk["mass_s"      ] = ms
#                 vtk["volume"      ] = Ω
#                 vtk["velocity_s"  ] = vs'
#                 vtk["nid"         ] = nid
#                 vtk["displacement"] = abs.(ξ .- ξ0)'
#                 args.coupling==:TS ? (
#                     vtk["strain_w"  ] = ϵijw[:, [1, 2, 4]]';
#                     vtk["pressure_w"] = σw                 ;
#                     vtk["velocity_w"] = vw'                ;
#                     vtk["porosity"  ] = n                  ;
#                 ) : nothing
#                 pvd[time] = vtk
#             end
#             next!(p)
#         end
#     end
#     # generate vtp files for nodes
#     vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:ni]
#     vtk_grid(nds_path, grid_ξ', vtp_cls, ascii=false) do vtk end
#     close(fid)
# end

"""
    animation(args::DeviceArgs3D{T1, T2})

Description:
---
Generates animation by using the data from HDF5 file (3D).
"""
@views function animation(args::DeviceArgs3D{T1, T2}) where {T1, T2}
    proj_path = joinpath(args.project_path, args.project_name)
    h5filedir = joinpath(proj_path, "$(args.project_name).h5")
    anim_path = mkpath(joinpath(proj_path, "animation"))
    nds_path  = joinpath(proj_path, "grid")
    mps_path  = joinpath(proj_path, args.project_name)
    fid       = h5open(h5filedir, "r")
    itr       = (read(fid["FILE_NUM"])-1) |> Int64
    nid       = fid["nid"        ] |> read
    ξ0        = fid["mp_coords0" ] |> read
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
# @views function animation(args::DeviceArgs3D{T1, T2}) where {T1, T2}
#     proj_path = joinpath(args.project_path, args.project_name)
#     h5filedir = joinpath(proj_path, "$(args.project_name).h5")
#     anim_path = mkpath(joinpath(proj_path, "animation"))
#     nds_path  = joinpath(proj_path, "grid")
#     mps_path  = joinpath(proj_path, args.project_name)
#     fid       = h5open(h5filedir, "r")
#     itr       = (read(fid["FILE_NUM"])-1) |> Int64
#     nid       = fid["nid"        ] |> read
#     ξ0        = fid["mp_coords0" ] |> read
#     grid_ξ    = fid["grid_coords"] |> read
#     ni        = size(grid_ξ, 1)
#     p         = Progress(length(1:1:itr) - 1; 
#         desc      = "\e[1;36m[ Info:\e[0m $(lpad("ani_vtp", 7))",
#         color     = :white,
#         barlen    = 12,
#         barglyphs = BarGlyphs(" ■■  ")
#     )
#     # generate files for particles
#     paraview_collection(mps_path) do pvd
#         @inbounds for i in 1:itr
#             # read data from HDF5 file
#             time = fid["group$(i)/time"      ] |> read
#             σij  = fid["group$(i)/stress"    ] |> read
#             ϵijs = fid["group$(i)/strain_s"  ] |> read
#             ϵq   = fid["group$(i)/eqstrain"  ] |> read
#             ϵk   = fid["group$(i)/ekstrain"  ] |> read
#             ϵv   = fid["group$(i)/eqrate"    ] |> read
#             vs   = fid["group$(i)/velocity_s"] |> read
#             ms   = fid["group$(i)/mass_s"    ] |> read
#             Ω    = fid["group$(i)/volume"    ] |> read
#             ξ    = fid["group$(i)/coords"    ] |> read
#             args.coupling==:TS ? (
#                 σw   = fid["group$(i)/pressure_w"] |> read;
#                 ϵijw = fid["group$(i)/strain_w"  ] |> read;
#                 vw   = fid["group$(i)/velocity_w"] |> read;
#                 n    = fid["group$(i)/porosity"  ] |> read;
#             ) : nothing            
#             # write data
#             vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:size(ξ0, 1)]
#             let vtk = vtk_grid(joinpath(anim_path, "iter_$(i)"), ξ', vtp_cls, ascii=false)
#                 vtk["stress"      ] = σij'
#                 vtk["strain_s"    ] = ϵijs'
#                 vtk["eqstrain"    ] = ϵq
#                 vtk["ekstrain"    ] = ϵk
#                 vtk["eqrate"      ] = ϵv
#                 vtk["mass_s"      ] = ms
#                 vtk["volume"      ] = Ω
#                 vtk["velocity_s"  ] = vs'
#                 vtk["nid"         ] = nid
#                 vtk["displacement"] = abs.(ξ .- ξ0)'
#                 args.coupling==:TS ? (
#                     vtk["strain_w"  ] = ϵijw';
#                     vtk["pressure_w"] = σw   ;
#                     vtk["velocity_w"] = vw'  ;
#                     vtk["porosity"  ] = n    ;
#                 ) : nothing
#                 pvd[time] = vtk
#             end
#             next!(p)
#         end
#     end
#     # generate vtp files for nodes
#     vtp_cls = [MeshCell(PolyData.Verts(), [i]) for i in 1:ni]
#     vtk_grid(nds_path, grid_ξ', vtp_cls, ascii=false) do vtk end
#     close(fid)
# end