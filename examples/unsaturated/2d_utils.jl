#==========================================================================================+
| Extension struct for the traction boundary                                               |
+==========================================================================================#
struct TractionBoundary{T<:AbstractArray} <: UserBoundaryExtra
    id::T
end

@user_struct TractionBoundary

struct NEWGrid{T<:AbstractArray} <: UserGridExtra
    nl::T
end

@user_struct NEWGrid

struct NEWProperty{T2} <: UserPropertyExtra
    S_min::T2 # min saturation for SWCC
    S_max::T2 # max saturation for SWCC
    P_ref::T2 # fitting parameter for SWCC
    λ    ::T2 # fitting parameter for SWCC
    α    ::T2 # fitting parameter for HCC
    β    ::T2 # fitting parameter for HCC
    Ŵy   ::T2 # infiltration velocity
end

@user_struct NEWGrid
#==========================================================================================#

using HDF5
using ProgressMeter
using WriteVTK

@views function Tanimation(args::DeviceArgs2D{T1, T2}) where {T1, T2}
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
            S    = fid["group$(i)/saturation"] |> read;
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
                vtk["saturation"] = S                  ;
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

function plotresult(hdf5_file)
    fid = h5open(hdf5_file, "r")
    timeset = [1, 2, 3, 4, 5, 6]
    init_ξ = fid["mp_coords0"] |> read
    np = size(init_ξ, 1)
    σw = zeros(np, 6)
    S  = zeros(np, 6)
    ξ  = zeros(np, 2, 6)
    for i in axes(timeset, 1)
        σw[:, i] .= vec(fid["group$(timeset[i])/pressure_w"] |> read)
        S[:, i]  .= vec(fid["group$(timeset[i])/saturation"] |> read)
        ξ[:, :, i] .= fid["group$(timeset[i])/coords"    ] |> read
    end
    close(fid)

    fig = Figure()
    ax1 = Axis(fig[1, 1], aspect=1)
    ax2 = Axis(fig[1, 2], aspect=DataAspect())
    ax3 = Axis(fig[2, 1], aspect=1)
    ax4 = Axis(fig[2, 2], aspect=DataAspect())

    

    display(fig)
    return nothing
end