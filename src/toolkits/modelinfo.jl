#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : modelinfo.jl                                                               |
|  Description: Export or import model in JSON format                                      |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 01. check_datasize                                                         |
|               02. @memcheck                                                              |
|               03. exportmodel                                                             |
|               04. importmodel                                                             |
+==========================================================================================#

export check_datasize
export @memcheck
export exportmodel
export importmodel

function check_datasize(args::     DeviceArgs{T1, T2}, 
                        grid::     DeviceGrid{T1, T2},
                        mp  :: DeviceParticle{T1, T2},
                        attr:: DeviceProperty{T1, T2},
                        bc  ::DeviceVBoundary{T1, T2}) where {T1, T2}
    args_mem  = Base.summarysize(args)
    grid_mem  = Base.summarysize(grid)
    mp_mem    = Base.summarysize(mp)
    attr_mem  = Base.summarysize(attr)
    bc_mem    = Base.summarysize(bc)
    mem_total = args_mem + grid_mem + mp_mem + attr_mem + bc_mem

    args_size = lpad(@sprintf("%.2f", args_mem / 1024 ^ 3), 5)
    grid_size = lpad(@sprintf("%.2f", grid_mem / 1024 ^ 3), 5)
    mp_size   = lpad(@sprintf("%.2f", mp_mem   / 1024 ^ 3), 5)
    attr_size = lpad(@sprintf("%.2f", attr_mem / 1024 ^ 3), 5)
    bc_size   = lpad(@sprintf("%.2f", bc_mem   / 1024 ^ 3), 5)
    
    argsp  = lpad(@sprintf("%.2f", args_mem / mem_total * 100), 5)
    gridp  = lpad(@sprintf("%.2f", grid_mem / mem_total * 100), 5)
    mpp    = lpad(@sprintf("%.2f", mp_mem   / mem_total * 100), 5)
    attr_p = lpad(@sprintf("%.2f", attr_mem / mem_total * 100), 5)
    bcp    = lpad(@sprintf("%.2f", bc_mem   / mem_total * 100), 5)

    tbar = string("─"^19, "┬", "─"^11, "┬", "─"^8)
    bbar = string("─"^19, "┴", "─"^11, "┴", "─"^8)
    @info """model data size
    $(tbar)
    MPM model args     │ $args_size GiB │ $argsp %
    background grid    │ $grid_size GiB │ $gridp %
    material particles │ $mp_size GiB │ $mpp %
    particle properties│ $attr_size GiB │ $attr_p %
    boundary conditions│ $bc_size GiB │ $bcp %
    $(bbar)
    """
    return nothing
end

macro memcheck(expr)
    return quote        
        data = Base.summarysize($(esc(expr))) / 1024 ^ 3
        data_str = @sprintf("%.2f GiB", data)
        @info "💾 $data_str" # print info
    end
end

function exportmodel(
    args      ::     DeviceArgs{T1, T2}, 
    grid      ::     DeviceGrid{T1, T2},
    mp        :: DeviceParticle{T1, T2},
    attr      :: DeviceProperty{T1, T2},
    bc        ::DeviceVBoundary{T1, T2}; 
    model_file::String = "temp.mpm",
    compress  ::Bool   = true
) where {T1, T2}
    # get model file path and timestamp
    model_file == "temp.mpm" && (model_file = joinpath(
        args.project_path, args.project_name, args.project_name * ".mpm"))
    ctime = now()
    timestamp = Int64(floor(datetime2unix(ctime)))
    
    # print info for the user
    model_status = compress ? "compressed model" : "uncompressed model"
    println("\e[1;33m[ Export:\e[0m", " $(model_status) data (.mpm)")
    println(" ", "\e[1;33m⦿\e[0m", " time stamp: $(ctime)")

    # save model data
    # 1) compress header, 2) timestamp, 3) args, 4) grid, 5) mp, 6) attr, 7) bc
    open(model_file, "w") do io
        # compress header（0x01 means compressed, 0x00 means uncompressed）
        compress ? write(io, UInt8(0x01)) : write(io, UInt8(0x00))
        # write timestamp        
        write(io, timestamp)
        cstream = compress ? ZlibCompressorStream(io) : io
        serialize(cstream, args)
        println(" ", "\e[1;33m⦿\e[0m", "       args: ", "\e[0;33mdone\e[0m")
        serialize(cstream, grid)
        println(" ", "\e[1;33m⦿\e[0m", "       grid: ", "\e[0;33mdone\e[0m")
        serialize(cstream, mp  )
        println(" ", "\e[1;33m⦿\e[0m", "         mp: ", "\e[0;33mdone\e[0m")
        serialize(cstream, attr)
        println(" ", "\e[1;33m⦿\e[0m", "       attr: ", "\e[0;33mdone\e[0m")
        serialize(cstream, bc  )
        println(" ", "\e[1;33m⦿\e[0m", "         bc: ", "\e[0;33mdone\e[0m")
        close(cstream)
    end

    # print info
    modelsize = round(filesize(model_file) / 1024^3, digits=2)
    println(" ", "\e[1;33m⦿\e[0m", "  file size: ", "$(modelsize) GiB")
    println(" ", "\e[1;33m⦿\e[0m", "  file path: ", model_file)
    return nothing
end

function importmodel(mpm_file::String)
    # input file check
    isfile(mpm_file) || throw(ArgumentError("file not found: $mpm_file"))
    splitext(mpm_file)[end] == ".mpm" || throw(ArgumentError(
        "invalid file format: $(splitext(mpm_file)[end])"))

    # read model data
    args, grid, mp, attr, bc = open(mpm_file, "r") do io
        # get comrression header
        compression_flag = read(io, UInt8)
        compress = compression_flag == 0x01 ? true : false
        # print info for the user
        model_status = compress ? "compressed model" : "uncompressed model"
        println("\e[1;33m[ Import:\e[0m", " $(model_status) data (.mpm)")
        # get the timestamp
        timestamp = read(io, Int64)
        creation_time = unix2datetime(Float64(timestamp))
        println(" ", "\e[1;33m⦿\e[0m", " time stamp: $(creation_time)")
        # read model data
        cstream = compress ? ZlibDecompressorStream(io) : io
        args = deserialize(cstream)
        println(" ", "\e[1;33m⦿\e[0m", "       args: ", "\e[0;33mdone\e[0m")
        grid = deserialize(cstream)
        println(" ", "\e[1;33m⦿\e[0m", "       grid: ", "\e[0;33mdone\e[0m")
        mp   = deserialize(cstream)
        println(" ", "\e[1;33m⦿\e[0m", "         mp: ", "\e[0;33mdone\e[0m")
        attr = deserialize(cstream)
        println(" ", "\e[1;33m⦿\e[0m", "       attr: ", "\e[0;33mdone\e[0m")
        bc   = deserialize(cstream)
        println(" ", "\e[1;33m⦿\e[0m", "         bc: ", "\e[0;33mdone\e[0m")
        close(cstream)
        println()
        args, grid, mp, attr, bc
    end

    return args, grid, mp, attr, bc
end