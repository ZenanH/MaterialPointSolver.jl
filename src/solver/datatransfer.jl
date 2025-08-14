export host2device, device2host!
export memorysize

function host2device(::CPU, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} 
    dev_grid = KAupload(Array, grid)
    dev_mpts = KAupload(Array, mpts)
    memsize = memorysize(dev_grid) + memorysize(dev_mpts)
    content = "uploaded $(@sprintf("%.2f", memsize)) GiB data to CPU"
    println("\e[1;32m[▲ I/O:\e[0m \e[0;32m$(content)\e[0m")
    return dev_grid, dev_mpts
end

host2device(::CPU, data) = KAupload(Array, data)

# Device to Host transfer (fully)
@inline device2host!(dst::AbstractArray, src::AbstractArray) = (copyto!(dst, src); dst)
@inline device2host!(dst::NamedTuple, src::NamedTuple) = (@inbounds foreach(n -> device2host!(getfield(dst,n), getfield(src,n)), fieldnames(typeof(dst))); dst)
@inline device2host!(::Any, ::Any) = nothing
@inline function device2host!(dst::DeviceParticle{T1, T2}, src::DeviceParticle{T1, T2}) where {T1, T2}
    @inbounds for n in fieldnames(typeof(dst))
        device2host!(getfield(dst,n), getfield(src,n))
    end
    content = "downloaded $(@sprintf("%.2f", memorysize(src))) GiB data to CPU"
    println("\e[1;31m[▼ I/O:\e[0m \e[0;31m$(content)\e[0m")
    dst
end
# Device to Host transfer (partially)
@inline _d2h!(d, s, fld::Symbol) = copyto!(getfield(d,fld), getfield(s,fld))
@inline _d2h!(d, s, fld::Tuple) = @inbounds foreach(f -> _d2h!(getfield(d,:ext), getfield(s,:ext), f), fld)
@inline _d2h!(::Any, ::Any, _) = nothing
@inline device2host!(dst, src, var::Symbol) = device2host!(dst, src, (var,))
@inline function device2host!(dst::DeviceParticle{T1, T2}, src::DeviceParticle{T1, T2}, vars::Tuple) where {T1, T2}
    @inbounds foreach(v -> _d2h!(dst, src, v), vars)
    dst
end

function memorysize(obj)
    bytes(x) = x isa AbstractArray ? sizeof(eltype(x)) * length(x) :
               x isa NamedTuple ? sum(bytes(v) for v in values(x)) :
               isbits(x) ? sizeof(x) :
               (fn = fieldnames(typeof(x)); isempty(fn) ? 0 : sum(bytes(getfield(x, f)) for f in fn))
    return bytes(obj) / 1024^3      # GiB
end