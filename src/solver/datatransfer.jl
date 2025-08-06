export host2device, device2host!
export memorysize

host2device(::CPU, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(Array, grid), KAupload(Array, mpts)

# Device to Host transfer (fully)
@inline device2host!(dst::AbstractArray, src::AbstractArray) = (copyto!(dst, src); dst)
@inline device2host!(dst::NamedTuple, src::NamedTuple) = (@inbounds foreach(n -> device2host!(getfield(dst,n), getfield(src,n)), fieldnames(typeof(dst))); dst)
@inline device2host!(::Any, ::Any) = nothing
@inline function device2host!(dst::DeviceParticle{T1, T2}, src::DeviceParticle{T1, T2}) where {T1, T2}
    @inbounds for n in fieldnames(typeof(dst))
        device2host!(getfield(dst,n), getfield(src,n))
    end
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
    bytes = obj isa AbstractArray                    ? sizeof(eltype(obj)) * length(obj)       :
            obj isa NamedTuple                       ? sum(memorysize(v) for v in values(obj)) :
            (f = fieldnames(typeof(obj)); isempty(f) ? sizeof(obj)                             : sum(memorysize(getfield(obj, n)) for n in f))
    return bytes / 1024^3
end