export host2device, device2host!
export memorysize

host2device(::CPU, grid::DeviceGrid{T1,T2}, mpts::DeviceParticle{T1,T2}) where {T1,T2} = KAupload(Array, grid), KAupload(Array, mpts)
host2device(::CPU, data) = KAupload(Array, data)

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
    bytes(x) = x isa AbstractArray ? sizeof(eltype(x)) * length(x) :
               x isa NamedTuple ? sum(bytes(v) for v in values(x)) :
               isbits(x) ? sizeof(x) :
               (fn = fieldnames(typeof(x)); isempty(fn) ? 0 : sum(bytes(getfield(x, f)) for f in fn))
    return bytes(obj) / 1024^3      # GiB
end