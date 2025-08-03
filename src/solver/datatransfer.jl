#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Manage how to transfer the data                                            |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

export dev_backend, host2device, device2host!
export @datasize, totalsize

dev_backend(sym::Symbol=:cpu) = dev_backend(Val(sym))
dev_backend(::Val{S}) where {S} = _missing_backend(S)
_missing_backend(S) = throw(ArgumentError(
    """
    Backend :$S is unavailable
    Please add and ensure compatibility with the corresponding packages:
    Nvidia GPU(s): CUDA.jl    │ backend name: :cuda (lowercase, julia symbol)
    AMD GPU(s)   : AMDGPU.jl  │ backend name: :rocm (lowercase, julia symbol)
    Intel GPU(s) : oneAPI.jl  │ backend name: :oneapi (lowercase, julia symbol)
    Apple GPU(s) : Metal.jl   │ backend name: :metal (lowercase, julia symbol)
    CPU          : by default │ backend name: :cpu (lowercase, julia symbol)
    in your environment, or confirm that you are using a compatible Julia version (>= 1.9) to enable extensions.
    """
))
dev_backend(::Val{:cpu}) = CPU()

@inline host2device(::CPU, host::NamedTuple) = KAupload(Array, host)
@inline host2device(::CPU, host::NamedTuple, hosts::NamedTuple...) = (KAupload(Array, host), map(nt -> KAupload(Array, nt), hosts)...)

@inline function device2host!(
    dst::NamedTuple{K},
    src::NamedTuple{K};
    var::NTuple{N, Symbol}=K
) where {K, N}
    @inbounds for j in 1:N
        idx = Base.fieldindex(NamedTuple{K}, var[j])
        h_v = getfield(dst, idx)
        d_v = getfield(src, idx)
        (h_v isa AbstractArray && d_v isa AbstractArray) && copyto!(h_v, d_v)
    end
end

_is_data(x) = x isa Real || (x isa AbstractArray{<:Real})
_bytes(x::Real) = sizeof(x)
_bytes(A::AbstractArray{T}) where {T<:Real} = sizeof(T) * length(A)
_bytes(::Nothing) = 0
_bytes(::Missing) = 0
_bytes(x) = 0
const _GiB = 2.0^30
_bytes_to_gib(b) = b / _GiB
_nt_bytes(nt::NamedTuple) = sum(v -> _is_data(v) ? _bytes(v) : 0, values(nt))

function _datasize(nts::Vararg{NamedTuple}; names)
    names = collect(String.(names))

    dots            = "."^10
    rows            = Vector{Vector{Any}}()
    body_hlines     = Int[]
    header_cells    = Tuple{Int,Int}[]
    total_rows      = Int[]
    section_heads   = Int[]
    row = 0

    for (i, nt) in enumerate(nts)
        # 表头行
        push!(rows, [names[i], "", "GiB"])
        row += 1
        push!(section_heads, row)
        append!(header_cells, [(row, 1), (row, 3)])
        push!(body_hlines, row)                      # 表头下划线

        # 字段行（只保留数据型）
        section_total = 0
        for (k, v) in pairs(nt)
            _is_data(v) || continue   # 直接跳过
            sz = _bytes(v)
            section_total += sz
            push!(rows, [String(k), dots, _bytes_to_gib(sz)])
            row += 1
        end

        push!(body_hlines, row)                      # Total 上划线

        # Total 行（字段之和）
        push!(rows, ["Total", dots, _bytes_to_gib(section_total)])
        row += 1
        push!(total_rows, row)

        # 段间空行
        if i < length(nts)
            push!(rows, ["", "", ""])
            row += 1
        end
    end

    # → n×3 矩阵
    data = Array{Any}(undef, length(rows), 3)
    for (i, r) in enumerate(rows)
        data[i, :] = r
    end

    cell_align = Dict((r, 1) => :l for r in section_heads)

    pretty_table(
        data;
        show_header        = false,
        tf                 = tf_borderless,
        formatters         = ft_printf("%8.2e", 3),
        body_hlines        = body_hlines,
        body_hlines_format = ('─','─','─','─'),    # 按你版本要求
        cell_alignment     = cell_align,
        highlighters       = (
            hl_cell(header_cells, crayon"bold"),
            hl_row(total_rows,   crayon"bold yellow"),
            hl_col(2,            crayon"dark_gray"),
        ),
    )
end

const _DEF_MOD = @__MODULE__

"""
    @datasize args...

Description:
===
A macro to compute and display the data size of the given arguments in a formatted table. Note
that the type of args is NamedTuple.

Usage:
===
```julia
julia> n = 100_000_000
julia> a = (x=rand(Float32, n, 2), y=rand(Float32, n, 3))
julia> b = (x=rand(Float32, n, 2), y=rand(Float32, n, 3))
julia> @datasize a b
  a                         GiB 
─────────────────────────────────
      x   ..........   7.45e-01
      y   ..........   1.12e+00
─────────────────────────────────
  Total   ..........   1.86e+00 
                     
  b                         GiB 
─────────────────────────────────
      x   ..........   7.45e-01
      y   ..........   1.12e+00
─────────────────────────────────
  Total   ..........   1.86e+00
```
"""
macro datasize(args...)
    names_vec = [string(a) for a in args]
    :( getfield($_DEF_MOD, :_datasize)($(esc.(args)...); names = $names_vec) )
end

"""
    totalsize(nts::Vararg{NamedTuple})

Description:
===
Calculates the total size of data in the provided NamedTuples and returns it in GiB.

Usage:
===
```julia
julia> n = 100_000_000
julia> a = (x=rand(Float32, n, 2), y=rand(Float32, n, 3))
julia> b = (x=rand(Float32, n, 2), y=rand(Float32, n, 3))
julia> totalsize(a, b)
3.725290298461914
julia> totalsize(a)
1.862645149230957
```
"""
function totalsize(nts::Vararg{NamedTuple})
    total_bytes = mapreduce(_nt_bytes, +, nts; init = 0)
    return _bytes_to_gib(total_bytes)
end