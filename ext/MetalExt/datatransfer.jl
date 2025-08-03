#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Datatransfer on Apple GPU                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

dev_backend(::Val{:metal}) = MetalBackend()

@inline host2device(::MetalBackend, host::NamedTuple) = KAupload(MtlArray, host)
@inline host2device(::MetalBackend, host::NamedTuple, hosts::NamedTuple...) = (KAupload(MtlArray, host), map(nt -> KAupload(MtlArray, nt), hosts)...)