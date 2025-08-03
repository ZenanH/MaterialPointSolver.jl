#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Datatransfer on Intel GPU                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

dev_backend(::Val{:oneapi}) = oneAPIBackend()

@inline host2device(::oneAPIBackend, host::NamedTuple) = KAupload(oneArray, host)
@inline host2device(::oneAPIBackend, host::NamedTuple, hosts::NamedTuple...) = (KAupload(oneArray, host), map(nt -> KAupload(oneArray, nt), hosts)...)