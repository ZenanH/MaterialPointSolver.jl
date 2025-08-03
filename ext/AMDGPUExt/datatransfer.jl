#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Datatransfer on AMD GPU                                                    |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

dev_backend(::Val{:cuda}) = ROCBackend()

@inline host2device(::ROCBackend, host::NamedTuple) = KAupload(ROCArray, host)
@inline host2device(::ROCBackend, host::NamedTuple, hosts::NamedTuple...) = (KAupload(ROCArray, host), map(nt -> KAupload(ROCArray, nt), hosts)...)