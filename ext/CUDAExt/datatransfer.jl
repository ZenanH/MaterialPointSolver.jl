#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  Description: Datatransfer on Nvidia GPU                                                 |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, ISTE, Université de Lausanne                                   |
|  Maintainer : Zenan Huo                                                                  |
+==========================================================================================#

dev_backend(::Val{:cuda}) = CUDABackend()

@inline host2device(::CUDABackend, host::NamedTuple) = KAupload(CuArray, host)
@inline host2device(::CUDABackend, host::NamedTuple, hosts::NamedTuple...) = (KAupload(CuArray, host), map(nt -> KAupload(CuArray, nt), hosts)...)