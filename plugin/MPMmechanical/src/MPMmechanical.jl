module MPMmechanical

using KernelAbstractions
using HDF5
using Printf
using MaterialPointSolver

export procedure!
export test

include(joinpath(@__DIR__, "kernels.jl"))
include(joinpath(@__DIR__, "solver.jl"))

test() = @info "sad"
end