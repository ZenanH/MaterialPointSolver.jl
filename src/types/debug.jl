#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : debug.jl                                                                   |
|  Description: Type system for debug configurations in MaterialPointSolver.jl             |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  License    : MIT License                                                                |
+==========================================================================================#

export DebugPlot
export DebugConfig

def_func(mp::DeviceParticle, x::Vector{Symbol}) = getproperty(mp, x[1])

"""
    DebugPlot

Description:
---
This struct is used to store the WGLMakie configurations for the in-situ visualization of 
the simulation.
"""
@kwdef struct DebugPlot
    axis      ::Bool   = false
    colorbar  ::Bool   = true
    colormap  ::Symbol = :redblue
    pointsize ::Real   = 1f-2
    attr      ::Vector{Symbol}
    cbname    ::String        
    cbfontsize::Real     = 24
    axfontsize::Real     = 24
    tickformat::String   = "{:9.2f}"
    calculate ::Function = def_func
end

def_debug = DebugPlot(attr=[:ρs], cbname="ρs")

"""
    DebugConfig

Description:
---
This struct is used to store the debug configurations for the in-situ visualization of the
simulation.
"""
@kwdef struct DebugConfig
    plot::DebugPlot=def_debug
    # here can put more debug configurations
end