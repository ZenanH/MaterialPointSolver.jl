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

"""
    DebugPlot

Description:
---
This struct is used to store the WGLMakie configurations for the in-situ visualization of 
the simulation.

- `axis::Bool=false: Whether to show the axis or not.
- `colorbar     ::Bool   = true`                       : Whether to show the colorbar or not
- `colormap     ::Symbol = :redblue`                   : The colormap used for the visualization
- `pointsize    ::Float32= 1e-2`                       : The size of the points
- `attr         ::Symbol = :ϵq`                        : The attribute used for the visualization
- `colorby      ::String = "equivalent plastic strain"`: The attribute used for the color
- `camoffset    ::Float32= 1.0`                        : The camera offset
- `cbfontsize   ::Int32  = 24`                         : The font size of the colorbar labels
- `axfontsize   ::Int32  = 24`                         : The font size of the axis labels
- `tickspace    ::Int32  = 9`                          : The space between the ticks
- `tickprecision::Int32  = 2`                          : The precision of the ticks
- `tickformat   ::String="f"`                          : The format of the ticks f/e
"""
@kwdef struct DebugPlot
    axis         ::Bool    = false
    colorbar     ::Bool    = true
    colormap     ::Symbol  = :redblue
    pointsize    ::Float32 = 1e-2
    attr         ::Symbol  = :ϵq
    colorby      ::String  = "equivalent plastic strain"
    camoffset    ::Float32 = 1.0
    cbfontsize   ::Int32   = 24
    axfontsize   ::Int32   = 24
    tickspace    ::Int32   = 9
    tickprecision::Int32   = 2
    tickformat   ::String  = "f"
end

def_debug = DebugPlot()

"""
    DebugConfig

Description:
---
This struct is used to store the debug configurations for the in-situ visualization of the
simulation.
"""
@kwdef struct DebugConfig
    plot::DebugPlot=def_debug
end