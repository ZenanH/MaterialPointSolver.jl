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

def_func(grid, mp, attr, bc, vid) = Array{Float32}(mp.ρs)

"""
    DebugPlot

Description:
---
This struct is used to store the WGLMakie configurations for the in-situ visualization of 
the simulation.

Users need to provice `attr` and `cbname` to specify the attribute to plot and the name of
the colorbar, even you choose `axis = false`.

- `axis::Bool = false` : Whether to show the axis
- `colorbar::Bool = true` : Whether to show the colorbar
- `colormap::Symbol = :readblue` : The colormap for the plot
- `pointsize::Real = 1f-2` : The size of the points
- `pointnum::Int = 500000` : The number of points to plot
- `tickformat::String = "{:9.2f}"` : The format of the ticks
- `calculate::Function = def_func` : The function to calculate the value to plot
- `axfontsize::Real = 24` : The font size of the axis
- `cbfontsize::Real = 24` : The font size of the colorbar
- `cbname::String` : The name of the colorbar
"""
@kwdef struct DebugPlot
    axis      ::Bool     = false
    colorbar  ::Bool     = true
    colormap  ::Symbol   = :redblue
    pointsize ::Real     = 1f-2
    pointnum  ::Int      = 500000
    tickformat::String   = "{:9.2f}"
    calculate ::Function = def_func
    axfontsize::Real     = 24
    cbfontsize::Real     = 24
    cbname    ::String
end


def_debug = DebugPlot(cbname="density")

"""
    DebugConfig

Description:
---
This struct is used to store the debug configurations for the in-situ visualization of the
simulation.

- `plot`::DebugPlot : The WGLMakie configurations for the debug mode
- ... : More debug configurations can be added in the future
"""
@kwdef struct DebugConfig
    plot::DebugPlot=def_debug
    # here can put more debug configurations
end