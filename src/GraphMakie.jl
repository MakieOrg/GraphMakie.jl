module GraphMakie

using NetworkLayout
using LightGraphs
using AbstractPlotting
using LinearAlgebra

using DocStringExtensions
import AbstractPlotting: DocThemer, ATTRIBUTES

include("recipes.jl")
include("interaction.jl")

end
