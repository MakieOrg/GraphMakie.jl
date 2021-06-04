module GraphMakie

using NetworkLayout
using LightGraphs
using Makie
using LinearAlgebra

using DocStringExtensions
import Makie: DocThemer, ATTRIBUTES, project, automatic

include("recipes.jl")
include("interaction.jl")

end
