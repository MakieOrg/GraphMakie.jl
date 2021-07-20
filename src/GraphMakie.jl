module GraphMakie

using NetworkLayout
using LightGraphs
using Makie
using LinearAlgebra

import Makie: DocThemer, ATTRIBUTES, project, automatic

include("beziercurves.jl")
include("recipes.jl")
include("interaction.jl")

end
