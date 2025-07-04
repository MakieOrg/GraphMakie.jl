module GraphMakie

using NetworkLayout
using Graphs
using Makie
using Makie: add_input!, add_constant!
using LinearAlgebra
using SimpleTraits

import Makie: DocThemer, ATTRIBUTES, project, automatic
import DataStructures: DefaultDict, DefaultOrderedDict

include("beziercurves.jl")
include("recipes.jl")
include("interaction.jl")
include("utils.jl")

end
