# # Examples
# ## Hello to the examples page!
using CairoMakie
CairoMakie.activate!() # hide
AbstractPlotting.inline!(true) # hide

plot([1,2,3],[1,2,3])

# maybe another one?
x = LinRange(0, 10, 100)
y = sin.(x)
lines(x, y)
