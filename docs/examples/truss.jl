#=
# Stress on Truss
In this example we'll plot an animation of the stress on some truss structure
using `GaphMakie.jl` and [`NetworkDynamics.jl`](https://github.com/PIK-ICoN/NetworkDynamics.jl)

![truss animation](truss.mp4)
=#
using NetworkDynamics
using OrdinaryDiffEqTsit5
using Graphs
using GraphMakie
using LinearAlgebra: norm
using Printf
using CairoMakie
CairoMakie.activate!()

#=
## Definition of the dynamical system
Define dynamic model for NetworkDynamics
For more information on the models see the corresponding [example in the NetworkDynamics.jl docs](https://juliadynamics.github.io/NetworkDynamics.jl/stable/generated/stress_on_truss/).
=#

function fixed_g(pos, x, p, t)
    pos .= p
end
vertex_fix = VertexModel(g=fixed_g, psym=[:xfix, :yfix], outsym=[:x, :y], ff=NoFeedForward())
function free_f(dx, x, Fsum, (M, γ, g), t)
    v = view(x, 1:2)
    dx[1:2] .= (Fsum .- γ .* v) ./ M
    dx[2] -= g
    dx[3:4] .= v
    nothing
end
vertex_free = VertexModel(f=free_f, g=3:4, sym=[:vx=>0, :vy=>0, :x, :y],
                             psym=[:M=>10, :γ=>200, :g=>9.81], insym=[:Fx, :Fy])
function edge_g!(F, pos_src, pos_dst, (K, L), t)
    dx = pos_dst[1] - pos_src[1]
    dy = pos_dst[2] - pos_src[2]
    d = sqrt(dx^2 + dy^2)
    Fabs = K * (L - d)
    F[1] = Fabs * dx / d
    F[2] = Fabs * dy / d
    nothing
end
function observedf(obsout, u, pos_src, pos_dst, (K, L), t)
    dx = pos_dst[1] .- pos_src[1]
    dy = pos_dst[2] .- pos_src[2]
    d = sqrt(dx^2 + dy^2)
    obsout[1] = K * (L - d)
    nothing
end
beam = EdgeModel(g=AntiSymmetric(edge_g!), psym=[:K=>0.5e6, :L], outsym=[:Fx, :Fy], obsf=observedf, obssym=[:Fabs])
nothing #hide

#=
Set up graph topology and initial positions.
=#
N = 5
dx = 1.0
shift = 0.2
g = SimpleGraph(2*N + 1)
for i in 1:N
    add_edge!(g, i, i+N); add_edge!(g, i, i+N)
    if i < N
        add_edge!(g, i+1, i+N); add_edge!(g, i, i+1); add_edge!(g, i+N, i+N+1)
    end
end
add_edge!(g, 2N, 2N+1)
pos0 = zeros(Point2f, 2N + 1)
pos0[1:N] = [Point((i-1)dx,0) for i in 1:N]
pos0[N+1:2*N] = [Point(i*dx + shift, 1) for i in 1:N]
pos0[2N+1] = Point(N*dx + 1, -1)
fixed = [1,4] # set fixed vertices
nothing #hide

#=
Now we can create:
- the `Network` object (rhs of the differential equation),
- the initial state `u0` and
- the parameters for the individual components.
=#
verts = VertexModel[vertex_free for i in 1:nv(g)]
for i in fixed
    verts[i] = vertex_fix # use the fixed vertex for the fixed points
end
nw = Network(g, verts, beam)
u0 = NWState(nw)
## set x/y initial conditions and xfix/yfix parameters
for i in eachindex(pos0, verts)
    if i in fixed
        u0.p.v[i, :xfix] = pos0[i][1]
        u0.p.v[i, :yfix] = pos0[i][2]
    else
        u0.v[i, :x] = pos0[i][1]
        u0.v[i, :y] = pos0[i][2]
    end
end
## set L for edges
for (i,e) in enumerate(edges(g))
    u0.p.e[i, :L] = norm(pos0[src(e)] - pos0[dst(e)])
end
## set damping and mass for "big mass" at the end
u0.p.v[11, :M] = 200
u0.p.v[11, :γ] = 100
nothing #hide

#=
With rhs, parameters and initial conditions constructed we can integrate the system.
=#
tspan = (0.0, 12.0)
prob = ODEProblem(nw, uflat(u0), tspan, pflat(u0))
sol  = solve(prob, Tsit5())
nothing #hide

#=
## Plot the solution
=#

fig = Figure(size=(1000,550));
fig[1,1] = title = Label(fig, "Stress on truss", fontsize=30)
title.tellwidth = false

fig[2,1] = ax = Axis(fig)
ax.aspect = DataAspect();
hidespines!(ax); # no borders
hidedecorations!(ax); # no grid, axis, ...
limits!(ax, -0.1, pos0[end][1]+0.3, pos0[end][2]-0.5, 1.15) # axis limits to show full plot

## get the maximum force during the simulation to get the color scale
## It is only possible to access `:Fabs` directly becaus we've define the observable function for it!
(fmin, fmax) = 0.3 .* extrema(Iterators.flatten(sol(sol.t, idxs=eidxs(nw, :, :Fabs))))

p = graphplot!(ax, g;
               edge_width = 4.0,
               node_size = 3*sqrt.(try u0.p.v[i, :M] catch; 10.0 end for i in 1:nv(g)),
               nlabels = [i in fixed ? "Δ" : "" for i in 1:nv(g)],
               nlabels_align = (:center,:top),
               nlabels_fontsize = 30,
               elabels = ["edge $i" for i in 1:ne(g)],
               elabels_side = Dict(ne(g)  => :right),
               edge_color = [0.0 for i in 1:ne(g)],
               edge_attr = (colorrange=(fmin,fmax),
                          colormap=:diverging_bkr_55_10_c35_n256))

## draw colorbar
fig[3,1] = cb = Colorbar(fig, get_edge_plot(p), label = "Axial force", vertical=false)

T = tspan[2]
fps = 30
trange = range(0.0, sol.t[end], length=Int(T * fps))
record(fig, "truss.mp4", trange; framerate=fps) do t
    title.text = @sprintf "Stress on truss (t = %.2f )" t
    s_at_t = NWState(sol, t)
    for i in eachindex(pos0)
        p[:node_pos][][i] = (s_at_t.v[i, :x], s_at_t.v[i, :y])
    end
    p[:node_pos][] = p[:node_pos][]
    load = s_at_t.e[:, :Fabs]
    p.edge_color[] = load
    p.elabels = [@sprintf("%.0f", l) for l in load]
    fig
end
nothing #hide
