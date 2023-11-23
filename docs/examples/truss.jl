#=
# Stress on Truss
In this example we'll plot an animation of the stress on some truss structure
using `GaphMakie.jl` and [`NetworkDynamics.jl`](https://github.com/PIK-ICoN/NetworkDynamics.jl)

![truss animation](truss.mp4)
=#
using NetworkDynamics
using OrdinaryDiffEq
using Graphs
using GraphMakie
using LinearAlgebra: norm, ⋅
using Printf
using CairoMakie
CairoMakie.activate!()

#=
## Definition of the dynamical system
Define dynamic model for NetworkDynamics
=#
function edge_f!(e, v_s, v_d, (K, L), _)
    Δd = view(v_d, 1:2) - view(v_s, 1:2)
    nd = norm(Δd)
    F = K * ( L - nd )
    e[1:2] = F .* view(Δd, 1:2) ./ nd
end

function vertex_free!(du, u, edges, (M, γ, g), _)
    F = sum(edges) - γ * view(u, 3:4) + [0, -M*g]
    du[3:4] .= F ./ M
    du[1:2] = view(u, 3:4)
end

function vertex_fixed!(du, u, edges, _, _)
    du[1:2] .= 0.0
end

edge = StaticEdge(f = edge_f!, dim=2, coupling=:antisymmetric)
vertex_free = ODEVertex(f = vertex_free!, dim=4, sym=[:x, :y, :v, :w])
vertex_fix  = ODEVertex(f = vertex_fixed!, dim=2, sym=[:x, :y])
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

# create networkdynamics object
verts = ODEVertex[vertex_free for i in 1:nv(g)]
for i in fixed
    verts[i] = vertex_fix # use the fixed vertex for the fixed
end
nd = network_dynamics(verts, edge, g);
nothing #hide

#=
write some auxiliary functions to map between state vector of DGL `u` and graph data
=#
x_idxs = [findfirst(isequal(Symbol("x_$i")), nd.syms) for i in 1:nv(g)]

"Set positions `pos` inside dgl state `u`"
function set_pos!(u, pos)
    for (i, idx) in enumerate(x_idxs)
        u[idx] = pos[i][1]
        u[idx+1] = pos[i][2]
    end
end

"Extract vector of Points from dgl state `u`"
function to_pos(u)
    pos = Vector{Point2f}(undef, nv(g))
    for (i, idx) in enumerate(x_idxs)
        pos[i] = Point(u[idx], u[idx+1])
    end
    return pos
end

"Calculate load on edge for given state."
function get_load(u, p, t=0.0)
    gd_nd = nd(u, p, t, GetGD) # exposes underlying graph data struct
    force = Vector{Float64}(undef, ne(g))
    pos = to_pos(u)
    for (i,e) in enumerate(edges(g))
        edgeval = get_edge(gd_nd, i)
        fvec = Point(edgeval[1], edgeval[2])
        dir = pos[dst(e)] .- pos[src(e)]
        force[i] = sign(fvec ⋅ dir) * norm(fvec)
    end
    return force
end
nothing #hide

#=
Set parameters.
=#
M = [10 for i in 1:nv(g)] # mass of the nodes
M[end] = 200 # heavy mass
gc = [9.81 for i in 1:nv(g)] # gravitational constant
γ = [200.0 for i in 1:nv(g)] # damping parameter
γ[end] = 100.0

L = [norm(pos0[src(e)] - pos0[dst(e)]) for e in edges(g)] # length of edges
K = [0.5e6 for i in 1:ne(g)] # spring constant of edges

## bundle parameters for NetworkDynamics
para = collect.((zip(M, γ, gc),
                 zip(K, L)))
nothing #hide

#=
Set initial state and solve the system
=#
u0 = zeros(length(nd.syms))
set_pos!(u0, pos0)

tspan = (0.0, 12.0)
prob = ODEProblem(nd, u0, tspan, para)
sol  = solve(prob, Tsit5());
nothing #hide

#=
## Plot the solution
=#

fig = Figure(size=(1000,550))
fig[1,1] = title = Label(fig, "Stress on truss", fontsize=30)
title.tellwidth = false

fig[2,1] = ax = Axis(fig)

## calculate some values for colorscaling
(fmin, fmax) = 0.3 .* extrema([(map(u->get_load(u, para), sol)...)...])
nlabels = [" " for i in 1:nv(g)]
nlabels[fixed] .= "Δ"
elabels = ["edge $i" for i in 1:ne(g)]
p = graphplot!(ax, g;
               edge_width=4.0,
               node_size=sqrt.(M)*3,
               nlabels=nlabels,
               nlabels_align=(:center,:top),
               nlabels_fontsize=30,
               elabels=elabels,
               elabels_side=Dict(ne(g) => :right),
               edge_color=[0.0 for i in 1:ne(g)],
               edge_attr=(colorrange=(fmin,fmax),
                          colormap=:diverging_bkr_55_10_c35_n256))
hidespines!(ax); hidedecorations!(ax); p[:node_pos][]=to_pos(u0); ax.aspect = DataAspect()
limits!(ax, -0.1, pos0[end][1]+0.3, pos0[end][2]-0.5, 1.15)

## draw colorbar
fig[3,1] = cb = Colorbar(fig, p.plots[1].plots[1], label = "Axial force", vertical=false)

T = tspan[2]
fps = 30
trange = range(0.0, sol.t[end], length=Int(T * fps))
record(fig, "truss.mp4", trange; framerate=fps) do t
    title.text = @sprintf "Stress on truss (t = %.2f )" t
    p[:node_pos][] = to_pos(sol(t))
    load = get_load(sol(t), para)
    p.elabels = [@sprintf("%.0f", l) for l in load]
    p.edge_color[] = load
end
nothing #hide
