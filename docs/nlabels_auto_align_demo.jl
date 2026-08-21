# Before/after figures for MakieOrg/GraphMakie.jl#257 (`nlabels_auto_align`).
#
# Honest comparison: package default align `(:left, :bottom)` vs auto-align,
# with the same `nlabels_distance`. Cases where a single fixed align sits on edges.
#
# Kept outside `docs/examples/` so Literate/Documenter do not execute it.
#
# From the GraphMakie.jl repo root:
#   julia --project=/path/to/env --threads=auto -e 'include("docs/nlabels_auto_align_demo.jl")'

using Pkg
const REPO = dirname(@__DIR__)
if !isdefined(Main, :GraphMakie) || !isdefined(Main, :CairoMakie)
    try
        @eval using GraphMakie
        @eval using CairoMakie
        @eval using Graphs
        @eval using GeometryBasics: Point2f
    catch
        Pkg.activate(temp = true)
        Pkg.develop(path = REPO)
        Pkg.add(["CairoMakie", "Graphs", "GeometryBasics"])
        @eval using GraphMakie
        @eval using CairoMakie
        @eval using Graphs
        @eval using GeometryBasics: Point2f
    end
end

CairoMakie.activate!(type = "png", px_per_unit = 2)

const OUTDIR = joinpath(@__DIR__, "src", "assets")
const TMPDIR = mktempdir()
mkpath(OUTDIR)

"""Save figure to a temp path, then copy into the package assets tree."""
function save_asset(name, fig)
    tmp = joinpath(TMPDIR, name)
    dest = joinpath(OUTDIR, name)
    save(tmp, fig)
    cp(tmp, dest; force = true)
    return dest
end

function set_padded_limits!(ax, positions; pad = 0.65)
    xs = float.((p[1] for p in positions))
    ys = float.((p[2] for p in positions))
    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)
    dx = max(xmax - xmin, 1.0)
    dy = max(ymax - ymin, 1.0)
    side = max(dx, dy)
    cx = (xmin + xmax) / 2
    cy = (ymin + ymax) / 2
    half = (0.5 + pad) * side
    limits!(ax, cx - half, cx + half, cy - half, cy + half)
    ax.aspect = DataAspect()
end

# Shared styling; only auto-align differs between panels.
const NLABELS_DISTANCE = 20
const GP_STYLE = (
    node_size = 28,
    node_color = :white,
    node_attr = (; strokecolor = :gray15, strokewidth = 2.0),
    edge_width = 2.0,
    edge_color = (:gray35, 0.95),
    arrow_size = 20,
    arrow_shift = :end,
    nlabels_fontsize = 16,
    nlabels_color = :gray10,
    nlabels_distance = NLABELS_DISTANCE,
)

function panel!(fig, row, col, g, labels, positions; title, auto)
    ax = Axis(fig[row, col]; title, titlesize = 14, titlegap = 8)
    if auto
        graphplot!(ax, g;
            layout = positions,
            nlabels = labels,
            nlabels_auto_align = true,
            GP_STYLE...,
        )
    else
        # Package defaults: nlabels_align = (:left, :bottom), auto off.
        # Same distance as the right panel — only the align strategy changes.
        graphplot!(ax, g;
            layout = positions,
            nlabels = labels,
            nlabels_auto_align = false,
            nlabels_align = (:left, :bottom),
            GP_STYLE...,
        )
    end
    hidedecorations!(ax)
    hidespines!(ax)
    set_padded_limits!(ax, positions; pad = 0.65)
    return ax
end

function before_after_row!(fig, row, g, labels, positions)
    panel!(fig, row, 1, g, labels, positions; title = "default", auto = false)
    panel!(fig, row, 2, g, labels, positions; title = "nlabels_auto_align = true", auto = true)
end

# Star: default NE offset places the hub label on a spoke; peripherals often look fine.
g_star = SimpleDiGraph(5)
for j in 2:5
    add_edge!(g_star, 1, j)
end
add_edge!(g_star, 3, 2)  # Mediator → Outcome
# Hub is Exposure with outgoing arrows: leaves must be effects of exposure
# (not Confounder / Instrument, which would require inward arrows).
labels_star = ["Exposure", "Outcome", "Mediator", "Biomarker", "Symptom"]
positions_star = Point2f[
    (0, 0),          # Exposure
    (1.55, 1.05),    # Outcome
    (-1.50, 1.15),   # Mediator
    (-1.45, -1.20),  # Biomarker
    (1.50, -1.10),   # Symptom
]

# Dense DAG: several nodes have edges into the NE quadrant where the default label sits.
g_dag = SimpleDiGraph(5)
for (i, j) in ((1, 2), (1, 5), (2, 3), (2, 5), (3, 5), (4, 3), (4, 5))
    add_edge!(g_dag, i, j)
end
labels_dag = ["U", "A", "M", "Z", "Y"]
positions_dag = Point2f[
    (-1.2,  1.15),  # U
    (-0.45, 0.15),  # A
    ( 0.55, 0.15),  # M
    ( 1.35, 0.95),  # Z
    ( 0.00, -1.15), # Y
]

fig = Figure(size = (1100, 900), fontsize = 12, backgroundcolor = :white, figure_padding = 18)
before_after_row!(fig, 1, g_star, labels_star, positions_star)
before_after_row!(fig, 2, g_dag, labels_dag, positions_dag)
rowgap!(fig.layout, 16)
colgap!(fig.layout, 16)
save_asset("nlabels_auto_align_before_after.png", fig)

fig_star = Figure(size = (1100, 460), fontsize = 12, backgroundcolor = :white, figure_padding = 18)
before_after_row!(fig_star, 1, g_star, labels_star, positions_star)
save_asset("nlabels_auto_align_star.png", fig_star)

fig_dag = Figure(size = (1100, 460), fontsize = 12, backgroundcolor = :white, figure_padding = 18)
before_after_row!(fig_dag, 1, g_dag, labels_dag, positions_dag)
save_asset("nlabels_auto_align_dag.png", fig_dag)

println("wrote assets under ", OUTDIR)
