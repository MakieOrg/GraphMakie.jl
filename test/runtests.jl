using CairoMakie
using GraphMakie
using GraphMakie.Graphs
using GraphMakie.NetworkLayout
using Makie.Colors
using StaticArrays
using Test
using StableRNGs

NetworkLayout.DEFAULT_RNG[] = StableRNG

include("beziercurves_test.jl")

@testset "GraphMakie.jl" begin
    g = Observable(wheel_digraph(10))
    f, ax, p = graphplot(g)
    f, ax, p = graphplot(g, node_attr=(visible=false,))
    f, ax, p = graphplot(g, node_attr=(;visible=true))

    # try to update graph
    add_edge!(g[], 2, 4)
    g[] = g[]
    # try to update network
    p.layout = SFDP()

    # update node observables
    p.node_color = :blue
    p.node_size = 30

    # update edge observables
    p.edge_width = 5.0
    p.edge_color = :green

    # it should be also possible to pass multiple values
    f, ax, p = graphplot(g, node_color=[rand([:blue,:red,:green]) for i in 1:nv(g[])])
end

@testset "graph without edges" begin
    g = SimpleDiGraph(10)
    graphplot(g)
    graphplot(g; edge_plottype=:linesegments)
    graphplot(g; edge_plottype=:beziersegments)
    g = SimpleGraph(10)
    graphplot(g)
    graphplot(g; edge_plottype=:linesegments)
    graphplot(g; edge_plottype=:beziersegments)
end

@testset "selfedge without neighbors" begin
    g = SimpleGraph(10)
    add_edge!(g, 1, 1)
    graphplot(g)
end

@testset "small graphs" begin
    g = SimpleGraph(1)
    fig, ax, p = graphplot(g)

    g = complete_graph(2)
    graphplot(g)
end

@testset "single line width per edge" begin
    g = complete_graph(3)
    graphplot(g)
    graphplot(g; edge_width=10)
    graphplot(g; edge_width=[5, 10, 15])
    graphplot(g; edge_width=[5, 10, 5, 10, 5, 10])
end

@testset "Hover, click and drag Interaction" begin
    g = wheel_graph(10)
    add_edge!(g, 1, 1)
    f, ax, p = graphplot(g,
                         edge_width = [3.0 for i in 1:ne(g)],
                         edge_color = [colorant"black" for i in 1:ne(g)],
                         node_size = [20 for i in 1:nv(g)],
                         node_color = [colorant"red" for i in 1:nv(g)])
    deregister_interaction!(ax, :rectanglezoom)

    function edge_hover_action(s, idx, event, axis)
        p.edge_width[][idx]= s ? 6.0 : 3.0
        p.edge_width[] = p.edge_width[] # trigger observable
    end
    ehover = EdgeHoverHandler(edge_hover_action)
    register_interaction!(ax, :ehover, ehover)

    function node_hover_action(s, idx, event, axis)
        p.node_size[][idx] = s ? 40 : 20
        p.node_size[] = p.node_size[] # trigger observable
    end
    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)

    function node_click_action(idx, args...)
        p.node_color[][idx] = rand(RGB)
        p.node_color[] = p.node_color[]
    end
    nclick = NodeClickHandler(node_click_action)
    register_interaction!(ax, :nclick, nclick)

    function edge_click_action(idx, args...)
        p.edge_color[][idx] = rand(RGB)
        p.edge_color[] = p.edge_color[]
    end
    eclick = EdgeClickHandler(edge_click_action)
    register_interaction!(ax, :eclick, eclick)

    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)

    mutable struct EdgeDragAction
        init::Union{Nothing, Point2f} # save click position
        src::Union{Nothing, Point2f}  # save src vertex position
        dst::Union{Nothing, Point2f}  # save dst vertex position
        EdgeDragAction() = new(nothing, nothing, nothing)
    end
    function (action::EdgeDragAction)(state, idx, event, axis)
        edge = collect(edges(g))[idx]
        if state == true
            if action.src===action.dst===action.init===nothing
                action.init = event.data
                action.src = p[:node_pos][][src(edge)]
                action.dst = p[:node_pos][][dst(edge)]
            end
            offset = event.data - action.init
            p[:node_pos][][src(edge)] = action.src + offset
            p[:node_pos][][dst(edge)] = action.dst + offset
            p[:node_pos][] = p[:node_pos][] # trigger change
        elseif state == false
            action.src = action.dst = action.init =  nothing
        end
    end
    edrag = EdgeDragHandler(EdgeDragAction())
    register_interaction!(ax, :edrag, edrag)
end

@testset "align_to_direction" begin
    using GraphMakie: align_to_dir
    using LinearAlgebra: normalize

    @test align_to_dir((:left, :center)) ≈ Point(1.0, 0)
    @test align_to_dir((:right, :center)) ≈ Point(-1.0, 0)
    @test align_to_dir((:center, :center)) ≈ Point(0.0, 0)

    @test align_to_dir((:left, :top)) ≈ normalize(Point(1.0, -1.0))
    @test align_to_dir((:right, :top)) ≈ normalize(Point(-1.0, -1))
    @test align_to_dir((:center, :top)) ≈ normalize(Point(0.0, -1))

    @test align_to_dir((:left, :bottom)) ≈ normalize(Point(1.0, 1.0))
    @test align_to_dir((:right, :bottom)) ≈ normalize(Point(-1.0, 1.0))
    @test align_to_dir((:center, :bottom)) ≈ normalize(Point(0.0, 1.0))

    # g = complete_graph(9)
    # nlabels_align = vec(collect(Iterators.product((:left,:center,:right),(:top,:center,:bottom))))
    # nlabels= repr.(nlabels_align)
    # graphplot(g; nlabels, nlabels_align, nlabels_distance=20)
end

@testset "test plot accessors" begin
    g = complete_graph(10)
    nlabels = repr.(1:nv(g))
    elabels = repr.(1:ne(g))
    fig, ax, p = graphplot(g)
    @test get_edge_plot(p) isa EdgePlot
    @test get_node_plot(p)[1][] == p[:node_pos][]
    @test get_arrow_plot(p).visible[] == false
    @test get_nlabel_plot(p) === nothing
    @test get_elabel_plot(p) === nothing

    fig, ax, p = graphplot(g; nlabels)
    @test get_nlabel_plot(p)[:text][] == nlabels
    @test get_elabel_plot(p) === nothing

    fig, ax, p = graphplot(g; elabels)
    @test get_nlabel_plot(p) === nothing
    @test get_elabel_plot(p)[:text][] == elabels

    fig, ax, p = graphplot(g; elabels, nlabels)
    @test get_nlabel_plot(p)[:text][] == nlabels
    @test get_elabel_plot(p)[:text][] == elabels
end

@testset "test Pointf" begin
    using GraphMakie: Pointf

    p = Point(0.0, 0.0)
    @test typeof(Pointf(p)) == Point2f
    @test Pointf(p) == Point2f(p)

    p = Point(1, 0)
    @test typeof(Pointf(p)) == Point2f
    @test Pointf(p) == Point2f(p)

    p = Point(0.0, 1.0, 2.0)
    @test typeof(Pointf(p)) == Point3f
    @test Pointf(p) == Point3f(p)

    @test Pointf(0.0, 0.0, 0.0) isa Point3f
    @test Pointf(1.0, 1.0) isa Point2f

    @test Pointf((0.0, 0.0, 0.0 )) isa Point3f
    @test Pointf((1.0, 1.0)) isa  Point2f

    @test Pointf([0.0, 0.0, 0.0]) isa Point3f
    @test Pointf([1.0, 1.0]) isa  Point2f

    @test Pointf(SA[0.0, 0.0, 0.0]) isa Point3f
    @test Pointf(SA[1.0, 1.0]) isa  Point2f

    g = complete_graph(3)
    pos1 = [(0,0),
            (1,1),
            (0,1)]
    pos2 = [[0,0],
            [1,1],
            [0,1]]
    pos3 = [SA[0,0],
            SA[1,1],
            SA[0,1]]
    @test isconcretetype(typeof(Pointf.(pos1)))
    @test isconcretetype(typeof(Pointf.(pos2)))
    @test isconcretetype(typeof(Pointf.(pos3)))
    graphplot(g; layout=(x)->pos1)
    graphplot(g; layout=(x)->pos2)
    graphplot(g; layout=(x)->pos3)
end

@testset "test edges distance parameters" begin
    g = SimpleDiGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 1)

    graphplot(g)

    add_edge!(g, 1, 3)

    graphplot(g; curve_distance=0.2)
    graphplot(g; curve_distance=0.2, curve_distance_usage=false)
    graphplot(g; curve_distance=0.2, curve_distance_usage=true)

    graphplot(g; curve_distance=collect(0.1:0.1:0.5))
    graphplot(g; curve_distance=collect(0.1:0.1:0.5), curve_distance_usage=true)
end

@testset "separate linestyle per edge" begin
    p = Point2f[(0,0), (0, 1), (0,0), (1,0)]
    @test_broken display(linesegments(p; linestyle = [:dot, :dash]))

    @test_throws ArgumentError graphplot(
        DiGraph([Edge(1 => 2), Edge(2 => 1)]),
        edge_attr = (; linestyle = [:dot, :dash]),
    )

    @test_throws ArgumentError graphplot(
        DiGraph([Edge(1 => 2), Edge(2 => 3)]),
        edge_attr = (; linestyle = [:dot, :dash]),
    )
end

@testset "pan in scene after removal of plot" begin
    g = wheel_graph(10)
    fig, ax, p = graphplot(g)
    delete!(ax, p)
    rem_edge!(g, 1, 2)
    ax.scene.camera.projectionview[] = ax.scene.camera.projectionview.val # Simlulate panning the screen
end

@testset "get correct axis type for arguments" begin
    g = smallgraph(:petersen)
    _, ax, _ = graphplot(g)
    @test ax isa Axis
    _, ax, _ = graphplot(g; layout=Stress(dim=3))
    @test ax isa LScene
    _, ax, _ = graphplot(g; layout=Stress(dim=2))
    @test ax isa Axis
    _, ax, _ = graphplot(g; layout=rand(Point2, nv(g)))
    @test ax isa Axis
    _, ax, _ = graphplot(g; layout=rand(Point3, nv(g)))
    @test ax isa LScene
    _, ax, _ = graphplot(g; layout=_->rand(Point2, nv(g)))
    @test ax isa Axis
    _, ax, _ = graphplot(g; layout=_->rand(Point3, nv(g)))
    @test ax isa LScene
end

include("referencetests.jl")
