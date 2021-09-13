using Test
using GraphMakie
using GeometryBasics

using GraphMakie: BezierPath, MoveTo, LineTo, CurveTo, interpolate, discretize, tangent, waypoints, Path, Line

@testset "interpolation and tangets" begin
    path = BezierPath([
        MoveTo(Point2f0(0,0)),
        LineTo(Point2f0(0.5,0.5)),
        CurveTo(Point2f0(1,0),
                Point2f0(1,0),
                Point2f0(1,1)),
        CurveTo(Point2f0(1,2),
                Point2f0(0,2),
                Point2f0(0,1))])

    lines(discretize(path))

    ts = 0:0.1:1.0
    pos = map(t->interpolate(path,t), ts)
    tan = map(t->tangent(path,t), ts)

    scatter!(pos)
    arrows!(pos, tan; lengthscale=0.1)

    path = BezierPath([
        MoveTo(Point3f0(0,0,0)),
        LineTo(Point3f0(1,0,0)),
        LineTo(Point3f0(1,1,1)),
        LineTo(Point3f0(0,1,1)),
        LineTo(Point3f0(0,0,1))])
    ts = 0:0.1:1.0
    pos = map(t->interpolate(path,t), ts)
    tan = map(t->tangent(path,t), ts)
    fig, ax, p = scatter(pos)
    arrows!(pos, tan; lengthscale=0.1)
end

@testset "natural spline constructor" begin
    path = Path(Point2f0(0.0,0.0), Point2f0(0.0,1.0))
    @test path isa GraphMakie.Line

    path = Path(Point2f0(0.0,0.0),
                Point2f0(0.0,1.0),
                Point2f0(1.0,1.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = Path(Point2f0(0.0,0.0),
                Point2f0(-0.5,1.0),
                Point2f0(0.5,1.0),
                Point2f0(0.0,0.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = Path(Point2f0(0.0,0.0),
                Point2f0(-0.5,1.0),
                Point2f0(1.5,1.0),
                Point2f0(2.0,0.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = Path(Point2f0(0.0,0.0),
                Point2f0(-0.5,1.0),
                Point2f0(1.5,1.0),
                Point2f0(2.0,0.0);
                tangents=(Point2f0(-1,0),
                          Point2f0(-1,0)))
    lines(discretize(path), linewidth=10)
    scatter!(waypoints(path))
end

@testset "two points and tangets" begin
    p1 = Point2f0(0,0)
    t1 = Point2f0(0,1)
    p2 = Point2f0(1,1)
    t2 = Point2f0(0,1)
    path = Path(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(1,0)*10
    p2 = Point2f0(1,1)
    t2 = Point2f0(0,1)*5
    path = Path(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(-1,0)
    p2 = Point2f0(1,1)
    t2 = Point2f0(0,1)
    path = Path(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(-.5,.5)
    p2 = Point2f0(0,0)
    t2 = Point2f0(-.5,-.5)
    path = Path(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    p2 = Point2f0(1,0)
    midp = Point2f0[]
    path = Path(p1, midp..., p2)
    lines(discretize(path))
end

@testset "straight lines with radi" begin
    using GraphMakie: plot_controlpoints!
    p1 = Point2f0(0,0)
    p2 = Point2f0(1,-.5)
    p3 = Point2f0(2,.5)
    p4 = Point2f0(3,0)
    path = Path(0.5, p1, p2, p3, p4)
    fig, ax, p = lines(discretize(path))
    plot_controlpoints!(ax, path)
    scatter!(ax, [p1,p2,p3,p4])

    path = Path(0.0, p1, p2, p3, p4)
    fig, ax, p = lines(discretize(path))
end

@testset "test beziersegments recipe" begin
    using GraphMakie: beziersegments
    paths = [Path(rand(Point2f0, 4)...) for _ in 1:4]
    fig, ax, p = beziersegments(paths; linewidth=[2,4,6,8], color=[1,2,3,4])
    p.attributes
end

@testset "selfloop test" begin
    using CairoMakie
    using LightGraphs, GraphMakie
    g = star_graph(10)
    add_edge!(g, 1, 1)
    add_edge!(g, 2, 2)
    fig, ax, p = graphplot(g)
    ax.aspect = DataAspect()

    g = star_graph(10)
    add_edge!(g, 1, 1)
    @test_throws ErrorException graphplot(g; layout=_->rand(Point3f0, 10))
end
