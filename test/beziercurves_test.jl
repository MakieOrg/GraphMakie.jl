using Test
using GraphMakie
using GraphMakie.Makie.GeometryBasics

using GraphMakie: BezierPath, MoveTo, LineTo, CurveTo, interpolate, discretize, tangent, waypoints

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
        MoveTo(Point2f0(0,0)),
        LineTo(Point2f0(1,0)),
        LineTo(Point2f0(1,1)),
        LineTo(Point2f0(0,1)),
        LineTo(Point2f0(0,0)),
                ])
    ts = 0:0.1:1.0
    pos = map(t->interpolate(path,t), ts)
    tan = map(t->tangent(path,t), ts)
    fig, ax, p = scatter(pos)
    arrows!(pos, tan; lengthscale=0.1)
end

@testset "natural spline constructor" begin
    path = BezierPath(Point2f0(0.0,0.0), Point2f0(0.0,1.0))
    @test length(path.commands) == 2
    @test path.commands[1] isa MoveTo
    @test path.commands[2] isa LineTo

    path = BezierPath(Point2f0(0.0,0.0),
                      Point2f0(0.0,1.0),
                      Point2f0(1.0,1.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = BezierPath(Point2f0(0.0,0.0),
                      Point2f0(-0.5,1.0),
                      Point2f0(0.5,1.0),
                      Point2f0(0.0,0.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = BezierPath(Point2f0(0.0,0.0),
                      Point2f0(-0.5,1.0),
                      Point2f0(1.5,1.0),
                      Point2f0(2.0,0.0))
    lines(discretize(path))
    scatter!(waypoints(path))

    path = BezierPath(Point2f0(0.0,0.0),
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
    path = BezierPath(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(1,0)*10
    p2 = Point2f0(1,1)
    t2 = Point2f0(0,1)*5
    path = BezierPath(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(-1,0)
    p2 = Point2f0(1,1)
    t2 = Point2f0(0,1)
    path = BezierPath(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))

    p1 = Point2f0(0,0)
    t1 = Point2f0(-.5,.5)
    p2 = Point2f0(0,0)
    t2 = Point2f0(-.5,-.5)
    path = BezierPath(p1, p2; tangents=(t1, t2))
    lines(discretize(path))
    scatter!(waypoints(path))
end

@testset "test beziersegments recipe" begin
    using GraphMakie: beziersegments
    paths = [BezierPath(rand(Point2f0, 4)...) for _ in 1:4]
    fig, ax, p = beziersegments(paths; linewidth=[2,4,6,8], color=[1,2,3,4])
    p.attributes
end
