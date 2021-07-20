using Test
using GraphMakie
using GraphMakie.Makie.GeometryBasics
using GLMakie

using GraphMakie: BezierPath, MoveTo, LineTo, CurveTo, interpolate, discretize, tangent


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
end

@testset "natural spline constructor" begin
    path = BezierPath(Point2f0(0.0,0.0), Point2f0(0.0,1.0))
    @test length(path.commands) == 2
    @test path.commands[1] isa MoveTo
    @test path.commands[2] isa LineTo

    path = BezierPath(Point2f0(0.0,0.0),
                      Point2f0(0.0,1.0),
                      Point2f0(0.2,1.0),
                      Point2f0(1.0,1.0))
    lines(discretize(path))

    path = BezierPath(Point2f0(0.0,0.0),
                      Point2f0(0.0,1.0),
                      Point2f0(0.0,0.0))

end
