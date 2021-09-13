using GeometryBasics
using StaticArrays
using LinearAlgebra: normalize, ⋅

####
#### Type definitions
####

# type defs from Julius Krumbiegels PR
# https://github.com/JuliaPlots/Makie.jl/pull/979
abstract type AbstractPath{PT<:AbstractPoint} end

# simple lines get there own type. This helps with type stability
# for graphs without curvy edges and improves performance
struct Line{PT} <: AbstractPath{PT}
    p0::PT
    p::PT
end

abstract type PathCommand{PT<:AbstractPoint} end

struct BezierPath{PT} <: AbstractPath{PT}
    commands::Vector{PathCommand{PT}}
end

struct MoveTo{PT} <: PathCommand{PT}
    p::PT
end

struct LineTo{PT} <: PathCommand{PT}
    p::PT
end

struct CurveTo{PT} <: PathCommand{PT}
    c1::PT
    c2::PT
    p::PT
end

ptype(::Union{AbstractPath{PT}, Type{<:AbstractPath{PT}}}) where {PT} = PT


####
#### Helper functions to work with bezier pathes
####

"""
    interpolate(p::AbstractPath, t)

Parametrize path `p` from `t ∈ [0, 1]`. Return postion at `t`.

TODO: Points are not necessarily evenly spaced!
"""
function interpolate(p::BezierPath{PT}, t) where PT
    @assert p.commands[begin] isa MoveTo
    N = length(p.commands) - 1

    tn = N*t
    seg = min(floor(Int, tn), N-1)
    tseg = tn - seg

    p0 = p.commands[seg+1].p
    return _interpolate(p.commands[seg+2], p0, tseg)
end

interpolate(l::Line{PT}, t) where PT = l.p0 + t*(l.p - l.p0) |> PT

"""
    tangent(p::AbstractPath, t)

Parametrize path `p` from `t ∈ [0, 1]`. Return tangent at `t`.
"""
function tangent(p::BezierPath, t)
    @assert p.commands[begin] isa MoveTo
    N = length(p.commands) - 1

    tn = N*t
    seg = min(floor(Int, tn), N-1)
    tseg = tn - seg

    p0 = p.commands[seg+1].p
    return _tangent(p.commands[seg+2], p0, tseg)
end
tangent(l::Line, _) = normalize(l.p - l.p0)

"""
    discretize(path::AbstractPath)

Return vector of points which represent the given `path`.
"""
function discretize(path::BezierPath{T}) where {T}
    v = Vector{T}()
    for c in path.commands
        _discretize!(v, c)
    end
    return v
end

discretize(l::Line) = [l.p0, l.p]

"""
    _interpolate(c::PathCommand, p0, t)

Returns positions along the path `c` starting from `p0` in range `t ∈ [0, 1]`.
"""
_interpolate(c::LineTo{PT}, p0, t) where {PT} = p0 + t*(c.p - p0) |> PT
function _interpolate(c::CurveTo{PT}, p0, t) where {PT}
    p1, p2, p3 = c.c1, c.c2, c.p
    (1 - t)^3 * p0 + 3(t - 2t^2 + t^3) * p1 + 3(t^2 -t^3) * p2 + t^3 * p3 |> PT
end

"""
    _tangent(c::PathCommand, p0, t)

Returns tangent vector along the path `c` starting from `p0` in range `t ∈ [0, 1]`.
"""
_tangent(c::LineTo, p0, _) = normalize(c.p - p0)
function _tangent(c::CurveTo{PT}, p0, t) where PT
    p1, p2, p3 = c.c1, c.c2, c.p
    normalize(-3(1 - t)^2 * p0 + 3(1 - 4t + 3t^2) * p1 + 3(2t -3t^2) * p2 + 3t^2 * p3) |> PT
end

"""
    _discretize!(v::Vector{AbstractPint}, c::PathCommand)

Append interpolated points of path `c` to pos vector `v`
"""
_discretize!(v::Vector{<:AbstractPoint}, c::Union{MoveTo, LineTo}) = push!(v, c.p)
function _discretize!(v::Vector{<:AbstractPoint}, c::CurveTo)
    N0 = length(v)
    p0 = v[end]
    N = 60 # TODO: magic number of points for discrtization
    resize!(v, N0 + N)
    dt = 1.0/N
    for (i, t) in enumerate(dt:dt:1.0)
        v[N0 + i] = _interpolate(c, p0, t)
    end
end

"""
    waypoints(p::BezierPath)

Returns all the characteristic points of the path. For debug reasons.
"""
function waypoints(p::BezierPath{PT}) where {PT}
    v = PT[]
    for c in p.commands
        if c isa CurveTo
            push!(v, c.c1)
            push!(v, c.c2)
        end
        push!(v, c.p)
    end
    return v
end

"""
    isline(p::AbstractPath)

True if the AbstractPath just represents a straight line.
"""
isline(p::Line) = true
isline(p::BezierPath) = length(p.commands)==2 && p.commands[1] isa MoveTo && p.commands[2] isa LineTo


####
#### Special constructors to create abstract pathes
####

"""
    Path(P::Vararg{PT, N}; tangents, tfactor=.5) where {PT<:AbstractPoint, N}

Create a bezier path by natural cubic spline interpolation of the points `P`.
If there are only two points and no tangents return a straight line.

The `tangets` kw allows you pass two vectors as tangents for the first and the
last point. The `tfactor` affects the curvature on the start and end given some
tangents.
"""
function Path(P::Vararg{PT, N}; tangents=nothing, tfactor=.5) where {PT<:AbstractPoint, N}
    @assert N>2

    # cubic_spline will work for each dimension separatly
    pxyz = _cubic_spline(map(p -> p[1], P)) # get first dimension
    for i in 2:length(PT) # append all other dims
        pxyz = hcat(pxyz, _cubic_spline(map(p -> p[i], P)))
    end

    # create waypoints from waypoints in separat dementions
    WP = SVector{length(P)-1, PT}(PT(p) for p in eachrow(pxyz))

    commands = Vector{PathCommand{PT}}(undef, N)
    commands[1] = MoveTo(P[1])

    # first command, recalculate WP if tangent is given
    first_wp = WP[1]
    if tangents !== nothing
        p1, p2, t = P[1], P[2], normalize(tangents[1])
        dir = p2 - p1
        d = tfactor * norm(dir ⋅ t)
        first_wp = PT(p1+d*t)
    end
    commands[2] = CurveTo(first_wp,
                          2*P[2] - WP[2],
                          P[2])
    # middle commands
    for i in 3:(N-1)
        commands[i] = CurveTo(WP[i-1],
                              2*P[i] - WP[i],
                              P[i])
    end
    # last command, recalculate last WP if tangent is given
    last_wp = (P[N] + WP[N-1])/2
    if tangents !== nothing
        p1, p2, t = P[N-1], P[N], normalize(tangents[2])
        dir = p2 - p1
        d = tfactor * norm(dir ⋅ t)
        last_wp = PT(p2-d*t)
    end
    commands[N] = CurveTo(WP[N-1],
                          last_wp,
                          P[N])

    BezierPath(commands)
end

# same function as above but for just 2 points specificially
function Path(P::Vararg{PT, 2}; tangents=nothing, tfactor=.5) where {PT<:AbstractPoint}
    if tfactor isa NTuple{2, <:Number}
        tf1, tf2 = tfactor
    else
        tf1 = tf2 = tfactor
    end

    p1, p2 = P
    if tangents === nothing
        return Line(p1, p2)
    else
        t1 = normalize(Pointf0(tangents[1]))
        t2 = normalize(Pointf0(tangents[2]))
        len = norm(p2 - p1)
        return BezierPath([MoveTo(p1),
                           CurveTo(PT(p1+len*tf1*t1),
                                   PT(p2-len*tf2*t2),
                                   p2)])
    end
end

"""
    Path(radius::Real, p::Vararg{PT, N}) where {PT<:AbstractPoint,N}

Draw straight lines through the points `p`. Within `radius` of each
point the line will be smoothly connected.
"""
function Path(radius::Real, p::Vararg{PT, N}) where {PT<:AbstractPoint,N}
    if iszero(radius)
        commands = PathCommand{PT}[MoveTo(PT(p[1]))]
        for pos in p[2:end]
            push!(commands, LineTo(PT(pos)))
        end
        return BezierPath(commands)
    else
        points = collect(reverse(p))
        pos = pop!(points)
        commands = PathCommand{PT}[MoveTo(PT(pos))]

        while length(points) > 1
            mid = pop!(points)
            dst = points[end] # due to reverse...
            dir1 = normalize(mid - pos)
            dir2 = normalize(dst - mid)
            r1 = mid - radius * dir1
            r2 = mid + radius * dir2
            c1 = mid - .25 * radius * dir1
            c2 = mid + .25 * radius * dir2
            push!(commands, LineTo(PT(r1)))
            push!(commands, CurveTo(PT(c1),
                                    PT(c2),
                                    PT(r2)))
            pos = r2
        end

        push!(commands, LineTo(points[end]))
        return BezierPath(commands)
    end
end


"""
    _cubic_spline(p)

Given a number of points in one dimension calculate waypoints between them.

    _cubic_spline([x1, x2, x3])

Will return the x coordinates of the waypoints `wp1` and `wp2`.
Those are the first waypoints between in the cubic bezier sense.
"""
function _cubic_spline(p)
    N = length(p) - 1

    M = SMatrix{N,N}(if i==j # diagonal
                         if i==1; 2; elseif i==N; 7; else 4 end
                     elseif i==j+1 # lower
                         if i==N; 2; else 1 end
                     elseif i==j-1 # upper
                         1
                     else
                         0
                     end for i in 1:N, j in 1:N)

    b = SVector{N}(if i == 1
                       p[i] + 2p[i+1]
                   elseif i == N
                       8p[i] + p[i+1]
                   else
                       4p[i] + 2p[i+1]
                   end for i in 1:N)
    return M \ b
end
