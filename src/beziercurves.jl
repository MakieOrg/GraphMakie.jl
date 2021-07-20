using Makie.GeometryBasics
using Makie.GeometryBasics.StaticArrays
using LinearAlgebra: normalize

####
#### Type definitions
####

abstract type PathCommand{PT<:AbstractPoint} end

struct BezierPath{PT}
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


####
#### Helper functions to work with bezier pathes
####

"""
    interpolate(c::PathCommand, p0, t)

Returns positions along the path `c` starting from `p0` in range `t ∈ [0, 1]`.
"""
interpolate(c::LineTo{PT}, p0, t) where {PT} = p0 + t*(c.p - p0) |> PT
function interpolate(c::CurveTo{PT}, p0, t) where {PT}
    p1, p2, p3 = c.c1, c.c2, c.p
    (1 - t)^3 * p0 + 3(t - 2t^2 + t^3) * p1 + 3(t^2 -t^3) * p2 + t^3 * p3 |> PT
end

"""
    tangent(c::PathCommand, p0, t)

Returns tanget vector along the path `c` starting from `p0` in range `t ∈ [0, 1]`.
"""
tangent(c::LineTo, p0, _) = normalize(c.p - p0)
function tangent(c::CurveTo{PT}, p0, t) where PT
    p1, p2, p3 = c.c1, c.c2, c.p
    normalize(-3(1 - t)^2 * p0 + 3(1 - 4t + 3t^2) * p1 + 3(2t -3t^2) * p2 + 3t^2 * p3) |> PT
end

"""
    discretize!(v::Vector{AbstractPint}, c::PathCommand)

Append interpolated points of path `c` to pos vector `v`
"""
discretize!(v::Vector{<:AbstractPoint}, c::Union{MoveTo, LineTo}) = push!(v, c.p)
function discretize!(v::Vector{<:AbstractPoint}, c::CurveTo)
    N0 = length(v)
    p0 = v[end]
    N = 100
    resize!(v, N0 + N)
    dt = 1.0/N
    for (i, t) in enumerate(dt:dt:1.0)
        v[N0 + i] = interpolate(c, p0, t)
    end
end

"""
    discretize(path::BezierPath)

Return vector of points which represent the given `path`.
"""
function discretize(path::BezierPath{T}) where {T}
    v = Vector{T}()
    for c in path.commands
        discretize!(v, c)
    end
    return v
end

"""
    interpolate(p::BezierPath, t)

Parametrize path `p` from `t ∈ [0, 1]`. Return postion at `t`.

TODO: Points are necessarily evenly spaced!
"""
function interpolate(p::BezierPath{PT}, t) where PT
    @assert p.commands[begin] isa MoveTo
    N = length(p.commands) - 1

    tn = N*t
    seg = min(floor(Int, tn), N-1)
    tseg = tn - seg

    p0 = p.commands[seg+1].p
    return interpolate(p.commands[seg+2], p0, tseg)
end

"""
    tangent(p::BezierPath, t)

Parametrize path `p` from `t ∈ [0, 1]`. Return tangent at `t`.
"""
function tangent(p::BezierPath, t)
    @assert p.commands[begin] isa MoveTo
    N = length(p.commands) - 1

    tn = N*t
    seg = min(floor(Int, tn), N-1)
    tseg = tn - seg

    p0 = p.commands[seg+1].p
    return tangent(p.commands[seg+2], p0, tseg)
end


####
#### Special constructors to create bezier pathes
####

"""
    BezierPath(P::Vararg{PT, N}) where {PT<:AbstractPoint, N}

Create a bezier path by natural cubic spline interpolation of the points `P`.
"""
function BezierPath(P::Vararg{PT, N}) where {PT<:AbstractPoint, N}
    if length(P) == 2
        return BezierPath([MoveTo(P[1]), LineTo(P[2])])
    end

    # cubic_spline will work for each dimension separatly
    pxyz = cubic_spline(map(p -> p[1], P)) # get first dimension
    for i in 2:length(PT) # append all other dims
        pxyz = hcat(pxyz, cubic_spline(map(p -> p[i], P)))
    end

    # create waypoints from waypoints in separat dementions
    WP = SVector{length(P)-1}(PT(p) for p in eachrow(pxyz))

    commands = Vector{PathCommand{PT}}(undef, N)
    commands[1] = MoveTo(P[1])
    for i in 2:(N-1)
        commands[i] = CurveTo(WP[i-1],
                              2*P[i] - WP[i],
                              P[i])
    end
    commands[N] = CurveTo(WP[N-1],
                          (P[N] + WP[N-1])/2,
                          P[N])
    BezierPath(commands)
end

function cubic_spline(p)
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
