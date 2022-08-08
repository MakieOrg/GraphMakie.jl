using Graphs
using CairoMakie
using GraphMakie
using ReferenceTests
using Literate
using FileIO

const ASSETS = joinpath(@__DIR__, "..", "assets")
const EXAMPLE_BASEPATH = joinpath(@__DIR__, "..", "docs", "examples")

const TMPDIR = joinpath(ASSETS, "tmp")
isdir(TMPDIR) && rm(TMPDIR; recursive=true)
mkdir(TMPDIR)

const IMAGE_COUNTERS = Dict{String, Int}()

"""
    @save_reference fig

Mark figs in example files as reference test image. Those images will
be saved as `filename-i.png` to the gobal temp dir defined in this script.
"""
macro save_reference(fig)
    f = splitpath(string(__source__.file))[end]
    postfix = if haskey(IMAGE_COUNTERS, f)
        IMAGE_COUNTERS[f] += 1
    else
        IMAGE_COUNTERS[f] = 1
    end
    quote
        path = joinpath(TMPDIR, $f*"-"*lpad($postfix, 2, "0")*".png")
        save(path, $(esc(fig)))
        println("   saved fig $path")
    end
end

@info "Generate reference images..."
# go through all examples in docs/examples
for exfile in filter(contains(r".jl$"), readdir(EXAMPLE_BASEPATH))
    expath = joinpath(EXAMPLE_BASEPATH, exfile)

    # check if `@save_reference` appears in the example
    # if not, don't evaluate file
    hasassets = false
    for l in eachline(expath)
        if contains(l, r"@save_reference")
            hasassets = true
            break
        end
    end

    if !hasassets
        @info "$exfile has no reference images!"
        continue
    end

    # create script, include script (saves files), remove exported script
    Literate.script(expath, TMPDIR)
    script = joinpath(TMPDIR, exfile)
    include(script)
    rm(script)
end

oldassets = filter(contains(r".png$"), readdir(ASSETS))
newassets = filter(contains(r".png$"), readdir(TMPDIR))

function compare(ref, x)
    if size(ref) != size(x)
        return 0
    end
    return ReferenceTests._psnr(ref, x)
end

# now test all the generated graphics in the TMPDIR and compare against files in assets dir
@testset "Reference Tests" begin
    for ass in oldassets
        # skip unresolved conflicts
        occursin(r"\+.png$", ass) && continue

        old = joinpath(ASSETS, ass)
        new = joinpath(TMPDIR, ass)

        if !isfile(new)
            @warn "New version for $ass missing! Delete file if not needed anymore."
            @test false
            continue
        end

        # equal = ReferenceTests.psnr_equality()(load(old), load(new))
        score = compare(load(old), load(new))
        MEH = 100
        GOOD = 200

        if score > GOOD
            printstyled(" ✓ [", repr(round(score, digits=1)), "] $ass\n"; color=:green)
            @test true
            rm(new)
        else
            if score > MEH
                printstyled(" ? [", repr(round(score, digits=1)), "] $ass\n"; color=:yellow)
                @test_broken false
            else
                printstyled(" × [", repr(round(score, digits=1)), "] $ass\n"; color=:red)
                @test false
            end
            parts = rsplit(ass, "."; limit=2)
            @assert length(parts) == 2
            newname = parts[1] * "+." *parts[2]
            mv(new, joinpath(ASSETS, newname), force=true)
            @warn "There is a difference in $(ass)! New version moved to $newname. Resolve manually!"
        end
    end

    for new in setdiff(newassets, oldassets)
        printstyled(" × Move new asset $(new)!\n"; color=:red)
        @test false
        mv(joinpath(TMPDIR, new), joinpath(ASSETS, new))
    end

    rm(TMPDIR)
end
