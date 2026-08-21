using Test
using GraphMakie
using GraphMakie.Graphs
using GraphMakie.NetworkLayout
using GraphMakie: compute_auto_label_aligns
using GeometryBasics

@testset "compute_auto_label_aligns" begin
    @testset "isolated nodes" begin
        # Single isolated node
        g = SimpleDiGraph(1)
        node_positions = [Point2f(0.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 1
        @test aligns[1] == (:right, :bottom)  # Default for isolated nodes
        
        # Multiple isolated nodes
        g = SimpleDiGraph(3)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0), Point2f(0.5, 1.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 3
        @test all(a -> a == (:right, :bottom), aligns)
    end
    
    @testset "single edge" begin
        # Chain: A → B
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 2
        # Node 1 (source): edge goes right (0°), label should be opposite (left or left-bottom)
        # Node 2 (destination): edge comes from left (180°), label should be opposite (right or right-bottom)
        # The exact alignment depends on the gap calculation, but should avoid the edge direction
        @test aligns[1] isa Tuple{Symbol, Symbol}
        @test aligns[2] isa Tuple{Symbol, Symbol}
        # Labels should not be in the edge direction
        @test aligns[1] != (:left, :center)  # Not in edge direction (right)
        @test aligns[2] != (:right, :center)  # Not in edge direction (left)
    end
    
    @testset "chain pattern" begin
        # Chain: A → B → C
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0), Point2f(2.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 3
        @test all(a -> a isa Tuple{Symbol, Symbol}, aligns)
        # Node 2 (middle) has edges from both sides, label should be above or below
        @test aligns[2] in [(:center, :bottom), (:center, :top)]
    end
    
    @testset "fork pattern" begin
        # Fork: A ← B → C
        g = SimpleDiGraph(3)
        add_edge!(g, 2, 1)  # B → A
        add_edge!(g, 2, 3)  # B → C
        node_positions = [Point2f(-1.0, 0.0), Point2f(0.0, 0.0), Point2f(1.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 3
        # Node 2 (center) has edges going left and right, label should be above or below
        @test aligns[2] in [(:center, :bottom), (:center, :top)]
    end
    
    @testset "collider pattern" begin
        # Collider: A → B ← C
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)  # A → B
        add_edge!(g, 3, 2)  # C → B
        node_positions = [Point2f(-1.0, 0.0), Point2f(0.0, 0.0), Point2f(1.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 3
        # Node 2 (collider) has edges coming from left and right, label should be above or below
        @test aligns[2] in [(:center, :bottom), (:center, :top)]
    end
    
    @testset "dense graph" begin
        # Complete graph with 4 nodes
        g = complete_graph(4)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0), Point2f(0.5, 1.0), Point2f(0.5, -1.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 4
        @test all(a -> a isa Tuple{Symbol, Symbol}, aligns)
        # All nodes should have valid alignments (not crash)
    end
    
    @testset "undirected graph" begin
        # Undirected chain
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0), Point2f(2.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 3
        @test all(a -> a isa Tuple{Symbol, Symbol}, aligns)
    end
    
    @testset "all edges in one direction" begin
        # Star graph: center node with all edges outgoing
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 1, 4)
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0), Point2f(0.0, 1.0), Point2f(-1.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 4
        # Node 1 (center) should have label opposite to where edges go
        @test aligns[1] isa Tuple{Symbol, Symbol}
    end
    
    @testset "return type consistency" begin
        g = SimpleDiGraph(5)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        add_edge!(g, 4, 5)
        node_positions = [Point2f(i, 0.0) for i in 0:4]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test aligns isa Vector{Tuple{Symbol, Symbol}}
        @test length(aligns) == 5
        @test all(a -> a isa Tuple{Symbol, Symbol} && length(a) == 2, aligns)
    end
end

@testset "nlabels_auto_align integration" begin
    @testset "basic usage" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # Test that nlabels_auto_align works without errors
        fig, ax, p = graphplot(g, nlabels = ["A", "B", "C"], nlabels_auto_align = true)
        @test p isa GraphMakie.GraphPlot
        
        # Test with auto_align = false (default behavior)
        fig, ax, p = graphplot(g, nlabels = ["A", "B", "C"], nlabels_auto_align = false)
        @test p isa GraphMakie.GraphPlot
    end
    
    @testset "isolated nodes use fallback" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        # Node 3 is isolated
        
        # With custom fallback alignment
        fig, ax, p = graphplot(g, 
            nlabels = ["A", "B", "C"], 
            nlabels_auto_align = true,
            nlabels_align = (:left, :top)  # Fallback for isolated nodes
        )
        @test p isa GraphMakie.GraphPlot
        # Isolated node should use fallback, nodes with edges should use auto-align
    end
    
    @testset "works with different layouts" begin
        g = complete_graph(5)
        layouts = [Spring(), Stress(), Shell()]
        
        for layout in layouts
            fig, ax, p = graphplot(g, 
                nlabels = ["A", "B", "C", "D", "E"],
                nlabels_auto_align = true,
                layout = layout
            )
            @test p isa GraphMakie.GraphPlot
        end
    end
    
    @testset "works with undirected graphs" begin
        g = SimpleGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 3, 4)
        
        fig, ax, p = graphplot(g, 
            nlabels = ["A", "B", "C", "D"],
            nlabels_auto_align = true
        )
        @test p isa GraphMakie.GraphPlot
    end
    
    @testset "backward compatibility" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # Default behavior (nlabels_auto_align = false) should work as before
        fig1, ax1, p1 = graphplot(g, nlabels = ["A", "B", "C"])
        fig2, ax2, p2 = graphplot(g, nlabels = ["A", "B", "C"], nlabels_auto_align = false)
        
        @test p1 isa GraphMakie.GraphPlot
        @test p2 isa GraphMakie.GraphPlot
    end
    
    @testset "works without labels" begin
        g = SimpleDiGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        
        # Should work fine even when nlabels is nothing
        fig, ax, p = graphplot(g, nlabels_auto_align = true)
        @test p isa GraphMakie.GraphPlot
    end
end

@testset "edge cases and robustness" begin
    @testset "empty graph" begin
        g = SimpleDiGraph(0)
        node_positions = Point2f[]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 0
    end
    
    @testset "single node graph" begin
        g = SimpleDiGraph(1)
        node_positions = [Point2f(0.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 1
        @test aligns[1] == (:right, :bottom)
    end
    
    @testset "self-loops" begin
        g = SimpleDiGraph(2)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 1)  # Self-loop
        node_positions = [Point2f(0.0, 0.0), Point2f(1.0, 0.0)]
        aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 2
        @test all(a -> a isa Tuple{Symbol, Symbol}, aligns)
    end
    
    @testset "large graph performance" begin
        # Test that it doesn't crash on larger graphs
        g = SimpleDiGraph(50)
        for i in 1:49
            add_edge!(g, i, i+1)
        end
        node_positions = [Point2f(i, 0.0) for i in 1:50]
        
        @time aligns = compute_auto_label_aligns(g, node_positions)
        @test length(aligns) == 50
        @test all(a -> a isa Tuple{Symbol, Symbol}, aligns)
    end
end
