# NEWS

## Breaking Changes

### Version 0.6.0 (Makie 0.24 Compatibility)

**Breaking Changes:**
- **Removed `edge_plottype` parameter**: The `edge_plottype` argument has been removed from `graphplot`. Edge plotting is now handled automatically with improved performance - curved edges (BezierPaths) are now rendered as efficiently as straight edges using a unified plotting system.

**New Features:**
- **Added `edge_linestyle` parameter**: New parameter to control line styles for individual edges. Supports vectors and dictionaries for per-edge customization. 
  *Performance note*: The new unified edge system maintains optimal performance for both straight and curved edges when using homogeneous linestyles. Only when using different linestyles for different edges (inhomogeneous) does the system fall back to creating separate line plots for each edge, which may reduce performance for graphs with many edges.

**Internal Changes:**
- Refactored reactive system to use `map!` instead of `@lift` macros for better performance and Makie 0.24 compatibility
- Unified edge plotting system replaces separate `linesegments` and `beziersegments` approaches
- Improved interaction system with updated selection handling
