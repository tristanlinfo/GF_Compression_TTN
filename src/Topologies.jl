using Graphs
using NamedGraphs

##############################################
# Classical topologies
##############################################

"""
    CTTN(R::Int)

Generate a NamedGraph representing a quantics Comb TTN topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).

Structure:
  x-chain: 1 — 2 — 3 — ... — R
  y-chain: (R+1) — (R+2) — ... — (2R)
  Link: 1 — (R+1)
"""
function CTTN(R::Int)
    @assert R ≥ 2 "R must be at least 2"
    g = NamedGraph(2R)

    # x chain
    for i in 1:(R-1)
        add_edge!(g, i, i + 1)
    end

    # x–y link
    add_edge!(g, 1, R + 1)

    # y chain
    for i in (R+1):(2R-1)
        add_edge!(g, i, i + 1)
    end

    return g

end

"""
    CTTN_Alt(R::Int)

Generate a NamedGraph representing a quantics Comb TTN topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices), using alternated ordering:
    x₁ — y₁ — x₂ — y₂ — ... — x_R — y_R
Structure:
    x-chain: 1 — 3 — 5 — ... — (2R-1)
    y-chain: 2 — 4 — 6 — ... — (2R)
    Link: 1 — 2
"""
function CTTN_Alt(R::Int)
    @assert R ≥ 2 "R must be at least 2"
    g = NamedGraph(2R)
    # x-chain (odd vertices)
    for i in 1:(R-1)
        add_edge!(g, 2i - 1, 2i + 1)
    end
    # y-chain (even vertices)
    for i in 1:(R-1)
        add_edge!(g, 2i, 2i + 2)
    end
    # x–y link
    add_edge!(g, 1, 2)
    return g
end

"""
    BTTN(R::Int)

Generate a NamedGraph representing a quantics Binary TTN topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).

Structure:
  - x-representation: binary tree of depth ≈ log2(R)
  - y-representation: same structure shifted by R
  - Link between roots: 1 — (R+1)
"""
function BTTN(R::Int)
    @assert R ≥ 2 "R must be at least 2"
    g = NamedGraph(2R)

    # ---- X binary tree ----
    for i in 1:R
        left = 2i
        right = 2i + 1
        if left ≤ R
            add_edge!(g, i, left)
        end
        if right ≤ R
            add_edge!(g, i, right)
        end
    end

    # Link between x-root and y-root
    add_edge!(g, 1, R + 1)

    # ---- Y binary tree ----
    for i in (R+1):(2R)
        left = R + 2 * (i - R)
        right = R + 2 * (i - R) + 1
        if left ≤ 2R
            add_edge!(g, i, left)
        end
        if right ≤ 2R
            add_edge!(g, i, right)
        end
    end

    return g
end

"""
    BTTN_3D(R::Int)

Generate a NamedGraph representing a quantics Binary TTN topology for 3D functions,
with `R` bits precision per dimension (total 3R vertices).

Structure:
  - x-representation: binary tree on vertices 1..R
  - y-representation: binary tree on vertices (R+1)..(2R)
  - z-representation: binary tree on vertices (2R+1)..(3R)
  - Links between roots: connect x-root (1) to y-root (R+1) and z-root (2R+1)
"""
function BTTN_3D(R::Int)
    @assert R ≥ 2 "R must be at least 2"
    g = NamedGraph(3R)

    # ---- X binary tree ----
    for i in 1:R
        left = 2i
        right = 2i + 1
        if left ≤ R
            add_edge!(g, i, left)
        end
        if right ≤ R
            add_edge!(g, i, right)
        end
    end

    # Links between roots: x-root connects to y-root and z-root
    add_edge!(g, 1, R + 1)
    add_edge!(g, 1, 2R + 1)

    # ---- Y binary tree ----
    for i in (R+1):(2R)
        left = R + 2 * (i - R)
        right = R + 2 * (i - R) + 1
        if left ≤ 2R
            add_edge!(g, i, left)
        end
        if right ≤ 2R
            add_edge!(g, i, right)
        end
    end

    # ---- Z binary tree ----
    for i in (2R+1):(3R)
        left = 2R + 2 * (i - 2R)
        right = 2R + 2 * (i - 2R) + 1
        if left ≤ 3R
            add_edge!(g, i, left)
        end
        if right ≤ 3R
            add_edge!(g, i, right)
        end
    end

    return g
end


"""
    BTTN_Alt(R::Int)

Generate a NamedGraph representing a quantics Binary TTN topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices), using alternated ordering:
    x₁ — y₁ — x₂ — y₂ — ... — x_R — y_R
Structure:
    - x-representation: binary tree on odd vertices
    - y-representation: binary tree on even vertices
    - Link between roots: 1 — 2
"""
function BTTN_Alt(R::Int)
    @assert R ≥ 2 "R must be at least 2"
    g = NamedGraph(2R)
    # X binary tree (odd vertices: 1, 3, 5, ...)
    for i in 1:R
        parent = 2i - 1
        left_idx = 2 * i
        right_idx = 2 * i + 1
        if left_idx <= R
            left = 2 * left_idx - 1
            add_edge!(g, parent, left)
        end
        if right_idx <= R
            right = 2 * right_idx - 1
            add_edge!(g, parent, right)
        end
    end
    # Y binary tree (even vertices: 2, 4, 6, ...)
    for i in 1:R
        parent = 2i
        left_idx = 2 * i
        right_idx = 2 * i + 1
        if left_idx <= R
            left = 2 * left_idx
            add_edge!(g, parent, left)
        end
        if right_idx <= R
            right = 2 * right_idx
            add_edge!(g, parent, right)
        end
    end
    # Link between x-root and y-root
    add_edge!(g, 1, 2)
    return g
end

"""
    QTT_Block(R::Int)

Generate a NamedGraph representing a tensor train topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).
All x coordinates first, then all y coordinates: x₁ — x₂ — ... — x_R — y₁ — y₂ — ... — y_R
"""
function QTT_Block(R::Int)
    @assert R ≥ 1 "R must be at least 1"
    g = NamedGraph(2R)
    for i in 1:(2R-1)
        add_edge!(g, i, i + 1)
    end
    return g
end

"""
    QTT_Block_Alt(R::Int)

    Generate a NamedGraph representing a tensor train topology for 2D functions,
    with `R` bits precision per dimension (total 2R vertices), using alternated ordering:
    Structure:
    x-chain: 1 — 3 — 5 — ... — (2R-1)
    y-chain: 2 — 4 — 6 — ... — (2R)
    Link: 2R-1 — 2
"""

function QTT_Block_Alt(R::Int)
    @assert R ≥ 1 "R must be at least 1"
    g = NamedGraph(2R)
    # x-chain (odd vertices)
    for i in 1:(R-1)
        add_edge!(g, 2i - 1, 2i + 1)
    end
    # y-chain (even vertices)
    for i in 1:(R-1)
        add_edge!(g, 2i, 2i + 2)
    end
    # link between x and y chains
    add_edge!(g, 2R - 1, 2)
    return g
end

"""
    QTT_Alt(R::Int)

Generate a NamedGraph representing a tensor train topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).
Vertices alternate between x and y coordinates: x₁ — y₁ — x₂ — y₂ — ... — x_R — y_R
"""
function QTT_Alt(R::Int)
    @assert R ≥ 1 "R must be at least 1"
    g = NamedGraph(2R)
    for i in 1:R
        add_edge!(g, i, R + i)
        if i < R
            add_edge!(g, R + i, i + 1)
        end
    end
    return g
end

"""
    QTT_Alt_Alt(R::Int)

Generate a NamedGraph representing a tensor train topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).
Vertices alternate between x and y coordinates: x₁ — y₁ — x₂ — y₂ — ... — x_R — y_R
Structure:
    x1-y1-x2-y2-...-xR-yR chain: 1 — 2 — 3 — 4 — ... — (2R-1) — (2R)
"""

function QTT_Alt_Alt(R::Int)
    @assert R ≥ 1 "R must be at least 1"
    g = NamedGraph(2R)
    for i in 1:(2R-1)
        add_edge!(g, i, i + 1)
    end
    return g
end

##############################################
# Improved topology between QTT and BTTN
##############################################

"""
    GF_topo(R::Int,r::Int)

Generate an intermidiate structure between Interleaved QTT and BTTN for 2D functions,
with 'R' bits precision per dimension (total 2R vertices). 'r' stands for the number of
bits that is used for the QTT represenation (for each dimension, so 2r in total).
The structure start with 'r' bits representation for the Interleaved QTT and then 'R-r'
bits for the BTTN.
Structure :
    x1-y1-x2-y2-...-xr-yr, then a BTTN starting from xr/yr for the x/y branch.

"""
function GF_topo(R::Int, r::Int)
    @assert R ≥ 1 "R must be at least 1"
    @assert 1 ≤ r ≤ R "r must satisfy 1 ≤ r ≤ R (number of QTT bits per dimension)"
    g = NamedGraph(2R)
    for i in 1:(2r-1)
        add_edge!(g, i, i + 1)
    end

    # ---- X binary tree ----
    N = R - r + 1  # number of odd vertices in the x-subtree (2r-1, 2r+1, ..., 2R-1)
    for m in 1:N
        parent = 2r - 1 + 2 * (m - 1)
        left_m = 2 * m
        right_m = 2 * m + 1
        if left_m <= N
            left = 2r - 1 + 2 * (left_m - 1)
            add_edge!(g, parent, left)
        end
        if right_m <= N
            right = 2r - 1 + 2 * (right_m - 1)
            add_edge!(g, parent, right)
        end
    end

    # ---- Y binary tree ----
    N = R - r + 1  # number of even vertices in the y-subtree (2r, 2r+2, ..., 2R)
    for m in 1:N
        parent = 2r + 2 * (m - 1)
        left_m = 2 * m
        right_m = 2 * m + 1
        if left_m <= N
            left = 2r + 2 * (left_m - 1)
            add_edge!(g, parent, left)
        end
        if right_m <= N
            right = 2r + 2 * (right_m - 1)
            add_edge!(g, parent, right)
        end
    end
    return g
end


##############################################
# Fused topology helpers
##############################################

"""
    fuse_graph_pairs(g::NamedGraph, groups::Vector{Vector{Int}})

Given an input graph `g` with vertices labeled 1..N and a partition `groups`
(vector of index vectors), return a new `NamedGraph` with one vertex per group
and an edge between group a and group b if any vertex in a is connected to any
vertex in b in `g`.
"""
function fuse_graph_pairs(g::NamedGraph, groups::Vector{Vector{Int}})
    m = length(groups)
    fg = NamedGraph(m)
    # For each pair of groups check if there's any edge between their members
    for a in 1:m
        for b in (a+1):m
            found = false
            for u in groups[a]
                for v in groups[b]
                    if has_edge(g, u, v) || has_edge(g, v, u)
                        add_edge!(fg, a, b)
                        found = true
                        break
                    end
                end
                if found
                    break
                end
            end
        end
    end
    return fg
end

"""
    fused_localdims(localdims::Vector{Int}, groups::Vector{Vector{Int}})

Compute the local dimensions for each fused group as the product of the
constituent local dims.
"""
function fused_localdims(localdims::Vector{Int}, groups::Vector{Vector{Int}})
    return [prod(localdims[idx] for idx in group) for group in groups]
end

"""
    compress_indexset(orig::Vector{Int}, groups::Vector{Vector{Int}}, localdims::Vector{Int})

Convert an original indexset (one integer per original site) to an indexset
for the fused groups. The mapping uses a mixed-radix ordering where the first
site in the group is the least significant digit.
"""
function compress_indexset(orig::Vector{Int}, groups::Vector{Vector{Int}}, localdims::Vector{Int})
    fused = Int[]
    for group in groups
        idx = 0
        radix = 1
        for site in group
            i = orig[site] - 1 # zero-based
            idx += i * radix
            radix *= localdims[site]
        end
        push!(fused, idx + 1)
    end
    return fused
end

"""
    expand_fused_indexset(fused::Vector{Int}, groups::Vector{Vector{Int}}, localdims::Vector{Int})

Inverse of `compress_indexset`: expands fused indices back to original per-site
indices.
"""
function expand_fused_indexset(fused::Vector{Int}, groups::Vector{Vector{Int}}, localdims::Vector{Int})
    N = sum(length(g) for g in groups)
    orig = Vector{Int}(undef, N)
    for (gi, group) in enumerate(groups)
        val = fused[gi] - 1
        for site in group
            d = localdims[site]
            orig[site] = (val % d) + 1
            val ÷= d
        end
    end
    return orig
end

"""
    make_fused_wrapper(f_orig, groups, localdims)

Return a function `f_fused` that accepts fused indexsets (one integer per
group) and calls `f_orig` with the expanded original indexset.
"""
function make_fused_wrapper(f_orig, groups::Vector{Vector{Int}}, localdims::Vector{Int})
    return function (fused_idx)
        orig = expand_fused_indexset(fused_idx, groups, localdims)
        return f_orig(orig)
    end
end