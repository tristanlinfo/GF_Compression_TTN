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
function CTTN(R::Int, d::Int=2)
    @assert R ≥ 2 "R must be at least 2"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    # build a linear chain for each dimension
    for dim in 0:(d-1)
        base = dim * R
        for i in 1:(R-1)
            add_edge!(g, base + i, base + i + 1)
        end
    end

    # link the root of the first chain to the roots of the other chains
    for dim in 1:(d-1)
        add_edge!(g, 1, dim * R + 1)
    end

    return g
end

"""
    qtt_int_2()

Return a NamedGraph for the specific 32-vertex TTN topology shown in `ttn32`.
This topology is a fixed 32-vertex graph (indices 1..32) useful for comparisons.
"""
function qtt_int_2()
    g = NamedGraph(32)

    # edges (from ttn32 topology)
    add_edge!(g, 1, 17)
    add_edge!(g, 2, 18)
    add_edge!(g, 3, 18)
    add_edge!(g, 3, 19)
    add_edge!(g, 3, 20)
    add_edge!(g, 4, 20)
    add_edge!(g, 5, 20)
    add_edge!(g, 5, 21)
    add_edge!(g, 5, 22)
    add_edge!(g, 6, 22)
    add_edge!(g, 7, 23)
    add_edge!(g, 8, 23)
    add_edge!(g, 8, 24)
    add_edge!(g, 8, 25)
    add_edge!(g, 9, 25)
    add_edge!(g, 10, 25)
    add_edge!(g, 10, 26)
    add_edge!(g, 10, 27)
    add_edge!(g, 11, 27)
    add_edge!(g, 12, 13)
    add_edge!(g, 12, 27)
    add_edge!(g, 12, 28)
    add_edge!(g, 13, 14)
    add_edge!(g, 13, 29)
    add_edge!(g, 14, 30)
    add_edge!(g, 14, 31)
    add_edge!(g, 15, 31)
    add_edge!(g, 16, 31)
    add_edge!(g, 16, 32)
    add_edge!(g, 17, 19)
    add_edge!(g, 22, 23)

    return g
end

"""
    qtt_int_3()

Return a NamedGraph for the specific 32-vertex TTN topology shown in `ttn32_b`.
This topology is a fixed 32-vertex graph (indices 1..32) useful for comparisons.
"""
function qtt_int_3()
    g = NamedGraph(32)

    # edges (from ttn32_b topology)
    add_edge!(g, 1, 17)
    add_edge!(g, 2, 18)
    add_edge!(g, 3, 17)
    add_edge!(g, 3, 18)
    add_edge!(g, 3, 19)
    add_edge!(g, 3, 20)
    add_edge!(g, 4, 20)
    add_edge!(g, 5, 20)
    add_edge!(g, 5, 21)
    add_edge!(g, 6, 22)
    add_edge!(g, 7, 22)
    add_edge!(g, 8, 23)
    add_edge!(g, 9, 25)
    add_edge!(g, 10, 25)
    add_edge!(g, 11, 26)
    add_edge!(g, 12, 28)
    add_edge!(g, 13, 28)
    add_edge!(g, 14, 29)
    add_edge!(g, 15, 31)
    add_edge!(g, 16, 31)
    add_edge!(g, 20, 22)
    add_edge!(g, 22, 23)
    add_edge!(g, 23, 24)
    add_edge!(g, 23, 25)
    add_edge!(g, 25, 26)
    add_edge!(g, 26, 27)
    add_edge!(g, 26, 28)
    add_edge!(g, 28, 29)
    add_edge!(g, 29, 30)
    add_edge!(g, 29, 31)
    add_edge!(g, 31, 32)

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
function CTTN_Alt(R::Int, d::Int=2)
    @assert R ≥ 2 "R must be at least 2"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    # vertices for dimension `dim` are at positions (dim-1)*R + 1 .. dim*R
    # connect within each dimension skipping according to interleaved layout
    for dim in 1:d
        for j in 1:(R-1)
            v1 = (dim - 1) * R + j
            v2 = (dim - 1) * R + j + 1
            add_edge!(g, v1, v2)
        end
    end

    # Link the first roots together (star from first root)
    for dim in 2:d
        add_edge!(g, 1, (dim - 1) * R + 1)
    end
    return g
end


"""
    BTTN2(R::Int, d::Int=2)

Construct a NamedGraph like the binary TTN shown in the attached screenshot:
two balanced binary trees (one per dimension) on the vertex blocks 1..R and (R+1)..(2R),
with the roots linked similarly to `BTTN`.

This implementation mirrors `BTTN` but is provided as a separate constructor
so you can call `BTTN2(R,d)` when you specifically want the screenshot-style
topology.
"""
function BTTN2(R::Int, d::Int=2)
    @assert R ≥ 2 "R must be at least 2"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    # For each dimension build a balanced binary tree on the block 1..R with offset
    for dim in 0:(d-1)
        base = dim * R
        for i in 1:R
            left = 2i
            right = 2i + 1
            if left ≤ R
                add_edge!(g, base + i, base + left)
            end
            if right ≤ R
                add_edge!(g, base + i, base + right)
            end
        end
    end

    # connect the roots of all dimensions to the first root (1)
    for dim in 1:(d-1)
        add_edge!(g, 1, dim * R + 1)
    end

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
function BTTN(R::Int, d::Int=2)
    @assert R ≥ 2 "R must be at least 2"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    # For each dimension build a binary tree on the block 1..R with offset
    for dim in 0:(d-1)
        base = dim * R
        for i in 1:R
            left = 2i
            right = 2i + 1
            if left ≤ R
                add_edge!(g, base + i, base + left)
            end
            if right ≤ R
                add_edge!(g, base + i, base + right)
            end
        end
    end

    # connect the roots of all dimensions to the first root (1)
    for dim in 1:(d-1)
        add_edge!(g, 1, dim * R + 1)
    end

    return g
end

function TTTN(R::Int, d::Int=2)
    @assert R ≥ 3 "R must be at least 3"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    for dim in 0:(d-1)
        base = dim * R
        # parents are base+1, base+3, base+5, ... (step 2)
        for p in (base+1):2:(base+R)
            if p + 1 <= base + R
                add_edge!(g, p, p + 1)
            end
            if p + 2 <= base + R
                add_edge!(g, p, p + 2)
            end
        end
    end

    # connect the roots of all dimensions to the first root (1)
    for dim in 1:(d-1)
        add_edge!(g, 1, dim * R + 1)
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
function BTTN_Alt(R::Int, d::Int=2)
    @assert R ≥ 2 "R must be at least 2"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)

    # For each dimension, treat positions (dim-1)*R + 1 .. dim*R as the binary tree nodes
    for dim in 0:(d-1)
        base = dim * R
        for i in 1:R
            parent = base + i
            left_idx = 2 * i
            right_idx = 2 * i + 1
            if left_idx <= R
                left = base + left_idx
                add_edge!(g, parent, left)
            end
            if right_idx <= R
                right = base + right_idx
                add_edge!(g, parent, right)
            end
        end
    end

    # Link roots together
    for dim in 1:(d-1)
        add_edge!(g, 1, dim * R + 1)
    end
    return g
end

"""
    QTT_Block(R::Int)

Generate a NamedGraph representing a tensor train topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).
All x coordinates first, then all y coordinates: x₁ — x₂ — ... — x_R — y₁ — y₂ — ... — y_R
"""
function QTT_Block(R::Int, d::Int=2)
    @assert R ≥ 1 "R must be at least 1"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)
    for i in 1:(N-1)
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

function QTT_Block_Alt(R::Int, d::Int=2)
    @assert R ≥ 1 "R must be at least 1"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)
    # interleaved ordering: for j in 1:R, dims 1..d
    seq = [(dim - 1) * R + j for j in 1:R for dim in 1:d]
    for i in 1:(length(seq)-1)
        add_edge!(g, seq[i], seq[i+1])
    end
    # optionally add a cross-link similar to the 2D version
    if d == 2 && R >= 2
        add_edge!(g, 2R - 1, 2)
    end
    return g
end

"""
    QTT_Alt(R::Int)

Generate a NamedGraph representing a tensor train topology for 2D functions,
with `R` bits precision per dimension (total 2R vertices).
Vertices alternate between x and y coordinates: x₁ — y₁ — x₂ — y₂ — ... — x_R — y_R
"""
function QTT_Alt(R::Int, d::Int=2)
    @assert R ≥ 1 "R must be at least 1"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)
    # connect dims across each position j: dim1_j - dim2_j - ... - dimd_j
    for j in 1:R
        for dim in 1:(d-1)
            add_edge!(g, (dim - 1) * R + j, dim * R + j)
        end
        # connect last dim at j to first dim at j+1 to create a chain along j
        if j < R
            add_edge!(g, (d - 1) * R + j, 1 + j)
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

function QTT_Alt_Alt(R::Int, d::Int=2)
    @assert R ≥ 1 "R must be at least 1"
    @assert d ≥ 1 "dimension d must be at least 1"
    N = d * R
    g = NamedGraph(N)
    for i in 1:(N-1)
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
