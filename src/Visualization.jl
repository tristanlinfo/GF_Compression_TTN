# --- Alternated BTTN topology visualization ---
function draw_BTTN_Alt_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 2 "R must be at least 2"
    n = 2R
    # Odd vertices: x-tree, Even vertices: y-tree
    xs = zeros(n)
    ys = zeros(n)
    for i in 1:R
        xs[2i-1] = i
        ys[2i-1] = 2.0
        xs[2i] = i
        ys[2i] = 1.0
    end
    # x-tree edges
    edge_x = []
    for i in 1:R
        parent = 2i - 1
        left_idx = 2 * i
        right_idx = 2 * i + 1
        if left_idx <= R
            left = 2 * left_idx - 1
            push!(edge_x, (parent, left))
        end
        if right_idx <= R
            right = 2 * right_idx - 1
            push!(edge_x, (parent, right))
        end
    end
    # y-tree edges
    edge_y = []
    for i in 1:R
        parent = 2i
        left_idx = 2 * i
        right_idx = 2 * i + 1
        if left_idx <= R
            left = 2 * left_idx
            push!(edge_y, (parent, left))
        end
        if right_idx <= R
            right = 2 * right_idx
            push!(edge_y, (parent, right))
        end
    end
    # Link
    edge_link = [(1, 2)]
    edges = vcat(edge_x, edge_y, edge_link)
    # Check bond_dims length
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, R + 1), ylim=(0.5, 2.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    bond_idx = 1
    for (src, dst) in edge_x
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    for (src, dst) in edge_y
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    src, dst = edge_link[1]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

# --- Alternated CTTN topology visualization ---
function draw_CTTN_Alt_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 2 "R must be at least 2"
    n = 2R
    xs = zeros(n)
    ys = zeros(n)
    for i in 1:R
        xs[2i-1] = i
        ys[2i-1] = 2.0
        xs[2i] = i
        ys[2i] = 1.0
    end
    edge_x = [(2i - 1, 2i + 1) for i in 1:(R-1)]
    edge_y = [(2i, 2i + 2) for i in 1:(R-1)]
    edge_link = [(1, 2)]
    edges = vcat(edge_x, edge_y, edge_link)
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, R + 1), ylim=(0.5, 2.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    bond_idx = 1
    for (src, dst) in edge_x
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    for (src, dst) in edge_y
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    src, dst = edge_link[1]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

# --- Alternated QTT_Block topology visualization ---
function draw_QTT_Block_Alt_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 1 "R must be at least 1"
    n = 2R
    xs = zeros(n)
    ys = zeros(n)
    for i in 1:R
        xs[2i-1] = i
        ys[2i-1] = 2.0
        xs[2i] = i
        ys[2i] = 1.0
    end
    edge_x = [(2i - 1, 2i + 1) for i in 1:(R-1)]
    edge_y = [(2i, 2i + 2) for i in 1:(R-1)]
    edge_link = [(2R - 1, 2)]
    edges = vcat(edge_x, edge_y, edge_link)
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, R + 1), ylim=(0.5, 2.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    bond_idx = 1
    for (src, dst) in edge_x
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    for (src, dst) in edge_y
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    src, dst = edge_link[1]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

# --- Alternated QTT_Alt topology visualization ---
function draw_QTT_Alt_Alt_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 1 "R must be at least 1"
    n = 2R
    xs = collect(1:n)
    ys = fill(1.5, n)
    edges = [(i, i + 1) for i in 1:(n-1)]
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, n + 1), ylim=(1, 2), aspect_ratio=2, xaxis=false, yaxis=false, grid=false)
    bond_idx = 1
    for (src, dst) in edges
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end
using Plots

"""
    draw_BTTN_topology(R::Int)

Draws the BTTN topology for a given R (bits per dimension) without constructing a graph object.
Visualizes:
  - x-binary tree: nodes 1 to R
  - y-binary tree: nodes R+1 to 2R
  - Link: edge from node 1 to node R+1
"""
function draw_BTTN_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 2 "R must be at least 2"
    n = 2R
    # Helper to compute tree layout positions
    function tree_layout(N, x0, y0, width, height)
        xs = zeros(N)
        ys = zeros(N)
        function place(i, x, y, w, h)
            xs[i] = x
            ys[i] = y
            left = 2i
            right = 2i + 1
            if left <= N
                place(left, x - w / 2, y - h, w / 2, h)
            end
            if right <= N
                place(right, x + w / 2, y - h, w / 2, h)
            end
        end
        place(1, x0, y0, width, height)
        return xs, ys
    end

    # x-tree layout (bottom, grows upward)
    xs_x, ys_x = tree_layout(R, 0.0, 1.0, 1.5, +0.5)
    # y-tree layout (top, grows downward)
    xs_y, ys_y = tree_layout(R, 0.0, 2.0, 1.5, -0.5)

    xs = vcat(xs_x, xs_y)
    ys = vcat(ys_x, ys_y)

    # Edges: x-binary tree
    edge_x = []
    for i in 1:R
        left = 2i
        right = 2i + 1
        if left ≤ R
            push!(edge_x, (i, left))
        end
        if right ≤ R
            push!(edge_x, (i, right))
        end
    end
    # Edges: y-binary tree
    edge_y = []
    for i in (R+1):n
        left = R + 2 * (i - R)
        right = R + 2 * (i - R) + 1
        if left ≤ n
            push!(edge_y, (i, left))
        end
        if right ≤ n
            push!(edge_y, (i, right))
        end
    end
    # Link
    edge_link = [(1, R + 1)]
    edges = vcat(edge_x, edge_y, edge_link)

    # Check bond_dims length
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end

    plt = plot(; legend=false, xlim=(-R, R), ylim=(-0.5, 3.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    # Map bond_dims: x-tree edges, link, y-tree edges
    bond_idx = 1
    # x-tree edges
    for (src, dst) in edge_x
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    # x-root/y-root link
    src, dst = edge_link[1]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
    end
    bond_idx += 1
    # y-tree edges
    for (src, dst) in edge_y
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

"""
    draw_QTT_Block_topology(R::Int)

Draws the QTT_Block topology for a given R (bits per dimension) without constructing a graph object.
Visualizes:
  - x-chain: nodes 1 to R
  - y-chain: nodes R+1 to 2R
  - All in a single chain: 1—2—...—R—R+1—...—2R
"""
function draw_QTT_Block_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 1 "R must be at least 1"
    n = 2R
    xs = collect(1:n)
    ys = fill(1.5, n)
    edges = [(i, i + 1) for i in 1:(n-1)]

    # Check bond_dims length
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, n + 1), ylim=(1, 2), aspect_ratio=2, xaxis=false, yaxis=false, grid=false)
    # Map bond_dims: x-chain edges, link, y-chain edges
    bond_idx = 1
    # x-chain
    for (src, dst) in edges[1:R-1]
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    # x-root/y-root link
    src, dst = edges[R]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[R]), :red, 8))
    end
    bond_idx += 1
    # y-chain
    for i in (R+1):(length(edges))
        src, dst = edges[i]
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

"""
    draw_QTT_Alt_topology(R::Int)

Draws the QTT_Alt topology for a given R (bits per dimension) without constructing a graph object.
Visualizes:
  - Alternating x and y nodes: x₁—y₁—x₂—y₂—...—x_R—y_R
  - Edges: i—R+i, R+i—i+1 (for i < R)
"""
function draw_QTT_Alt_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 1 "R must be at least 1"
    n = 2R
    xs = zeros(n)
    ys = zeros(n)
    for i in 1:R
        xs[i] = 2i - 1
        ys[i] = 2.0
        xs[R+i] = 2i
        ys[R+i] = 1.0
    end
    edges = []
    for i in 1:R
        push!(edges, (i, R + i))
        if i < R
            push!(edges, (R + i, i + 1))
        end
    end

    # Check bond_dims length
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end
    plt = plot(; legend=false, xlim=(0, 2R + 1), ylim=(0.5, 2.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    # Map bond_dims: x-y edges, link, y-x edges
    bond_idx = 1
    # x-y edges
    for i in 1:R
        src, dst = edges[bond_idx]
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    # y-x edges
    for i in bond_idx:length(edges)
        src, dst = edges[i]
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[i]), :red, 8))
        end
    end
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    for i in 1:n
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end

"""
	draw_CTTN_topology(R::Int)

Draws the CTTN topology for a given R (bits per dimension) without constructing a graph object.
Visualizes:
  - x-chain: nodes 1 to R in a horizontal line
  - y-chain: nodes R+1 to 2R in a parallel line below
  - Link: edge from node 1 to node R+1
"""
function draw_CTTN_topology(R::Int, bond_dims::Vector{<:Any}=nothing, savepath::Union{Nothing,String}=nothing)
    @assert R ≥ 2 "R must be at least 2"
    # Node positions
    x_x = collect(1:R)
    y_x = fill(2.0, R)
    x_y = collect(1:R)
    y_y = fill(1.0, R)

    # Combine positions
    xs = vcat(x_x, x_y)
    ys = vcat(y_x, y_y)

    # Edges: x-chain
    edge_x = [(i, i + 1) for i in 1:(R-1)]
    # Edges: y-chain
    edge_y = [(R + i, R + i + 1) for i in 1:(R-1)]
    # Link
    edge_link = [(1, R + 1)]
    edges = vcat(edge_x, edge_y, edge_link)

    # Check bond_dims length
    if bond_dims !== nothing
        if length(bond_dims) != length(edges)
            error("bond_dims length $(length(bond_dims)) does not match number of edges $(length(edges))")
        end
    end

    # Plot
    plt = plot(; legend=false, xlim=(0, R + 1), ylim=(0.5, 2.5), aspect_ratio=1, xaxis=false, yaxis=false, grid=false)
    # Draw edges
    # Map bond_dims: x-chain edges, link, y-chain edges
    bond_idx = 1
    # x-chain
    for (src, dst) in edge_x
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    # x-root/y-root link
    src, dst = edge_link[1]
    plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
    if bond_dims !== nothing
        mx = (xs[src] + xs[dst]) / 2
        my = (ys[src] + ys[dst]) / 2 + 0.12
        annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
    end
    bond_idx += 1
    # y-chain
    for (src, dst) in edge_y
        plot!([xs[src], xs[dst]], [ys[src], ys[dst]], lw=2, color=:black)
        if bond_dims !== nothing
            mx = (xs[src] + xs[dst]) / 2
            my = (ys[src] + ys[dst]) / 2 + 0.12
            annotate!(mx, my, text(string(bond_dims[bond_idx]), :red, 8))
        end
        bond_idx += 1
    end
    # Draw nodes
    scatter!(xs, ys, marker=:circle, ms=8, color=:blue)
    # Label nodes
    for i in 1:(2R)
        annotate!(xs[i], ys[i], text(string(i), :white, 5))
    end
    if savepath !== nothing
        savefig(plt, savepath)
    else
        display(plt)
    end
end