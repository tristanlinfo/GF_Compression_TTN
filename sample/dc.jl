using Random
using Random
using Plots
gr()

include(joinpath(@__DIR__, "..", "src", "topologies.jl"))
include(joinpath(@__DIR__, "..", "src", "utils.jl"))
using .TTNUtils: bits2decimal, fused_preparations, compress_indexset
using TreeTCI

function sampled_error(f, ttn, nsamples, bits, d)
    """ Compute sampled errors between function f and ttn approximation over nsamples random inputs of length 2*bits."""
    eval_ttn = if ttn isa TreeTCI.SimpleTCI
        sitetensors = TreeTCI.fillsitetensors(ttn, f)
        TreeTCI.TreeTensorNetwork(ttn.g, sitetensors)
    else
        ttn
    end
    error_l1 = 0.0
    for _ in 1:nsamples
        # Generate a random 3R sequence of 1s and 2s
        x = rand(1:2, d * bits)
        # Evaluate the concrete TreeTensorNetwork (it provides evaluate/call)
        approx = TreeTCI.evaluate(eval_ttn, x)
        err = abs(f(x) - approx)
        error_l1 += err
    end
    return error_l1 / nsamples
end

function dcgf(v)
    """ Fictious Green's function undergoing a DC field in 2D. """
    x = bits2decimal(v[1:div(length(v), 2)])
    y = bits2decimal(v[div(length(v), 2)+1:end])
    A = 10.0 # amplitude of the dc field
    B = 10.0 # damping factor
    return 1im * exp(-1im * A * (sin(10x) + cos(10x))) * exp(1im * A * (sin(10y) + cos(10y))) * exp(-B * abs(x - y))
end

function plotdcgf()
    A = 1.0
    B = 5.0
    f(x, y) = 1im * exp(-1im * A * (sin(10x) + cos(10x))) *
              exp(1im * A * (sin(10y) + cos(10y))) *
              exp(-B * abs(x - y))

    xs = ys = range(0, 1, length=400)
    zs = [f(x, y) for x in xs, y in ys]
    realzs = real.(zs)
    imagzs = imag.(zs)

    p1 = heatmap(xs, ys, realzs, title="Real part", xlabel="x", ylabel="y")
    p2 = heatmap(xs, ys, imagzs, title="Imag part", xlabel="x", ylabel="y")

    plt = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    display(plt)
    savefig(plt, "SVG/GF_DC/dcgf_A=$(Int(A))_B=$(Int(B)).svg")

    return plt
end


function main()
    # Parameters for TCI
    R = 16 # number of bits per dimension
    d = 2 # spatial dimension
    localdims = fill(2, d * R) # local dimensions for d dimensions with R bits each

    groups = [[i, R + i] for i in 1:R]
    fused_g, fused_localdims, f_fused = fused_preparations(dcgf, groups, localdims)
    # generate 10 random global pivot of length dR
    global_pivots = Vector{Vector{Int}}()
    fused_global_pivots = Vector{Vector{Int}}()
    for _ in 1:10
        push!(global_pivots, rand(1:2, length(localdims)))
        push!(fused_global_pivots, compress_indexset(rand(1:2, length(localdims)), groups, localdims))
    end
    topo = Dict(
        "BTTN" => BTTN(R, d),
        "QTT_Seq" => QTT_Block(R, d),
        "QTT_Int" => QTT_Alt(R, d),
        "CTTN" => CTTN(R, d),
        "Fused_BTTN" => fused_g,
        "TTTN" => TTTN(R, d)
    )

    ntopos = length(topo)
    nsteps = 6
    step = 10
    maxit = 5
    nsamples = 1000

    # Output storage
    error_l1 = zeros(ntopos * 2, nsteps)
    error_pivot = zeros(ntopos * 2, nsteps)
    rankendlist = zeros(ntopos * 2, nsteps)
    ranklist = zeros(ntopos * 2, nsteps)



    #-----------------------------
    # Main loop over maxbonddim
    #-----------------------------
    for i in 1:nsteps
        maxbd = step * i
        println("Max bond dimension: $maxbd")
        tstart = @elapsed begin

            # -----------------------------
            # Construct initial TCIs
            # -----------------------------
            for (j, (toponame, topology)) in enumerate(topo)
                println("Topology: $toponame")
                for k in 1:2

                    if toponame == "Fused_BTTN"
                        # Build Simple TCIs so each run starts from identical initial pivots
                        ttn = TreeTCI.SimpleTCI{ComplexF64}(f_fused, fused_localdims, fused_g)
                        if k == 2
                            TreeTCI.addglobalpivots!(ttn, fused_global_pivots)
                        end
                        # Optimize TCI
                        ranks, errors, bonds = TreeTCI.optimize!(ttn, f_fused; tolerance=1e-16, maxiter=maxit, maxbonddim=maxbd, sweepstrategy=TreeTCI.LocalAdjacentSweep2sitePathProposer())
                        errl1_fused = 0.0
                        sitetensors = TreeTCI.fillsitetensors(ttn, f_fused)
                        eval_ttn = TreeTCI.TreeTensorNetwork(ttn.g, sitetensors)
                        for k in 1:nsamples
                            orig = rand(1:2, 2R)
                            fused_idx = compress_indexset(orig, groups, localdims)
                            approx = TreeTCI.evaluate(eval_ttn, fused_idx)
                            # compare the fused-network approximation to the original function evaluated
                            # on the full indexset `orig` (not to `f_fused(orig)` — that expects fused indices)
                            errl1_fused += abs(dcgf(orig) - approx)
                        end
                        if k == 1
                            jidx = j
                        else
                            jidx = j + ntopos
                        end
                        error_l1[jidx, i] = errl1_fused / nsamples
                        error_pivot[jidx, i] = errors[end]
                        rankendlist[jidx, i] = ranks[end]
                        ranklist[jidx, i] = maxbd
                    else
                        ttn = TreeTCI.SimpleTCI{ComplexF64}(dcgf, localdims, topology)
                        if k == 2
                            TreeTCI.addglobalpivots!(ttn, global_pivots)
                        end
                        # Optimize TCI
                        ranks, errors, bonds = TreeTCI.optimize!(ttn, dcgf; tolerance=1e-16, maxiter=maxit, maxbonddim=maxbd, sweepstrategy=TreeTCI.LocalAdjacentSweep2sitePathProposer())
                        # Compute sampled L1 error
                        errl1 = sampled_error(dcgf, ttn, nsamples, R, d)
                        # Store results

                        if k == 1
                            jidx = j
                        else
                            jidx = j + ntopos
                        end
                        error_l1[jidx, i] = errl1
                        error_pivot[jidx, i] = errors[end]
                        rankendlist[jidx, i] = ranks[end]
                        ranklist[jidx, i] = maxbd
                    end
                end
            end
        end # tstart elapsed
        println("Elapsed time for step $i: $tstart seconds")
    end

    topo_names = collect(keys(topo))

    p1 = plot(title="Sampled L1 Error (DC GF)",
        xlabel="Max Bond Dimension",
        ylabel="Sampled L1 Error",
        yscale=:log10,
        size=(1200, 675),
        legend=:outertopright)
    x = [i * step for i in 1:nsteps]
    colours = ["blue", "red", "green", "orange", "purple", "brown", "pink", "gray"]
    for j in 1:ntopos
        plot!(p1, x, error_l1[j, :], label=topo_names[j], color=colours[j], marker=:o, lw=2, markersize=6)
        plot!(p1, x, error_l1[j+ntopos, :], label=topo_names[j] * "_with_globalpivots", color=colours[j], marker=:star, lw=2, markersize=6)
    end

    savefig(p1, "SVG/GF_DC/GF_DC_A=10_B=10_piv.svg")
    display(p1)
end