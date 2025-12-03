using Random
using Random
using Plots
gr()

include(joinpath(@__DIR__, "..", "src", "topologies.jl"))
include(joinpath(@__DIR__, "..", "src", "utils.jl"))
using .TTNUtils: bits2decimal, fused_preparations, compress_indexset
using TreeTCI
using Integrals

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

function integrand(t)
    # Intx = \int_0^t dz cos(-A * sin(omega*z)/omega) + sin(-A * sin(omega*z)/omega)
    A = 50.0 # amplitude of the ac field
    ω = 5.0  # frequency of the ac field
    return cos(-A * sin(ω * t) / ω) + sin(-A * sin(ω * t) / ω)
end

function acgf(v)
    """ Fictious Green's function undergoing a ac field in 2D. """
    x = bits2decimal(v[1:div(length(v), 2)])
    y = bits2decimal(v[div(length(v), 2)+1:end])
    B = 1.0 # damping factor 

    Intx = Integrals.quadgk(integrand, 0, x)[1]
    Inty = Integrals.quadgk(integrand, 0, y)[1]
    return 1im * exp(-1im * Intx) * exp(1im * Inty) * exp(-B * abs(x - y))
end

function plotacgf()
    B = 5.0
    f(x, y) = 1im * exp(-1im * Integrals.quadgk(integrand, 0, x)[1]) *
              exp(1im * Integrals.quadgk(integrand, 0, y)[1]) *
              exp(-B * abs(x - y))

    xs = ys = range(0, 1, length=400)
    zs = [f(x, y) for x in xs, y in ys]
    realzs = real.(zs)
    imagzs = imag.(zs)

    p1 = heatmap(xs, ys, realzs, title="Real part", xlabel="x", ylabel="y")
    p2 = heatmap(xs, ys, imagzs, title="Imag part", xlabel="x", ylabel="y")

    plt = plot(p1, p2, layout=(1, 2), size=(1200, 500))
    display(plt)
    return plt
end

function main()
    # Parameters for TCI
    R = 16 # number of bits per dimension
    d = 2 # spatial dimension
    localdims = fill(2, d * R) # local dimensions for d dimensions with R bits each

    groups = [[i, R + i] for i in 1:R]
    fused_g, fused_localdims, f_fused = fused_preparations(acgf, groups, localdims)
    topo = Dict(
        "BTTN" => BTTN(R, d),
        "QTT_Seq" => QTT_Block(R, d),
        "QTT_Int" => QTT_Alt(R, d),
        "CTTN" => CTTN(R, d),
        "Fused_BTTN" => fused_g,
        "TTTN" => TTTN(R, d)
    )

    ntopos = length(topo)
    nsteps = 10
    step = 3
    maxit = 5
    nsamples = 1000

    # Output storage
    error_l1 = zeros(ntopos, nsteps)
    error_pivot = zeros(ntopos, nsteps)
    rankendlist = zeros(ntopos, nsteps)
    ranklist = zeros(ntopos, nsteps)

    # generate 10 random global pivot of length 3R
    global_pivots = Vector{Vector{Int}}()
    for _ in 1:10
        push!(global_pivots, rand(1:2, length(localdims)))
    end

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
                if toponame == "Fused_BTTN"
                    # Build Simple TCIs so each run starts from identical initial pivots
                    ttn = TreeTCI.SimpleTCI{ComplexF64}(f_fused, fused_localdims, fused_g)
                    #TreeTCI.addglobalpivots!(ttn,global_pivots)
                    # Optimize TCI
                    ranks, errors = TreeTCI.optimize!(ttn, f_fused; tolerance=1e-16, maxiter=maxit, maxbonddim=maxbd, sweepstrategy=TreeTCI.LocalAdjacentSweep2sitePathProposer())
                    errl1_fused = 0.0
                    sitetensors = TreeTCI.fillsitetensors(ttn, f_fused)
                    eval_ttn = TreeTCI.TreeTensorNetwork(ttn.g, sitetensors)
                    for k in 1:nsamples
                        orig = rand(1:2, 2R)
                        fused_idx = compress_indexset(orig, groups, localdims)
                        approx = TreeTCI.evaluate(eval_ttn, fused_idx)
                        # compare the fused-network approximation to the original function evaluated
                        # on the full indexset `orig` (not to `f_fused(orig)` — that expects fused indices)
                        errl1_fused += abs(acgf(orig) - approx)
                    end
                    error_l1[j, i] = errl1_fused / nsamples
                    error_pivot[j, i] = errors[end]
                    rankendlist[j, i] = ranks[end]
                    ranklist[j, i] = maxbd
                else
                    ttn = TreeTCI.SimpleTCI{ComplexF64}(acgf, localdims, topology)
                    #TreeTCI.addglobalpivots!(ttn,global_pivots)
                    # Optimize TCI
                    ranks, errors = TreeTCI.optimize!(ttn, acgf; tolerance=1e-16, maxiter=maxit, maxbonddim=maxbd, sweepstrategy=TreeTCI.LocalAdjacentSweep2sitePathProposer())
                    # Compute sampled L1 error
                    errl1 = sampled_error(acgf, ttn, nsamples, R, d)
                    # Store results
                    error_l1[j, i] = errl1
                    error_pivot[j, i] = errors[end]
                    rankendlist[j, i] = ranks[end]
                    ranklist[j, i] = maxbd
                end
            end
        end # tstart elapsed
        println("Elapsed time for step $step: $tstart seconds")
    end

    p1 = plot(title="Sampled L1 Error vs Max Bond Dimension", xlabel="Max Bond Dimension", ylabel="Sampled L1 Error", yscale=:log10)
    topo_names = collect(keys(topo))
    x = [i * step for i in 1:nsteps]
    for j in 1:ntopos
        plot!(p1, x, error_l1[j, :], label=topo_names[j], marker=:o)
    end

    savefig(p1, "SVG/GF_AC/GF_AC_A=50_B=1.svg")
    display(p1)
end