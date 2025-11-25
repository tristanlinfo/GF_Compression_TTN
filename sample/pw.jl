using Random
using Plots
using Logging
include(joinpath(@__DIR__, "..", "src", "functions.jl"))
import .SampleTreeTCI: pw, sampled_error
include(joinpath(@__DIR__, "..", "src", "topologies.jl"))
using TreeTCI

function main()
    # Parameters for TCI
    R = 16 # number of bits per dimension
    d = 3 # spatial dimension
    localdims = fill(2, d * R) # local dimensions for d dimensions with R bits each

    topo = Dict(
        "BTTN" => BTTN(R, d),
        "QTT_Seq" => QTT_Block(R, d),
        "QTT_Int" => QTT_Alt(R, d),
        "CTTN" => CTTN(R, d)
    )

    ntopos = length(topo)
    nsteps = 20
    step = 3
    maxit = 5

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
                # Build Simple TCIs so each run starts from identical initial pivots
                ttn = TreeTCI.SimpleTCI{Float64}(pw, localdims, topology)
                #TreeTCI.addglobalpivots!(ttn,global_pivots)
                # Optimize TCI
                ranks, errors = TreeTCI.optimize!(ttn, pw; maxiter=maxit, maxbonddim=maxbd, sweepstrategy=TreeTCI.LocalAdjacentSweep2sitePathProposer())

                # Compute sampled L1 error
                errl1 = sampled_error(pw, ttn, 1000, R, d)


                # Store results
                error_l1[j, i] = errl1
                error_pivot[j, i] = errors[end]
                rankendlist[j, i] = ranks[end]
                ranklist[j, i] = maxbd
            end

        end # tstart elapsed
        println("Elapsed time for step $step: $tstart seconds")
    end

    p1 = plot(xlabel="Max bond dimension", ylabel="Sampled L1 Error", yscale=:log10)
    topo_names = collect(keys(topo))
    for j in 1:ntopos
        plot!(p1, 5 * collect(1:nsteps), error_l1[j, :], label=topo_names[j], marker=:o)
    end
    display(p1)
    savefig(p1, "pw_tci_sampled_l1_error.png")
end
