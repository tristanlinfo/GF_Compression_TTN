module optiTTN

using Random
using Logging
using TreeTCI
using NamedGraphs: NamedGraph, add_edge!, vertices, NamedEdge, has_edge

include(joinpath(@__DIR__, "utils.jl"))
import .TTNUtils: bits2decimal

include(joinpath(@__DIR__, "treetci_src.jl"))
import .SampleTreeTCI: make_gf, set_gf_params, set_k, l1_error, run_experiment_gf, export
optimize_ttn_gf

"""
    optimize_ttn_gf(R::Int, nstep::Int, maxbonddim_step::Int, nsamples::Int,
                    topology::NamedGraph, n::Int, eps::Vector{Float64},
                    lambda::Vector{Float64})
Optimize a TTN representation of the fictitious Green's function
with given parameters and return the optimized TTN, final ranks, and errors.
"""




end