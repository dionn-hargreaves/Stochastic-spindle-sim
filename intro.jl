using DrWatson
@quickactivate "Stochastic-spindle-sim"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
""")

@info "Loading test parameters"
include("scripts/TestParameters.jl")

@info "Precompiling project"

push!(LOAD_PATH,"src")
using Revise
using BenchmarkTools
using StochasticSpindleSim
