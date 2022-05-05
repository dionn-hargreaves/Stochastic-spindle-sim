
#push!(LOAD_PATH,"src")

using StochasticSpindleSim
if contains(pwd(),"/Users/")
    include("scripts/TestParameters.jl")
    stochasticSpindleSim(workingFolder,Notes,NumGenerators,NumStates,finalTime, burnTime, maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)
else
#    import Pkg
#    Pkg.add("CircularArrayBuffers")
    include("scripts/TestParameters.jl")

    stochasticSpindleSim(remoteFolder,Notes,NumGenerators,NumStates,finalTime, burnTime, maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)
end
