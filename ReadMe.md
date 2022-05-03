# Stochastic-spindle-sim

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Stochastic-spindle-sim

It is authored by dionn-hargreaves.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.



Make sure `src` is in your `LOAD_PATH`
Start REPL
Run `includet("scripts/TestParameters.jl")`
Run `using StochasticSpindleSim`
Run `stochasticSpindleSim(workingFolder,Notes,NumGenerators,NumStates,finalTime, burnTime, maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)`
