


module StochasticSpindleSim

using DrWatson
using Dates
using Simulate
using Presets


function stochasticSpindleSim(Notes,NumGenerators,NumStates,finalTime,burnTime,maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)

    # cd(workingFolder)

    # folderNotFound = 0
    # folderCounter = 0
    # #folderName = "$workingFolder/output/GillespieRunData_$Notes$folderCounter"
    # while folderNotFound==0
    #     if isdir(datadir("sims",Notes,"GillespieRunData_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))"))
    #         folderCounter += 1
    #     else
    #         folderName = datadir("sims",Notes,"GillespieRunData_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))")
    #         folderNotFound = 1
    #     end
    # end
    folderName = datadir("sims",Notes,"GillespieRunData_$(Dates.format(Dates.now(),"yy-mm-dd-HH-MM-SS"))")
    mkpath(folderName)

    ExtList = LinRange(0,maxExt,NumStates) # full list  of available extensions
    dExt = ExtList[2]-ExtList[1] # dL
    v = 1 .- ExtList # velocity term, need a ±z_t term eventually

    p = (folderName, NumGenerators, NumStates, burnTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on)

    ## saving parameters in a text file
    save_params = open("$folderName/run_parameters.txt","w")
    write(save_params, "NumGenerators NumStates maxExt finalTime α β Γ dExt w_on w_0   \n")
    write(save_params, string(NumGenerators, ", ", NumStates,", ", maxExt,", ", finalTime,", ", α,", ", β,", ", Γ,", ", dExt, ", ", ω_on, ", ", ω_0))
    close(save_params)

    initialStates = preset(p)

    q = (folderName, NumGenerators, NumStates, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on)
    # println("Time to simulate :)")

    simulate(Notes, q, initialStates)

    return 0

end

export stochasticSpindleSim

end
