


module StochasticSpindleSim

using Simulate
using Presets


function stochasticSpindleSim(dataDirectory,Notes,NumGenerators,NumStates,finalTime,maxExt,α,β,Γ,γ,z,μ,K,ω_0,ω_on)

    folderNotFound = 0
    folderCounter = 0
    folderRoot = "$dataDirectory/GillespieRunData_$Notes"
    folderName = ""
    while folderNotFound==0
        if isdir("$folderRoot$folderCounter")
            folderCounter += 1
        else
            folderName = "$folderRoot$folderCounter"
            folderNotFound = 1
        end
    end

    ExtList = LinRange(0,maxExt,NumStates) # full list  of available extensions
    dExt = ExtList[2]-ExtList[1] # dL
    v = 1 .- ExtList # velocity term, need a ±z_t term eventually

    p = (folderName, NumGenerators, NumStates, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on)

    ## saving parameters in a text file
    mkpath(folderName)
    save_params = open("$folderName/run_parameters.txt","w")
    write(save_params, "NumGenerators NumStates maxExt finalTime α β Γ dExt \n")
    write(save_params, string(NumGenerators, ", ", NumStates,", ", maxExt,", ", finalTime,", ", α,", ", β,", ", Γ,", ", dExt))
    close(save_params)

    @time initialStates = preset(p)

    q = (folderName, NumGenerators, NumStates, 100000000, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on)
    println("Time to simulate :)")

    @time simulate(Notes, q, initialStates)

    return 0

end

export stochasticSpindleSim

end
