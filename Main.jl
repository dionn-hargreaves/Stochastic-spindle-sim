using Simulate
using Presets

workingFolder = "/Volumes/DH_Simulations/"
cd(workingFolder)


Notes = "N30_Alpha_Beta_Larger_bigger_onrate_smallerK_newfiles"


folderNotFound = 0
folderCounter = 0
folderName = "$workingFolder/output/GillespieRunData_$Notes$folderCounter"
while folderNotFound==0
    if isdir("$workingFolder/output/GillespieRunData_$Notes$folderCounter")
        global folderCounter += 1
    else
        folderName = "$workingFolder/output/GillespieRunData_$Notes$folderCounter"
        folderNotFound = 1
    end
end



# set up seed for random numbers, set parameters
NumGenerators = 30 # per cortex
NumStates = 6000 # number of states for extension length
finalTime = 700000

maxExt = 6 # maximum extension available
ExtList = LinRange(0,maxExt,NumStates) # full list  of available extensions
α = 0.0771604938271605
β = 0.0412615740740741
Γ = 1/50
dExt = ExtList[2]-ExtList[1] # dL
v = 1 .- ExtList # velocity term, need a ±z_t term eventually
γ = 2.0
z = 0
μ = 50
K = 5e-3
ω_0 = 0.001
ω_on = 0.03

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
