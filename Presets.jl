module Presets

using GillespieTransitions
using Random

@inline function preset(p)
    folderName, NumGenerators, NumStates, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on = p

    GenList = zeros(2,NumGenerators*2)
    MastParams = zeros((NumGenerators*3)*2)

    Random.seed!(1234) # random number seed
    # parameters: bound generators
    paramB = zeros(3,NumStates)
    paramB[1,2:NumStates] .= α/(dExt^2) .- v[2:NumStates]./(2*dExt); # backward
    paramB[2,1:NumStates-1] .= α/(dExt^2) .+ v[1:NumStates-1]./(2*dExt); #forward
    paramB[3, :] .= ω_0*exp.(γ.*ExtList) #unbind

    UpparamB = zeros(3,NumStates)
    UpparamB[3, :] .= ω_0*exp.(γ.*ExtList)
    DownparamB = zeros(3,NumStates)
    DownparamB[3, :] .= ω_0*exp.(γ.*ExtList)


    # parameters: unbound generators
    paramU = zeros(3,NumStates)
    paramU[1,2:NumStates] .=  Γ.*(β/(dExt^2) .+ ExtList[2:NumStates]./(2*dExt)); # backward
    paramU[2,1:NumStates-1] .= Γ.*(β/(dExt^2) .- ExtList[1:NumStates-1]./(2*dExt)); # forward
    paramU[3, :] .= ω_on; # bind

    UpparamU = zeros(3,NumStates)
    UpparamU[3, :] .= ω_on; # binding
    DownparamU = zeros(3,NumStates)
    DownparamU[3, :] .= ω_on; # binding


    # Generate list of generators, holds extension position and bind state
    GenList[1, :] = rand(1:NumStates, NumGenerators*2) # randomoly choose a beginning state for all generators
    rand_index =  rand(1:NumGenerators) # choose which generators are bound
    GenList[end,1:rand_index] .= 1 ## set x0[end] = 1 for bound, = -1 for unbound
    GenList[end,NumGenerators+1:NumGenerators+rand_index] .= 1
    GenList[end,rand_index+1:NumGenerators] .= -1
    GenList[end,NumGenerators+rand_index+1:end] .= -1

    #OVERRIDE: ALL GENERATORS BEGIN BOUND
    GenList[end,:] .= 1

    GenList = convert(Array{Int64, 2}, GenList) # so that we can use values as indices later
    tPassed = zeros(2) # array to hold time passed


    for i in 1:NumGenerators*2 # fill master list of all possible state change probabilities for each generator
        if GenList[2,i] == 1
            MastParams[3*i-2:3*i] .= paramB[:, GenList[1,i]]
        else
            MastParams[3*i-2:3*i] .= paramU[:, GenList[1,i]]
        end
    end



    j = 1 # time counter
    while j <= finalTime
        j+=1 # tick time
        genInd, chState, tPassed[2] = gillespieTran!(MastParams, tPassed[1])
        ## returns genInd, inidex of changed generator, chState, how generator is affected, and new sim time
        if chState == 1.0 # retraction
            GenList[1, genInd] -= 1
        elseif chState == 2.0 # extension
            GenList[1, genInd] += 1
        else # bind/unbind
            GenList[2, genInd] = GenList[2, genInd]*-1
        end

#=
### block to catch negative extensions
        if GenList[1, genInd] < 0
            println("Something wrong with index ")
            println(genInd)
            println("Doing method")
            println(chState)
        end
=#
        if tPassed[2]<tPassed[1] # backwards in time flag
            println("Backwards in time! We've broken the biscuits!")
            println(findall(x->x<0, MastParams))
            #return @save "$folderName/$Notes+_results.bson" BStateX Hold_index z tPassed
        end

        ## update parameters based on new system
        BoundUp = findall(x->x>0, GenList[2,1:NumGenerators])
        BoundDown = findall(x->x>0, GenList[2,NumGenerators+1:end]).+NumGenerators
        DzDt = (1/0.625).*( -K*z[1]-(sum(ExtList[GenList[1,BoundDown]])-sum(ExtList[GenList[1,BoundUp]]))) # new spindle velocity
        z[2] = z[1]+(tPassed[2]-tPassed[1])*DzDt # new spindle position, forward Euler


        upV = 1.0 .- ExtList .+ DzDt # new v+ for parameters
        downV = 1.0 .- ExtList .- DzDt # new v- for parameters
        UpparamB[1,2:NumStates] .= α/(dExt^2) .- upV[2:NumStates]/(2*dExt); # updating parameter
        UpparamB[2,1:NumStates-1] .= α/(dExt^2) .+ upV[1:NumStates-1]/(2*dExt); # updating parameter
        DownparamB[1,2:NumStates] .= α/(dExt^2) .- downV[2:NumStates]/(2*dExt); # updating parameter
        DownparamB[2,1:NumStates-1] .= α/(dExt^2) .+ downV[1:NumStates-1]/(2*dExt); # updating parameter


        UpparamU[1,2:NumStates] .= Γ.*(β/(dExt^2) .+ (ExtList[2:NumStates]./(2*dExt))); # backward
        UpparamU[2,1:NumStates-1] .= Γ.*(β/(dExt^2) .- (ExtList[2:NumStates]./(2*dExt))); # forward
        DownparamU[1,2:NumStates] .= Γ.*(β/(dExt^2) .+ (ExtList[2:NumStates]./(2*dExt))); # backward
        DownparamU[2,1:NumStates-1] .= Γ.*(β/(dExt^2) .- (ExtList[2:NumStates]./(2*dExt))); # forward

        # update master parameters based on current states
        for i in 1:NumGenerators # upper cortex
            if GenList[2,i] == 1
                MastParams[3*i-2:3*i] .= UpparamB[:, GenList[1,i]]
            else
                MastParams[3*i-2:3*i] .= UpparamU[:, GenList[1,i]]
            end
        end

        for i in NumGenerators+1:2*NumGenerators # lower cortex
            if GenList[2,i] == 1
                MastParams[3*i-2:3*i] .= DownparamB[:, GenList[1,i]]
            else
                MastParams[3*i-2:3*i] .= DownparamU[:, GenList[1,i]]
            end
        end

        tPassed[1] = tPassed[2]
        z[1] = z[2]
    end

initialStates = (z[2], MastParams, tPassed[2], GenList)

return initialStates
end

export preset

end
