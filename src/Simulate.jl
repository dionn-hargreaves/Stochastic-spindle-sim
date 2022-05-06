#
# Simulate.jl
# Stochastic model
#
# Created by Dionn Hargreaves 02/07/2021

module Simulate

# import Julia packages
# using Base.Threads
using Random
using DelimitedFiles
using Statistics: mean
using CircularArrayBuffers
using FastBroadcast
# import local modules
using GillespieTransitions


@views function simulate(Notes, p, )

    # import parameters
    folderName, NumGenerators, NumStates, burnTime, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on = p

# **************** Setup arrays

    # GenList[i,1] gives the extention of the ith generator; GenList[i,2] is its on/off binding state
    GenList = zeros(NumGenerators*2,2)
    MastParams = zeros((NumGenerators*3)*2)

    Random.seed!(1234) # random number seed
    # parameters: bound generators
    paramB = zeros(3,NumStates)
    paramB[1,2:NumStates] .= α/(dExt^2) .- v[2:NumStates]./(2*dExt) # backward
    paramB[2,1:NumStates-1] .= α/(dExt^2) .+ v[1:NumStates-1]./(2*dExt) #forward
    paramB[3, :] .= ω_0*exp.(γ.*ExtList) #unbind

    # set forward, backward and switching parameters (upper and lower cortex)
    UpparamB = zeros(NumStates,3)
    UpparamB[:,3] .= ω_0*exp.(γ.*ExtList)
    DownparamB = zeros(NumStates,3)
    DownparamB[:,3] .= ω_0*exp.(γ.*ExtList)

    # parameters: unbound generators
    paramU = zeros(3,NumStates)
    paramU[1,2:NumStates] .=  Γ.*(β/(dExt^2) .+ ExtList[2:NumStates]./(2*dExt)) # backward
    paramU[2,1:NumStates-1] .= Γ.*(β/(dExt^2) .- ExtList[1:NumStates-1]./(2*dExt)) # forward
    paramU[3, :] .= ω_on # bind

    UpparamU = zeros(NumStates,3)
    UpparamU[:,3] .= ω_on # binding
    DownparamU = zeros(NumStates,3)
    DownparamU[:,3] .= ω_on # binding

    # Generate list of generators, holds extension position and bind state
    GenList[:,1] = rand(1:NumStates,NumGenerators*2) # randomoly choose a beginning state for all generators
    rand_index =  rand(1:NumGenerators) # choose which generators are bound
    GenList[1:rand_index,end] .= 1 ## set x0[end] = 1 for bound, = -1 for unbound
    GenList[NumGenerators+1:NumGenerators+rand_index,end] .= 1
    GenList[rand_index+1:NumGenerators,end] .= -1
    GenList[NumGenerators+rand_index+1:end,end] .= -1

    #OVERRIDE: ALL GENERATORS BEGIN BOUND
    GenList[:,end] .= 1
    GenList = convert(Array{Int64, 2}, GenList) # so that we can use values as indices later

    for i in 1:NumGenerators*2 # fill master list of all possible state change probabilities for each generator
        if GenList[i,2] == 1
            MastParams[3*i-2:3*i] .= paramB[:, GenList[i,1]]
        else
            MastParams[3*i-2:3*i] .= paramU[:, GenList[i,1]]
        end
    end

# ****************

    #open files
    upperFile = open("$folderName/$Notes+_upper.txt","a")
    lowerFile = open("$folderName/$Notes+_lower.txt","a")
    poleFile = open("$folderName/$Notes+_pole.txt","a")
    timeFile = open("$folderName/$Notes+_timestamp.txt","a")

    tPassed = 0.0
    newtPassed = 0.0
    z = 0.0
    noBoundVec = zeros(Int64,2)# [noBound_Up, noBound_Down]
    avgYVec = zeros(Float64,4) # [avgYBound_Up, avgYUnbound_Up, avgYBound_Down, avgYUnbound_Down]

    stateIndVec = zeros(Int64,2) # [chState, genInd]

    j = 1 # time counter
    FileCount = 1 # file counter
    while j != finalTime*10

        j+=1 # tick time

        newtPassed = gillespieTran!(stateIndVec, MastParams, tPassed)
        ## returns stateIndVec[2], inidex of changed generator, stateIndVec[1], how generator is affected, and new sim time

        if stateIndVec[1] == 1 # retraction
            GenList[stateIndVec[2],1] -= 1
        elseif stateIndVec[1] == 2 # extension
            GenList[stateIndVec[2],1] += 1
        else # bind/unbind
            GenList[stateIndVec[2],2] = GenList[stateIndVec[2],2]*-1
        end

        if newtPassed<tPassed # backwards in time flag
            println("Backwards in time!")
            println(findall(x->x<0, MastParams))
        end

        ## update parameters based on new system
        BoundUp = findall(x->x>0, GenList[1:NumGenerators,2])
        BoundDown = findall(x->x>0, GenList[NumGenerators+1:end,2]).+NumGenerators
        DzDt = (1/0.625).*( -K*z[end]-(sum(ExtList[GenList[BoundDown,1]])-sum(ExtList[GenList[BoundUp,1]]))) # new spindle velocity
        z = z + (newtPassed-tPassed)*DzDt # new spindle position, forward Euler
        tPassed = newtPassed

        upV = 1.0 .- ExtList .- DzDt        # new v+ for parameters
        downV = 1.0 .- ExtList .+ DzDt      # new v- for parameters
        @.. thread=true UpparamB[2:NumStates,1] = α/(dExt^2) - upV[2:NumStates]/(2*dExt)          # updating parameter
        @.. thread=true UpparamB[1:NumStates-1,2] = α/(dExt^2) + upV[1:NumStates-1]/(2*dExt)      # updating parameter
        @.. thread=true DownparamB[2:NumStates,1] = α/(dExt^2) - downV[2:NumStates]/(2*dExt)      # updating parameter
        @.. thread=true DownparamB[1:NumStates-1,2] = α/(dExt^2) + downV[1:NumStates-1]/(2*dExt)  # updating parameter


        @.. thread=true UpparamU[2:NumStates,1] = Γ*(β/(dExt^2) + (ExtList[2:NumStates]/(2*dExt))) # backward
        @.. thread=true UpparamU[1:NumStates-1,2] = Γ*(β/(dExt^2) - (ExtList[2:NumStates]/(2*dExt))) # forward
        @.. thread=true DownparamU[2:NumStates,1] = Γ*(β/(dExt^2) + (ExtList[2:NumStates]/(2*dExt))) # backward
        @.. thread=true DownparamU[1:NumStates-1,2] = Γ*(β/(dExt^2) - (ExtList[2:NumStates]/(2*dExt))) # forward

        # update master parameters based on current states
        for i in 1:NumGenerators # upper cortex
            if GenList[i,2] == 1
                MastParams[3*i-2:3*i] .= UpparamB[GenList[i,1],:]
            else
                MastParams[3*i-2:3*i] .= UpparamU[GenList[i,1],:]
            end
        end

        for i in NumGenerators+1:2*NumGenerators # lower cortex
            if GenList[i,2] == 1
                MastParams[3*i-2:3*i] .= DownparamB[GenList[i,1],:]
            else
                MastParams[3*i-2:3*i] .= DownparamU[GenList[i,1],:]
            end
        end

# *** Wrap all of these parameters into an array and just update each component of the array
        noBoundVec[1] = count(x->x>0, GenList[1:NumGenerators,2])
        noBoundVec[2] = count(x->x>0, GenList[(1+NumGenerators):2*NumGenerators,2])
        avgYVec[1]    = mean(ExtList[convert(Array{Int64,1},GenList[findall(x -> x > 0, GenList[1:NumGenerators,2]),1])])
        avgYVec[3]    = mean(ExtList[convert(Array{Int64,1},GenList[findall(x -> x > 0, GenList[(1+NumGenerators):2*NumGenerators,2]).+NumGenerators,1])])
        avgYVec[2]    = mean(ExtList[convert(Array{Int64,1},GenList[findall(x -> x < 0, GenList[1:NumGenerators,2]),1])])
        avgYVec[4]    = mean(ExtList[convert(Array{Int64,1},GenList[findall(x -> x < 0, GenList[(1+NumGenerators):2*NumGenerators,2]).+NumGenerators,1])])
# ***************

        if (mod(j,10000) == 0) && (tPassed > burnTime)
# *** Write a slice of the variables array to file; don't create a new array within the function call
            writedlm(upperFile, [ noBoundVec[1] avgYVec[1] avgYVec[2] ])
            writedlm(lowerFile, [ noBoundVec[2] avgYVec[3] avgYVec[4] ])
            writedlm(poleFile, [ z ])
            writedlm(timeFile, [ tPassed ])
# ***************
            if mod(j, 100000000) == 0
                flush(upperFile)
                flush(lowerFile)
                flush(poleFile)
                flush(timeFile)
            end
        end
    end

    close(upperFile)
    close(lowerFile)
    close(poleFile)
    close(timeFile)
    return 1
end

export simulate

end
