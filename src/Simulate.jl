#
# Simulate.jl
# Stochastic model
#
# Created by Dionn Hargreaves 02/07/2021

module Simulate

# import Julia packages
using Base.Threads
using Plots
using Random
using DelimitedFiles
using LaTeXStrings
using Statistics: mean
using CircularArrayBuffers
# import local modules
using GillespieTransitions


@inline function simulate(Notes, p, initialStates)

    # import parameters
    folderName, NumGenerators, NumStates, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on = p
    initZ, MastParams, inittPassed, GenList = initialStates
    #open files
    upperFile = open("$folderName/$Notes+_upper.txt","a")
    lowerFile = open("$folderName/$Notes+_lower.txt","a")
    poleFile = open("$folderName/$Notes+_pole.txt","a")
    timeFile = open("$folderName/$Notes+_timestamp.txt","a")

    # initialise arrays
    #GenList = zeros(2,NumGenerators*2)
    tPassed = CircularArrayBuffer{Float64}(2) #zeros(finalTime*10) # array to hold time passed
    z = CircularArrayBuffer{Float64}(2) #zeros(finalTime*10)

    push!(z, initZ)
    push!(tPassed, inittPassed)

    # set forward, backward and switching parameters (upper and lower cortex)
    UpparamB = zeros(3,NumStates)
    UpparamB[3, :] .= ω_0*exp.(γ.*ExtList)
    DownparamB = zeros(3,NumStates)
    DownparamB[3, :] .= ω_0*exp.(γ.*ExtList)

    UpparamU = zeros(3,NumStates)
    UpparamU[3, :] .= ω_on; # binding
    DownparamU = zeros(3,NumStates)
    DownparamU[3, :] .= ω_on; # binding


    noBound_Up = CircularArrayBuffer{Int64}(2) #   zeros(finalTime*10)
    noBound_Down = CircularArrayBuffer{Int64}(2) #zeros(finalTime*10)
    avgYBound_Up = CircularArrayBuffer{Float64}(2) #zeros(finalTime*10)
    avgYUnbound_Up = CircularArrayBuffer{Float64}(2)#zeros(finalTime*10)
    avgYBound_Down = CircularArrayBuffer{Float64}(2)#zeros(finalTime*10)
    avgYUnbound_Down = CircularArrayBuffer{Float64}(2)#zeros(finalTime*10)

    j = 1 # time counter
    FileCount = 1 # file counter
    while j != finalTime*10
        j+=1 # tick time
        genInd, chState, newtPassed = gillespieTran!(MastParams, tPassed[end])
        ## returns genInd, inidex of changed generator, chState, how generator is affected, and new sim time
        push!(tPassed, newtPassed)
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
        end




        ## update parameters based on new system
        BoundUp = findall(x->x>0, GenList[2,1:NumGenerators])
        BoundDown = findall(x->x>0, GenList[2,NumGenerators+1:end]).+NumGenerators
        DzDt = (1/0.625).*( -K*z[end]-(sum(ExtList[GenList[1,BoundDown]])-sum(ExtList[GenList[1,BoundUp]]))) # new spindle velocity
        newZ = z[end]+(tPassed[end]-tPassed[1])*DzDt # new spindle position, forward Euler
        push!(z, newZ)

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

        push!(noBound_Up, length(findall(x -> x > 0, GenList[2, 1:NumGenerators])))
        push!(noBound_Down, length(findall(x -> x > 0, GenList[2, (1+NumGenerators):2*NumGenerators])))
        push!(avgYBound_Up, mean(ExtList[convert(Array{Int64,1},GenList[1, findall(x -> x > 0, GenList[2, 1:NumGenerators])])]))
        push!(avgYBound_Down, mean(ExtList[convert(Array{Int64,1},GenList[1, findall(x -> x > 0, GenList[2, (1+NumGenerators):2*NumGenerators])])]))
        push!(avgYUnbound_Up, mean(ExtList[convert(Array{Int64,1},GenList[1, findall(x -> x < 0, GenList[2, 1:NumGenerators])])]))
        push!(avgYUnbound_Down, mean(ExtList[convert(Array{Int64,1},GenList[1, findall(x -> x < 0, GenList[2, (1+NumGenerators):2*NumGenerators])])]))



        if mod(j,1000) == 0
            #sBStateX = BStateX[:,(j+1-100000):j]
            #sHold_index = Hold_index[:,(j+1-100000):j]
            writedlm(upperFile, [ noBound_Up[end] avgYBound_Up[end] avgYUnbound_Up[end]])
            writedlm(lowerFile, [ noBound_Down[end] avgYBound_Down[end] avgYUnbound_Down[end]])
            writedlm(poleFile, [ z[end]])
            writedlm(timeFile, [ tPassed[end]])
            if mod(j, 10000) == 0
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
