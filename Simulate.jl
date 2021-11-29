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
using BSON: @save
using LaTeXStrings
using Statistics: mean
# import local modules
using GillespieTransitions

@inline function simulate(Notes, p, initialStates)


    folderName, NumGenerators, NumStates, finalTime, maxExt, ExtList, α, β, Γ, dExt, v, γ, z, μ, K, ω_0, ω_on = p

    GenList = zeros(2,NumGenerators*2)
    tPassed = zeros(finalTime*10) # array to hold time passed
    z = zeros(finalTime*10)

    z[1], MastParams, tPassed[1], GenList = initialStates

    UpparamB = zeros(3,NumStates)
    UpparamB[3, :] .= ω_0*exp.(γ.*ExtList)
    DownparamB = zeros(3,NumStates)
    DownparamB[3, :] .= ω_0*exp.(γ.*ExtList)

    UpparamU = zeros(3,NumStates)
    UpparamU[3, :] .= ω_on; # binding
    DownparamU = zeros(3,NumStates)
    DownparamU[3, :] .= ω_on; # binding


    Hold_index = zeros(NumGenerators*2,finalTime*10) # saves extension states for all time
    #println("here")
    Hold_index[:,1] .= GenList[1,:]
    # create a vector to hold number of bound generators
    noBound_Up = zeros(finalTime*10)
    noBound_Down = zeros(finalTime*10)
    avgYBound_Up = zeros(finalTime*10)
    avgYUnbound_Up = zeros(finalTime*10)
    avgYBound_Down = zeros(finalTime*10)
    avgYUnbound_Down = zeros(finalTime*10)

    BStateX = zeros(NumGenerators*2,finalTime*10) # saves binding states for all time
    BStateX[:,1] .= GenList[2,:]

    j = 1 # time counter
    FileCount = 1 # file counter
    while j != finalTime*10
        j+=1 # tick time
        genInd, chState, tPassed[j] = gillespieTran!(MastParams, tPassed[j-1])
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
        if tPassed[j]<tPassed[j-1] # backwards in time flag
            println("Backwards in time! We've broken the biscuits!")
            println(findall(x->x<0, MastParams))
            return @save "$folderName/$Notes+_results.bson" BStateX Hold_index z tPassed
        end


        Hold_index[:, j] .= GenList[1,:] # update index list
        BStateX[:,j] .= GenList[2,:] # update bind state list

        ## update parameters based on new system
        BoundUp = findall(x->x>0, GenList[2,1:NumGenerators])
        BoundDown = findall(x->x>0, GenList[2,NumGenerators+1:end]).+NumGenerators
        DzDt = (1/0.625).*( -K*z[j-1]-(sum(ExtList[GenList[1,BoundDown]])-sum(ExtList[GenList[1,BoundUp]]))) # new spindle velocity
        z[j] = z[j-1]+(tPassed[j]-tPassed[j-1])*DzDt # new spindle position, forward Euler


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

        noBound_Up[j] = length(findall(x -> x > 0, BStateX[1:NumGenerators,j]))
        noBound_Down[j] = length(findall(x -> x > 0, BStateX[(1+NumGenerators):2*NumGenerators,j]))
        avgYBound_Up[j] = mean(ExtList[convert(Array{Int64,1},Hold_index[findall(x -> x > 0, BStateX[1:NumGenerators,j]),j])])
        avgYBound_Down[j] = mean(ExtList[convert(Array{Int64,1},Hold_index[findall(x -> x > 0, BStateX[(1+NumGenerators):2*NumGenerators,j]),j])])
        avgYUnbound_Up[j] = mean(ExtList[convert(Array{Int64,1},Hold_index[findall(x -> x < 0, BStateX[1:NumGenerators,j]),j])])
        avgYUnbound_Down[j] = mean(ExtList[convert(Array{Int64,1},Hold_index[findall(x -> x < 0, BStateX[(1+NumGenerators):2*NumGenerators,j]),j])])

        if mod(j,100000) == 0
            #sBStateX = BStateX[:,(j+1-100000):j]
            #sHold_index = Hold_index[:,(j+1-100000):j]
            snoBound_Up = noBound_Up[(j+1-100000):j]
            snoBound_Down = noBound_Down[(j+1-100000):j]
            savgYBound_Up = avgYBound_Up[(j+1-100000):j]
            savgYBound_Down = avgYBound_Down[(j+1-100000):j]
            savgYUnbound_Up = avgYUnbound_Up[(j+1-100000):j]
            savgYUnbound_Down = avgYUnbound_Down[(j+1-100000):j]
            sz = z[(j+1-100000):j]
            st = tPassed[(j+1-100000):j]
            #@save "$folderName/$Notes+$FileCount+_results.bson" sBStateX sHold_index sz st
            @save "$folderName/$Notes+$FileCount+_results.bson" snoBound_Up snoBound_Down savgYBound_Up savgYBound_Down savgYUnbound_Up savgYUnbound_Down sz st
            FileCount+=1
        end
    end


    #### Plootting sections
    p1 = plot(legend = false)
    p2 = plot(legend = false)
    for i in 1:NumGenerators
        #println("in scatter loop")
        UpBXinds = findall(x -> x > 0, BStateX[i,:])
        DownBXinds = findall(x -> x > 0, BStateX[NumGenerators+i,:])
        # scatter bound indices
        p1 = scatter!(p1, tPassed[UpBXinds], ExtList[convert(Array{Int64,1},Hold_index[i,UpBXinds])], markercolor = :green, markersize = 0.5)
        p2 = scatter!(p2, tPassed[DownBXinds], ExtList[convert(Array{Int64,1},Hold_index[NumGenerators+i,DownBXinds])], markercolor = :green, markersize = 0.5)

        #p1 = plot!(p1, tPassed[:], ExtList[convert(Array{Int64,1},Hold_index[i,:])], xlabel = "Time", ylabel = "extension state")
        #p2 = plot!(p2, tPassed[:], ExtList[convert(Array{Int64,1},Hold_index[NumGenerators+i,:])], xlabel = "Time", ylabel = "extension state")


        UpUXinds = findall(x -> x < 0, BStateX[i,:])
        DownUXinds = findall(x -> x < 0, BStateX[NumGenerators+i,:])
        # scatter unbound indices
        p1 = scatter!(p1, tPassed[UpUXinds], ExtList[convert(Array{Int64,1},Hold_index[i,UpUXinds])], markercolor = :red, markersize = 0.5)
        p1 = plot!( ylabel = L"\bar{y}^+", yguidefontsize = 12)
        p2 = scatter!(p2, tPassed[DownUXinds], ExtList[convert(Array{Int64,1},Hold_index[NumGenerators+i,DownUXinds])], markercolor = :red, markersize = 0.5)
        p2 = plot!(ylabel = L"\bar{y}^-", yguidefontsize = 12, xlabelfontsize = 12, xlabel = L"\bar{t}")
    end

     #p = plot!(ylim = (5e-8, 6e-8))
    p = plot(p1, p2, layout = (2,1))
    display(p)
    savefig(p,"$folderName/$Notes")

    plot(tPassed, z, xlabel = L"\bar{t}", ylabel = L"\bar{z}", yguidefontsize = 12, xguidefontsize = 12)
    savefig("$folderName/PolePosition")
return 1
end

export simulate

end
