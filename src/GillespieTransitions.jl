#
# GillespieTransitions.jl
# Stochastic model
#
# Created by Dionn Hargreaves on 01/07/2021
#
#
# Run Gillespie Algorithm


module GillespieTransitions

using Random



@inline @views function gillespieTran!(stateIndVec, param, tPassed)
    # choose a random variable
    choiceParam, timeParam = rand(2)


    # find time to next event
    τ = (1/sum(param)).*log(1/timeParam)
    # save time stamp
    tPassed = tPassed+τ

    #param .= param.*τ
    # accumulate probabilities
    param_line = cumsum(param)./sum(param)
    # Find where along the accumulation line the second random variable lies
    choosing = param_line .<= choiceParam

    # find value to use (first item of choosing that == 0)
    if all(x-> x == 1,  choosing)
        chooseFunc = findall(x->x>0, choosing)[end]
    else
        chooseFunc = findall(x->x==0, choosing)[1]
    end

    #if tPassed - τ == 1.7522253429036234
    #    println(chooseFunc)
    #    println(choosing)
    #end
    #if chooseFunc != length(param_line) && param_line[chooseFunc] - param_line[chooseFunc-1] < 1e-6
    #    chooseFunc += 1
    #end

    # which state is it?
    stateIndVec[1] = mod(chooseFunc, 3)

    # which generator is it?
    stateIndVec[2] = ceil(Int64, chooseFunc/3)


    return tPassed
end


export gillespieTran!

end
