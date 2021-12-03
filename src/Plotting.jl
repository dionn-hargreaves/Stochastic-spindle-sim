using BSON: @load
using Plots
using LaTeXStrings

chRes = 1000
saveSize = 100000
numberGens = 10
Size = convert(Int64,ceil(saveSize/chRes))

NumberOfFiles = 430
BStateX = zeros(numberGens*2, NumberOfFiles*Size)
Hold_index = zeros(numberGens*2, NumberOfFiles*Size)
z = zeros(NumberOfFiles*Size)
tPassed = zeros(NumberOfFiles*Size)
cd("/Users/mbcx4dh9/Documents/Woolner,S_Jensen,O/PhD/Julia_code/Stochastics_bioParams_precomp/")
Path="output/GillespieRunData_N10_Alpha_Beta_Larger0"

for i in 1:NumberOfFiles
    println(i)
    @load "$Path/N10_Alpha_Beta_Larger+$i+_results.bson" sBStateX sHold_index sz st
    BStateX[:,((i-1)*Size+1):(i*Size)] = sBStateX[:,1:chRes:end]
    Hold_index[:,((i-1)*Size+1):(i*Size)] = sHold_index[:,1:chRes:end]
    z[((i-1)*Size+1):(i*Size)] = sz[1:chRes:end]
    tPassed[((i-1)*Size+1):(i*Size)] = st[1:chRes:end]
end

Bound_Holdindex = similar(Hold_index) .- 1
UnBound_Holdindex = similar(Hold_index) .- 1

for i in 1:length(Hold_index[:,1])
      BXinds = findall(x -> x >0, BStateX[i,:])
      Bound_Holdindex[i, BXinds] = Hold_index[i, BXinds]
      UXinds = findall(x -> x<0, BStateX[i,:])
      UnBound_Holdindex[i, UXinds] = Hold_index[i, UXinds]
end


 #=
@gif for i in 1:10:NumberOfFiles*Size
      histogram(Bound_Holdindex[:,i],bins = 100, xlim = (0.0,15000), ylim = (0.0,6.0), alpha = 0.5, label = "bound")
      histogram!(UnBound_Holdindex[:,i], bins = 100, alpha=0.5, label = "unbound")
end
=#
#### Plootting sections
p1 = plot(legend = false)
p2 = plot(legend = false)

NumGenerators = convert(Int64,length(Hold_index[:,1])/2)
for i in 1:NumGenerators
    #println("in scatter loop")
    UpBXinds = findall(x -> x > 0, BStateX[i,:])
    DownBXinds = findall(x -> x > 0, BStateX[NumGenerators+i,:])
    # scatter bound indices
    global p1 = scatter!(p1, tPassed[UpBXinds], Hold_index[i,UpBXinds], markercolor = :green, markersize = 0.5)
    global p2 = scatter!(p2, tPassed[DownBXinds], Hold_index[NumGenerators+i,DownBXinds], markercolor = :green, markersize = 0.5)
  
    #p1 = plot!(p1, tPassed[:], ExtList[convert(Array{Int64,1},Hold_index[i,:])], xlabel = "Time", ylabel = "extension state")
    #p2 = plot!(p2, tPassed[:], ExtList[convert(Array{Int64,1},Hold_index[NumGenerators+i,:])], xlabel = "Time", ylabel = "extension state")


    UpUXinds = findall(x -> x < 0, BStateX[i,:])
    DownUXinds = findall(x -> x < 0, BStateX[NumGenerators+i,:])
    # scatter unbound indices
    global p1 = scatter!(p1, tPassed[UpUXinds], Hold_index[i,UpUXinds], markercolor = :red, markersize = 0.5)
    global p1 = plot!( ylabel = L"\bar{y}^+", yguidefontsize = 12)
    global p2 = scatter!(p2, tPassed[DownUXinds], Hold_index[NumGenerators+i,DownUXinds], markercolor = :red, markersize = 0.5)
    global p2 = plot!(ylabel = L"\bar{y}^-", yguidefontsize = 12, xlabelfontsize = 12, xlabel = L"\bar{t}")
end

 #p = plot!(ylim = (5e-8, 6e-8))
p = plot(p1, p2, layout = (2,1), dpi = 300)
display(p)
savefig(p,"$Path/Pops")


plot(tPassed, z, legend = false, ylabel = L"\bar{z}", xlabel = L"\bar{t}", dpi=300, linewidth = 2)
savefig("$Path/z")
