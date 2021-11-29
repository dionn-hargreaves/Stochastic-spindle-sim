using BSON: @load
using Plots
using LaTeXStrings

chRes = 1000
saveSize = 100000
numberGens = 30
Size = convert(Int64,ceil(saveSize/chRes))

NumberOfFiles = 500

noBound_Up = zeros(NumberOfFiles*Size)
noBound_Down = zeros(NumberOfFiles*Size)
avgYBound_Up = zeros(NumberOfFiles*Size)
avgYBound_Down = zeros(NumberOfFiles*Size)
avgYUnbound_Up = zeros(NumberOfFiles*Size)
avgYUnbound_Down = zeros(NumberOfFiles*Size)

z = zeros(NumberOfFiles*Size)
tPassed = zeros(NumberOfFiles*Size)
Notes = "N30_Alpha_Beta_Larger_bigger_onrate_smallerK"
cd("/Users/mbcx4dh9/Documents/Woolner,S_Jensen,O/PhD/Julia_code/Stochastics_bioParams_precomp/")
Path="output/GillespieRunData_N30_Alpha_Beta_Larger_bigger_onrate_smallerK5"

for i in 1:NumberOfFiles
    println(i)
    @load "$Path/$Notes+$i+_results.bson" snoBound_Up snoBound_Down savgYBound_Up savgYBound_Down savgYUnbound_Up savgYUnbound_Down sz st

    noBound_Up[((i-1)*Size+1):(i*Size)] = snoBound_Up[1:chRes:end]
    noBound_Down[((i-1)*Size+1):(i*Size)] = snoBound_Down[1:chRes:end]
    avgYBound_Up[((i-1)*Size+1):(i*Size)] = savgYBound_Up[1:chRes:end]
    avgYBound_Down[((i-1)*Size+1):(i*Size)] = savgYBound_Down[1:chRes:end]
    avgYUnbound_Up[((i-1)*Size+1):(i*Size)] = savgYUnbound_Up[1:chRes:end]
    avgYUnbound_Down[((i-1)*Size+1):(i*Size)] = savgYUnbound_Down[1:chRes:end]

    z[((i-1)*Size+1):(i*Size)] = sz[1:chRes:end]
    tPassed[((i-1)*Size+1):(i*Size)] = st[1:chRes:end]
end


#### Plootting sections
p1 = plot(legend = false, ylim=(0.0,6.0), xlabel=L"\bar{t}", ylabel=L"\bar{y}^+_{\textrm{avg}}")
p2 = plot(legend = false, ylim=(0.0,6.0), xlabel=L"\bar{t}", ylabel=L"\bar{y}^+_{\textrm{avg}}")
p1 = scatter!(p1, tPassed, avgYBound_Up, markercolor = :green, markersize = 0.1)
p1 = scatter!(p1, tPassed, avgYUnbound_Up, markercolor = :red, markersize = 0.1)
p2 = scatter!(p2, tPassed, avgYBound_Down, markercolor = :green, markersize = 0.1)
p2 = scatter!(p2, tPassed, avgYUnbound_Down, markercolor = :red, markersize = 0.1)

 #p = plot!(ylim = (5e-8, 6e-8))
p = plot(p1, p2, layout = (2,1), dpi = 300)
#display(p)
savefig(p,"$Path/Pops")

p1 = plot(legend = false, xlabel=L"\bar{t}", ylabel=L"n_b^+")
p2 = plot(legend = false, xlabel=L"\bar{t}", ylabel=L"n_b^-")
p1 = scatter!(p1, tPassed, noBound_Up, markercolor = :green, markersize = 1)
p2 = scatter!(p2, tPassed, noBound_Down, markercolor = :green, markersize = 1)


 #p = plot!(ylim = (5e-8, 6e-8))
p = plot(p1, p2, layout = (2,1), dpi = 300)

#display(p)
savefig(p,"$Path/Pops_nums")


plot(tPassed, z, legend = false, ylabel = L"\bar{z}", xlabel = L"\bar{t}", dpi=300, linewidth = 2)
savefig("$Path/z")
