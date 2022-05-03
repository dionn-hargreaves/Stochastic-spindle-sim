using DelimitedFiles
using Plots
using LaTeXStrings

numberGens = 15

Notes = "N15_test"
#cd("/Volumes/DH_Simulations/")
Path="output/GillespieRunData_N15_test0"

upperCortex = readdlm("$Path/$Notes+_upper.txt")
lowerCortex = readdlm("$Path/$Notes+_lower.txt")
pole = readdlm("$Path/$Notes+_pole.txt")
useTime = readdlm("$Path/$Notes+_timestamp.txt")

noBound_Up = upperCortex[:,1]
noBound_Down = lowerCortex[:,1]
avgYBound_Up = upperCortex[:,2]
avgYBound_Down = lowerCortex[:,2]
avgYUnbound_Up = upperCortex[:,3]
avgYUnbound_Down = lowerCortex[:,3]

#### Plootting sections
p1 = plot(legend = false, ylim=(0.0,6.0), xlabel=L"\bar{t}", ylabel=L"\bar{y}^+_{\textrm{avg}}")
p2 = plot(legend = false, ylim=(0.0,6.0), xlabel=L"\bar{t}", ylabel=L"\bar{y}^-_{\textrm{avg}}")
p1 = scatter!(p1, useTime, avgYBound_Up, markercolor = :green, markersize = 0.1)
p1 = scatter!(p1, useTime, avgYUnbound_Up, markercolor = :red, markersize = 0.1)
p2 = scatter!(p2, useTime, avgYBound_Down, markercolor = :green, markersize = 0.1)
p2 = scatter!(p2, useTime, avgYUnbound_Down, markercolor = :red, markersize = 0.1)

 #p = plot!(ylim = (5e-8, 6e-8))
p = plot(p1, p2, layout = (2,1), dpi = 300)
#display(p)
savefig(p,"$Path/Pops")

p1 = plot(legend = false, xlabel=L"\bar{t}", ylabel=L"n_b^+", ylim=(0.0,14.0))
p2 = plot(legend = false, xlabel=L"\bar{t}", ylabel=L"n_b^-", ylim=(0.0,14.0))
p1 = scatter!(p1, useTime, noBound_Up, markercolor = :green, markersize = 2)
p2 = scatter!(p2, useTime, noBound_Down, markercolor = :green, markersize = 2)


 #p = plot!(ylim = (5e-8, 6e-8))
p = plot(p1, p2, layout = (2,1), dpi = 300)

#display(p)
savefig(p,"$Path/Pops_nums")


plot(useTime, pole, legend = false, ylabel = L"\bar{z}", xlabel = L"\bar{t}", dpi=300, linewidth = 1)
savefig("$Path/z")
