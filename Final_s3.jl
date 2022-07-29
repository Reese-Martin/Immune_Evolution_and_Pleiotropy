using FileIO
using Statistics 
using Fontconfig
using DataFrames
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using KernelDensity
using Statistics
using Trapz

const V = 2
const UseCoef = .01
const Infs = 3
const Runs = 100

#region Load Data 
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"

FixRplotDir = string(topDir,"/Data/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/Draft_2_Data/")
slplotDir = string(topDir,"/Data/SlowEvo/Draft_2_Data/")
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> !contains(x, "Network"), FixRFiles)
cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> !contains(x, "Network"), FixUDFiles)
cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> !contains(x, "Network"), slFiles)

#group data from saved files
uncDataSet = filter(x -> contains(x, "Uncon"), FixRFiles);
FixRDataSet = filter(x -> contains(x, "Con"), FixRFiles);
FixUDataSet = filter(x -> contains(x, "upreg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "downreg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

uncHostFit = []
FixRHostFit = []
FixUHostFit = []
FixDHostFit = []
sl100HostFit = []

uncParFit = []
FixRParFit = []
FixUParFit = []
FixDParFit = []
sl100ParFit = []
#load data
for i in 1:Infs 
    push!(uncHostFit,FileIO.load(string(FixRplotDir,uncDataSet[i]))["HostFit"])
    push!(FixRHostFit,FileIO.load(string(FixRplotDir,FixRDataSet[i]))["HostFit"])
    push!(FixUHostFit,FileIO.load(string(FixUDplotDir,FixUDataSet[i]))["HostFit"])
    push!(FixDHostFit,FileIO.load(string(FixUDplotDir,FixDDataSet[i]))["HostFit"])
    push!(sl100HostFit,FileIO.load(string(slplotDir,sl100DataSet[i]))["HostFit"])

    push!(uncParFit,FileIO.load(string(FixRplotDir,uncDataSet[i]))["ParFit"])
    push!(FixRParFit,FileIO.load(string(FixRplotDir,FixRDataSet[i]))["ParFit"])
    push!(FixUParFit,FileIO.load(string(FixUDplotDir,FixUDataSet[i]))["ParFit"])
    push!(FixDParFit,FileIO.load(string(FixUDplotDir,FixDDataSet[i]))["ParFit"])
    push!(sl100ParFit,FileIO.load(string(slplotDir,sl100DataSet[i]))["ParFit"])
end
#endregion

#region supplemental figure 3
colSchem = cgrad(:haline, 100, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1750, 1100))
lims = (0,50,nothing,1)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = (0,50,0,1), xticks = ([1,25,50],["1","250","500"]),xlabel = "Generations", yticks = ([0,1],["0", "1"]), ylabel = "Absolute Fitness", title = "10% Chance of Infection")
gb= f[1, 2] = GridLayout()
gb[1,1] = Axis(f, limits = (0,50,0,.5), xticks = ([1,25,50],["1","250","500"]),xlabel = "Generations", yticks = ([0,.5],["0", ".5"]), title = "50% Chance of Infection")
gc= f[1, 3] = GridLayout()
gc[1,1] = Axis(f, limits = (0,50,0,.5), xticks = ([1,25,50],["1","250","500"]),xlabel = "Generations", yticks = ([0,.5],["0", ".5"]), title = "90% Chance of Infection")
sizes = vcat([20,18,16,14],fill(21,46))
Spacer = 10
i = 1
unc = lines!(ga[1,1],vec(mean(uncHostFit[i], dims = 1)[1:Spacer:end]), linewidth = 3, color = colSchem[90])
FixR = scatter!(ga[1,1],vec(mean(FixRHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :rect, color = colSchem[1])
FixU = scatter!(ga[1,1],vec(mean(FixUHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :utriangle, color = colSchem[30])
FixD = scatter!(ga[1,1],vec(mean(FixDHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :diamond, color = colSchem[50])
sl100 = scatter!(ga[1,1],vec(mean(sl100HostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :circle, color = colSchem[70])
i = 2
unc = lines!(gb[1,1],vec(mean(uncHostFit[i], dims = 1)[1:Spacer:end]), linewidth = 3, color = colSchem[90])
FixR = scatter!(gb[1,1],vec(mean(FixRHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :rect, color = colSchem[1])
FixU = scatter!(gb[1,1],vec(mean(FixUHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :utriangle, color = colSchem[30])
FixD = scatter!(gb[1,1],vec(mean(FixDHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :diamond, color = colSchem[50])
sl100 = scatter!(gb[1,1],vec(mean(sl100HostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :circle, color = colSchem[70])
i = 3
unc = lines!(gc[1,1],vec(mean(uncHostFit[i], dims = 1)[1:Spacer:end]), linewidth = 3, color = colSchem[90])
FixR = scatter!(gc[1,1],vec(mean(FixRHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :rect, color = colSchem[1])
FixU = scatter!(gc[1,1],vec(mean(FixUHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :utriangle, color = colSchem[30])
FixD = scatter!(gc[1,1],vec(mean(FixDHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :diamond, color = colSchem[50])
sl100 = scatter!(gc[1,1],vec(mean(sl100HostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :circle, color = colSchem[70])

e1 = MarkerElement( marker = :rect, color = colSchem[1], markersize = 20)
e2 = MarkerElement(marker = :utriangle, color = colSchem[30], markersize = 20)
e3 = MarkerElement(marker = :diamond, color = colSchem[50], markersize = 20)
e4 = MarkerElement(marker = :circle, color = colSchem[70], markersize = 20)

Legend(f[1,end+1],[unc,e1,e2,e3,e4],["Non-Pleiotropic","Fixed Random", "Fixed Up","Fixed Down","Slow"])
cd(string(topDir,"/Images"))
CairoMakie.save("s3_Host_Fit.svg", f, pt_per_unit = 1)

#region supplemental figure 5
colSchem = cgrad(:haline, 100, categorical = true);
fig = Figure(resolution = (1200,1200))
ax1 = Axis(fig[1,1], xticklabelsvisible = false, xticksvisible = false, limits = (-2,50,0,1), title = "Host fitness- 10%")
ax2 = Axis(fig[1,2], xticklabelsvisible = false, xticksvisible = false, limits = (-2,50,0,1), yticklabelsvisible = false, yticksvisible = false, title = "Parasite fitness- 10%")

ax3 = Axis(fig[2,1], xticklabelsvisible = false, xticksvisible = false, limits = (-2,50,0,1), title = "Host Fitness- 50%")
ax4 = Axis(fig[2,2], xticklabelsvisible = false, xticksvisible = false, limits = (-2,50,0,1), yticklabelsvisible = false, yticksvisible = false, title = "Parasite fitness- 50%")

ax5 = Axis(fig[3,1], xticks = ([1,25,50],["1","250","500"]), limits = (-2,50,0,1), title = "Host fitness- 90%")
ax6 = Axis(fig[3,2], xticks = ([1,25,50],["1","250","500"]), limits = (-2,50,0,1), yticklabelsvisible = false, yticksvisible = false, title = "Parasite fitness- 90%")
sizes = vcat([16,14,12,10],fill(8,46))

for i in 1:Infs
    Spacer = 10
    unc = lines!(fig[i,1],vec(mean(uncHostFit[i], dims = 1)[1:Spacer:end]), linewidth = 3, color = colSchem[90])
    FixR = scatter!(fig[i,1],vec(mean(FixRHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :rect,color = colSchem[1])
    FixU = scatter!(fig[i,1],vec(mean(FixUHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :utriangle, color = colSchem[30])
    FixD = scatter!(fig[i,1],vec(mean(FixDHostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :diamond, color = colSchem[50])
    sl100 = scatter!(fig[i,1],vec(mean(sl100HostFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :circle, color = colSchem[70])
    
    unc = lines!(fig[i,2],vec(mean(uncParFit[i], dims = 1)[1:Spacer:end]), linewidth = 3, color = colSchem[90])
    FixR = scatter!(fig[i,2],vec(mean(FixRParFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :rect, color = colSchem[1])
    FixU = scatter!(fig[i,2],vec(mean(FixUParFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :utriangle, color = colSchem[30])
    FixD = scatter!(fig[i,2],vec(mean(FixDParFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :diamond, color = colSchem[50])
    sl100 = scatter!(fig[i,2],vec(mean(sl100ParFit[i],dims = 1)[1:Spacer:end]), markersize = sizes, marker = :circle, color = colSchem[70])
    Legend(fig[2,3],[unc,FixR,FixU,FixD,sl100],["Non-Pleio","Fixed Random", "Fixed Up","Fixed Down","100x slowed"])
end
Label(fig[:,0], "Fitness", rotation = pi/2)
Label(fig[end+1,:], "Generations")
cd(string(topDir,"/Images"))
GLMakie.save("s5_Host_Par_Fit.png", fig, pt_per_unit = 1)
#endregion
