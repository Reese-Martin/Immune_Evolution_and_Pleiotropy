#v2 cleans up old functions and unused plots
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using FileIO
using Statistics
using HypothesisTests
using KernelDensity
using DataFrames
using StatsBase

ToSave = true
Gens = 500
Runs = 100
Infs = 3
Hosts = 500
ImmuneTypes = 4 #0- no response, 1- constitutive, 2- Inducible, 3- mixed indu/const
#region Functions
function MakeFigStackDens(lims = (nothing, nothing, nothing, nothing))
    labels = ["None" "Constitutive" "Mixed" "Inducible"]
    fig = Figure()

    ax1 = Axis(fig[1,1],xticklabelsvisible = false,xticksvisible = false, title = "Unc v FixR", ylabel = "10%",limits = lims)
    ax2 = Axis(fig[2,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "50%",limits = lims)
    ax3 = Axis(fig[3,1], ylabel = "90%",limits = lims,xticks = ([1,128],["Cosnt", "Indu"]))

    ax4 = Axis(fig[1,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,title = "Unc v FixU",limits = lims)
    ax5 = Axis(fig[2,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,limits = lims)
    ax6 = Axis(fig[3,2],yticklabelsvisible = true, yticksvisible = false,limits = lims,xticks = ([1,128],["Cosnt", "Indu"]))

    ax7 = Axis(fig[1,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,title = "Unc v FixD",limits = lims)
    ax8 = Axis(fig[2,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,limits = lims)
    ax9 = Axis(fig[3,3],yticklabelsvisible = true, yticksvisible = false,limits = lims,xticks = ([1,128],["Cosnt", "Indu"]))

    ax10 = Axis(fig[1,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,title = "Unc v SLx100",limits = lims)
    ax11 = Axis(fig[2,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = true, yticksvisible = false,limits = lims)
    ax12 = Axis(fig[3,4],yticklabelsvisible = true, yticksvisible = false,limits = lims,xticks = ([1,128],["Cosnt", "Indu"]))
    return fig
end
function DetPropIndu(x) #x should be a single infection percentage
    tmp = x.+.1*10^-10 #add a small value to each element to prevent 0/0 events
    LenRealData = sum((!).(isnan.(x[1,:,1])))
    perConst = round.(tmp[:,1:LenRealData,1]./tmp[:,1:LenRealData,2],digits = 2) #round removes the error introduced by epsilon addition, gives percent of peak response that was constitutive
    perIndu = 1 .- perConst
    return perIndu
end
#endregion
#region processing
# topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Co-evo model/"
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Co-evo model/Julia Code/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/IMMUNE_PLEIO_MODEL/Julia Code"

FixRplotDir = string(topDir,"/Data/FixedRand/Draft_2_Data/")#string(topDir,"/Data/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/Draft_2_Data/")#string(topDir,"/Data/FixedUpDown/Draft_2_Data/")
slplotDir = string(topDir,"/Data/SlowEvo/Draft_2_Data/")#string(topDir,"/Data/SlowEvo/Draft_2_Data/")
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> !contains(x, "Network"), FixRFiles)
FixRFiles = filter(x -> !contains(x, "short"), FixRFiles)

cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> !contains(x, "Network"), FixUDFiles)
FixUDFiles = filter(x -> !contains(x, "short"), FixUDFiles)

cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> !contains(x, "Network"), slFiles)
slFiles = filter(x -> !contains(x, "short"), slFiles)

#group data from saved files
uncDataSet = filter(x -> contains(x, "Uncon"), FixRFiles);
FixRDataSet = filter(x -> contains(x, "Con"), FixRFiles);
FixUDataSet = filter(x -> contains(x, "upreg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "downreg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

uncImmMag = zeros(Infs,Runs,Gens,2); #first :,:,1 is pre inf eq, :,:,2 is peak of immune response
FixRImmMag = zeros(Infs,Runs,Gens,2);
FixUImmMag = zeros(Infs,Runs,Gens,2);
FixDImmMag = zeros(Infs,Runs,Gens,2);
sl100ImmMag = zeros(Infs,Runs,Gens,2);

for i = 1:Infs
    unctmpDict = FileIO.load(string(FixRplotDir,uncDataSet[i]))
    FixRtmpDict = FileIO.load(string(FixRplotDir,FixRDataSet[i]))
    FixUtmpDict = FileIO.load(string(FixUDplotDir,FixUDataSet[i]))
    FixDtmpDict = FileIO.load(string(FixUDplotDir,FixDDataSet[i]))
    sl100tmpDict = FileIO.load(string(slplotDir,sl100DataSet[i]))

    uncImmMag[i,:,:,:] = unctmpDict["ImmMagInf"]
    FixRImmMag[i,:,:,:] = FixRtmpDict["ImmMagInf"]
    FixUImmMag[i,:,:,:] = FixUtmpDict["ImmMagInf"]
    FixDImmMag[i,:,:,:] = FixDtmpDict["ImmMagInf"]
    sl100ImmMag[i,:,:,:] = sl100tmpDict["ImmMagInf"]
end

uncInduDens = []
FixRInduDens = []
FixUInduDens = []
FixDInduDens = []
sl100InduDens = []
for i in 1:Infs
    tmp1 = DetPropIndu(uncImmMag[i,:,:,:]) #returns the proportion of the total immune response that is inducible
    tmp2 = DetPropIndu(FixRImmMag[i,:,:,:])
    tmp3 = DetPropIndu(FixUImmMag[i,:,:,:])
    tmp4 = DetPropIndu(FixDImmMag[i,:,:,:])
    tmp5 = DetPropIndu(sl100ImmMag[i,:,:,:])

    push!(uncInduDens, tmp1)
    push!(FixRInduDens, tmp2)
    push!(FixUInduDens, tmp3)
    push!(FixDInduDens, tmp4)
    push!(sl100InduDens, tmp5)
end
#endregion
#region Figure 2
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1500))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood", xlabel = "Immune Response")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]),xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]

colSchem = cgrad(:haline, 100, categorical = true);
for i in 1:Infs
    tmp1 = kde(vec(uncInduDens[i]), npoints = 128).density
    tmp2 = kde(vec(FixRInduDens[i]), npoints = 128).density
    tmp3 = kde(vec(FixUInduDens[i]), npoints = 128).density
    tmp4 = kde(vec(FixDInduDens[i]), npoints = 128).density
    tmp5 = kde(vec(sl100InduDens[i]), npoints = 128).density

    lines!(ga[i,1],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])
    lines!(ga[i,1],tmp2./maximum(tmp2),linewidth = 3,color = colSchem[1])
    text!(ga[i,1],string(round(cor(tmp1./maximum(tmp1),tmp2./maximum(tmp2)),digits = 2)),position = (50,.5))
    lines!(ga[i,2],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])
    lines!(ga[i,2],tmp3./maximum(tmp3),linewidth = 3,color = colSchem[1])
    text!(ga[i,2],string(round(cor(tmp1./maximum(tmp1),tmp3./maximum(tmp3)),digits = 2)),position = (50,.5))
    lines!(ga[i,3],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])
    lines!(ga[i,3],tmp4./maximum(tmp4),linewidth = 3,color = colSchem[1])
    text!(ga[i,3],string(round(cor(tmp1./maximum(tmp1),tmp4./maximum(tmp4)),digits = 2)),position = (50,.5))
    lines!(ga[i,4],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])
    lines!(ga[i,4],tmp5./maximum(tmp5),linewidth = 3,color = colSchem[1])
    text!(ga[i,4],string(round(cor(tmp1./maximum(tmp1),tmp5./maximum(tmp5)),digits = 2)),position = (50,.5))
end

#density calculated separately then added together
gb = f[2,1]
Axis(gb,xticks = ([0,1],["Constitutive", "Inducible"]), yticks = ([0,1], ["Min.", "Max."]), ylabel = "Peak Immune Effector Activity", xlabel = "Immune Response")
dens = []
numpars = [50,250,450]
colSchem = cgrad(:haline,rev = true);
numpoints = 128
for i in 1:Infs
    np = numpars[i]
    Uncx = vec(uncInduDens[i])
    Uncy = vec(uncImmMag[i,:,1:np,2])
    UncDens = kde((Uncx,Uncy),npoints = (numpoints,numpoints)).density

    FixRx = vec(FixRInduDens[i])
    FixRy = vec(FixRImmMag[i,:,1:np,2])
    FixRDens = kde((FixRx,FixRy),npoints = (numpoints,numpoints)).density

    FixUx = vec(FixUInduDens[i])
    FixUy = vec(FixUImmMag[i,:,1:np,2])
    FixUDens = kde((FixUx,FixUy),npoints = (numpoints,numpoints)).density

    FixDx = vec(FixDInduDens[i])
    FixDy = vec(FixDImmMag[i,:,1:np,2])
    FixDDens = kde((FixDx,FixDy),npoints = (numpoints,numpoints)).density

    sl100x = vec(sl100InduDens[i])
    sl100y = vec(sl100ImmMag[i,:,1:np,2])
    sl100Dens = kde((sl100x,sl100y),npoints = (numpoints,numpoints)).density

    push!(dens,UncDens./maximum(UncDens) + FixRDens./maximum(FixRDens) + FixUDens./maximum(FixUDens) + FixDDens./maximum(FixDDens) +sl100Dens./maximum(sl100Dens))

end
hm = heatmap!(gb,0:1/numpoints:1,0:1/numpoints:1,dens[1]+dens[2]+dens[3], colormap = Reverse(:haline))
Colorbar(f[2,end+1], limits = (0,5), colormap = Reverse(:haline), ticks = ([0,5], ["Low Density", "High Density"]))
# Label(ga[1,0], "10% Chance of Infection", rotation = pi/2)
# Label(ga[2,0], "50% Chance of Infection", rotation = pi/2)
# Label(ga[3,0], "90% Chance of Infection", rotation = pi/2)
# Label(f[2,0], "Effector Activity vs Induciblity of Response", rotation = pi/2)
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)

cd(string(topDir,"/Images"))
CairoMakie.save("f2_kde_heatmap.svg", f, pt_per_unit = 1)
#endregion
#region s7
NumPars = [50,250,450]
thresh = .1:.1:1
colSchem = cgrad(:haline, 100, categorical = true);
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,10,nothing,1)
#ylabels = ["50%","90%"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "Proportion of Inducible Hosts",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "Proportion of Inducible Hosts")
ga[3, 1] = Axis(f, limits = lims, xticks = ([1,5,9],[".1",".5",".9"]), xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "Proportion of Inducible Hosts", xlabel = "Percent of Response Induced")
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticklabelsvisible = false,yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticks =  ([1,5,9],[".1",".5",".9"]), xgridvisible = false, ygridvisible = false, yticklabelsvisible = false, yticksvisible = false, xlabel = "Percent of Response Induced") for j in 2:5]

for i in 1:Infs
    uncPerofIndivsIndu = zeros(length(thresh))
    FixRPerofIndivsIndu = zeros(length(thresh))
    FixUPerofIndivsIndu = zeros(length(thresh))
    FixDPerofIndivsIndu = zeros(length(thresh))
    sl100PerofIndivsIndu = zeros(length(thresh))
    count = 1
    for j in thresh
        uncRowsWithIndu = unique(getindex.(findall(x -> x >j, uncInduDens[i]),1))
        uncPerofIndivsIndu[count] = mean(sum(uncInduDens[i][uncRowsWithIndu,:].>j,dims = 2)./NumPars[i])

        FixRRowsWithIndu = unique(getindex.(findall(x -> x >j, FixRInduDens[i]),1))
        FixRPerofIndivsIndu[count] = mean(sum(FixRInduDens[i][FixRRowsWithIndu,:].>j,dims = 2)./NumPars[i])

        FixURowsWithIndu = unique(getindex.(findall(x -> x >j, FixUInduDens[i]),1))
        FixUPerofIndivsIndu[count] = mean(sum(FixUInduDens[i][FixURowsWithIndu,:].>j,dims = 2)./NumPars[i])

        FixDRowsWithIndu = unique(getindex.(findall(x -> x >j, FixDInduDens[i]),1))
        FixDPerofIndivsIndu[count] = mean(sum(FixDInduDens[i][FixDRowsWithIndu,:].>j,dims = 2)./NumPars[i])

        sl100RowsWithIndu = unique(getindex.(findall(x -> x >j, sl100InduDens[i]),1))
        sl100PerofIndivsIndu[count] = mean(sum(sl100InduDens[i][sl100RowsWithIndu,:].>j,dims = 2)./NumPars[i])
        count += 1
    end

    barplot!(ga[i,1],uncPerofIndivsIndu)#, color = [colSchem[1],colSchem[40],colSchem[80]])
    barplot!(ga[i,2],FixRPerofIndivsIndu)#, color =  [colSchem[1],colSchem[40],colSchem[80]])
    barplot!(ga[i,3],FixUPerofIndivsIndu)#, color =  [colSchem[1],colSchem[40],colSchem[80]])
    barplot!(ga[i,4],FixDPerofIndivsIndu)#, color =  [colSchem[1],colSchem[40],colSchem[80]])
    barplot!(ga[i,5],sl100PerofIndivsIndu)#, color =  [colSchem[1],colSchem[40],colSchem[80]])
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
# Label(fig[end+1, :], text = "Percent of Hosts inducible", textsize = 25)
# Label(fig[:, 0], text = "Chance of Infection", textsize = 25,
#     rotation = pi/2)
cd(string(topDir,"/Images"))
CairoMakie.save("s7_InduFixation_Bar.svg", f, pt_per_unit = 1)
#endregion