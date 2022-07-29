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
points = 11
#region Functions
function MakeFigInduFix(lims = (nothing, nothing, nothing, nothing))
    labels = ["None" "Constitutive" "Mixed" "Inducible"]
    fig = Figure()

    ax1 = Axis(fig[1,1],xticklabelsvisible = false,xticksvisible = false, title = "Unconstrained", ylabel = "10%",limits = lims)
    ax2 = Axis(fig[2,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "50%",limits = lims)
    ax3 = Axis(fig[3,1], ylabel = "90%",limits = lims, xticks = (1:4:9,[".1",".5",".9"]))

    ax4 = Axis(fig[1,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Random",limits = lims)
    ax5 = Axis(fig[2,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax6 = Axis(fig[3,2],yticklabelsvisible = false, yticksvisible = false,limits = lims, xticks = (1:4:9,[".1",".5",".9"]))

    ax7 = Axis(fig[1,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Up",limits = lims)
    ax8 = Axis(fig[2,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax9 = Axis(fig[3,3],yticklabelsvisible = false, yticksvisible = false,limits = lims, xticks = (1:4:9,[".1",".5",".9"]))

    ax10 = Axis(fig[1,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Down",limits = lims)
    ax11 = Axis(fig[2,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax12 = Axis(fig[3,4],yticklabelsvisible = false, yticksvisible = false,limits = lims, xticks = (1:4:9,[".1",".5",".9"]))

    ax13 = Axis(fig[1,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "100x Slower",limits = lims)
    ax14 = Axis(fig[2,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax15 = Axis(fig[3,5],yticklabelsvisible = false, yticksvisible = false,limits = lims, xticks = (1:4:9,[".1",".5",".9"]))
    return fig
end
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
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Co-evo model/Julia Code/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/IMMUNE_PLEIO_MODEL/Julia Code"

FixRplotDir = string(topDir,"/Data/FixedRand/")#string(topDir,"/Data/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/")#string(topDir,"/Data/FixedUpDown/Draft_2_Data/")
slplotDir = string(topDir,"/Data/SlowEvo/")#string(topDir,"/Data/SlowEvo/Draft_2_Data/")
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> !contains(x, "Network"), FixRFiles)
FixRFiles = filter(x -> contains(x, "short"), FixRFiles)
cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> !contains(x, "Network"), FixUDFiles)
FixUDFiles = filter(x -> contains(x, "short"), FixUDFiles)
cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> !contains(x, "Network"), slFiles)
slFiles = filter(x -> contains(x, "short"), slFiles)

#group data from saved files
uncDataSet = filter(x -> contains(x, "Uncon"), FixRFiles);
FixRDataSet = filter(x -> contains(x, "Con"), FixRFiles);
FixUDataSet = filter(x -> contains(x, "upreg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "downreg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

uncImmMag = zeros(Infs,Runs,points,Gens,2); #first :,:,1 is pre inf eq, :,:,2 is peak of immune response
FixRImmMag = zeros(Infs,Runs,points,Gens,2);
FixUImmMag = zeros(Infs,Runs,points,Gens,2);
FixDImmMag = zeros(Infs,Runs,points,Gens,2);
sl100ImmMag = zeros(Infs,Runs,points,Gens,2);

for i = 1:Infs
    unctmpDict = FileIO.load(string(FixRplotDir,uncDataSet[i]))
    FixRtmpDict = FileIO.load(string(FixRplotDir,FixRDataSet[i]))
    FixUtmpDict = FileIO.load(string(FixUDplotDir,FixUDataSet[i]))
    FixDtmpDict = FileIO.load(string(FixUDplotDir,FixDDataSet[i]))
    sl100tmpDict = FileIO.load(string(slplotDir,sl100DataSet[i]))

    uncImmMag[i,:,:,:,:] = unctmpDict["ImmMagInfSnaps"]
    FixRImmMag[i,:,:,:,:] = FixRtmpDict["ImmMagInfSnaps"]
    FixUImmMag[i,:,:,:,:] = FixUtmpDict["ImmMagInfSnaps"]
    FixDImmMag[i,:,:,:,:] = FixDtmpDict["ImmMagInfSnaps"]
    sl100ImmMag[i,:,:,:,:] = sl100tmpDict["ImmMagInfSnaps"]

    unctmpDict = nothing
    FixRtmpDict = nothing
    FixUtmpDict = nothing
    FixDtmpDict = nothing
    sl100tmpDict = nothing
end

uncInduDens = []
FixRInduDens = []
FixUInduDens = []
FixDInduDens = []
sl100InduDens = []
for i in 1:Infs
    tmp11 = []
    tmp21 = []
    tmp31 = []
    tmp41 = []
    tmp51 = []
    for j in 1:points
        tmp1 = DetPropIndu(uncImmMag[i,:,j,:,:]) #returns the proportion of the total immune response that is inducible
        tmp2 = DetPropIndu(FixRImmMag[i,:,j,:,:])
        tmp3 = DetPropIndu(FixUImmMag[i,:,j,:,:])
        tmp4 = DetPropIndu(FixDImmMag[i,:,j,:,:])
        tmp5 = DetPropIndu(sl100ImmMag[i,:,j,:,:])
        push!(tmp11,tmp1)
        push!(tmp21,tmp2)
        push!(tmp31,tmp3)
        push!(tmp41,tmp4)
        push!(tmp51,tmp5)
    end
    push!(uncInduDens, tmp11)
    push!(FixRInduDens, tmp21)
    push!(FixUInduDens, tmp31)
    push!(FixDInduDens, tmp41)
    push!(sl100InduDens, tmp51)
end
#endregion

#region 10% infection immune response over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 1000))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP","RF","FU","FD","Sl"]
ga = f[1, 1] = GridLayout()
ga[1:5, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i]) for i in 1:5]
ga[1:5, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:5 for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "0" )
Label(ga[1,2, Top()], "50" )
Label(ga[1,3, Top()], "100" )
Label(ga[1,4, Top()], "150" )
Label(ga[1,5, Top()], "200" )
Label(ga[1,6, Top()], "250" )
Label(ga[1,7, Top()], "300" )
Label(ga[1,8, Top()], "350" )
Label(ga[1,9, Top()], "400" )
Label(ga[1,10, Top()], "450" )
Label(ga[1,11, Top()], "500" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
for i in 1:points
    tmp1 = kde(vec(uncInduDens[1][i]), npoints = 128).density
    tmp2 = kde(vec(FixRInduDens[1][i]), npoints = 128).density
    tmp3 = kde(vec(FixUInduDens[1][i]), npoints = 128).density
    tmp4 = kde(vec(FixDInduDens[1][i]), npoints = 128).density
    tmp5 = kde(vec(sl100InduDens[1][i]), npoints = 128).density
    
    lines!(ga[1,i],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])
    lines!(ga[2,i],tmp2./maximum(tmp2),linewidth = 3,color = colSchem[1]) 
    lines!(ga[3,i],tmp3./maximum(tmp3),linewidth = 3,color = colSchem[1])
    lines!(ga[4,i],tmp4./maximum(tmp4),linewidth = 3,color = colSchem[1])
    lines!(ga[5,i],tmp5./maximum(tmp5),linewidth = 3,color = colSchem[1])
end
cd(string(topDir,"/Images"))
GLMakie.save("IoT_kde_10_500.png", f, pt_per_unit = 1)
#endregion

#region 50% infection immune response over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 1000))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP","RF","FU","FD","Sl"]
ga = f[1, 1] = GridLayout()
ga[1:5, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i]) for i in 1:5]
ga[1:5, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:5 for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "0" )
Label(ga[1,2, Top()], "50" )
Label(ga[1,3, Top()], "100" )
Label(ga[1,4, Top()], "150" )
Label(ga[1,5, Top()], "200" )
Label(ga[1,6, Top()], "250" )
Label(ga[1,7, Top()], "300" )
Label(ga[1,8, Top()], "350" )
Label(ga[1,9, Top()], "400" )
Label(ga[1,10, Top()], "450" )
Label(ga[1,11, Top()], "500" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
for i in 1:points
    tmp1 = kde(vec(uncInduDens[2][i]), npoints = 128).density
    tmp2 = kde(vec(FixRInduDens[2][i]), npoints = 128).density
    tmp3 = kde(vec(FixUInduDens[2][i]), npoints = 128).density
    tmp4 = kde(vec(FixDInduDens[2][i]), npoints = 128).density
    tmp5 = kde(vec(sl100InduDens[2][i]), npoints = 128).density
    
    lines!(ga[1,i],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])

    lines!(ga[2,i],tmp2./maximum(tmp2),linewidth = 3,color = colSchem[1]) 

    lines!(ga[3,i],tmp3./maximum(tmp3),linewidth = 3,color = colSchem[1])

    lines!(ga[4,i],tmp4./maximum(tmp4),linewidth = 3,color = colSchem[1])
   
    lines!(ga[5,i],tmp5./maximum(tmp5),linewidth = 3,color = colSchem[1])
end
cd(string(topDir,"/Images"))
GLMakie.save("IoT_kde_50_500.png", f, pt_per_unit = 1)
#endregion

#region 90% infection immune response over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1000, 1000))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP","RF","FU","FD","Sl"]
ga = f[1, 1] = GridLayout()
ga[1:5, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i]) for i in 1:5]
ga[1:5, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:5 for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "0" )
Label(ga[1,2, Top()], "50" )
Label(ga[1,3, Top()], "100" )
Label(ga[1,4, Top()], "150" )
Label(ga[1,5, Top()], "200" )
Label(ga[1,6, Top()], "250" )
Label(ga[1,7, Top()], "300" )
Label(ga[1,8, Top()], "350" )
Label(ga[1,9, Top()], "400" )
Label(ga[1,10, Top()], "450" )
Label(ga[1,11, Top()], "500" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
for i in 1:points
    tmp1 = kde(vec(uncInduDens[3][i]), npoints = 128).density
    tmp2 = kde(vec(FixRInduDens[3][i]), npoints = 128).density
    tmp3 = kde(vec(FixUInduDens[3][i]), npoints = 128).density
    tmp4 = kde(vec(FixDInduDens[3][i]), npoints = 128).density
    tmp5 = kde(vec(sl100InduDens[3][i]), npoints = 128).density
    
    lines!(ga[1,i],tmp1./maximum(tmp1),linewidth = 3,color = colSchem[70])

    lines!(ga[2,i],tmp2./maximum(tmp2),linewidth = 3,color = colSchem[1]) 

    lines!(ga[3,i],tmp3./maximum(tmp3),linewidth = 3,color = colSchem[1])

    lines!(ga[4,i],tmp4./maximum(tmp4),linewidth = 3,color = colSchem[1])
   
    lines!(ga[5,i],tmp5./maximum(tmp5),linewidth = 3,color = colSchem[1])
end
cd(string(topDir,"/Images"))
GLMakie.save("IoT_kde_90_500.png", f, pt_per_unit = 1)
#endregion

#region 10% infection immune response x magnitude over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1800, 1500))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP \n Peak Immune Effector Activity","RF \n Peak Immune Effector Activity","FU \n Peak Immune Effector Activity","FD \n Peak Immune Effector Activity","Sl \n Peak Immune Effector Activity"]
ga = f[1, 1] = GridLayout()
ga[1:4, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[i]) for i in 1:4]
ga[5,1] = Axis(f,xticks = ([0,1],["C", "I"]), yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[5], xlabel = "Immune Response") 
ga[1:4, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:4 for j in 2:11]
ga[5, 2:11] = [Axis(f,xticks = ([0,1],["C", "I"]),yticklabelsvisible = false,yticksvisible = false, xlabel = "Immune Response") for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "Gen: 0" )
Label(ga[1,2, Top()], "Gen: 5" )
Label(ga[1,3, Top()], "Gen: 10" )
Label(ga[1,4, Top()], "Gen: 15" )
Label(ga[1,5, Top()], "Gen: 20" )
Label(ga[1,6, Top()], "Gen: 25" )
Label(ga[1,7, Top()], "Gen: 30" )
Label(ga[1,8, Top()], "Gen: 35" )
Label(ga[1,9, Top()], "Gen: 40" )
Label(ga[1,10, Top()], "Gen: 45" )
Label(ga[1,11, Top()], "Gen: 50" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
dens = []
numpars = [50,250,450]
colSchem = cgrad(:haline,rev = true);
numpoints = 128
for i in 1:points
    np = numpars[1]
    Uncx = vec(uncInduDens[1][i])
    Uncy = vec(uncImmMag[1,:,i,1:np,2])
    UncDens = kde((Uncx,Uncy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[1,i],0:1/numpoints:1,0:1/numpoints:1,UncDens./maximum(UncDens), colormap = Reverse(:haline))

    FixRx = vec(FixRInduDens[1][i])
    FixRy = vec(FixRImmMag[1,:,i,1:np,2])
    FixRDens = kde((FixRx,FixRy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[2,i],0:1/numpoints:1,0:1/numpoints:1,FixRDens./maximum(FixRDens), colormap = Reverse(:haline))

    FixUx = vec(FixUInduDens[1][i])
    FixUy = vec(FixUImmMag[1,:,i,1:np,2])
    FixUDens = kde((FixUx,FixUy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[3,i],0:1/numpoints:1,0:1/numpoints:1,FixUDens./maximum(FixUDens), colormap = Reverse(:haline))

    FixDx = vec(FixDInduDens[1][i])
    FixDy = vec(FixDImmMag[1,:,i,1:np,2])
    FixDDens = kde((FixDx,FixDy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[4,i],0:1/numpoints:1,0:1/numpoints:1,FixDDens./maximum(FixDDens), colormap = Reverse(:haline))

    sl100x = vec(sl100InduDens[1][i])
    sl100y = vec(sl100ImmMag[1,:,i,1:np,2])
    sl100Dens = kde((sl100x,sl100y),npoints = (numpoints,numpoints)).density
    heatmap!(ga[5,i],0:1/numpoints:1,0:1/numpoints:1,sl100Dens./maximum(sl100Dens), colormap = Reverse(:haline))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s9_10.eps", f, pt_per_unit = 1)
#endregion

#region 50% infection immune response x magnitude over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1800, 1500))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP \n Peak Immune Effector Activity","RF \n Peak Immune Effector Activity","FU \n Peak Immune Effector Activity","FD \n Peak Immune Effector Activity","Sl \n Peak Immune Effector Activity"]
ga = f[1, 1] = GridLayout()
ga[1:4, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[i]) for i in 1:4]
ga[5,1] = Axis(f,xticks = ([0,1],["C", "I"]), yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[5], xlabel = "Immune Response") 
ga[1:4, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:4 for j in 2:11]
ga[5, 2:11] = [Axis(f,xticks = ([0,1],["C", "I"]),yticklabelsvisible = false,yticksvisible = false, xlabel = "Immune Response") for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "Gen: 0" )
Label(ga[1,2, Top()], "Gen: 5" )
Label(ga[1,3, Top()], "Gen: 10" )
Label(ga[1,4, Top()], "Gen: 15" )
Label(ga[1,5, Top()], "Gen: 20" )
Label(ga[1,6, Top()], "Gen: 25" )
Label(ga[1,7, Top()], "Gen: 30" )
Label(ga[1,8, Top()], "Gen: 35" )
Label(ga[1,9, Top()], "Gen: 40" )
Label(ga[1,10, Top()], "Gen: 45" )
Label(ga[1,11, Top()], "Gen: 50" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
dens = []
numpars = [50,250,450]
colSchem = cgrad(:haline,rev = true);
numpoints = 128
for i in 1:points
    np = numpars[2]
    Uncx = vec(uncInduDens[2][i])
    Uncy = vec(uncImmMag[2,:,i,1:np,2])
    UncDens = kde((Uncx,Uncy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[1,i],0:1/numpoints:1,0:1/numpoints:1,UncDens./maximum(UncDens), colormap = Reverse(:haline))

    FixRx = vec(FixRInduDens[2][i])
    FixRy = vec(FixRImmMag[2,:,i,1:np,2])
    FixRDens = kde((FixRx,FixRy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[2,i],0:1/numpoints:1,0:1/numpoints:1,FixRDens./maximum(FixRDens), colormap = Reverse(:haline))

    FixUx = vec(FixUInduDens[2][i])
    FixUy = vec(FixUImmMag[2,:,i,1:np,2])
    FixUDens = kde((FixUx,FixUy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[3,i],0:1/numpoints:1,0:1/numpoints:1,FixUDens./maximum(FixUDens), colormap = Reverse(:haline))

    FixDx = vec(FixDInduDens[2][i])
    FixDy = vec(FixDImmMag[2,:,i,1:np,2])
    FixDDens = kde((FixDx,FixDy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[4,i],0:1/numpoints:1,0:1/numpoints:1,FixDDens./maximum(FixDDens), colormap = Reverse(:haline))

    sl100x = vec(sl100InduDens[2][i])
    sl100y = vec(sl100ImmMag[2,:,i,1:np,2])
    sl100Dens = kde((sl100x,sl100y),npoints = (numpoints,numpoints)).density
    heatmap!(ga[5,i],0:1/numpoints:1,0:1/numpoints:1,sl100Dens./maximum(sl100Dens), colormap = Reverse(:haline))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s9_50.eps", f, pt_per_unit = 1)
#endregion

#region 90% infection immune response x magnitude over time
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1800, 1500))
lims = (nothing,nothing,nothing,1.2)
ylabels = ["NP \n Peak Immune Effector Activity","RF \n Peak Immune Effector Activity","FU \n Peak Immune Effector Activity","FD \n Peak Immune Effector Activity","Sl \n Peak Immune Effector Activity"]
ga = f[1, 1] = GridLayout()
ga[1:4, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[i]) for i in 1:4]
ga[5,1] = Axis(f,xticks = ([0,1],["C", "I"]), yticks = ([0,1], ["Min.", "Max."]), ylabel = ylabels[5], xlabel = "Immune Response") 
ga[1:4, 2:11] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for i in 1:4 for j in 2:11]
ga[5, 2:11] = [Axis(f,xticks = ([0,1],["C", "I"]),yticklabelsvisible = false,yticksvisible = false, xlabel = "Immune Response") for j in 2:11]
#region label plots
Label(ga[1,1, Top()], "Gen: 0" )
Label(ga[1,2, Top()], "Gen: 5" )
Label(ga[1,3, Top()], "Gen: 10" )
Label(ga[1,4, Top()], "Gen: 15" )
Label(ga[1,5, Top()], "Gen: 20" )
Label(ga[1,6, Top()], "Gen: 25" )
Label(ga[1,7, Top()], "Gen: 30" )
Label(ga[1,8, Top()], "Gen: 35" )
Label(ga[1,9, Top()], "Gen: 40" )
Label(ga[1,10, Top()], "Gen: 45" )
Label(ga[1,11, Top()], "Gen: 50" )
#endregion
colSchem = cgrad(:haline, 100, categorical = true);
dens = []
numpars = [50,250,450]
colSchem = cgrad(:haline,rev = true);
numpoints = 128
for i in 1:points
    np = numpars[3]
    Uncx = vec(uncInduDens[3][i])
    Uncy = vec(uncImmMag[3,:,i,1:np,2])
    UncDens = kde((Uncx,Uncy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[1,i],0:1/numpoints:1,0:1/numpoints:1,UncDens./maximum(UncDens), colormap = Reverse(:haline))

    FixRx = vec(FixRInduDens[3][i])
    FixRy = vec(FixRImmMag[3,:,i,1:np,2])
    FixRDens = kde((FixRx,FixRy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[2,i],0:1/numpoints:1,0:1/numpoints:1,FixRDens./maximum(FixRDens), colormap = Reverse(:haline))

    FixUx = vec(FixUInduDens[3][i])
    FixUy = vec(FixUImmMag[3,:,i,1:np,2])
    FixUDens = kde((FixUx,FixUy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[3,i],0:1/numpoints:1,0:1/numpoints:1,FixUDens./maximum(FixUDens), colormap = Reverse(:haline))

    FixDx = vec(FixDInduDens[3][i])
    FixDy = vec(FixDImmMag[3,:,i,1:np,2])
    FixDDens = kde((FixDx,FixDy),npoints = (numpoints,numpoints)).density
    heatmap!(ga[4,i],0:1/numpoints:1,0:1/numpoints:1,FixDDens./maximum(FixDDens), colormap = Reverse(:haline))

    sl100x = vec(sl100InduDens[3][i])
    sl100y = vec(sl100ImmMag[3,:,i,1:np,2])
    sl100Dens = kde((sl100x,sl100y),npoints = (numpoints,numpoints)).density
    heatmap!(ga[5,i],0:1/numpoints:1,0:1/numpoints:1,sl100Dens./maximum(sl100Dens), colormap = Reverse(:haline))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s9_90.eps", f, pt_per_unit = 1)
#endregion