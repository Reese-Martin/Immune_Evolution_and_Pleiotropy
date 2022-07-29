#packages
using FileIO
using Statistics 
using DataFrames
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using KernelDensity
using Trapz
using HypothesisTests

#variables
const V = 2
const UseCoef = .01
const Infs = 3
const infs = [.1,.5,.9]
const NumHosts = 500
const Runs = 100
const TimePoints = 19

#functions 
function DetPropIndu(x) #x should be a single infection percentage
    tmp = x.+.1*10^-10 #add a small value to each element to prevent 0/0 events
    perConst = round.(tmp[1]/tmp[2],digits = 2) #round removes the error introduced by epsilon addition, gives percent of peak response that was constitutive
    perIndu = 1 .- perConst
    return perIndu
end
function DetLinImmRes(x,Runs,TimePoints)
    RunLinImmRes = []
    for j in 1:Runs
        GenerationalLineage = x[j]
        GenLinImmRes = []
        for k in 1:TimePoints
            UnqHosts = collect(values(GenerationalLineage[k]))
            UnqLinImmRes = []
            for ii in UnqHosts
                push!(UnqLinImmRes,[ii[1],DetPropIndu(ii[2])])
            end
            push!(GenLinImmRes,UnqLinImmRes)
        end
        push!(RunLinImmRes,GenLinImmRes)
    end
    return RunLinImmRes
end
function DetSplitLin(x,Thresh)
    AllDiffs = []
    tmp = x #get lin and immres for last generation of the jth run
    tmp = tmp[(!).(isnan.(getindex.(tmp,2))),:]#filter out hosts with nan values for their perIndu
    Unqs = unique(tmp,dims = 1) #unique lin x ImmRes combinations
    UnqLins = unique(collect(getindex.(Unqs,1))) #just unq lineages
    for i in UnqLins
        SameLins = getindex.(Unqs,1) .== i
        SameLinsImmRes = getindex.(Unqs[SameLins],2)
        Diffs = SameLinsImmRes .- SameLinsImmRes'
        push!(AllDiffs,maximum(Diffs))
    end
    if maximum(AllDiffs) > Thresh #if more lin x ImRes than Lin, then at least one lineage is split
        return 1, maximum(AllDiffs)
    else
        return 0, 0
    end
end
function jitter(a::Array, factor=1.0)
    tmpA = copy(a)
    for itr = 1:length(a)
        tmp = rand()
        if tmp < .5
            tmpA[itr] -= rand() .* factor
        else
            tmpA[itr] += rand() .* factor
        end
    end
    return tmpA
end
function AddMedianLines(x,Data,Runs, row,col)
    tmp = median(Data,dims = 1)
    lines!(ga[row,col],x[1],[tmp[1],tmp[1]], linewidth = 5, color = :Black)
end
function DetMaxSplitLinDistance(x) #for each for each lineage present in the given trajectory snap shot, identify how far apart the representatives of a give lineage are from each other in PorpIndu x ImmMag space
    LinsPres = getindex.(collect(values(x)),1) #find all lineages present
    UnqLins = unique(LinsPres)
    MaxDistances = zeros(size(UnqLins)[1])
    for i in 1:size(UnqLins)[1]
        DiffDescs = LinsPres .== UnqLins[i] #find the descendants of the lineage i
        IR = getindex.(collect(values(x)),2)[DiffDescs] #Immune response of selected lineage i
        Dyn = DetPropIndu.(IR) #percent indu of lineage i
        LinCoords = hcat(Dyn,getindex.(IR,2)) #prop indu x immune magnitude coordinates
        DistMat =  zeros(size(LinCoords)[1],size(LinCoords)[1])
        for j in 1:size(LinCoords)[1]-1
            for k in j+1:size(LinCoords)[1]
                DistMat[j,k] = sqrt((LinCoords[j,1]-LinCoords[k,1])^2+(LinCoords[j,2]-LinCoords[k,2])^2)
            end
        end
        MaxDistances[i] = maximum(DistMat[(!).(isnan.(DistMat))])
    end
    return MaxDistances
end
function DetMeanSplitLinDistance(x) #for each for each lineage present in the given trajectory snap shot, identify how far apart the representatives of a give lineage are from each other in PorpIndu x ImmMag space
    LinsPres = getindex.(collect(values(x)),1) #find all lineages present
    UnqLins = unique(LinsPres)
    MeanDistances = zeros(size(UnqLins)[1])
    for i in 1:size(UnqLins)[1]
        DiffDescs = LinsPres .== UnqLins[i] #find the descendants of the lineage i
        IR = getindex.(collect(values(x)),2)[DiffDescs] #Immune response of selected lineage i
        Dyn = DetPropIndu.(IR) #percent indu of lineage i
        LinCoords = hcat(Dyn,getindex.(IR,2)) #prop indu x immune magnitude coordinates
        DistMat =  zeros(size(LinCoords)[1],size(LinCoords)[1])
        DistMat .= NaN
        for j in 1:size(LinCoords)[1]-1
            for k in j+1:size(LinCoords)[1]
                DistMat[j,k] = sqrt((LinCoords[j,1]-LinCoords[k,1])^2+(LinCoords[j,2]-LinCoords[k,2])^2)
            end
        end
        MeanDistances[i] = mean(DistMat[(!).(isnan.(DistMat))])
    end
    return mean(MeanDistances[(!).(isnan.(MeanDistances))])
end
function nanmean(x)
    mn = zeros(size(x)[2])
    for i in 1:size(x)[2]
        mn[i] = mean(x[(!).(isnan.(x[:,i])),i])
    end
    return mn
end

#region Processing
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"

FixRplotDir = string(topDir,"/Data/FixedRand/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/")
slplotDir = string(topDir,"/Data/SlowEvo/")
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

uncTraj = []
FixRTraj = []
FixUTraj = []
FixDTraj = []
sl100Traj = []

#load data
for i in 1:3 
    push!(uncTraj,FileIO.load(string(FixRplotDir,uncDataSet[i]))["Trajectories"])
    push!(FixRTraj,FileIO.load(string(FixRplotDir,FixRDataSet[i]))["Trajectories"])
    push!(FixUTraj,FileIO.load(string(FixUDplotDir,FixUDataSet[i]))["Trajectories"])
    push!(FixDTraj,FileIO.load(string(FixUDplotDir,FixDDataSet[i]))["Trajectories"])
    push!(sl100Traj,FileIO.load(string(slplotDir,sl100DataSet[i]))["Trajectories"])
end

#determine the euc distance in propIndu x immMag space for all hosts from the same initial lineage, return maximum
UncMaxDist = zeros(Infs, Runs, TimePoints)
FixRMaxDist = zeros(Infs, Runs, TimePoints)
FixUMaxDist = zeros(Infs, Runs, TimePoints)
FixDMaxDist = zeros(Infs, Runs, TimePoints)
sl100MaxDist = zeros(Infs, Runs, TimePoints)
for ii in 1:TimePoints
    for i in 1:Infs
        for j in 1:Runs
            UncMaxDist[i,j,ii] = maximum(DetMaxSplitLinDistance(uncTraj[i][j][ii]))
            FixRMaxDist[i,j,ii] = maximum(DetMaxSplitLinDistance(FixRTraj[i][j][ii]))
            FixUMaxDist[i,j,ii] = maximum(DetMaxSplitLinDistance(FixUTraj[i][j][ii]))
            FixDMaxDist[i,j,ii] = maximum(DetMaxSplitLinDistance(FixDTraj[i][j][ii]))
            sl100MaxDist[i,j,ii] = maximum(DetMaxSplitLinDistance(sl100Traj[i][j][ii]))
        end
    end
end

#determine the euc distance in propIndu x immMag space for all hosts from the same initial lineage, return mean
UncMeanDist = zeros(Infs, Runs, TimePoints)
FixRMeanDist = zeros(Infs, Runs, TimePoints)
FixUMeanDist = zeros(Infs, Runs, TimePoints)
FixDMeanDist = zeros(Infs, Runs, TimePoints)
sl100MeanDist = zeros(Infs, Runs, TimePoints)
for ii in 1:TimePoints
    for i in 1:Infs
        for j in 1:Runs
            UncMeanDist[i,j,ii] = DetMeanSplitLinDistance(uncTraj[i][j][ii])
            FixRMeanDist[i,j,ii] = DetMeanSplitLinDistance(FixRTraj[i][j][ii])
            FixUMeanDist[i,j,ii] = DetMeanSplitLinDistance(FixUTraj[i][j][ii])
            FixDMeanDist[i,j,ii] = DetMeanSplitLinDistance(FixDTraj[i][j][ii])
            sl100MeanDist[i,j,ii] = DetMeanSplitLinDistance(sl100Traj[i][j][ii])
        end
    end
end

#calculate perIndu for each host and extract the lineage of that host
UncLinImmRes = []
FixRLinImmRes = []
FixULinImmRes = []
FixDLinImmRes = []
sl100LinImmRes = []
for i in 1:Infs
    push!(UncLinImmRes,DetLinImmRes(uncTraj[i],Runs,TimePoints))
    push!(FixRLinImmRes,DetLinImmRes(FixRTraj[i],Runs,TimePoints))
    push!(FixULinImmRes,DetLinImmRes(FixUTraj[i],Runs,TimePoints))
    push!(FixDLinImmRes,DetLinImmRes(FixDTraj[i],Runs,TimePoints))
    push!(sl100LinImmRes,DetLinImmRes(sl100Traj[i],Runs,TimePoints))
end

UncSplitLins = zeros(Infs,Runs)
FixRSplitLins = zeros(Infs,Runs)
FixUSplitLins = zeros(Infs,Runs)
FixDSplitLins = zeros(Infs,Runs)
sl100SplitLins = zeros(Infs,Runs)

UncDegSplit = zeros(Infs,Runs)
FixRDegSplit = zeros(Infs,Runs)
FixUDegSplit = zeros(Infs,Runs)
FixDDegSplit = zeros(Infs,Runs)
sl100DegSplit = zeros(Infs,Runs)

for i in 1:Infs
    for j in 1:Runs
        UncSplitLins[i,j], UncDegSplit[i,j] = DetSplitLin(UncLinImmRes[i][j][19],.1)
        FixRSplitLins[i,j], FixRDegSplit[i,j] = DetSplitLin(FixRLinImmRes[i][j][19],.1)
        FixUSplitLins[i,j], FixUDegSplit[i,j] = DetSplitLin(FixULinImmRes[i][j][19],.1)
        FixDSplitLins[i,j], FixDDegSplit[i,j] = DetSplitLin(FixDLinImmRes[i][j][19],.1)
        sl100SplitLins[i,j], sl100DegSplit[i,j] = DetSplitLin(sl100LinImmRes[i][j][19],.1)
    end
end
#endregion

#region Supplemental Figure 10: Figure Showing number of runs with a split lineage
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,100)
ylabels = ["50%","90%"] 
xlabels = ["Non-Pleio.","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10% Chance of Infection \n Number of Split Lineages",title = xlabels[1],limits = lims)
ga[1, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:5]
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "50% Chance of Infection \n Number of Split Lineages",limits = lims)
ga[3, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "90% Chance of Infection \n Number of Split Lineages",limits = lims)
ga[2:3, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:5]
for i in 1:Infs
    barplot!(ga[i,1],[1],[sum(UncSplitLins[i,:])], color = colSchem[1])
    barplot!(ga[i,2],[1],[sum(FixRSplitLins[i,:])], color = colSchem[1])
    barplot!(ga[i,3],[1],[sum(FixUSplitLins[i,:])], color = colSchem[1])
    barplot!(ga[i,4],[1],[sum(FixDSplitLins[i,:])], color = colSchem[1])
    barplot!(ga[i,5],[1],[sum(sl100SplitLins[i,:])], color = colSchem[1])
end
cd(string(topDir,"/Images"))
CairoMakie.save("S10_Split_Lin_Runs.svg", f, pt_per_unit = 1)
#endregion

#region figure showing the difference in degree of inducibility for each split lineage (ie 444 has .1 and .9 hosts, would show as .8)
a = fill(1.0,Runs)
Jit = jitter(a,.3)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ylabels = ["50%","90%"]
xlabels = ["Non-Pleio.","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10%",title = xlabels[1],limits = lims)
ga[1, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:5]
ga[2:3, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i],limits = lims) for i in 1:2]
ga[2:3, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:5]
for i in 1:Infs
    scatter!(ga[i,1],Jit,UncDegSplit[i,:], color = colSchem[1])
    AddMedianLines(xLocs,UncDegSplit[i,:],Runs,i,1)
    scatter!(ga[i,2],Jit,FixRDegSplit[i,:], color = colSchem[1])
    AddMedianLines(xLocs,FixRDegSplit[i,:],Runs,i,2)
    scatter!(ga[i,3],Jit,FixUDegSplit[i,:], color = colSchem[1])
    AddMedianLines(xLocs,FixUDegSplit[i,:],Runs,i,3)
    scatter!(ga[i,4],Jit,FixDDegSplit[i,:], color = colSchem[1])
    AddMedianLines(xLocs,FixDDegSplit[i,:],Runs,i,4)
    scatter!(ga[i,5],Jit,sl100DegSplit[i,:], color = colSchem[1])
    AddMedianLines(xLocs,sl100DegSplit[i,:],Runs,i,5)
end
Label(ga[0,:], "Difference in Inducibility between split lineages")
cd(string(topDir,"/Images/SplitLin"))
GLMakie.save("Split_Lin_Diffs.png", f, pt_per_unit = 1)
#endregion

#region Figure showing the maximum distance between hosts that share a lineage in a given run in the last generation
a = fill(1.0,Runs)
Jit = jitter(a,.3)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ylabels = ["50%","90%"]
xlabels = ["Non-Pleio.","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10%",title = xlabels[1],limits = lims)
ga[1, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:5]
ga[2:3, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i],limits = lims) for i in 1:2]
ga[2:3, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:5]
for i in 1:Infs
    scatter!(ga[i,1],Jit,UncMaxDist[i,:,TimePoints], color = colSchem[1])
    AddMedianLines(xLocs,UncMaxDist[i,:,TimePoints],Runs,i,1)
    scatter!(ga[i,2],Jit,FixRMaxDist[i,:,TimePoints], color = colSchem[1])
    AddMedianLines(xLocs,FixRMaxDist[i,:,TimePoints],Runs,i,2)
    scatter!(ga[i,3],Jit,FixUMaxDist[i,:,TimePoints], color = colSchem[1])
    AddMedianLines(xLocs,FixUMaxDist[i,:,TimePoints],Runs,i,3)
    scatter!(ga[i,4],Jit,FixDMaxDist[i,:,TimePoints], color = colSchem[1])
    AddMedianLines(xLocs,FixDMaxDist[i,:,TimePoints],Runs,i,4)
    scatter!(ga[i,5],Jit,sl100MaxDist[i,:,TimePoints], color = colSchem[1])
    AddMedianLines(xLocs,sl100MaxDist[i,:,TimePoints],Runs,i,5)
end
Label(ga[0,:], "Max distance in split lineage")
cd(string(topDir,"/Images/SplitLin"))
GLMakie.save("Split_Lin_MaxDist.png", f, pt_per_unit = 1)
#endregion

#region Figure showing how the maximum distance between shared lineages changes over time
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1.45)
ylabels = ["50%","90%"]
xlabels = ["Non-Pleio.","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10%",title = xlabels[1],limits = lims)
ga[1, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:5]
ga[2:3, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i],limits = lims) for i in 1:2]
ga[2:3, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:5]
for i in 1:Infs
    lines!(ga[i,1],vec(mean(UncMaxDist[i,:,:], dims = 1)), color = colSchem[1])
    lines!(ga[i,2],vec(mean(FixRMaxDist[i,:,:], dims = 1)), color = colSchem[1])
    lines!(ga[i,3],vec(mean(FixUMaxDist[i,:,:], dims = 1)), color = colSchem[1])
    lines!(ga[i,4],vec(mean(FixDMaxDist[i,:,:], dims = 1)), color = colSchem[1])
    lines!(ga[i,5],vec(mean(sl100MaxDist[i,:,:], dims = 1)), color = colSchem[1])
end
Label(ga[0,:], "Max distance over time")
cd(string(topDir,"/Images/SplitLin"))
GLMakie.save("Split_Lin_MaxDist_overtime.png", f, pt_per_unit = 1)
#endregion

#region Figure showing how the mean distance between shared lineages changes over time
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1.45)
ylabels = ["50%","90%"]
xlabels = ["Non-Pleio.","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10%",title = xlabels[1],limits = lims)
ga[1, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:5]
ga[2:3, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i],limits = lims) for i in 1:2]
ga[2:3, 2:5] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:5]
for i in 1:Infs
    lines!(ga[i,1],vec(nanmean(UncMeanDist[i,:,:])), color = colSchem[1])
    lines!(ga[i,2],vec(nanmean(FixRMeanDist[i,:,:])), color = colSchem[1])
    lines!(ga[i,3],vec(nanmean(FixUMeanDist[i,:,:])), color = colSchem[1])
    lines!(ga[i,4],vec(nanmean(FixDMeanDist[i,:,:])), color = colSchem[1])
    lines!(ga[i,5],vec(nanmean(sl100MeanDist[i,:,:])), color = colSchem[1])
end
Label(ga[0,:], "Mean distance over time")
cd(string(topDir,"/Images/SplitLin"))
GLMakie.save("Split_Lin_MeanDist_overtime.png", f, pt_per_unit = 1)
#endregion

#find all lineages represented across a single run's generations
tmp = []
for i in 3:19
    push!(tmp,unique(getindex.(collect(values(uncTraj[1][1][i])),1)))
end
