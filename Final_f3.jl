#purpose of this script is to determine the resilience (ability to produce the consistent output during disruption) of the most dominant networks from the evolution simulations (non-competition)
#based on work by soyer and salathe with some modification 
#v2 trims down on unneccessary plots and functions
using FileIO
using Statistics 
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using KernelDensity
using Trapz
using HypothesisTests

const V = 2
const UseCoef = .01
const Infs = 3
const Runs = 100

#region Functions
function AcheiveEq(Net, Conc, RemProt, StepLim, InfEq) #modified from original acheive eq to always run until the limit is reached 
    Equilib = false
    Step = 1
    tConc = deepcopy(Conc)
    if RemProt != 0
        tConc[end,RemProt] = 0
    end
    while ~Equilib 
    
        tmpConcDelta = zeros(1,length(Net[1,:]))

        for prot in 1:length(tmpConcDelta)

            InitProtConc = tConc[end,prot]
            ActingProts = Net[:,prot]
            UsedIn = sum(Net[prot,:].!= 0)

            UpregCoefs = ActingProts[findall(ActingProts.>0)]
            UpRegConc = tConc[end,findall(ActingProts.>0)]

            DownRegCoefs = ActingProts[findall(ActingProts.<0)]
            DownProtConc = tConc[end,findall(ActingProts.<0)]

            tmpUp = UpregCoefs.*UpRegConc
            tmpDown = DownRegCoefs.*DownProtConc

            if isempty(tmpUp)
                tmpUp = 0
            end

            Inc = (1-InitProtConc)*sum(tmpUp)

            if isempty(tmpDown)
                tmpDown = 0
            end

            Dec = (InitProtConc)*sum(abs.(tmpDown))

            tmpConcDelta[prot] = Inc - Dec - (UseCoef*UsedIn)

        end

        if RemProt != 0
            tmpConcDelta[RemProt] = 0
        end

        tConc = vcat(tConc,tConc[[end],:]+tmpConcDelta)

        if any(tConc[[end],:].>1)
            tConc[[end],findall(>(1),tConc[end,:])].=1
        end

        if any(tConc[[end],:].<1e-2)
            tConc[[end],findall(<(1e-2),tConc[end,:])].=0
        end
        tConc
        Step += 1

        if InfEq
            if Step >= StepLim
                Equilib = true
            end
        else
            tmpDif = tConc[[end],[end]]-tConc[[end-1],[end]]
            if abs(tmpDif[1])<1e-2
                Equilib = true
            end

            if Step> StepLim
                Equilib = true
            end
        end
    end 
    return tConc
end
function InfectHosts(ParNetworks, HostNetworks, InfHosts)
    #this has been modified to suit single use infected hosts for the Paper_resilience
    HostNet = HostNetworks
    ParNet = ParNetworks
    
    tmp = vcat(HostNet,zeros(length(HostNet[1,:]))')
    tmp[end,1] = 1
    tmp[end,2:end-1] = ParNet[2:end-1]
    parCol = zeros(length(tmp[:,1]))
    parCol[end] = .8
    parCol[end-1] = -1
    tmp = hcat(tmp,parCol) 
    push!(InfHosts, tmp)
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
#resilience as measured by producing the same out put while losing nodes
function DetermineEffectorChange(Net,Limit = 20) #measure of effector production suppression using absolute diffecence between points in effector dynamics 
    Net = Net[1] #this line assumes that the networks come in have a extra array wrapper
    
    tmp = zeros(1,size(Net,2))
    Conc = vcat(tmp.+0,tmp.+.5)
    tmp[1]= 1
    tmp[end] = .8
    Par = tmp
    
    tmpEffScore= zeros(size(Net,1)-2) #each column corresponds to e^-(area of effector with deletion i) - e^-(area of effector with no deletion) should have two fewer columns than total proteins
    tmpConc = AcheiveEq(Net,Conc,0,Limit,false)
    blah = []
    InfectHosts(Par,Net,blah)

    tmpPostInfConc = AcheiveEq(blah[1],hcat(tmpConc[end-1:end,:],[.5;.5]),0,Limit,true)
    IntactEffDynamics = tmpPostInfConc[:,end-1]

    for j in 1:size(Net,1)-2
        tmpInfConc = hcat(tmpConc[end-1:end,:],[.5;.5])
        tmpPostInfConc = AcheiveEq(blah[1],tmpInfConc,j+1,Limit,true)
        tmpDelEffDynamics = tmpPostInfConc[:,end-1]
        Diff = abs.(tmpDelEffDynamics - IntactEffDynamics)
        tmpEffScore[j] = mean(Diff)
    end
    return tmpEffScore #0 is identical behavior, gets exponentially worse as number increases maxes out at .63 ( that is one area = 1 and the other = 0)
end
function AddMedianLines(x,Data,Runs, row,col)
    tmp = median(Data,dims = 1)
    lines!(ga[row,col],x[1],[tmp[1],tmp[1]], linewidth = 5, color = :Black)
end
function RemNans(a)
    return a[(!isnan).(a)]
end
#endregion
#region Load Data 
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"

FixRplotDir = string(topDir,"/Data/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/Draft_2_Data/")
slplotDir = string(topDir,"/Data/SlowEvo/Draft_2_Data/")
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> contains(x, "Network"), FixRFiles)
cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> contains(x, "Network"), FixUDFiles)
cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> contains(x, "Network"), slFiles)

#group data from saved files
uncDataSet = filter(x -> contains(x, "Uncon"), FixRFiles);
FixRDataSet = filter(x -> contains(x, "Con"), FixRFiles);
FixUDataSet = filter(x -> contains(x, "upreg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "downreg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

uncNets = []
FixRNets = []
FixUNets = []
FixDNets = []
sl100Nets = []
#load data
for i in 1:Infs 
    push!(uncNets,FileIO.load(string(FixRplotDir,uncDataSet[i]))["Networks"])
    push!(FixRNets,FileIO.load(string(FixRplotDir,FixRDataSet[i]))["Networks"])
    push!(FixUNets,FileIO.load(string(FixUDplotDir,FixUDataSet[i]))["Networks"])
    push!(FixDNets,FileIO.load(string(FixUDplotDir,FixDDataSet[i]))["Networks"])
    push!(sl100Nets,FileIO.load(string(slplotDir,sl100DataSet[i]))["Networks"])
end
#endregion
#region process data
uncEC = []
FixREC = []
FixUEC = []
FixDEC = []
sl100EC = []
#mean absolute difference following deletion
for i in 1:Infs
    push!(uncEC,DetermineEffectorChange.(uncNets[i]))
    push!(FixREC,DetermineEffectorChange.(FixRNets[i]))
    push!(FixUEC,DetermineEffectorChange.(FixUNets[i]))
    push!(FixDEC,DetermineEffectorChange.(FixDNets[i]))
    push!(sl100EC,DetermineEffectorChange.(sl100Nets[i]))
end
#endregion

#region figure 3
a = fill(1.0,Runs)
Jit = jitter(a,.3)
colSchem = cgrad(:haline, 100, categorical = true);

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1500))
lims = (nothing,nothing,0,1)
#ylabels = ["50%","90%"]
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "10% Chance of Infection \n Effector Mean Absolute Difference",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "50% Chance of Infection \n Effector Mean Absolute Difference")
ga[3, 1] = Axis(f, limits = lims, xticks = ([1,2],["Plei.","NP."]), xgridvisible = false, ygridvisible = false, yticks = ([0,1],["0%", "100%"]),ylabel = "10% Chance of Infection \n Effector Mean Absolute Difference")
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, xgridvisible = false, ygridvisible = false, yticklabelsvisible = false,yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticks =  ([1,2],["Plei.","NP"]), xgridvisible = false, ygridvisible = false, yticklabelsvisible = false, yticksvisible = false, xlabel = "Percent of Response Induced") for j in 2:5]
Variances = zeros(3,5,2)
Ttests = zeros(3,5)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
xLocs2 = []
for i in 1:5
    push!(xLocs2,[1.5,2.5])
end
for i in 1:Infs
    tmp1P = getindex.(uncEC[i],1)
    tmp2P = getindex.(FixREC[i],1)
    tmp3P = getindex.(FixUEC[i],1)
    tmp4P = getindex.(FixDEC[i],1)
    tmp5P = getindex.(sl100EC[i],1)

    tmp1NP = Vector{Float64}()
    tmp2NP = Vector{Float64}()
    tmp3NP = Vector{Float64}()
    tmp4NP = Vector{Float64}()
    tmp5NP = Vector{Float64}()

    for j in 1:Runs
        append!(tmp1NP, uncEC[i][j][2:end])
        append!(tmp2NP, FixREC[i][j][2:end])
        append!(tmp3NP, FixUEC[i][j][2:end])
        append!(tmp4NP, FixDEC[i][j][2:end])
        append!(tmp5NP, sl100EC[i][j][2:end])
    end

    Variances[i,1,1] = var(tmp1P)
    Variances[i,2,1] = var(tmp2P)
    Variances[i,3,1] = var(tmp3P)
    Variances[i,4,1] = var(tmp4P)
    Variances[i,5,1] = var(tmp5P)

    Variances[i,1,2] = var(tmp1NP)
    Variances[i,2,2] = var(tmp2NP)
    Variances[i,3,2] = var(tmp3NP)
    Variances[i,4,2] = var(tmp4NP)
    Variances[i,5,2] = var(tmp5NP)

    Ttests[i,1] = pvalue(EqualVarianceTTest(RemNans(tmp1P),RemNans(tmp1NP)))
    Ttests[i,2] = pvalue(EqualVarianceTTest(RemNans(tmp2P),RemNans(tmp2NP)))
    Ttests[i,3] = pvalue(EqualVarianceTTest(RemNans(tmp3P),RemNans(tmp3NP)))
    Ttests[i,4] = pvalue(EqualVarianceTTest(RemNans(tmp4P),RemNans(tmp4NP)))
    Ttests[i,5] = pvalue(EqualVarianceTTest(RemNans(tmp5P),RemNans(tmp5NP)))

    scatter!(ga[i,1],Jit,tmp1P, color = colSchem[1])
    scatter!(ga[i,1],jitter(fill(2.0,length(tmp1NP)),.3),tmp1NP, color = colSchem[70])
    AddMedianLines(xLocs,tmp1P,Runs,i,1)
    AddMedianLines(xLocs2,tmp1NP,length(tmp1NP),i,1)

    scatter!(ga[i,2],Jit,tmp2P, color = colSchem[1])
    scatter!(ga[i,2],jitter(fill(2.0,length(tmp2NP)),.3),tmp2NP, color = colSchem[70])
    AddMedianLines(xLocs,tmp2P,Runs,i,2)
    AddMedianLines(xLocs2,tmp2NP,length(tmp2NP),i,2)

    scatter!(ga[i,3],Jit,tmp3P, color = colSchem[1])
    scatter!(ga[i,3],jitter(fill(2.0,length(tmp3NP)),.3),tmp3NP, color = colSchem[70])
    AddMedianLines(xLocs,tmp3P,Runs,i,3)
    AddMedianLines(xLocs2,tmp3NP,length(tmp3NP),i,3)

    scatter!(ga[i,4],Jit,tmp4P, color = colSchem[1])
    scatter!(ga[i,4],jitter(fill(2.0,length(tmp4NP)),.3),tmp4NP, color = colSchem[70])
    AddMedianLines(xLocs,tmp4P,Runs,i,4)
    AddMedianLines(xLocs2,tmp4NP,length(tmp4NP),i,4)

    scatter!(ga[i,5],Jit,tmp5P, color = colSchem[1])
    scatter!(ga[i,5],jitter(fill(2.0,length(tmp5NP)),.3),tmp5NP, color = colSchem[70])
    AddMedianLines(xLocs,tmp5P,Runs,i,5)
    AddMedianLines(xLocs2,tmp5NP,length(tmp5NP),i,5)
end 
Sigs = findall(x -> x<.05/15,Ttests)
for i in Sigs
    text!(ga[i[1],i[2]],"*",position = (1.4,.8),textsize = 50)
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("f3_Pleio_NonPleioKODiff.svg", f, pt_per_unit = 1)
#endregion
