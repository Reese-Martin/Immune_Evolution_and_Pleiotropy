using FileIO
using Graphs
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using GraphMakie #
Gens = 500
Runs = 100
Infs = 3
Hosts = 500
ImmuneTypes = 4 #0- no response, 1- constitutive, 2- Inducible, 3- mixed indu/const
UseCoef = .01

#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Co-evo model/Julia Code/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"

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
function GenAvgHost(x)
    MaxProts = maximum(size.(x,1))
    Det = zeros(1,MaxProts)
    Sigs = zeros(MaxProts-2,MaxProts)
    Eff = zeros(1,MaxProts)

    for i in x
        Det[2:size(i,1)-1] += i[1,2:end-1]
        Eff[2:size(i,1)-1] += i[end,2:end-1]
        Sigs[1:size(i,1)-2,1:size(i,1)-1] += i[2:end-1,1:end-1] 
        Sigs[1:size(i,1)-2,end] += i[2:end-1,end]
    end
    AllNets = vcat(Det,Sigs,Eff)
    ScaledNet = AllNets./(maximum(abs.(AllNets)))
    return ScaledNet
end
function PlotAvgHost(Loc,x)
    NumProts = size(x,1)
    NumSigs = size(x,1)-2
    DetPosit = [2,ceil(Int64,NumSigs/3)+2]'
    SigPosits = zeros(NumSigs,2)
    c1 = ceil(Int64,NumSigs/3)+1
    c2 = 1
    for i in 1:NumSigs
        SigPosits[i,1] = c2
        SigPosits[i,2] = c1
        c2 += 1
        if c2 > 3
            c2 = 1
            c1 -= 1
        end
    end
    EffPosit = [2,1]'
    jitteredSigs = jitter(SigPosits,0)
    Points = vcat(DetPosit,jitteredSigs,EffPosit)

    Plah(_) = [(Points[i,1],Points[i,2]) for i in 1:NumProts]
    g = SimpleDiGraph(NumProts)
    Cons = findall(x -> x != 0, x)
    edgecolors = []

    for i in Cons
        add_edge!(g,i[1],i[2])
    end
    for i in collect(edges(g))
        sor = i.src 
        dest = i.dst
        Reg = findall(x -> x[1]==sor , Cons)
        SubCons = Cons[Reg]
        Reg = findall(x -> x[2] == dest, SubCons)
        Connection = SubCons[Reg][1]

        if x[Connection]>0
            push!(edgecolors,(:red,x[Connection]))
        else
            push!(edgecolors,(:blue,abs(x[Connection])))
        end
    end
    Markers = vcat(:rect, fill(:circle,NumSigs), :utriangle)
    labels =  repr.(1:ne(g))
    p = graphplot!(ga[Loc[1],Loc[2]],g, edge_color = edgecolors, node_marker = Markers,
    node_size = 30)
    p.layout = Plah
    p.edge_width = 4
end
function CalcPerIndu(Net, Limit = 20) #determine how inducible a given network is with non-interfering parasite
    #Net = Net[1] #this line assumes that the networks come in have a extra array wrapper
    
    tmp = zeros(1,size(Net,2))
    Conc = vcat(tmp.+0,tmp.+.5)
    tmp[1]= 1
    tmp[end] = .8
    Par = tmp
    
    tmpEffScore= zeros(size(Net,1)-2) #each column corresponds to e^-(area of effector with deletion i) - e^-(area of effector with no deletion) should have two fewer columns than total proteins
    tmpConc = AcheiveEq(Net,Conc,0,Limit,false)
    Pre_Eq = tmpConc[end,end]
    blah = []
    InfectHosts(Par,Net,blah)
    tmpPostInfConc = AcheiveEq(blah[1],hcat(tmpConc[end-1:end,:],[.5;.5]),0,Limit,true)
    MaxEff = maximum(tmpPostInfConc[:,end-1])

    return round(1-(Pre_Eq+1e-10)/MaxEff,digits = 2) #0 is identical behavior, gets exponentially worse as number increases maxes out at .63 ( that is one area = 1 and the other = 0)
end
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
#region Processing
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
 
#Direct Connections to Effector 
uncNets = []
FixRNets = []
FixUNets = []
FixDNets = []
sl100Nets = []

for i = 1:Infs
    unctmp = getindex.(FileIO.load(string(FixRplotDir,uncDataSet[i]))["Networks"],1)
    FixRtmp = getindex.(FileIO.load(string(FixRplotDir,FixRDataSet[i]))["Networks"],1)
    FixUtmp = getindex.(FileIO.load(string(FixUDplotDir,FixUDataSet[i]))["Networks"],1)
    FixDtmp = getindex.(FileIO.load(string(FixUDplotDir,FixDDataSet[i]))["Networks"],1)
    sl100tmp = getindex.(FileIO.load(string(slplotDir,sl100DataSet[i]))["Networks"],1)

    push!(uncNets, unctmp)
    push!(FixRNets, FixRtmp)
    push!(FixUNets, FixUtmp)
    push!(FixDNets, FixDtmp)
    push!(sl100Nets, sl100tmp)
    
    unctmpDict = nothing
    FixRtmpDict = nothing
    FixUtmpDict = nothing
    FixDtmpDict = nothing
    sl10tmpDict = nothing
    sl100tmpDict = nothing
end

AvgUncHost = []
AvgFixRHost = []
AvgFixUHost = []
AvgFixDHost = []
AvgSl100Host = []
for i in 1:Infs
    push!(AvgUncHost,GenAvgHost(uncNets[i]))
    push!(AvgFixRHost,GenAvgHost(FixRNets[i]))
    push!(AvgFixUHost,GenAvgHost(FixUNets[i]))
    push!(AvgFixDHost,GenAvgHost(FixDNets[i]))
    push!(AvgSl100Host,GenAvgHost(sl100Nets[i]))
end

uncConstHosts = []
uncInduHosts = []
FixDConstHosts = []
FixDInduHosts = []

uncPerIndu = collect(CalcPerIndu.(uncNets[3]))
FixDPerIndu = collect(CalcPerIndu.(FixDNets[3]))

for i in 1:Runs
    if uncPerIndu[i] > .9
        push!(uncInduHosts, uncNets[3][i])
    elseif uncPerIndu[i]<.1
        push!(uncConstHosts, uncNets[3][i])
    end

    if FixDPerIndu[i] > .9
        push!(FixDInduHosts, FixDNets[3][i])
    elseif FixDPerIndu[i]<.1
        push!(FixDConstHosts, FixDNets[3][i])
    end
end
AvgUncConstHost = []
AvgUncInduHost = []

AvgFixDConstHost = []
AvgFixDInduHost = []

push!(AvgUncConstHost,GenAvgHost(uncConstHosts))
push!(AvgUncInduHost,GenAvgHost(uncInduHosts))
push!(AvgFixDConstHost,GenAvgHost(FixDConstHosts))
push!(AvgFixDInduHost,GenAvgHost(FixDInduHosts))

#endregion
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,4,0,6)
#ylabels = ["50%","90%"]
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]

ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, ylabel = "10% Chance of Infection",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, ylabel = "50% Chance of Infection",)
ga[3, 1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, ylabel = "90% Chance of Infection",)
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
ga[2, 6] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, title = "Non-Pleiotropic \n Constitutive")
ga[2, 7] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, title = "Non-Pleiotropic \n Inducible")
ga[3, 6] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, title = "Downregulatory \n Constitutive")
ga[3, 7] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false, title = "Downregulatory \n Inducible")
for i in 1:Infs
    PlotAvgHost([i,1], AvgUncHost[i])
    PlotAvgHost([i,2], AvgFixRHost[i])
    PlotAvgHost([i,3], AvgFixUHost[i])
    PlotAvgHost([i,4], AvgFixDHost[i])
    PlotAvgHost([i,5], AvgSl100Host[i])
end

PlotAvgHost([2,6], AvgUncConstHost[1])
PlotAvgHost([2,7], AvgUncInduHost[1])
PlotAvgHost([3,6], AvgFixDConstHost[1])
PlotAvgHost([3,7], AvgFixDInduHost[1])
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s8_avg_Host_net.svg", f, pt_per_unit = 1)