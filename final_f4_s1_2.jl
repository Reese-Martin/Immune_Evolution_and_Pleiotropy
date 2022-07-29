using FileIO
using Statistics
using Fontconfig
using DataFrames
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using KernelDensity
ToSave = true    
FixGens = 1000
Runs = 100
Infs = 3
Hosts = 500
ImmuneTypes = 4 #0- no response, 1- constitutive, 2- Inducible, 3- mixed indu/const
UseCoef = .01
EQGens = 1

#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Immune_Pleio_Model/Julia Code"
FixRplotDir = string(topDir,"/Data/Competition/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/Competition/FixedUD/Draft_2_Data/")
slplotDir = string(topDir,"/Data/Competition/SlowEvo/Draft_2_Data/")

#region Functions for network analysis
function AcheiveEq(Net, Conc, StepLim, InfEq)
    Equilib = false
    Step = 1
    
    while ~Equilib 
    
        tmpConcDelta = zeros(1,length(Net[1,:]))

        for prot in 1:length(tmpConcDelta)

            InitProtConc = Conc[end,prot]
            ActingProts = Net[:,prot]
            UsedIn = sum(Net[prot,:].!= 0)

            UpregCoefs = ActingProts[findall(ActingProts.>0)]
            UpRegConc = Conc[end,findall(ActingProts.>0)]

            DownProtConc = ActingProts[findall(ActingProts.<0)]
            DownRegCoefs = Conc[end,findall(ActingProts.<0)]

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

        Conc = vcat(Conc,Conc[[end],:]+tmpConcDelta)
        #Conc = vcat(Conc,round.(Conc[[end],:]+tmpConcDelta,digits = 2))

        if any(Conc[[end],:].>1)
            Conc[[end],findall(>(1),Conc[end,:])].=1
        end

        if any(Conc[[end],:].<1e-2)
            Conc[[end],findall(<(1e-2),Conc[end,:])].=0
        end

        Step += 1

        if InfEq
            if Conc[end,end] <= 1e-2
                Equilib = true
            elseif Step == StepLim
                Equilib = true
            end
        else
            tmpDif = Conc[[end],[end]]-Conc[[end-1],[end]]
            if abs(tmpDif[1])<1e-2
                Equilib = true
            end

            if Step> StepLim
                Equilib = true
            end
        end
    end 
    return Conc
end
function InfectHosts(ParNetworks, HostNetworks, InfHosts)
    par = 1
    HostNet = HostNetworks
    ParNet = ParNetworks[par]
    ParTarg = findall(x->x != 0,ParNet)
    if length(ParTarg) == 1
        ParTarg = ParTarg[1]
    end

    if ParTarg in 2:length(HostNet[1,:])-1
        #if length(ParNet) < length(HostNet[1,:])
            #   zertoadd = length(HostNet[1,:]) - length(ParNet)
            #   for q in 1:zertoadd
        #        ParNet = vcat(ParNet[1:ParTarg],0,ParNet[ParTarg+1:end])
        #    end
        # elseif length(ParNet) > length(HostNet[1,:])

        #end
        tmp = vcat(HostNet,zeros(length(HostNet[1,:]))')
        tmp[end,ParTarg] = ParNet[ParTarg]
        tmp[end,1] = 1

        parCol = zeros(length(tmp[:,1]))
        parCol[end] = .8
        parCol[end-1] = -1
        tmp = hcat(tmp,parCol) 
    else    
        ParTarg = 2
        parRow = zeros(length(HostNet[1,:]))
        parRow[ParTarg] = 0
        parRow[1] = 1
        tmp = vcat(HostNet,parRow')

        parCol = zeros(length(tmp[:,1]))
        parCol[end] = .8
        parCol[end-1] = -1
        tmp = hcat(tmp,parCol)
    end

    push!(InfHosts, tmp)
end
function MakeFigCompBarPlot(size,lims = (nothing, nothing, nothing, nothing))
    fig = Figure(resolution = size)

    ax1 = Axis(fig[1,1],xticklabelsvisible = false,xticksvisible = false, title = "Fixed Random", ylabel = "10%",limits = lims)
    ax2 = Axis(fig[2,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "50%",limits = lims)
    ax3 = Axis(fig[3,1], ylabel = "90%",xticklabelsvisible = false,xticksvisible = false, xlabel = "Unevolved",limits = lims)

    ax4 = Axis(fig[1,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Up",limits = lims)
    ax5 = Axis(fig[2,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax6 = Axis(fig[3,2],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Unevolved",limits = lims)

    ax7 = Axis(fig[1,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Down",limits = lims)
    ax8 = Axis(fig[2,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax9 = Axis(fig[3,3],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Unevolved",limits = lims)

    ax10 = Axis(fig[1,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "100x Slower",limits = lims)
    ax11 = Axis(fig[2,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax12 = Axis(fig[3,4],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Unevolved",limits = lims)

    ax13 = Axis(fig[1,5],xticklabelsvisible = false,xticksvisible = false, title = "Fixed Random",yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax14 = Axis(fig[2,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax15 = Axis(fig[3,5],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Evolved",limits = lims)

    ax16 = Axis(fig[1,6],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Up",limits = lims)
    ax17 = Axis(fig[2,6],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax18 = Axis(fig[3,6],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Evolved",limits = lims)

    ax19 = Axis(fig[1,7],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Down",limits = lims)
    ax20 = Axis(fig[2,7],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax21 = Axis(fig[3,7],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Evolved",limits = lims)

    ax22 = Axis(fig[1,8],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "100x Slower",limits = lims)
    ax23 = Axis(fig[2,8],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax24 = Axis(fig[3,8],yticklabelsvisible = false, yticksvisible = false,xticklabelsvisible = false,xticksvisible = false, xlabel = "Evolved",limits = lims)

    return fig
end
function SortNets(LoadDir,PWinNets,UWinNets,PLoseNets,ULoseNets)
    tmp = FileIO.load(LoadDir)
    Winners = tmp["WinInfRun"]
    Winners = Winners[Winners[:,1] .!=0,:]
    PWinners = tmp["WinNets"][Winners[:,3] .== 1, :]
    push!(PWinNets,PWinners)
    UWinners = tmp["WinNets"][Winners[:,3] .== 0, :]
    push!(UWinNets,UWinners)

    Losers = tmp["LoseInfRun"]
    Losers = Losers[Losers[:,1] .!=0,:]
    PLosers = tmp["LastLoser"][Losers[:,3] .==1, :]
    push!(PLoseNets,PLosers)
    ULosers = tmp["LastLoser"][Losers[:,3] .==0, :]
    push!(ULoseNets,ULosers)
end
function DetermineImmuneResponse(Nets,pars,Concs,Limit)
    ImmMag = zeros(length(Nets),2)
    for i in 1:length(Nets)
        tmpConc = AcheiveEq(Nets[i],Concs[i],Limit,false)
        blah = []
        InfectHosts(pars[i],Nets[i],blah)
        tmpInfConc = hcat(tmpConc[end-1:end,:],[.5;.5])
        tmpPostInfConc = AcheiveEq(blah[1],tmpInfConc,Limit,true)
        ImmMag[i,1] = tmpPostInfConc[2,end-1]
        ImmMag[i,2] = maximum(tmpPostInfConc[:,end-1])
    end
    return ImmMag
end
function DensityofImmuneResponses(networks)
    tmpPars = []
    tmpConcs = []
    for j in networks
        tmp = zeros(1,size(j,2))
        push!(tmpConcs,vcat(tmp.+0,tmp.+.5))
        tmp[1]= 1
        push!(tmpPars,tmp)
    end
    ImmMags = DetermineImmuneResponse(networks,tmpPars,tmpConcs,20)
    ImmMags.+= 1*10^-6
    PropIndu = 1 .- (ImmMags[:,1]./ImmMags[:,2])
    return kde(PropIndu, npoints = 128).density
end
function MakeFigStackDens(size,lims = (nothing, nothing, nothing, nothing))
    labels = ["None" "Constitutive" "Mixed" "Inducible"]
    fig = Figure(resolution = size)

    ax1 = Axis(fig[1,1],xticklabelsvisible = false,xticksvisible = false, title = "FixR", ylabel = "10%",limits = lims)
    ax2 = Axis(fig[2,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "50%",limits = lims)
    ax3 = Axis(fig[3,1], ylabel = "90%",limits = lims,xticks = ([1,128],["Const", "Indu"]))

    ax4 = Axis(fig[1,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "FixU",limits = lims)
    ax5 = Axis(fig[2,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax6 = Axis(fig[3,2],yticklabelsvisible = false, yticksvisible = false,limits = lims,xticks = ([1,128],["Const", "Indu"]))

    ax7 = Axis(fig[1,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "FixD",limits = lims)
    ax8 = Axis(fig[2,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax9 = Axis(fig[3,3],yticklabelsvisible = false, yticksvisible = false,limits = lims,xticks = ([1,128],["Const", "Indu"]))

    ax10 = Axis(fig[1,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "SLx100",limits = lims)
    ax11 = Axis(fig[2,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax12 = Axis(fig[3,4],yticklabelsvisible = false, yticksvisible = false,limits = lims,xticks = ([1,128],["Const", "Indu"]))
    return fig
end
#endregion

#region processing
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> !contains(x, "Network"), FixRFiles)

cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> !contains(x, "Network"), FixUDFiles)

cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> !contains(x, "Network"), slFiles)

FixRFiles250 = filter(x -> contains(x, "250"), FixRFiles)
FixUDFiles250 = filter(x -> contains(x, "250"), FixUDFiles)
slFiles250 = filter(x -> contains(x, "250"), slFiles)

FixRFiles = filter(x -> !contains(x, "250"), FixRFiles)
FixUDFiles = filter(x -> !contains(x, "250"), FixUDFiles)
slFiles = filter(x -> !contains(x, "250"), slFiles)

#group data from saved files
FixRDataSet = FixRFiles;
FixUDataSet = filter(x -> contains(x, "UpReg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "DownReg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

FixRDataSet250 = FixRFiles250;
FixUDataSet250 = filter(x -> contains(x, "UpReg"), FixUDFiles250);
FixDDataSet250 = filter(x -> contains(x, "DownReg"), FixUDFiles250);
sl100DataSet250 = filter(x -> contains(x, "100X"), slFiles250);

FixRPrev = zeros(Infs,Runs, FixGens, 2); #first 'sheet' is constrained, second is unconstrained 
FixUPrev = zeros(Infs,Runs, FixGens, 2);
FixDPrev = zeros(Infs,Runs, FixGens, 2);
sl100Prev = zeros(Infs,Runs, FixGens, 2);

FixRPrev250 = zeros(Infs,Runs, FixGens, 2); #first 'sheet' is constrained, second is unconstrained 
FixUPrev250 = zeros(Infs,Runs, FixGens, 2);
FixDPrev250 = zeros(Infs,Runs, FixGens, 2);
sl100Prev250 = zeros(Infs,Runs, FixGens, 2);

for i = 1:Infs
    FixRtmpDict = FileIO.load(string(FixRplotDir,FixRDataSet[i]))
    FixUtmpDict = FileIO.load(string(FixUDplotDir,FixUDataSet[i]))
    FixDtmpDict = FileIO.load(string(FixUDplotDir,FixDDataSet[i]))
    sl100tmpDict = FileIO.load(string(slplotDir,sl100DataSet[i]))

    FixRPrev[i,:,:,:] = FixRtmpDict["PopulationDistribution"]
    FixUPrev[i,:,:,:] = FixUtmpDict["PopulationDistribution"]
    FixDPrev[i,:,:,:] = FixDtmpDict["PopulationDistribution"]
    sl100Prev[i,:,:,:] = sl100tmpDict["PopulationDistribution"]

    FixRtmpDict = nothing
    FixUtmpDict = nothing
    FixDtmpDict = nothing
    sl100tmpDict = nothing

    FixRtmpDict250 = FileIO.load(string(FixRplotDir,FixRDataSet250[i]))
    FixUtmpDict250 = FileIO.load(string(FixUDplotDir,FixUDataSet250[i]))
    FixDtmpDict250 = FileIO.load(string(FixUDplotDir,FixDDataSet250[i]))
    sl100tmpDict250 = FileIO.load(string(slplotDir,sl100DataSet250[i]))

    FixRPrev250[i,:,:,:] = FixRtmpDict250["PopulationDistribution"]
    FixUPrev250[i,:,:,:] = FixUtmpDict250["PopulationDistribution"]
    FixDPrev250[i,:,:,:] = FixDtmpDict250["PopulationDistribution"]
    sl100Prev250[i,:,:,:] = sl100tmpDict250["PopulationDistribution"]

    FixRtmpDict250 = nothing
    FixUtmpDict250 = nothing
    FixDtmpDict250 = nothing
    sl100tmpDict250 = nothing
end

#set pleio/conc prevalence for all Generations
for j in 1:Infs
    for i in 1:Runs
        FixREndInd = findall(x -> x == Hosts,FixRPrev[j,i,:,:])
        FixUEndInd = findall(x -> x == Hosts,FixUPrev[j,i,:,:])
        FixDEndInd = findall(x -> x == Hosts,FixDPrev[j,i,:,:])
        sl100EndInd = findall(x -> x == Hosts,sl100Prev[j,i,:,:])
        
        if ~isempty(FixREndInd)
            FixRPrev[j,i,FixREndInd[1][1]:end,FixREndInd[1][2]].=Hosts 
        end
        if ~isempty(FixUEndInd)
            FixUPrev[j,i,FixUEndInd[1][1]:end,FixUEndInd[1][2]].=Hosts 
        end
        if ~isempty(FixDEndInd)
            FixDPrev[j,i,FixDEndInd[1][1]:end,FixDEndInd[1][2]].=Hosts 
        end
        if ~isempty(sl100EndInd)
            sl100Prev[j,i,sl100EndInd[1][1]:end,sl100EndInd[1][2]].=Hosts 
        end

        FixREndInd250 = findall(x -> x == Hosts,FixRPrev250[j,i,:,:])
        FixUEndInd250 = findall(x -> x == Hosts,FixUPrev250[j,i,:,:])
        FixDEndInd250 = findall(x -> x == Hosts,FixDPrev250[j,i,:,:])
        sl100EndInd250 = findall(x -> x == Hosts,sl100Prev250[j,i,:,:])
        
        if ~isempty(FixREndInd250)
            FixRPrev250[j,i,FixREndInd250[1][1]:end,FixREndInd250[1][2]].=Hosts 
        end
        if ~isempty(FixUEndInd250)
            FixUPrev250[j,i,FixUEndInd250[1][1]:end,FixUEndInd250[1][2]].=Hosts 
        end
        if ~isempty(FixDEndInd250)
            FixDPrev250[j,i,FixDEndInd250[1][1]:end,FixDEndInd250[1][2]].=Hosts 
        end
        if ~isempty(sl100EndInd250)
            sl100Prev250[j,i,sl100EndInd250[1][1]:end,sl100EndInd250[1][2]].=Hosts 
        end
    end
end

FixRFixedPleio = zeros(3,3)
FixUFixedPleio = zeros(3,3)
FixDFixedPleio = zeros(3,3)
sl100FixedPleio = zeros(3,3)

FixRFixedPleio250 = zeros(3,3)
FixUFixedPleio250 = zeros(3,3)
FixDFixedPleio250 = zeros(3,3)
sl100FixedPleio250 = zeros(3,3)

for i in 1:Infs
    FixRPcount = findall(x -> x == Hosts,FixRPrev[i,:,end,1])
    FixUPcount = findall(x -> x == Hosts,FixUPrev[i,:,end,1])
    FixDPcount = findall(x -> x == Hosts,FixDPrev[i,:,end,1])
    sl100Pcount = findall(x -> x == Hosts,sl100Prev[i,:,end,1])

    FixRUcount = findall(x -> x == Hosts,FixRPrev[i,:,end,2])
    FixUUcount = findall(x -> x == Hosts,FixUPrev[i,:,end,2])
    FixDUcount = findall(x -> x == Hosts,FixDPrev[i,:,end,2])
    sl100Ucount = findall(x -> x == Hosts,sl100Prev[i,:,end,2])

    FixRFixedPleio[i,:] = [size(FixRPcount)[1],size(FixRUcount)[1],(Runs-(size(FixRPcount)[1]+size(FixRUcount)[1]))]./Runs
    FixUFixedPleio[i,:] = [size(FixUPcount)[1],size(FixUUcount)[1],(Runs-(size(FixRPcount)[1]+size(FixRUcount)[1]))]./Runs
    FixDFixedPleio[i,:]= [size(FixDPcount)[1],size(FixDUcount)[1],(Runs-(size(FixDPcount)[1]+size(FixDUcount)[1]))]./Runs
    sl100FixedPleio[i,:] = [size(sl100Pcount)[1],size(sl100Ucount)[1],(Runs-(size(sl100Pcount)[1]+size(sl100Ucount)[1]))]./Runs

    FixRPcount250 = findall(x -> x == Hosts,FixRPrev250[i,:,end,1])
    FixUPcount250 = findall(x -> x == Hosts,FixUPrev250[i,:,end,1])
    FixDPcount250 = findall(x -> x == Hosts,FixDPrev250[i,:,end,1])
    sl100Pcount250 = findall(x -> x == Hosts,sl100Prev250[i,:,end,1])

    FixRUcount250 = findall(x -> x == Hosts,FixRPrev250[i,:,end,2])
    FixUUcount250 = findall(x -> x == Hosts,FixUPrev250[i,:,end,2])
    FixDUcount250 = findall(x -> x == Hosts,FixDPrev250[i,:,end,2])
    sl100Ucount250 = findall(x -> x == Hosts,sl100Prev250[i,:,end,2])

    FixRFixedPleio250[i,:] = [size(FixRPcount250)[1],size(FixRUcount250)[1],(Runs-(size(FixRPcount250)[1]+size(FixRUcount250)[1]))]./Runs
    FixUFixedPleio250[i,:] = [size(FixUPcount250)[1],size(FixUUcount250)[1],(Runs-(size(FixRPcount250)[1]+size(FixRUcount250)[1]))]./Runs
    FixDFixedPleio250[i,:]= [size(FixDPcount250)[1],size(FixDUcount250)[1],(Runs-(size(FixDPcount250)[1]+size(FixDUcount250)[1]))]./Runs
    sl100FixedPleio250[i,:] = [size(sl100Pcount250)[1],size(sl100Ucount250)[1],(Runs-(size(sl100Pcount250)[1]+size(sl100Ucount250)[1]))]./Runs
end

#endregion

#region Supplemental Figure 1
colSchem = cgrad(:haline, 10, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ylabels = ["50% Chance of Infection ","90% Chance of Infection"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = "10% Chance of Infection",title = xlabels[1],limits = lims)
ga[1, 2:8] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i],limits = lims) for i in 2:8]
ga[2:3, 1] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,ylabel = ylabels[i],limits = lims) for i in 1:2]
ga[2:3, 2:8] = [Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,limits = lims) for i in 2:3 for j in 2:8]
for i in 1:Infs
    barplot!(ga[i,1],[1],FixRFixedPleio[i,:],stack = [1,2,3], color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,2],[1],FixUFixedPleio[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,3],[1],FixDFixedPleio[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,4],[1],sl100FixedPleio[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])

    barplot!(ga[i,5],[1],FixRFixedPleio250[i,:],stack = [1,2,3], color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,6],[1],FixUFixedPleio250[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,7],[1],FixDFixedPleio250[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])
    barplot!(ga[i,8],[1],sl100FixedPleio250[i,:],stack = [1,2,3],color = [colSchem[1],colSchem[5],colSchem[10]])
end
for i in 1:Infs
    for j in 1:8
        lines!(ga[i,j],[.65,1.42],[.5,.5],linestyle = :dash, color = :white, linewidth = 2)
    end
end
# Label(ga[end+1,2], "Unevolved", textsize = 25)
# Label(ga[end,7], "Evolved", textsize = 25)
e1 = PolyElement(color = colSchem[1])
e2 = PolyElement(color = colSchem[5])
e3 = PolyElement(color = colSchem[10])
Legend(ga[2,9],[e1, e2, e3], ["Pleio.","Non-Pleio.","Draw"])

cd(string(topDir,"/Images"))
CairoMakie.save(string("s1_fixed_runs_V_constraint_Bar.svg"), f, pt_per_unit = 1)
#endregion

#region Network data from saved files
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> contains(x, "Network"), FixRFiles)

cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> contains(x, "Network"), FixUDFiles)

cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> contains(x, "Network"), slFiles)

FixRFiles250 = filter(x -> contains(x, "250"), FixRFiles)
FixUDFiles250 = filter(x -> contains(x, "250"), FixUDFiles)
slFiles250 = filter(x -> contains(x, "250"), slFiles)

FixRFiles = filter(x -> !contains(x, "250"), FixRFiles)
FixUDFiles = filter(x -> !contains(x, "250"), FixUDFiles)
slFiles = filter(x -> !contains(x, "250"), slFiles)

#group data from saved files
FixRDataSet = FixRFiles;
FixUDataSet = filter(x -> contains(x, "UpReg"), FixUDFiles);
FixDDataSet = filter(x -> contains(x, "DownReg"), FixUDFiles);
sl100DataSet = filter(x -> contains(x, "100X"), slFiles);

FixRDataSet250 = FixRFiles250;
FixUDataSet250 = filter(x -> contains(x, "UpReg"), FixUDFiles250);
FixDDataSet250 = filter(x -> contains(x, "DownReg"), FixUDFiles250);
sl100DataSet250 = filter(x -> contains(x, "100X"), slFiles250);

#region assign data to appropriate arrays
FixRPWinners = []
FixRULosers = []

FixRUWinners = []
FixRPLosers = []

FixUPWinners = []
FixUULosers = []

FixUUWinners = []
FixUPLosers = []

FixDPWinners = []
FixDULosers = []

FixDUWinners = []
FixDPLosers = []

sl100PWinners = []
sl100ULosers = []

sl100UWinners = []
sl100PLosers = []
for i = 1:Infs
    SortNets(string(FixRplotDir,FixRDataSet250[i]),FixRPWinners,FixRUWinners,FixRPLosers,FixRULosers)
    SortNets(string(FixUDplotDir,FixUDataSet250[i]),FixUPWinners,FixUUWinners,FixUPLosers,FixUULosers)
    SortNets(string(FixUDplotDir,FixDDataSet250[i]),FixDPWinners,FixDUWinners,FixDPLosers,FixDULosers)
    SortNets(string(slplotDir,sl100DataSet250[i]),sl100PWinners,sl100UWinners,sl100PLosers,sl100ULosers)
end
#endregion
#endregion

#region Figure 4
colSchem = cgrad(:haline, 100, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRPWinDens = DensityofImmuneResponses(FixRPWinners[i])
    FixUPWinDens = DensityofImmuneResponses(FixUPWinners[i])
    FixDPWinDens = DensityofImmuneResponses(FixDPWinners[i])
    sl100PWinDens = DensityofImmuneResponses(sl100PWinners[i])
    
    FixRULoseDens = DensityofImmuneResponses(FixRULosers[i])
    FixUULoseDens = DensityofImmuneResponses(FixUULosers[i])
    FixDULoseDens = DensityofImmuneResponses(FixDULosers[i])
    sl100ULoseDens = DensityofImmuneResponses(sl100ULosers[i])

    lines!(ga[i,1],FixRPWinDens./maximum(FixRPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,1],FixRULoseDens./maximum(FixRULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,1],string(round(cor(FixRPWinDens./maximum(FixRPWinDens),FixRULoseDens./maximum(FixRULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUPWinDens./maximum(FixUPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,2],FixUULoseDens./maximum(FixUULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,2],string(round(cor(FixUPWinDens./maximum(FixUPWinDens),FixUULoseDens./maximum(FixUULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDPWinDens./maximum(FixDPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,3],FixDULoseDens./maximum(FixDULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,3],string(round(cor(FixDPWinDens./maximum(FixDPWinDens),FixDULoseDens./maximum(FixDULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100PWinDens./maximum(sl100PWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,4],sl100ULoseDens./maximum(sl100ULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,4],string(round(cor(sl100PWinDens./maximum(sl100PWinDens),sl100ULoseDens./maximum(sl100ULoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
# e1 = PolyElement(color = colSchem[1])
# e2 = PolyElement(color = colSchem[70])
# Legend(ga[2,5],[e1, e2], ["Pleio.","Non-Pleio."])
cd(string(topDir,"/Images"))
CairoMakie.save("f4_WinLoseDens_Pleio_V_NP.svg", f, pt_per_unit = 1)
#endregion

#region supplemental figure 2
#sup 2a
colSchem = cgrad(:haline, 100, categorical = true);
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRPWinDens = DensityofImmuneResponses(FixRPWinners[i])
    FixUPWinDens = DensityofImmuneResponses(FixUPWinners[i])
    FixDPWinDens = DensityofImmuneResponses(FixDPWinners[i])
    sl100PWinDens = DensityofImmuneResponses(sl100PWinners[i])
    
    FixRULoseDens = DensityofImmuneResponses(FixRULosers[i])
    FixUULoseDens = DensityofImmuneResponses(FixUULosers[i])
    FixDULoseDens = DensityofImmuneResponses(FixDULosers[i])
    sl100ULoseDens = DensityofImmuneResponses(sl100ULosers[i])

    lines!(ga[i,1],FixRPWinDens./maximum(FixRPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,1],FixRULoseDens./maximum(FixRULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,1],string(round(cor(FixRPWinDens./maximum(FixRPWinDens),FixRULoseDens./maximum(FixRULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUPWinDens./maximum(FixUPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,2],FixUULoseDens./maximum(FixUULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,2],string(round(cor(FixUPWinDens./maximum(FixUPWinDens),FixUULoseDens./maximum(FixUULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDPWinDens./maximum(FixDPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,3],FixDULoseDens./maximum(FixDULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,3],string(round(cor(FixDPWinDens./maximum(FixDPWinDens),FixDULoseDens./maximum(FixDULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100PWinDens./maximum(sl100PWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,4],sl100ULoseDens./maximum(sl100ULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,4],string(round(cor(sl100PWinDens./maximum(sl100PWinDens),sl100ULoseDens./maximum(sl100ULoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = :Blue)
e2 = PolyElement(color = :Red)
Legend(ga[2,5],[e1, e2], ["Pleiotropic Winners"," Non-Pleiotropic Losers"])
cd(string(topDir,"/Images")) 
CairoMakie.save("s2a_PWinNPLose.svg", f, pt_per_unit = 1)

#sup 2b
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRPWinDens = DensityofImmuneResponses(FixRPWinners[i])
    FixUPWinDens = DensityofImmuneResponses(FixUPWinners[i])
    FixDPWinDens = DensityofImmuneResponses(FixDPWinners[i])
    sl100PWinDens = DensityofImmuneResponses(sl100PWinners[i])
    
    FixRPLoseDens = DensityofImmuneResponses(FixRPLosers[i])
    FixUPLoseDens = DensityofImmuneResponses(FixUPLosers[i])
    FixDPLoseDens = DensityofImmuneResponses(FixDPLosers[i])
    sl100PLoseDens = DensityofImmuneResponses(sl100PLosers[i])

    lines!(ga[i,1],FixRPWinDens./maximum(FixRPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,1],FixRPLoseDens./maximum(FixRPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,1],string(round(cor(FixRPWinDens./maximum(FixRPWinDens),FixRPLoseDens./maximum(FixRPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUPWinDens./maximum(FixUPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,2],FixUPLoseDens./maximum(FixUPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,2],string(round(cor(FixUPWinDens./maximum(FixUPWinDens),FixUPLoseDens./maximum(FixUPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDPWinDens./maximum(FixDPWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,3],FixDPLoseDens./maximum(FixDPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,3],string(round(cor(FixDPWinDens./maximum(FixDPWinDens),FixDPLoseDens./maximum(FixDPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100PWinDens./maximum(sl100PWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,4],sl100PLoseDens./maximum(sl100PLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,4],string(round(cor(sl100PWinDens./maximum(sl100PWinDens),sl100PLoseDens./maximum(sl100PLoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = :Blue)
e2 = PolyElement(color = :Red)
Legend(ga[2,5],[e1, e2], ["Pleiotropic Winners","Pleiotropic Losers"])
cd(string(topDir,"/Images"))
CairoMakie.save("s2b_PWinPLose.svg", f, pt_per_unit = 1)

#sup 2c
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRUWinDens = DensityofImmuneResponses(FixRUWinners[i])
    FixUUWinDens = DensityofImmuneResponses(FixUUWinners[i])
    FixDUWinDens = DensityofImmuneResponses(FixDUWinners[i])
    sl100UWinDens = DensityofImmuneResponses(sl100UWinners[i])
    
    FixRPLoseDens = DensityofImmuneResponses(FixRPLosers[i])
    FixUPLoseDens = DensityofImmuneResponses(FixUPLosers[i])
    FixDPLoseDens = DensityofImmuneResponses(FixDPLosers[i])
    sl100PLoseDens = DensityofImmuneResponses(sl100PLosers[i])

    lines!(ga[i,1],FixRUWinDens./maximum(FixRUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,1],FixRPLoseDens./maximum(FixRPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,1],string(round(cor(FixRUWinDens./maximum(FixRUWinDens),FixRPLoseDens./maximum(FixRPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUUWinDens./maximum(FixUUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,2],FixUPLoseDens./maximum(FixUPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,2],string(round(cor(FixUUWinDens./maximum(FixUUWinDens),FixUPLoseDens./maximum(FixUPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDUWinDens./maximum(FixDUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,3],FixDPLoseDens./maximum(FixDPLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,3],string(round(cor(FixDUWinDens./maximum(FixDUWinDens),FixDPLoseDens./maximum(FixDPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100UWinDens./maximum(sl100UWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,4],sl100PLoseDens./maximum(sl100PLoseDens),linewidth = 3, color = :Red)
    text!(ga[i,4],string(round(cor(sl100UWinDens./maximum(sl100UWinDens),sl100PLoseDens./maximum(sl100PLoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = :Blue)
e2 = PolyElement(color = :Red)
Legend(ga[2,5],[e1, e2], ["Non-Pleiotropic Winners","Pleiotropic Losers"])
cd(string(topDir,"/Images"))
CairoMakie.save("s2c_NPWinPLose.svg", f, pt_per_unit = 1)

#sup 2d

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRUWinDens = DensityofImmuneResponses(FixRUWinners[i])
    FixUUWinDens = DensityofImmuneResponses(FixUUWinners[i])
    FixDUWinDens = DensityofImmuneResponses(FixDUWinners[i])
    sl100UWinDens = DensityofImmuneResponses(sl100UWinners[i])
    
    FixRULoseDens = DensityofImmuneResponses(FixRULosers[i])
    FixUULoseDens = DensityofImmuneResponses(FixUULosers[i])
    FixDULoseDens = DensityofImmuneResponses(FixDULosers[i])
    sl100ULoseDens = DensityofImmuneResponses(sl100ULosers[i])

    lines!(ga[i,1],FixRUWinDens./maximum(FixRUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,1],FixRULoseDens./maximum(FixRULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,1],string(round(cor(FixRUWinDens./maximum(FixRUWinDens),FixRULoseDens./maximum(FixRULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUUWinDens./maximum(FixUUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,2],FixUULoseDens./maximum(FixUULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,2],string(round(cor(FixUUWinDens./maximum(FixUUWinDens),FixUULoseDens./maximum(FixUULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDUWinDens./maximum(FixDUWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,3],FixDULoseDens./maximum(FixDULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,3],string(round(cor(FixDUWinDens./maximum(FixDUWinDens),FixDULoseDens./maximum(FixDULoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100UWinDens./maximum(sl100UWinDens),linewidth = 3, color = :Blue)
    lines!(ga[i,4],sl100ULoseDens./maximum(sl100ULoseDens),linewidth = 3, color = :Red)
    text!(ga[i,4],string(round(cor(sl100UWinDens./maximum(sl100UWinDens),sl100ULoseDens./maximum(sl100ULoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = :Blue)
e2 = PolyElement(color = :Red)
Legend(ga[2,5],[e1, e2], ["Non-pleiotropic Winners","Non-pleiotropic Losers"])
cd(string(topDir,"/Images"))
CairoMakie.save("s2d_NPWinNPLose.svg", f, pt_per_unit = 1)

#sup 2e
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRUWinDens = DensityofImmuneResponses(FixRUWinners[i])
    FixUUWinDens = DensityofImmuneResponses(FixUUWinners[i])
    FixDUWinDens = DensityofImmuneResponses(FixDUWinners[i])
    sl100UWinDens = DensityofImmuneResponses(sl100UWinners[i])
    
    FixRPWinDens = DensityofImmuneResponses(FixRPWinners[i])
    FixUPWinDens = DensityofImmuneResponses(FixUPWinners[i])
    FixDPWinDens = DensityofImmuneResponses(FixDPWinners[i])
    sl100PWinDens = DensityofImmuneResponses(sl100PWinners[i])

    lines!(ga[i,1],FixRPWinDens./maximum(FixRPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,1],FixRUWinDens./maximum(FixRUWinDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,1],string(round(cor(FixRUWinDens./maximum(FixRUWinDens),FixRPWinDens./maximum(FixRPWinDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUPWinDens./maximum(FixUPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,2],FixUUWinDens./maximum(FixUUWinDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,2],string(round(cor(FixUUWinDens./maximum(FixUUWinDens),FixUPWinDens./maximum(FixUPWinDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDPWinDens./maximum(FixDPWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,3],FixDUWinDens./maximum(FixDUWinDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,3],string(round(cor(FixDUWinDens./maximum(FixDUWinDens),FixDPWinDens./maximum(FixDPWinDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100PWinDens./maximum(sl100PWinDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,4],sl100UWinDens./maximum(sl100UWinDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,4],string(round(cor(sl100UWinDens./maximum(sl100UWinDens),sl100PWinDens./maximum(sl100PWinDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = colSchem[1])
e2 = PolyElement(color = colSchem[70])
Legend(ga[2,5],[e2, e1], ["Pleiotropic Winners","Non-pleiotropic Winners"])
cd(string(topDir,"/Images"))
CairoMakie.save("s2e_PWinNPWin.svg", f, pt_per_unit = 1)

#2up 2f
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (0,128,nothing,1.2)
#ylabels = ["50%","90%"]
xlabels = ["Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "10% Chance of Infection \n Response Likelihood",title = xlabels[1])
ga[1, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = xlabels[i]) for i in 2:4]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false, yticks = ([0,1],["Min.", "Max."]),ylabel = "50% Chance of Infection \n Response Likelihood")
ga[3, 1] = Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response", yticks = ([0,1],["Min.", "Max."]),ylabel = "90% Chance of Infection \n Response Likelihood")
ga[2, 2:4] = [Axis(f, limits = lims, xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
ga[3, 2:4] = [Axis(f, limits = lims, xticks = ([0,128], ["C","I"]), xlabel = "Immune Response",yticklabelsvisible = false,yticksvisible = false) for j in 2:4]
for i in 1:Infs
    FixRULoseDens = DensityofImmuneResponses(FixRULosers[i])
    FixUULoseDens = DensityofImmuneResponses(FixUULosers[i])
    FixDULoseDens = DensityofImmuneResponses(FixDULosers[i])
    sl100ULoseDens = DensityofImmuneResponses(sl100ULosers[i])
    
    FixRPLoseDens = DensityofImmuneResponses(FixRPLosers[i])
    FixUPLoseDens = DensityofImmuneResponses(FixUPLosers[i])
    FixDPLoseDens = DensityofImmuneResponses(FixDPLosers[i])
    sl100PLoseDens = DensityofImmuneResponses(sl100PLosers[i])

    lines!(ga[i,1],FixRPLoseDens./maximum(FixRPLoseDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,1],FixRULoseDens./maximum(FixRULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,1],string(round(cor(FixRULoseDens./maximum(FixRULoseDens),FixRPLoseDens./maximum(FixRPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,2],FixUPLoseDens./maximum(FixUPLoseDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,2],FixUULoseDens./maximum(FixUULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,2],string(round(cor(FixUULoseDens./maximum(FixUULoseDens),FixUPLoseDens./maximum(FixUPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,3],FixDPLoseDens./maximum(FixDPLoseDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,3],FixDULoseDens./maximum(FixDULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,3],string(round(cor(FixDULoseDens./maximum(FixDULoseDens),FixDPLoseDens./maximum(FixDPLoseDens)),digits = 2)),position = (80,.89))

    lines!(ga[i,4],sl100PLoseDens./maximum(sl100PLoseDens),linewidth = 3, color = colSchem[1])
    lines!(ga[i,4],sl100ULoseDens./maximum(sl100ULoseDens),linewidth = 3, color = colSchem[70])
    text!(ga[i,4],string(round(cor(sl100ULoseDens./maximum(sl100ULoseDens),sl100PLoseDens./maximum(sl100PLoseDens)),digits = 2)),position = (80,.89))
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
e1 = PolyElement(color = colSchem[1])
e2 = PolyElement(color = colSchem[70])
Legend(ga[2,5],[e2, e1], ["Pleiotropic Losers.","Non-pleiotropic Losers"])
cd(string(topDir,"/Images"))
CairoMakie.save("s2f_PLoseNPLose.svg", f, pt_per_unit = 1)
#endregion