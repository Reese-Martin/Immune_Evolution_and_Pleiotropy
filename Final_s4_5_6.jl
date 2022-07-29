using FileIO
using Statistics 
#using DataFrames
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using KernelDensity
using HypothesisTests

Gens = 500
Runs = 100
Infs = 3
Hosts = 500
ImmuneTypes = 4 #0- no response, 1- constitutive, 2- Inducible, 3- mixed indu/const

#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/IMMUNE_PLEIO_MODEL/Julia Code"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/IMMUNE_PLEIO_MODEL/Julia Code"

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
function issubsequence(A,B)
    for i in 1:length(B)-length(A)
      if B[i:i+length(A)-1] == A
        return true
      end
    end
    return false
end
function DEPaths(Mat)
    #take in a directed matrix and find all paths between v1 and vE
    # tmp = copy(Mat) # make a copy then get rid of self regulatory edges
    # for j in 1:length(Mat[1,end])
    #     tmp[j,j] = 0
    # end
    edges = findall(x -> x!=0,Mat)
    d = Dict()
    for k in edges
        if in(k[1],keys(d))
            push!(d[k[1]],k[2])
        else 
            d[k[1]] = [k[2]]
        end
    end

    chains = []
    if in(1,keys(d))
        for x in d[1]
            push!(chains,[1,x])
        end
    end
    storage = [] 

    while ~isempty(chains)
        tmp = copy(chains)
        DeadEndCheck = []
        for ii in 1:length(tmp)
            if ~in(chains[ii][end],keys(d))
                push!(DeadEndCheck,ii)
            end 
        end
        reverse!(DeadEndCheck)
        for X in DeadEndCheck
            deleteat!(tmp,X)
        end

        for i in 1:length(chains)
            if in(chains[i][end],keys(d))
                for k in d[chains[i][end]]
                    push!(tmp, vcat(chains[i],k))
                end
            end
        end
        if ~isempty(tmp)
            subCheck = []
            LongestPath = maximum(length.(tmp))
            for j in 1:length(tmp)
                if length(tmp[j]) < LongestPath
                    for k in j+1:length(tmp)
                        if issubsequence(tmp[j],tmp[k])
                            push!(subCheck,j)
                            break
                        end
                    end
                end
            end
            reverse!(subCheck)
            for X in subCheck
                deleteat!(tmp,X)
            end
        end

        cycleCheck = []
        for j in 1:length(tmp) 
            if  in(tmp[j][end],tmp[j][1:end-1])
                push!(cycleCheck,j)
            end
        end
        reverse!(cycleCheck)
        for X in cycleCheck
            deleteat!(tmp,X)
        end

        chains = copy(tmp)
        endCheck = []
        for k in 1:length(chains)
            if chains[k][end] == length(Mat[1,:]) 
                push!(storage,chains[k])
                push!(endCheck,k)
            end
        end
        reverse!(endCheck)
        for X in endCheck
            deleteat!(chains,X)
        end
    end
    return(storage)
end
function DistinctPaths(a)
    tmp = trim.(a)
    tmp2 = []
    toRem = []
    for i in 1:length(tmp)
        tmp3 = []
        for j in 1:length(tmp[i])
            if j == 1
                push!(tmp3,tmp[i][j])# start here when fixing, need to cycle through each entry in tmp[i]
            else 
                tmp4 = 0
                for k in tmp3
                    if !isempty(intersect(k,tmp[i][j])) 
                        tmp4 += 1
                    end
                end
                if tmp4 == 0
                    push!(tmp3,tmp[i][j])
                end   
            end
        end
        push!(tmp2,tmp3)
    end
    return tmp2
end
function trim(x)
    tmp = []
    for i in x
        push!(tmp,i[2:end-1])
    end
    return tmp
end
function MakeFig(size,lims = (nothing, nothing, nothing, nothing))
    
    fig = Figure(resolution = size)

    ax1 = Axis(fig[1,1],xticklabelsvisible = false,xticksvisible = false, title = "Unconstrained", ylabel = "10%",limits = lims)
    ax2 = Axis(fig[2,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "50%",limits = lims)
    ax3 = Axis(fig[3,1],xticklabelsvisible = false,xticksvisible = false, ylabel = "90%",limits = lims)

    ax4 = Axis(fig[1,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Random",limits = lims)
    ax5 = Axis(fig[2,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax6 = Axis(fig[3,2],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)

    ax7 = Axis(fig[1,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Up",limits = lims)
    ax8 = Axis(fig[2,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax9 = Axis(fig[3,3],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)

    ax10 = Axis(fig[1,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Fixed Down",limits = lims)
    ax11 = Axis(fig[2,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax12 = Axis(fig[3,4],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)

    ax13 = Axis(fig[1,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,title = "Slow",limits = lims)
    ax14 = Axis(fig[2,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    ax15 = Axis(fig[3,5],xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false, yticksvisible = false,limits = lims)
    return fig
end
function AddMedianLines(x,Data,Runs, row,col)
    Rows = Runs
    Cols = Int64(length(Data)/Rows)
    tmp = median(reshape(Data,(Rows,Cols)),dims = 1)
    lines!(ga[row,col],x[1],[tmp[1],tmp[1]], linewidth = 5, color = :Black)
end
function DetConnect(net)
    TotalPossibleCons = size(net,1)^2 - 2
    TotalCons = size(findall(x -> x!= 0, net),1)
    return TotalCons/TotalPossibleCons
end

#region Processing
FixRplotDir = string(topDir,"/Data/FixedRand/Draft_2_Data/")
FixUDplotDir = string(topDir,"/Data/FixedUpDown/Draft_2_Data/")
slplotDir = string(topDir,"/Data/SlowEvo/Draft_2_Data/")
cd(FixRplotDir)
FixRFiles = filter(x->endswith(x, ".jld2"), readdir())
FixRFiles = filter(x -> contains(x, "Network"), FixRFiles)
#FixRFiles = filter(x -> contains(x, "test_"), FixRFiles)

cd(FixUDplotDir)
FixUDFiles = filter(x->endswith(x, ".jld2"), readdir())
FixUDFiles = filter(x -> contains(x, "Network"), FixUDFiles)
#FixUDFiles = filter(x -> contains(x, "test_"), FixUDFiles)

cd(slplotDir)
slFiles = filter(x->endswith(x, ".jld2"), readdir())
slFiles = filter(x -> contains(x, "Network"), slFiles)
#slFiles = filter(x -> contains(x, "test_"), slFiles)

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

#Final generation Connections to Effector
uncDEPaths = []
FixRDEPaths = []
FixUDEPaths = []
FixDDEPaths = []
sl10DEPaths = []
sl100DEPaths = []

#Final generation Connections to Effector
uncNetSize = []
FixRNetSize = []
FixUNetSize = []
FixDNetSize = []
sl100NetSize = []

#Mean data for plotting 
for i in 1:Infs
    push!(uncDEPaths, DEPaths.(uncNets[i]))
    push!(FixRDEPaths, DEPaths.(FixRNets[i]))
    push!(FixUDEPaths, DEPaths.(FixUNets[i]))
    push!(FixDDEPaths, DEPaths.(FixDNets[i]))
    push!(sl100DEPaths, DEPaths.(sl100Nets[i]))

    push!(uncNetSize, size.(uncNets[i],1))
    push!(FixRNetSize, size.(FixRNets[i],1))
    push!(FixUNetSize, size.(FixUNets[i],1))
    push!(FixDNetSize, size.(FixDNets[i],1))
    push!(sl100NetSize, size.(sl100Nets[i],1))
end
uncDistDE = []
FixRDistDE = []
FixUDistDE = []
FixDDistDE = []
sl100DistDE = []
for i in 1:Infs
    push!(uncDistDE, DistinctPaths(uncDEPaths[i]))
    push!(FixRDistDE, DistinctPaths(FixRDEPaths[i]))
    push!(FixUDistDE, DistinctPaths(FixUDEPaths[i]))
    push!(FixDDistDE, DistinctPaths(FixDDEPaths[i]))
    push!(sl100DistDE, DistinctPaths(sl100DEPaths[i]))
end
uncCon = []
FixRCon = []
FixUCon = []
FixDCon = []
sl100Con = []
for i in 1:Infs
    push!(uncCon,DetConnect.(uncNets[i]))
    push!(FixRCon,DetConnect.(FixRNets[i]))
    push!(FixUCon,DetConnect.(FixUNets[i]))
    push!(FixDCon,DetConnect.(FixDNets[i]))
    push!(sl100Con,DetConnect.(sl100Nets[i]))
end
#endregion
colSchem = cgrad(:haline, 100, categorical = true);

#region supplemental 6

a = fill(1.0,Runs)
Jit = jitter(a,.3)
colSchem = cgrad(:haline, 100, categorical = true);

f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (nothing,nothing,-1,10)
#ylabels = ["50%","90%"]
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "10% Chance of Infection \n Distinct Connections to Effector",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "50% Chance of Infection \n Distinct Connections to Effector")
ga[3, 1] = Axis(f, limits = lims, ygridvisible = false, xticklabelsvisible = false, xticksvisible = false, ylabel = "90% Chance of Infection \n Distinct Connections to Effector")
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
Ttests = zeros(3,4)
ImTypes = fill(1.0,100)
xJitter = jitter(ImTypes,.3)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
for i in 1:Infs
    tmp1 = size.(uncDistDE[i],1)
    tmp2 = size.(FixRDistDE[i],1)
    tmp3 = size.(FixUDistDE[i],1)
    tmp4 = size.(FixDDistDE[i],1)
    tmp5 = size.(sl100DistDE[i],1)

    Ttests[i,1] = pvalue(EqualVarianceTTest(tmp1,tmp2))
    Ttests[i,2] = pvalue(EqualVarianceTTest(tmp1,tmp3))
    Ttests[i,3] = pvalue(EqualVarianceTTest(tmp1,tmp4))
    Ttests[i,4] = pvalue(EqualVarianceTTest(tmp1,tmp5))

    scatter!(ga[i,1], xJitter, tmp1, color = colSchem[70])
    AddMedianLines(xLocs,tmp1,Runs, i,1)
    scatter!(ga[i,2], xJitter, tmp2, color = colSchem[1])
    AddMedianLines(xLocs,tmp2,Runs, i,2)
    scatter!(ga[i,3], xJitter, tmp3, color = colSchem[1])
    AddMedianLines(xLocs,tmp3,Runs, i,3)
    scatter!(ga[i,4], xJitter, tmp4, color = colSchem[1])
    AddMedianLines(xLocs,tmp4,Runs, i,4)
    scatter!(ga[i,5], xJitter, tmp5, color = colSchem[1])
    AddMedianLines(xLocs,tmp5,Runs, i,5) 
end
Sigs = findall(x -> x<.05/12,Ttests)
for i in Sigs
    text!(ga[i[1],i[2]+1],"*",position = (.25,7.5),textsize = 65)
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s6_DistConDE.svg", f, pt_per_unit = 1)
#endregion

#region supplemental 4
Ttests = zeros(3,4)
ImTypes = fill(1.0,100)
xJitter = jitter(ImTypes,.1)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (.5,1.5,3, 11)
#ylabels = ["50%","90%"]
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "10% Chance of Infection \n Network Size",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "50% Chance of Infection \n Network Size")
ga[3, 1] = Axis(f, limits = lims, ygridvisible = false, xticklabelsvisible = false, xticksvisible = false, ylabel = "90% Chance of Infection \n Network Size")
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
for i in 1:Infs
    tmp1 = uncNetSize[i]
    tmp2 = FixRNetSize[i]
    tmp3 = FixUNetSize[i]
    tmp4 = FixDNetSize[i]
    tmp5 = sl100NetSize[i]

    Ttests[i,1] = pvalue(EqualVarianceTTest(tmp1,tmp2))
    Ttests[i,2] = pvalue(EqualVarianceTTest(tmp1,tmp3))
    Ttests[i,3] = pvalue(EqualVarianceTTest(tmp1,tmp4))
    Ttests[i,4] = pvalue(EqualVarianceTTest(tmp1,tmp5))
    
    scatter!(ga[i,1], xJitter, tmp1, color = colSchem[70])
    AddMedianLines(xLocs,tmp1,Runs, i,1)
    scatter!(ga[i,2], xJitter, tmp2, color = colSchem[1])
    AddMedianLines(xLocs,tmp2,Runs, i,2)
    scatter!(ga[i,3], xJitter, tmp3, color = colSchem[1])
    AddMedianLines(xLocs,tmp3,Runs, i,3)
    scatter!(ga[i,4], xJitter, tmp4, color = colSchem[1])
    AddMedianLines(xLocs,tmp4,Runs, i,4)
    scatter!(ga[i,5], xJitter, tmp5, color = colSchem[1])
    AddMedianLines(xLocs,tmp5,Runs, i,5)
end
Sigs = findall(x -> x<.05/12,Ttests)
for i in Sigs
    text!(ga[i[1],i[2]+1],"*",position = (.6,9),textsize = 65)
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s4_NetSize.svg", f, pt_per_unit = 1)
#endregion

#region supplemental 5
Ttests = zeros(3,4)
ImTypes = fill(1.0,100)
xJitter = jitter(ImTypes,.1)
xLocs = []
for i in 1:5
    push!(xLocs,[.5,1.5])
end
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1100))
lims = (.5,1.5,0,1)
#ylabels = ["50%","90%"]
xlabels = ["Non-pleiotropic","Fixed Random", "Fixed Up", "Fixed Down", "Slow"]
ga = f[1, 1] = GridLayout()
ga[1,1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "10% Chance of Infection \n Connectivity",title = xlabels[1])
ga[1, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false,title = xlabels[i]) for i in 2:5]
ga[2, 1] = Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, ylabel = "50% Chance of Infection \n Connectivity")
ga[3, 1] = Axis(f, limits = lims, ygridvisible = false, xticklabelsvisible = false, xticksvisible = false, ylabel = "90% Chance of Infection \n Connectivity")
ga[2, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
ga[3, 2:5] = [Axis(f, limits = lims, xticklabelsvisible = false, xticksvisible = false, yticklabelsvisible = false, yticksvisible = false) for j in 2:5]
for i in 1:Infs
    tmp1 = uncCon[i]
    tmp2 = FixRCon[i]
    tmp3 = FixUCon[i]
    tmp4 = FixDCon[i]
    tmp5 = sl100Con[i]

    Ttests[i,1] = pvalue(EqualVarianceTTest(tmp1,tmp2))
    Ttests[i,2] = pvalue(EqualVarianceTTest(tmp1,tmp3))
    Ttests[i,3] = pvalue(EqualVarianceTTest(tmp1,tmp4))
    Ttests[i,4] = pvalue(EqualVarianceTTest(tmp1,tmp5))

    scatter!(ga[i,1], xJitter, tmp1, color = colSchem[70])
    AddMedianLines(xLocs,tmp1,Runs, i,1)
    scatter!(ga[i,2], xJitter, tmp2, color = colSchem[1])
    AddMedianLines(xLocs,tmp2,Runs, i,2)
    scatter!(ga[i,3], xJitter, tmp3, color = colSchem[1])
    AddMedianLines(xLocs,tmp3,Runs, i,3)
    scatter!(ga[i,4], xJitter, tmp4, color = colSchem[1])
    AddMedianLines(xLocs,tmp4,Runs, i,4)
    scatter!(ga[i,5], xJitter, tmp5, color = colSchem[1])
    AddMedianLines(xLocs,tmp5,Runs, i,5)
end
Sigs = findall(x -> x<.05/12,Ttests)
for i in Sigs
    text!(ga[i[1],i[2]+1],"*",position = (.6,.75),textsize = 65)
end
Font_Theme = Theme(font = :Arial)
set_theme!(Font_Theme)
cd(string(topDir,"/Images"))
CairoMakie.save("s5_MeanNetCon.svg", f, pt_per_unit = 1)
#endregion