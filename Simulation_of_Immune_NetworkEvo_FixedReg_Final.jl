#region Documentation
# this is the script necessary to simulate fixed up and downregulation
# this was written with flexability in mind, so changes to infection chance, the number of runs, and other such parameters should be tolerated with minimal work 
# to simulate fixed upregulation set Fixedup to true, to simulate fixed downregulation set fixed up  to false. For non-pleiotropic evolution set protnode to false 
# set the save directory  (FileName) to the desired location for saved data 
#endregion
using Core: GeneratedFunctionStub
using Base: RangeStepIrregular, end_base_include, Float64, StatusActive

#Import packages
using Random
using Trapz
using Dates
using JLD2
using FileIO
using StatsBase
using Statistics

#region Declare variables
const NumHosts = 500 #declares number of hosts to be used in each generation
const InfRatio = [.1,.5,.9] #how much of the population is infected per generation
const workDir = pwd() #gets current directory for saving data
const ProtNode = true #if true protects first node from evolutio
const Fixedup = false #if true fixes upregulatory behavior of first signaling protein in constrained case

#endregion

#region Funtions
function Simulation(NumHosts,NumPars,Generations,Start_Time,Run,Fixedup,ProtNode) #function that executes a single simulation consisting of a given number of generations 
    DeathCoef = .3 #amount of parasites and hosts that die each generation
    DeathThreshold = .9 #parasites that exceed this threshold kill their host
    V = 2 #damage coefficient of infection
    gen = 1
    HostNetworks = []
    HostConcs = []
    ParNetworks = []
    HostLins = []

    for i in 1:NumHosts
        tmp = rand(0:1,(5,5))
        tmp[[1,5,21,25]].= 0
        tmp1 = (rand(Float64,(5,5))*2).-1
        tmpNet = tmp.*tmp1
        tmpInitConc = [0 0 0 0 0; .5 .5 .5 .5 .5]

        if ProtNode & Fixedup
            tmpNet[2,:] .= 0
            tmpNet[2,5] = .5
        elseif ProtNode & !Fixedup
            tmpNet[2,:] .= 0
            tmpNet[2,5] = -.5
        end

        push!(HostNetworks,tmpNet)
        push!(HostConcs,tmpInitConc)
        push!(HostLins, i)
        if i <= NumPars
            partmp = zeros(5)
            partmp[rand(2:4)] = rand()*2-1
            if sum(partmp) == 0
                partmp[2] = .5
            end
            push!(ParNetworks,partmp)
        end
    end
    ImmuneMagnitudeInf = zeros(Float64, Generations, NumHosts,2);
    ImmuneMagnitudeInf.= NaN;
    HostFit = zeros(Float64, Generations, NumHosts);
    ParFit = zeros(Float64, Generations, NumPars);
    
    GenLineageVector = []
    while gen <= Generations
        HostDead = 0
        #region Phase One: Hosts Achieve pre-infection Equilib, create blank infection concs
        for i in 1:NumHosts
            HostConcs[i] = AchieveEq(HostNetworks[i], HostConcs[i], 20, false)
        end 
        for i in 1:length(HostConcs)
            tmp = HostConcs[i]
            if i <= NumPars
                ImmuneMagnitudeInf[ gen, i,1] = round(tmp[end,end], digits = 2)
            end
        end      
        #endregion
        #region Phase Two: Infect Hosts
        InfHosts = []
        InfConcs = []
        InfectHosts(ParNetworks, HostNetworks, InfHosts)
        for i in 1:NumPars
            push!(InfConcs, hcat(HostConcs[i][end,:]', .5))
        end
        #endregion
        #region Phase Three: Infected Equilib record infected immune responses
        for i in 1:length(InfHosts)
            InfConcs[i] = AchieveEq(InfHosts[i], InfConcs[i], 20, true)
            ImmuneMagnitudeInf[ gen, i,2] = round(maximum(InfConcs[i][:,end-1]), digits = 2)
        end
        #endregion 
        #region Phase Four: calculate post infection Equilib
        PostInfConcs = []
        for i in 1:length(InfConcs)
            push!(PostInfConcs,AchieveEq(HostNetworks[i], InfConcs[i][end,1:end-1]', 5, false))
        end
        #endregion
        #region Phase Five: Calculate Fitness
        HostToPop = []
        HostDead = 0
        for i in 1:NumHosts
            if i <= NumPars
                a = InfConcs[i][:,end]
                NormArea = trapz(1:length(a),a)/trapz(1:length(a),ones(length(a)))
                if NormArea > DeathThreshold
                    push!(HostToPop,i)
                    HostDead = HostDead + 1
                end

                ParFit[gen,i] = NormArea
                HostFit[gen,i] = OrgFitness(HostConcs[i][end,end], V, NormArea, InfConcs[i][end,end-1], length(HostNetworks[i][1,:]))
            else
                #no infection so norm area becomes 0, and post inf eq is the same as the pre-inf
                HostFit[gen,i] = OrgFitness(HostConcs[i][end,end], V, 0, HostConcs[i][end,end], length(HostNetworks[i][1,:]))
            end
        end
        while HostDead > NumHosts*DeathCoef
            rem = rand(1:HostDead)
            deleteat!(HostToPop,rem)
            HostDead -= 1 
        end
        if (gen%5 == 0) & (gen <= 50)
            tmp = Dict()
            UnqNetCM = countmap(HostNetworks)
            UnqNets = collect(keys(UnqNetCM))
            for k in UnqNets
                HostNum = findfirst(x -> x == k, HostNetworks)
                lin = HostLins[HostNum]
                ImmRes = ImmuneMagnitudeInf[gen,HostNum,:]
                tmp[k] = [lin,ImmRes,UnqNetCM[k],gen]
            end
            push!(GenLineageVector,tmp)
        end
        if (gen%50 == 0) & (gen > 50) & (gen <= 500)
            tmp = Dict()
            UnqNetCM = countmap(HostNetworks)
            UnqNets = collect(keys(UnqNetCM))
            for k in UnqNets
                HostNum = findfirst(x -> x == k, HostNetworks)
                lin = HostLins[HostNum]
                ImmRes = ImmuneMagnitudeInf[gen,HostNum,:]
                tmp[k] = [lin,ImmRes,UnqNetCM[k],gen]
            end
            push!(GenLineageVector,tmp)
        end
        #endregion
        HostDeaths(gen,HostDead,HostLins,HostFit,HostToPop,NumHosts,DeathCoef,HostNetworks)

        tmpParFit = ParFit[gen,:]
        ProgNum = zeros(NumPars)
        for i in 1:NumPars
            if tmpParFit[i] < .33
                ProgNum[i] = 1
            elseif .34 < tmpParFit[i] < .66
                ProgNum[i] = 2
            else
                ProgNum[i] = 3
            end 
        end

        p = sortperm(vec(tmpParFit), rev = true)
        ParNetworks = ParNetworks[p]
        ProgNum = ProgNum[p]
        ParNetworks = ParNetworks[1:NumPars - Int64(NumPars*DeathCoef)]                
        DeadHosts = NumHosts-length(HostNetworks)     
        #region Phase Seven: repopulate Hosts and parasitse allowing for evolution
        OrgEvolution(HostNetworks,HostLins,HostFit[gen,:],ParFit[gen,:],DeathThreshold,ProtNode,DeadHosts,NumHosts)
        MeanSigs = floor(Statistics.mean([length(x[1,:]) for x in HostNetworks]))-2
        ParEvolution(ParNetworks, ProgNum, NumPars*DeathCoef, Int64(MeanSigs))
        
        HostShuffle = Random.randperm(NumHosts)  
        HostNetworks = HostNetworks[HostShuffle]
        HostLins = HostLins[HostShuffle]
        for i in 1:length(HostConcs)
            HostConcs[i] = zeros(2, length(HostNetworks[i][1,:]))
            HostConcs[i][2,:] .+= .5
        end
        
        ParShuffle = Random.randperm(NumPars)
        ParNetworks = ParNetworks[ParShuffle]
        #endregion
        if gen%100 == 0
            print('\n',Run, ' ', gen,' ',(Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000), ' ', length(HostNetworks))
        end
        
        gen = gen+1
    end
    return HostNetworks, HostFit, ParNetworks, ParFit, ImmuneMagnitudeInf, GenLineageVector
end
function HostDeaths(gen,HostDead,HostLins,HostFit,HostToPop,NumHosts,DeathCoef,HostNetworks) #kills specified number of hosts in a population
    culled = false
    Host=1
    while !culled
        if HostDead >= DeathCoef*NumHosts
            culled = true
        elseif Host == NumHosts
            culled = true
        end
        if !culled 
            Culpoint = count(x->x < HostFit[gen,Host],HostFit[gen,:])/NumHosts
            
            if HostFit[gen,Host] < 1e-3
                push!(HostToPop,Host)
                HostDead += 1  
            elseif rand() > Culpoint
                push!(HostToPop,Host)
                HostDead += 1  
            end
            
            Host += 1
        end
    end

    while HostDead > NumHosts*DeathCoef
        rem = rand(1:HostDead)
        deleteat!(HostToPop,rem)
        HostDead -= 1 
    end
    HostToPop = sort(HostToPop, rev = true)
    for rem in HostToPop
        deleteat!(HostNetworks, rem)
        deleteat!(HostLins,rem)
    end
end
function AchieveEq(Net, Conc, StepLim, InfEq) #calculates changes in system of DEs that control network interactions
    Equilib = false
    Step = 1
    UseCoef = .01 #amount of protein deactivated for simply interacting with other proteins
    NumProts = length(Net[1,:])
    tmpConc = zeros(StepLim+2,NumProts)
    tmpConc[2,:] = Conc[end,:]
    counter = 2
    AllAct = Vector{Vector{Float64}}(undef,NumProts)
    AllUse = Vector{Int64}(undef,NumProts)
    for i in 1:NumProts
        AllAct[i] = Net[:,i]
        AllUse[i] = sum(Net[i,:].!= 0)
    end
    tmp = zeros(NumProts)
    while ~Equilib 
        for prot in 1:NumProts #collect up and down regulating actions for each protein, sum them, and have that as the conDelta for each protein
            InitProtConc = tmpConc[counter,prot]
            ActingProts = AllAct[prot]
            UsedIn = AllUse[prot]
            for j in 1:length(tmp)
                tmp[j] = ActingProts[j]*tmpConc[counter,j]
            end
            for i in 1:length(ActingProts)
                 if ActingProts[i] > 0 
                    tmp[i] = tmp[i]*(1-InitProtConc)
                 elseif ActingProts[i] < 0
                    tmp[i] = tmp[i]*(InitProtConc)
                 end
            end
            tmpConc[counter+1,prot] = tmpConc[counter,prot] + (sum(tmp) - (UseCoef*UsedIn))
        end
        for j in 1:NumProts
            if tmpConc[counter+1,j] > 1
                tmpConc[counter+1,j] = 1
            elseif tmpConc[counter+1,j] < 0
                tmpConc[counter+1,j] = 0
            end
        end

        Step += 1

        if InfEq #escape functions differently in infected cases, need to go a set number of steps with the parasite
            if tmpConc[counter+1,end] <= 1e-2
                Equilib = true
            elseif Step == StepLim
                Equilib = true
            end
        else
            tmpDif = tmpConc[counter+1,[end]]-tmpConc[counter,[end]]
            if abs(tmpDif[1])<1e-2
                Equilib = true
            end

            if Step> StepLim
                Equilib = true
            end
        end
        counter += 1
    end 
    return tmpConc[1:counter,:]
end
function InfectHosts(ParNetworks, HostNetworks, InfHosts) #combines hosts and parasites to create infected hosts
    for par in 1:length(ParNetworks)
        HostNet = HostNetworks[par]
        ParNet = ParNetworks[par]
        ParTarg = findall(x->x != 0,ParNet)
        if length(ParTarg) == 1
            ParTarg = ParTarg[1]
        end

        if ParTarg in 2:length(HostNet[1,:])-1
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
end
function OrgEvolution(HostNetworks,HostLins,HostFit,ParFit,DeathThreshold,ProtNode,ToRep,NumHosts)
    NewOrgs = []
    NewLins = []
    while length(NewOrgs) < ToRep
        RepHost = rand(1:length(HostNetworks))
        RepChance = count(x->x <= HostFit[RepHost],HostFit)/NumHosts
        RepCheck = rand()
        if (RepChance > RepCheck) 
            Dupe = copy(HostNetworks[RepHost])
            ToMute = rand()
            if ToMute < 5e-3 #if the mut threshold is passed, go on to mutations
                Mutation = rand()

                if Mutation <= .25 #add edge to graph Lineage code: 1
                    #print(' ', 1)
                    if ~ProtNode
                        len = length(Dupe[1,:])
                        Zers = findall(x->x==0,Dupe)

                        deleteat!(Zers,findall(x->x == CartesianIndex(1,1), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(1,len), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(len,1), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(len,len), Zers)[1])

                    elseif ProtNode
                        len = length(Dupe[1,:])
                        Zers = findall(x->x==0,Dupe)
                        deleteat!(Zers,findall(x->x == CartesianIndex(1,1), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(1,len), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(len,1), Zers)[1])
                        deleteat!(Zers,findall(x->x == CartesianIndex(len,len), Zers)[1])

                        Count = 1
                        store = []
                        for X in Zers
                            if X[1] == 2
                                push!(store,Count)
                            end
                            Count +=1
                        end

                        reverse!(store)
                        for X in store
                            deleteat!(Zers,X)
                        end
                    end

                    if length(Zers) >0
                        AddEdge = Zers[rand(1:length(Zers))]
                        Dupe[AddEdge[1],AddEdge[2]] = (rand()*2)-1
                    end

                elseif (Mutation > .25) & (Mutation <= .5) #delete edge from graph

                    #print(' ', 2)
                    if ~ProtNode
                        Ones = findall(x -> x!=0,Dupe)

                    elseif ProtNode
                        Ones = findall(x -> x!=0,Dupe)

                        Count = 1
                        store = []
                        for X in Ones
                            if X[1] == 2
                                push!(store,Count)
                            end
                            Count +=1
                        end
                        reverse!(store)
                        for X in store
                            deleteat!(Ones,X)
                        end
                    end

                    if length(Ones) >1
                        DelEdge = Ones[rand(1:length(Ones))]
                        Dupe[DelEdge[1],DelEdge[2]] = 0
                    end
                elseif (Mutation > .5) & (Mutation <= .8) #change coefficient by 10% randomly up or down
                    if ~ProtNode
                        Coefs = findall(x -> x!=0, Dupe)

                    elseif ProtNode
                        Coefs = findall(x -> x!=0,Dupe)

                        Count = 1
                        store = []
                        for X in Coefs
                            if X[1] == 2
                                push!(store,Count)
                            end
                            Count +=1
                        end
                        reverse!(store)
                        for X in store
                            deleteat!(Coefs,X)
                        end
                    end

                    if length(Coefs) > 0
                        ChanCoef = Coefs[rand(1:length(Coefs))]
                        delta = Dupe[ChanCoef[1],ChanCoef[2]]*.1
                        tmp = rand()
                        if tmp > .5
                            Dupe[ChanCoef[1],ChanCoef[2]] = Dupe[ChanCoef[1],ChanCoef[2]]+delta
                        else
                            Dupe[ChanCoef[1],ChanCoef[2]] = Dupe[ChanCoef[1],ChanCoef[2]]-delta
                        end
                        if Dupe[ChanCoef[1],ChanCoef[2]] > 1
                            Dupe[ChanCoef[1],ChanCoef[2]] = 1
                        elseif Dupe[ChanCoef[1],ChanCoef[2]] < -1
                            Dupe[ChanCoef[1],ChanCoef[2]] = -1
                        end
                    end
                elseif (Mutation >.8) & (Mutation <= .9) #delete a protein from the network
                    #print(' ', 4)
                    if ~ProtNode
                        if length(Dupe[1,:]) > 3
                            ToDel = rand(2:length(Dupe[1,:])-1)
                        else
                            ToDel = NaN
                        end
                    elseif ProtNode
                        if length(Dupe[1,:]) > 3
                            ToDel = rand(3:length(Dupe[1,:])-1)
                        else
                            ToDel = NaN
                        end
                    end

                    if ~isnan(ToDel)
                        Dupe = Dupe[1:end .!= ToDel,1:end .!= ToDel ]
                    end
                else #duplicate a protein in the network
                    toDup = rand(2:length(Dupe[1,:])-1)
                    tmpRow = Dupe[toDup,:]
                    tmpRowStack = vcat(Dupe[1:toDup-1,:],tmpRow',Dupe[toDup:end,:])
                    tmpCol = tmpRowStack[:,toDup]
                    tmpColStack = hcat(tmpRowStack[:,1:toDup-1],tmpCol,tmpRowStack[:,toDup:end])
                    Dupe = tmpColStack
                end
            end
            push!(NewOrgs,Dupe)
            push!(NewLins,HostLins[RepHost])
        end
    end
    for k in NewOrgs
        push!(HostNetworks,k)
    end
    for k in NewLins
        push!(HostLins,k)
    end
end
function OrgFitness(PEP, V, Area, PEPO ,NumProts) #calculates organismal fitness
    if NumProts < 11
        ProtCost = 0
    else 
        ProtCost = 1.1^(NumProts-10)
    end
    exp(-(PEP + (V*Area) + PEPO + ProtCost))
end
function ParEvolution(ParNetworks, ProgNum, ToRep, MeanSigs) #repopulates parasites with chance for mutation
    NewPars = []
    Rep = false
    par = 1
    while !Rep
        prog = ProgNum[par]
        for Surv in 1:prog
            Dupe = copy(ParNetworks[par])
            if length(Dupe) < 2+MeanSigs
                tmp = zeros(2+MeanSigs)
                tmp[1:length(Dupe)] = Dupe
                Dupe = tmp
            end
        
            if rand() < 1e-2
                Mutation = rand()
                if Mutation < .5
                    if MeanSigs == 1
                        prevTarg = findall(x->x !=0,Dupe)[1]
                        Dupe[2] = Dupe[prevTarg]
                        if prevTarg != 2
                            Dupe[prevTarg] = 0
                        end
                    else
                        prevTarg = findall(x->x !=0,Dupe)[1]
                        NewTarg = Int64(rand(2:MeanSigs+1))
                        Dupe[NewTarg] = Dupe[prevTarg]
                        if NewTarg != prevTarg
                            Dupe[prevTarg] = 0
                        end
                    end
                else
                    prevTarg = findall(x->x !=0,Dupe)[1]
                    Dupe[prevTarg] = rand()*2-1
                end
            end
            push!(NewPars,Dupe)
        end
        par += 1
        if length(NewPars) >= ToRep
            Rep = true
        end
    end
    for i in 1:Int64(ToRep)
        push!(ParNetworks, NewPars[i])
    end
end
#endregion
for inf in InfRatio #cycles through specified infections chances, conducts (Runs) simulations at a given level of infection 
    #collection of arrays for the saving of data
    Runs = 100 #number of runs the sim goes through
    Generations = 500 #number of generations the sim runs for
    NumPars = Int64(NumHosts*inf); #how many parasites to create in each run
    HostFit = zeros(Float64, Runs,  Generations, NumHosts); #tracks host fitness across Runs, and each generation within a run
    ParFit = zeros(Float64, Runs,  Generations, NumPars); #tracks parasites fitness across Runs, and each generation within a run
    
    ImmuneMagnitudeInf = zeros(Float64, Runs,  Generations, NumHosts,2); #pre infection immune effector levels and peri infection peak immune effector
    ImmuneMagnitudeInf.= NaN;

    MostCommonNet =  [] #array to save the most common host in each simulation
    RunLineageVector = []

    ParTargets = zeros(Float64, Runs,  Generations, NumPars); #what signaling protein the parasites are acting on

    Run = 1
    Start_Time = Dates.Time(Dates.now())
    while Run <= Runs #conducts Runs number of simulations
        #call simulation function
        HostNetworks, HostFit[Run,:,:], ParNetworks, ParFit[Run,:,:], ImmuneMagnitudeInf[Run,:,:,:], Blah = Simulation(NumHosts,NumPars,Generations,Start_Time,Run,Fixedup,ProtNode)
        push!(RunLineageVector,Blah)
        #identify most common Host and save
        UnqNetCM = countmap(HostNetworks)
        Vals = collect(values(UnqNetCM))
        UnqNets = collect(keys(UnqNetCM))
        MostCom = findall(x -> x == maximum(Vals),Vals)
        push!(MostCommonNet, UnqNets[MostCom])
        Run += 1
        if Fixedup
            print('\n',"******** just did fixed up regulation ", inf, "**********  ", (Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000))
        else
            print('\n',"******** just did fixed down regulation ", inf, "********** ", (Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000))
        end
    end

    #saves data where each generation has been averaged (no individual host response, but significantly smaller file sizes)
    FinalImmuneMagnitudeInf = ImmuneMagnitudeInf[:,Generations,:,:]
    #saves snapshots of immune activity for change over time plots
    Short = collect(1:Generations/100:Generations/10)
    push!(Short,Generations/10)
    Short = convert.(Int64,Short)
    ShortImmuneMagSnaps = ImmuneMagnitudeInf[:,Short,:,:]

    Long = collect(1:Generations/10:Generations)
    push!(Long,Generations)
    Long = convert.(Int64,Long)
    LongImmuneMagSnaps = ImmuneMagnitudeInf[:,Long,:,:]
    
    HostFit = reshape(mean(HostFit, dims = 3),(Runs,Generations))
    ParFit = reshape(mean(ParFit, dims = 3),(Runs,Generations))
    if ProtNode & Fixedup
        FileName = string(workDir,"/Data/FixedupDown/",string(inf*100)[1:end-2],"_PercentInf_Fixed_upreg_Schrom_Selection_Martin_Deaths.jld2")
    elseif ProtNode &! Fixedup
        FileName = string(workDir,"/Data/FixedUpDown/",string(inf*100)[1:end-2],"_PercentInf_Fixed_downreg_Schrom_Selection_Martin_Deaths.jld2")
    else
        FileName = string(workDir,"/Data/FixedUpDown/",string(inf*100)[1:end-2],"_PercentInf.jld2")
    end
    
    save(FileName,"ImmMagInf",FinalImmuneMagnitudeInf,"ShortImmMagInfSnaps",ShortImmuneMagSnaps,"LongImmMagInfSnaps",LongImmuneMagSnaps,"Trajectories",RunLineageVector,"ParTarg",ParTargets,"HostFit",HostFit,"ParFit",ParFit)
    save(string(FileName[1:end-5],"_Networks.jld2"),"Networks",MostCommonNet)
end
