#region Documentation
# this is the script necessary to simulate competion between non-pleiotropic and Fixed up or downregulation
# this was written with flexability in mind, so changes to infection chance, the number of runs, and other such parameters should be tolerated with minimal work 
# set the save directory  (FileName) to the desired location for saved data 
#endregion
#Import packages
using Base: RangeStepIrregular, end_base_include, Float64, StatusActive
using Random
using Trapz
using Dates
using JLD2
using FileIO
using StatsBase
import Statistics
using NaNStatistics


#region Declare variables
NumHostsCon = NumHostsUnc = 250 #declares number of unconstrained and constrained hosts. if these numbers differ the code breaks (as it is currently written)
EqGens = 1 #number of generations the populations equilibrate for
Runs = 100 #number of runs the sim goes through
FixationGens = 1000

InfRatio = [.1, .5, .9] #how much of the population is infected

StartingParCon = .5 #initial parasite ammount
StartConc = [0 0 0 0 0; .5 .5 .5 .5 .5] #initial network amounts
UseCoef = .01 #amount of protein deactivated for simply interacting with other proteins

V = 2 #damage coefficient of infection

DeathCoef = .3 #amount of parasites and hosts that die each generation
DeathThreshold = .9 #parasites that exceed this threshold kill their host

workDir = pwd() #gets current directory for saving
#endregion

#region Funtions
function AcheiveEq(Net, Conc, StepLim, InfEq) #step by step calculates the system of DEs that control network interactions
    Equilib = false
    Step = 1
    
    while ~Equilib 
    
        tmpConcDelta = zeros(1,length(Net[1,:]))

        for prot in 1:length(tmpConcDelta) #collect up and down regulating actions for each protein, sum them, and have that as the conDelta for each protein

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

        if any(Conc[[end],:].>1) #keeps concentration <= 1
            Conc[[end],findall(>(1),Conc[end,:])].=1
        end
        if any(Conc[[end],:].<1e-2) #keeps concentration >=0
            Conc[[end],findall(<(1e-2),Conc[end,:])].=0
        end

        Step += 1

        if InfEq #escape functions differently in infected cases, need to go a set number of steps with the parasite
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
function DetImmRes(ImmMag, Run, gen, Hosts)
    tmpImmRes = zeros(Hosts)
    for i in 1:Hosts
        if ImmMag[Run,gen,i,1] >= .1
            if (ImmMag[Run,gen,i,1] <= ImmMag[Run,gen,i,2]+.05) & (ImmMag[Run,gen,i,1] >= ImmMag[Run,gen,i,2]-.05)
                tmpImmRes[i] = 1 #constitutive response
            elseif  (ImmMag[Run,gen,i,1] > ImmMag[Run,gen,i,2]-.05)
                tmpImmRes[i] = 1 #constitutive response
            else
                tmpImmRes[i] = 3 #mixed response
            end
        else
            if ImmMag[Run,gen,i,2] < .1
                tmpImmRes[i] = 0 # No response
            else
                tmpImmRes[i] = 2 #pure Inducible
            end
        end
    end
    return tmpImmRes
end
function InfectHosts(ParNetworks, HostNetworks, InfHosts)
    for par in 1:length(ParNetworks)
        HostNet = HostNetworks[par]
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
end
function OrgEvolution(HostNetworks,HostCons, ToRep)
    NewOrgs = []
    NewCons = []
    for Surv in 1:ToRep #run through survivors and potentially mutate
        RepHost = rand(1:length(HostNetworks))
        Dupe = copy(HostNetworks[RepHost])
        ToMute = rand()
        HostCon = HostCons[RepHost]
        if ToMute < 5e-3 #if the mut threshold is passed, go on to mutations
            Mutation = rand()

            if Mutation <= .25 #add edge to graph Lineage code: 1
                #print(' ', 1)
                if HostCon
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
                else            
                    len = length(Dupe[1,:])
                    Zers = findall(x->x==0,Dupe)

                    deleteat!(Zers,findall(x->x == CartesianIndex(1,1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(1,len), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(len,1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(len,len), Zers)[1])
                end

                if length(Zers) >0
                    AddEdge = Zers[rand(1:length(Zers))]
                    Dupe[AddEdge[1],AddEdge[2]] = rand()
                end
                
            elseif (Mutation > .25) & (Mutation <= .5) #delete edge from graph
                
                #print(' ', 2)
                if HostCon
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
                else
                    Ones = findall(x -> x!=0,Dupe)
                end

                if length(Ones) >1
                    DelEdge = Ones[rand(1:length(Ones))]
                    Dupe[DelEdge[1],DelEdge[2]] = 0
                end 
            elseif (Mutation > .5) & (Mutation <= .8) #change coefficient by 10% randomly up or down
               # print(' ', 3)
                if HostCon
                    Coefs = findall(x -> x!=0, Dupe)
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
                else
                    Coefs = findall(x -> x!=0,Dupe)
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
                if HostCon 
                    if length(Dupe[1,:]) > 3
                        ToDel = rand(2:length(Dupe[1,:])-1)
                    else 
                        ToDel = NaN
                    end
                else
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
        push!(NewCons, HostCon)
    end
    for k in NewOrgs
        push!(HostNetworks,k)
    end
    for k in NewCons
        push!(HostCons,k)
    end
end
function OrgFitness(PEP, V, Area, PEPO ,NumProts)
    if NumProts < 11
        ProtCost = 0
    else 
        ProtCost = 1.1^(NumProts-10)
    end
    exp(-(PEP + (V*Area) + PEPO + ProtCost))
end
function ParEvolution(ParNetworks, ProgNum, ToRep, MeanSigs)
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
    for i in 1:Int64(round(ToRep))
        push!(ParNetworks, NewPars[i])
    end
end

#endregion
for inf in InfRatio
    #collection of arrays for the saving of data

    PopulationDistribution = zeros(Int64, Runs, FixationGens, 2)
    WinningNets = [] 
    WNInfRun = zeros(Float64,Runs,3) #Inf%, Run, constrained/unconstrained (1 = constrained)
    LastLoser = []
    LLInfRun = zeros(Float64,Runs,3)
    Run = 1
    Start_Time = Dates.Time(Dates.now())

    while Run <= Runs
        gen = 1
        NumPars = Int64(NumHostsUnc*inf);
        UncHostNetworks = []
        UncHostConcs = []
        UncHostCons = []

        ConHostNetworks = []
        ConHostConcs = []
        ConHostCons = []

        UncParNetworks = []
        ConParNetworks = []
        for i in 1:NumHostsUnc
            tmp = rand(0:1,(5,5))
            tmp[[1,5,21,25]].= 0
            tmp1 = (rand(Float64,(5,5))*2).-1
            tmpNet = tmp.*tmp1
            tmpInitConc = StartConc
            push!(UncHostNetworks,tmpNet)
            push!(UncHostConcs,tmpInitConc)
            push!(UncHostCons,false)

            if i <= NumPars
                partmp = zeros(5)
                partmp[rand(2:4)] = rand()*2-1
                if sum(partmp) == 0
                    partmp[2] = .5
                end
                push!(UncParNetworks,partmp)
            end
        end 
             
        ConHostNetworks = deepcopy(UncHostNetworks)
        for i in 1:NumHostsCon

            # tmp = rand(0:1,(5,5))
            # tmp[[1,5,21,25]].= 0
            # tmp1 = (rand(Float64,(5,5))*2).-1
            # tmpNet = tmp.*tmp1
             tmpInitConc = StartConc
            # push!(ConHostNetworks,tmpNet)
            push!(ConHostConcs,tmpInitConc)
            push!(ConHostCons,true)

            if i <= NumPars
                partmp = zeros(5)
                partmp[rand(2:4)] = rand()*2-1
                if sum(partmp) == 0
                    partmp[2] = .5
                end
                push!(ConParNetworks,partmp)
            end
        end

        while gen <= EqGens # this defines the first step of the evolutionary competition, prior to the two populations meeting. No data will be recorded here other than networks.
            UncHostFit = zeros(Float64, NumHostsUnc);
            ConHostFit = zeros(Float64, NumHostsCon)
            UncParFit = zeros(Float64,NumPars);
            ConParFit = zeros(Float64,NumPars);

            UncHostDead = 0
            UncHostToRem = []
            UncParDead = 0
            UncParToRem = []
            
            ConHostDead = 0
            ConHostToRem = []
            ConParDead = 0
            ConParToRem = [] 

            #region Phase One: Hosts Acheive pre-infection Equilib
            for i in 1:NumHostsUnc
                UncHostConcs[i] = AcheiveEq(UncHostNetworks[i], UncHostConcs[i], 200, false)
            end
            for i in 1:NumHostsCon
                ConHostConcs[i] = AcheiveEq(ConHostNetworks[i], ConHostConcs[i], 200, false)
            end     
            #endregion
            #region Phase Two: Infect Hosts
            UncInfHosts = []
            UncInfConcs = []
            ConInfHosts = []
            ConInfConcs = []
            InfectHosts(UncParNetworks, UncHostNetworks, UncInfHosts)
            for i in 1:NumPars
                push!(UncInfConcs, hcat(UncHostConcs[i][end,:]', .5))
            end
            InfectHosts(ConParNetworks, ConHostNetworks, ConInfHosts)
            for i in 1:NumPars
                push!(ConInfConcs, hcat(ConHostConcs[i][end,:]', .5))
            end
            #endregion
            #region Phase Three: Infected Equilib record infected immune responses
            for i in 1:length(UncInfHosts)
                UncInfConcs[i] = AcheiveEq(UncInfHosts[i], UncInfConcs[i], 20, true)
            end
            for i in 1:length(ConInfHosts)
                ConInfConcs[i] = AcheiveEq(ConInfHosts[i], ConInfConcs[i], 20, true)
            end
            #endregion 
            #region Phase Four: calculate post infection Equilib
            UncPostInfConcs = []
            ConPostInfConcs = []
            for i in 1:length(UncInfConcs)
                push!(UncPostInfConcs,AcheiveEq(UncHostNetworks[i], UncInfConcs[i][end,1:end-1]', 5, false))
                push!(ConPostInfConcs,AcheiveEq(ConHostNetworks[i], ConInfConcs[i][end,1:end-1]', 5, false))
            end
            #endregion
            #region Phase Five: Calculate Fitness
            UncHostToPop = []
            UncHostDead = 0
            ConHostToPop = []
            ConHostDead = 0
            for i in 1:NumHostsUnc
                if i <= NumPars
                    Unca = UncInfConcs[i][:,end]
                    Cona = ConInfConcs[i][:,end]
                    UncNormArea = trapz(1:length(Unca),Unca)/trapz(1:length(Unca),ones(length(Unca)))
                    ConNormArea = trapz(1:length(Cona),Cona)/trapz(1:length(Cona),ones(length(Cona)))
                    if UncNormArea > DeathThreshold
                        push!(UncHostToPop,i)
                        UncHostDead = UncHostDead + 1
                    end
                    if ConNormArea > DeathThreshold
                        push!(ConHostToPop,i)
                        ConHostDead = ConHostDead + 1
                    end

                    UncParFit[i] = UncNormArea
                    ConParFit[i] = ConNormArea
                    UncHostFit[i] = OrgFitness(UncHostConcs[i][end,end], V, UncNormArea, UncPostInfConcs[i][end,end-1], length(UncHostNetworks[i][1,:]))
                    ConHostFit[i] = OrgFitness(ConHostConcs[i][end,end], V, ConNormArea, ConPostInfConcs[i][end,end-1], length(ConHostNetworks[i][1,:]))
                else
                    #no infection so norm area becomes 0, and post inf eq is the same as the pre-inf
                    UncHostFit[i] = OrgFitness(UncHostConcs[i][end,end], V, 0, UncHostConcs[i][end,end], length(UncHostNetworks[i][1,:]))
                    ConHostFit[i] = OrgFitness(ConHostConcs[i][end,end], V, 0, ConHostConcs[i][end,end], length(ConHostNetworks[i][1,:]))
                end
            end
            while UncHostDead > NumHostsUnc*DeathCoef
                rem = rand(1:UncHostDead)
                deleteat!(UncHostToPop,rem)
                UncHostDead -= 1 
            end
            while ConHostDead > NumHostsCon*DeathCoef
                rem = rand(1:ConHostDead)
                deleteat!(ConHostToPop,rem)
                ConHostDead -= 1 
            end
            #endregion
            #region Phase Six: Cull the population in a fitness weighted manner
            Uncculled = false
            UncHost=1
            while !Uncculled
                if UncHostDead >= DeathCoef*NumHostsUnc
                    Uncculled = true
                elseif UncHost == NumHostsUnc
                    Uncculled = true
                end
                if !Uncculled 
                    Culpoint = length(findall(x->x < UncHostFit[UncHost], UncHostFit[:]))/NumHostsUnc
                    
                    if UncHostFit[UncHost] < 1e-3
                        push!(UncHostToPop,UncHost)
                        UncHostDead += 1  
                    elseif rand() > Culpoint
                        push!(UncHostToPop,UncHost)
                        UncHostDead += 1  
                    end
                    UncHost += 1
                end
            end

            Conculled = false
            ConHost=1
            while !Conculled
                if ConHostDead >= DeathCoef*NumHostsCon
                    Conculled = true
                elseif ConHost == NumHostsCon
                    Conculled = true
                end
                if !Conculled 
                    Culpoint = length(findall(x->x < ConHostFit[ConHost],ConHostFit[:]))/NumHostsCon
                    if ConHostFit[ConHost] < 1e-3
                        push!(ConHostToPop,ConHost)
                        ConHostDead += 1  
                    elseif rand() > Culpoint
                        push!(ConHostToPop,ConHost)
                        ConHostDead += 1  
                    end
                    ConHost += 1
                end
            end

            while UncHostDead > NumHostsUnc*DeathCoef
                rem = rand(1:UncHostDead)
                deleteat!(UncHostToPop,rem)
                UncHostDead -= 1 
            end
            while ConHostDead > NumHostsCon*DeathCoef
                rem = rand(1:ConHostDead)
                deleteat!(ConHostToPop,rem)
                ConHostDead -= 1 
            end

            UncHostToPop = sort(UncHostToPop, rev = true)
            for rem in UncHostToPop
                deleteat!(UncHostNetworks, rem)
                deleteat!(UncHostCons, rem)
            end
            ConHostToPop = sort(ConHostToPop, rev = true)
            for rem in ConHostToPop
                deleteat!(ConHostNetworks, rem)
                deleteat!(ConHostCons, rem)
            end
            
            UncProgNum = zeros(NumPars)
            ConProgNum = zeros(NumPars)
            for i in 1:NumPars
                if UncParFit[i] < .33
                    UncProgNum[i] = 1
                elseif .34 < UncParFit[i] < .66
                    UncProgNum[i] = 2
                else
                    UncProgNum[i] = 3
                end 
                if ConParFit[i] < .33
                    ConProgNum[i] = 1
                elseif .34 < ConParFit[i] < .66
                    ConProgNum[i] = 2
                else
                    ConProgNum[i] = 3
                end 
            end

            p = sortperm(vec(UncParFit), rev = true)
            UncParNetworks = UncParNetworks[p]
            UncProgNum = UncProgNum[p]
            UncParNetworks = UncParNetworks[1:NumPars - Int64(round(NumPars*DeathCoef))]
            
            p = sortperm(vec(UncParFit), rev = true)
            ConParNetworks = ConParNetworks[p]
            ConProgNum = ConProgNum[p]
            ConParNetworks = ConParNetworks[1:NumPars - Int64(round(NumPars*DeathCoef))]   
            #endregion
            #region Phase Seven: repopulate Hosts and parasitse allowing for evolution
            OrgEvolution(UncHostNetworks, UncHostCons, NumHostsUnc*DeathCoef);
            OrgEvolution(ConHostNetworks, ConHostCons, NumHostsCon*DeathCoef);
            UncMeanSigs = floor(Statistics.mean([length(x[1,:]) for x in UncHostNetworks]))-2;
            ConMeanSigs = floor(Statistics.mean([length(x[1,:]) for x in ConHostNetworks]))-2;
            ParEvolution(UncParNetworks, UncProgNum, NumPars*DeathCoef, Int64(UncMeanSigs));
            ParEvolution(ConParNetworks, ConProgNum, NumPars*DeathCoef, Int64(ConMeanSigs));
            
            UncHostShuffle = Random.randperm(NumHostsUnc);  
            UncHostNetworks = UncHostNetworks[UncHostShuffle];
            ConHostShuffle = Random.randperm(NumHostsCon);  
            ConHostNetworks = ConHostNetworks[ConHostShuffle];

            for i in 1:length(UncHostConcs)
                UncHostConcs[i] = zeros(2, length(UncHostNetworks[i][1,:]));
                UncHostConcs[i][2,:] .+= .5;
            end
            for i in 1:length(ConHostConcs)
                ConHostConcs[i] = zeros(2, length(ConHostNetworks[i][1,:]));
                ConHostConcs[i][2,:] .+= .5;
            end

            UncParShuffle = Random.randperm(NumPars);
            UncParNetworks = UncParNetworks[UncParShuffle];
            ConParShuffle = Random.randperm(NumPars);
            ConParNetworks = ConParNetworks[ConParShuffle];
            # if gen%50 == 0
            #     print('\n',Run, ' ', gen,' ',(Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000))
            # end
            
            gen = gen+1;
            #endregion
        end
        
        #region Combine host populations
        NumPars = NumPars*2
        HostNetworks = vcat(ConHostNetworks,UncHostNetworks)
        NumHosts = NumHostsUnc+NumHostsUnc
        HostConcs = vcat(ConHostConcs, UncHostConcs)
        HostCons = vcat(ConHostCons, UncHostCons)
        ParNetworks = vcat(ConParNetworks, UncParNetworks)

        HostShuffle = Random.randperm(NumHostsCon+NumHostsUnc)
        HostNetworks = HostNetworks[HostShuffle]
        HostConcs = HostConcs[HostShuffle]
        HostCons = HostCons[HostShuffle]
        ParShuffle = Random.randperm(NumPars)
        ParNetworks = ParNetworks[ParShuffle]
        #endregion

        Fixed = false
        CompGen = 1
        MostCommonConNet = Array{Float64,2}
        MostCommonunConNet = Array{Float64,2}

        while !Fixed
            HostDead = 0
            HostToRem = []
            ParDead = 0
            ParToRem = []
            ParFit = zeros(NumPars)
            HostFit = zeros(NumHostsUnc+NumHostsCon)

            PopDist = countmap(HostCons)
            if in(true, keys(PopDist)) & in(false, keys(PopDist))
                PopulationDistribution[Run,CompGen,1] = PopDist[true] #constrained population
                PopulationDistribution[Run,CompGen,2] = PopDist[false] #unconstrained population 

                conHosts = findall(x -> x == true, HostCons)
                UnqConNetCM = countmap(HostNetworks[conHosts])
                UnqConNets = collect(keys(UnqConNetCM))
                UnqConValues = collect(values(UnqConNetCM))
                MostCommonConNet = UnqConNets[findall(x -> x == maximum(UnqConValues), UnqConValues)[1]]

                unconHosts = findall(x -> x == false, HostCons)
                UnqunConNetCM = countmap(HostNetworks[unconHosts])
                UnqunConNets = collect(keys(UnqunConNetCM))
                UnqunConValues = collect(values(UnqunConNetCM))
                MostCommonunConNet = UnqunConNets[findall(x -> x == maximum(UnqunConValues), UnqunConValues)[1]]
            else 
                if in(true, keys(PopDist))
                    PopulationDistribution[Run,CompGen,1] = PopDist[true]
                    WNInfRun[Run,:] = [inf,Run,1] #1 indicates winner was constrained
                    LLInfRun[Run,:] = [inf,Run,0] 
                    push!(WinningNets,MostCommonConNet)
                    push!(LastLoser,MostCommonunConNet)
                else
                    PopulationDistribution[Run,CompGen,2] = PopDist[false]
                    WNInfRun[Run,:] = [inf,Run,0] #0 indicates winner was constrained
                    LLInfRun[Run,:] = [inf,Run,1] 
                    push!(LastLoser,MostCommonConNet)
                    push!(WinningNets,MostCommonunConNet)
                end
                Fixed = true
            end
            if CompGen == FixationGens
                Fixed = true
                PopulationDistribution[Run,CompGen,1] = PopDist[true]
                PopulationDistribution[Run,CompGen,2] = PopDist[false] 
            end

            #region Phase One: Hosts Acheive pre-infection Equilib
            for i in 1:NumHosts
                HostConcs[i] = AcheiveEq(HostNetworks[i], HostConcs[i], 20, false)
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
            #region Phase Three: Infected Equilib 
            for i in 1:length(InfHosts)
                InfConcs[i] = AcheiveEq(InfHosts[i], InfConcs[i], 20, true)
            end
            #endregion 
            #region Phase Four: calculate post infection Equilib
            PostInfConcs = []
            for i in 1:length(InfConcs)
                push!(PostInfConcs,AcheiveEq(HostNetworks[i], InfConcs[i][end,1:end-1]', 5, false))
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
                    ParFit[i] = NormArea
                    HostFit[i] = OrgFitness(HostConcs[i][end,end], V, NormArea, PostInfConcs[i][end,end-1], length(HostNetworks[i][1,:]))
                else
                    #no infection so norm area becomes 0, and post inf eq is the same as the pre-inf
                    HostFit[i] = OrgFitness(HostConcs[i][end,end], V, 0, HostConcs[i][end,end], length(HostNetworks[i][1,:]))
                end
            end
            while HostDead > NumHosts*DeathCoef
                rem = rand(1:HostDead)
                deleteat!(HostToPop,rem)
                HostDead -= 1 
            end
            #endregion
            #region Phase Six: Cull the population in a fitness weighted manner
            culled = false
            Host=1
            while !culled
                if HostDead >= DeathCoef*NumHosts
                    culled = true
                elseif Host == NumHosts
                    culled = true
                end
                if !culled 
                    Culpoint = length(findall(x->x < HostFit[Host],HostFit[:]))/NumHosts
                    
                    if HostFit[Host] < 1e-3
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
                deleteat!(HostCons, rem)
            end
            
            tmpParFit = ParFit[:]
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
            ParNetworks = ParNetworks[1:NumPars - Int64(round(NumPars*DeathCoef))]          
            #endregion
            #region Phase Seven: repopulate Hosts and parasitse allowing for evolution
            OrgEvolution(HostNetworks, HostCons, NumHosts*DeathCoef)
            MeanSigs = floor(Statistics.mean([length(x[1,:]) for x in HostNetworks]))-2
            ParEvolution(ParNetworks, ProgNum, NumPars*DeathCoef, Int64(MeanSigs))
            
            HostShuffle = Random.randperm(NumHosts)  
            HostNetworks = HostNetworks[HostShuffle]
            HostCons = HostCons[HostShuffle]
            for i in 1:length(HostConcs)
                HostConcs[i] = zeros(2, length(HostNetworks[i][1,:]))
                HostConcs[i][2,:] .+= .5
            end
            
            ParShuffle = Random.randperm(NumPars)
            ParNetworks = ParNetworks[ParShuffle]
            CompGen +=1
            #endregion
        end

        print('\n',Run," ******** just did ", inf, "***** With Constraint ", " Random Fixed ",EqGens, " Equilibrating Generations **********")
        Run += 1 
    end

    FileName = string(workDir,"/Data/Competition/FixedRand/",string(inf*100)[1:end-2],"_PercentInf","_",EqGens,"_EqGens_Fixed_Random_competition.jld2")

    FileIO.save(FileName,"PopulationDistribution",PopulationDistribution)
    FileIO.save(string(FileName[1:end-5],"_Networks.jld2"),"WinNets",WinningNets,"LastLoser",LastLoser,"WinInfRun",WNInfRun,"LoseInfRun",LLInfRun)
end