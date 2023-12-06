# # Step 0: Activate environment
# using Pkg
# # Pkg.activate(@__DIR__)
# # Pkg.instantiate()
# Pkg.add("Ipopt")
# Pkg.add("PowerModels")
# Pkg.add("JuMP")
# Pkg.add("Gurobi")
# Pkg.add("Juniper")
# Pkg.add("Printf")
# Pkg.add("DataFrames")
# Pkg.add("XLSX")
# using PowerModels, Ipopt, JuMP
# using Gurobi, Juniper
# using Printf, DataFrames
# using XLSX

# Define solver

## Step 1: input data
case_5 = "case24.m"
## Here we are using the parser of Powermodels for convenience 
data = PowerModels.parse_file(case_5)

GSmin = [2, 5, 8, 12, 17]
GSmax = [4, 6, 10, 14, 19]


for t in 1:1

    ## Step 2: create the JuMP model & pass data to model
    m = Model(optimizer_with_attributes(Gurobi.Optimizer))

    # Step 2a: create sets
    function define_sets!(m::Model, data::Dict)
        m.ext[:sets] = Dict()
        # nodes and branches
        N = m.ext[:sets][:N] = 1:length(data["bus"])
        NT = m.ext[:sets][:NT] = [bt for bt = 1:4]
        ND = m.ext[:sets][:ND] = [bd for bd = 5:9]

        B = m.ext[:sets][:B] = 1:length(data["branch"]) 
        BT = m.ext[:sets][:BT] = [bt for bt = 1:4]
        BD = m.ext[:sets][:BD] = [bd for bd = 5:8]

        # ac topology from side (i->j) and to side (j->i)
        B_ac_fr = m.ext[:sets][:B_ac_fr] = [(b, data["branch"][string(b)]["f_bus"], data["branch"][string(b)]["t_bus"]) for b in B] 
        B_ac_frd = m.ext[:sets][:B_ac_frd] = [(bd, data["branch"][string(bd)]["f_bus"], data["branch"][string(bd)]["t_bus"]) for bd in BD] 
        B_ac_frt = m.ext[:sets][:B_ac_frt] = [(bt, data["branch"][string(bt)]["f_bus"], data["branch"][string(bt)]["t_bus"]) for bt in BT] 
        B_ac_to = m.ext[:sets][:B_ac_to] = [(b, data["branch"][string(b)]["t_bus"], data["branch"][string(b)]["f_bus"]) for b in B]
        B_ac_tod = m.ext[:sets][:B_ac_tod] = [(bd, data["branch"][string(bd)]["t_bus"], data["branch"][string(bd)]["f_bus"]) for bd in BD]
        B_ac_tot = m.ext[:sets][:B_ac_tot] = [(bt, data["branch"][string(bt)]["t_bus"], data["branch"][string(bt)]["f_bus"]) for bt in BT]
        # Build union set of both sets above
        # B_ac = m.ext[:sets][:B_ac] = [B_ac_fr; B_ac_to]
        # Branch connectivity to buses, e.g. which branches are connected to a certain node, used in nodal power balance equations
        bus_arcs_fr = Dict((i, Tuple{Int,Int,Int}[]) for i in N)
        for (l,i,j) in B_ac_fr
            push!(bus_arcs_fr[i], (l,i,j))
        end
        B_arcs_fr = m.ext[:sets][:B_arcs_fr] = bus_arcs_fr

        bus_arcs_frd = Dict((i, Tuple{Int,Int,Int}[]) for i in ND)
        for (l,i,j) in B_ac_frd
            if i > 3
                push!(bus_arcs_frd[i], (l,i,j))  
            end
        end
        B_arcs_frd = m.ext[:sets][:B_arcs_frd] = bus_arcs_frd

        bus_arcs_frt = Dict((i, Tuple{Int,Int,Int}[]) for i in NT)
        for (l,i,j) in B_ac_frt
            push!(bus_arcs_frt[i], (l,i,j))
        end
        B_arcs_frt = m.ext[:sets][:B_arcs_frt] = bus_arcs_frt

        bus_arcs_to = Dict((i, Tuple{Int,Int,Int}[]) for i in N)
        for (l,i,j) in B_ac_to
            push!(bus_arcs_to[i], (l,j,i))
        end
        B_arcs_to = m.ext[:sets][:B_arcs_to] = bus_arcs_to

        bus_arcs_tod = Dict((i, Tuple{Int,Int,Int}[]) for i in ND)
        for (l,i,j) in B_ac_tod
            push!(bus_arcs_tod[i], (l,j,i))
        end
        B_arcs_tod = m.ext[:sets][:B_arcs_tod] = bus_arcs_tod
        
        bus_arcs_tot = Dict((i, Tuple{Int,Int,Int}[]) for i in NT)
        for (l,i,j) in B_ac_tot
            if i < 4
                push!(bus_arcs_tot[i], (l,j,i))
            end
        end
        B_arcs_tot = m.ext[:sets][:B_arcs_tot] = bus_arcs_tot

        gsmin = GSmin[t]
        gsmax = GSmax[t]

        # Set of generators
        G = m.ext[:sets][:G] = 1:length(data["gen"])
        # GS = m.ext[:sets][:GS] = []
        GS = m.ext[:sets][:GS] = [gs for gs = gsmin:gsmax]
        GC = m.ext[:sets][:GC] = symdiff(G, GS)
        GD = m.ext[:sets][:GD] = [gd for gd = 21:50]
        GT = m.ext[:sets][:GT] = [gt for gt = 1:45]
        GDOD = m.ext[:sets][:GDOD] = [gdod for gdod = 46:49]
        GUPD = m.ext[:sets][:GUPD] = symdiff(GD,GDOD)
        GUP = m.ext[:sets][:GUP] = symdiff(G,GDOD)
        GUPT = m.ext[:sets][:GUPT] = [gupt for gupt = 1:20]
        GDOT = m.ext[:sets][:GDOT] = symdiff(GT,GUPT)

        GDODS = m.ext[:sets][:GDODS] = intersect(GDOD,GS)
        GUPDS = m.ext[:sets][:GUPDS] = intersect(GUPD,GS)
        GUPS = m.ext[:sets][:GUPS] = intersect(GUP,GS)
        GUPC = m.ext[:sets][:GUPC] = intersect(GUP,GC)
        GDODC = m.ext[:sets][:GDODC] = intersect(GDOD,GC)
        GUPDC = m.ext[:sets][:GUPDC] = intersect(GUPD,GC)
        GDTS = m.ext[:sets][:GDTS] = intersect(GUPT,GS)
        GDTC = m.ext[:sets][:GDTC] = intersect(GUPT,GC)
        GUPTS = m.ext[:sets][:GUPTS] = intersect(GUPT,GS)
        GUPTC = m.ext[:sets][:GUPTC] = intersect(GUPT,GC)
        GDOTS = m.ext[:sets][:GDOTS] = intersect(GDOT,GS)
        GDOTC = m.ext[:sets][:GDOTC] = intersect(GDOT,GC)

        # new_set = Dict{Int, Array{Int}}

        # Generator connectivity, e.g. which generators are connected to a certain node, used in nodal power balance equations
        bus_gens = Dict((i, Int[]) for i in N)
        for (g, gen) in data["gen"]
            push!(bus_gens[gen["gen_bus"]], parse(Int, g))
        end
        G_ac = m.ext[:sets][:G_ac] = bus_gens
        # which nodes are connected to a certain generator

        new_setup = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsup = filter(x -> x in GUPD, generators)
            if !isempty(new_generatorsup)
                new_setup[bus] = new_generatorsup
            end
        end
        G_acupd = m.ext[:sets][:G_acupd] = new_setup

        new_setdo = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsdo = filter(x -> x in GDOD, generators)
            if !isempty(new_generatorsdo)
                new_setdo[bus] = new_generatorsdo
            end
        end
        G_acdod = m.ext[:sets][:G_acdod] = new_setdo

        new_setupt = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsupt = filter(x -> x in GUPT, generators)
            if !isempty(new_generatorsupt)
                new_setupt[bus] = new_generatorsupt
            end
        end
        G_acupt = m.ext[:sets][:G_acupt] = new_setupt

        new_setuptc = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsuptc = filter(x -> x in GUPTC, generators)
            if !isempty(new_generatorsuptc)
                new_setuptc[bus] = new_generatorsuptc
            end
        end
        G_acuptc = m.ext[:sets][:G_acuptc] = new_setuptc

        new_setupts = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsupts = filter(x -> x in GUPTS, generators)
            if !isempty(new_generatorsupts)
                new_setupts[bus] = new_generatorsupts
            end
        end
        G_acupts = m.ext[:sets][:G_acupts] = new_setupts

        new_setdot = Dict((i, Int[]) for i in N)
        for (bus, generators) in m.ext[:sets][:G_ac]
            new_generatorsdot = filter(x -> x in GDOT, generators)
            if !isempty(new_generatorsdot)
                new_setdot[bus] = new_generatorsdot
            end
        end
        G_acdot = m.ext[:sets][:G_acdot] = new_setdot

        gens_bus = Dict((g, Int[]) for g in G)
        for g in 1:length(data["gen"])
            push!(gens_bus[g], data["gen"][string(g)]["gen_bus"])
        end
        N_ac = m.ext[:sets][:N_ac] = gens_bus

        gens_buss = Dict((gs, Int[]) for gs in GT)
        for gs in GT
            push!(gens_buss[gs], data["gen"][string(gs)]["gen_bus"])
        end
        N_act = m.ext[:sets][:N_act] = gens_buss

        gens_busd = Dict((gs, Int[]) for gs in GD)
        for gs in GD
            push!(gens_busd[gs], data["gen"][string(gs)]["gen_bus"])
        end
        N_acd= m.ext[:sets][:N_acd] = gens_busd

        # Set of loads
        L = m.ext[:sets][:L] = 1:length(data["load"])
        # Load connectivity, e.g. which loads are connected to a certain node, used in nodal power balance equations
        bus_loads = Dict((i, Int[]) for i in N)
        for (l, load) in data["load"]
            push!(bus_loads[load["load_bus"]], parse(Int, l))
        end
        L_ac = m.ext[:sets][:L_ac] = bus_loads 

        return m
    end
    define_sets!(m, data)

    # step 2b: process input parameters
    function process_parameters!(m::Model, data::Dict)
        # extract sets
        N = m.ext[:sets][:N]
        B = m.ext[:sets][:B]
        G = m.ext[:sets][:G]
        L = m.ext[:sets][:L]

        # Create parameter dictionary
        m.ext[:parameters] = Dict()

        # Bus parameters
        vmmin = m.ext[:parameters][:vmmin] = Dict(i => data["bus"][string(i)]["vmin"] for i in N) # minimum voltage magnitude
        vmmax = m.ext[:parameters][:vmmax] = Dict(i => data["bus"][string(i)]["vmax"] for i in N) # maximum voltage magnitude

        # Branch parameters
        br = m.ext[:parameters][:br] = Dict(b => data["branch"][string(b)]["br_r"] for b in B) # branch resistance
        bx = m.ext[:parameters][:bx] = Dict(b => data["branch"][string(b)]["br_x"] for b in B) # branch reactance
        pbmax = m.ext[:parameters][:pbmax] = [30000, 2000, 26000, 15000, 100,40,100,100]/100
        qbmax = m.ext[:parameters][:qbmax] = [30000, 2000, 26000, 21000, 75,30,75,75]/100

        # load parameters: Assuming a fixed demand!
        pd = m.ext[:parameters][:pd] = Dict(l => data["load"][string(l)]["pd"] for l in L)  # active power demand in pu
        qd = m.ext[:parameters][:qd] = Dict(l => data["load"][string(l)]["qd"] for l in L)  # reactive power demand in pu
        
        # generator parameters:
        pmax = m.ext[:parameters][:pmax] = Dict(g => data["gen"][string(g)]["pmax"] for g in G)  # maximum active power in pu
        qmax = m.ext[:parameters][:qmax] = Dict(g => data["gen"][string(g)]["qmax"] for g in G)  # maximum reactive power in pu
        cg = m.ext[:parameters][:cg] = Dict(g => data["gen"][string(g)]["mbase"] for g in G) # linear generation cost

        M = m.ext[:parameters][:M] = 600

        k1 = m.ext[:parameters][:k1] = 20
        M1 = m.ext[:parameters][:M1] = 2^k1
        gamma = m.ext[:parameters][:gamma] = 10/M1
        GG = m.ext[:parameters][:GG] = 1000
        K = m.ext[:sets][:K] = [k for k = 1:k1]

        return m
    end

    # call functions
    process_parameters!(m, data)

    ## Step 3: construct your model
    function build_model!(m::Model)
        # Create m.ext entries "variables", "expressions" and "constraints"
        m.ext[:variables] = Dict()
        m.ext[:expressions] = Dict()
        m.ext[:constraints] = Dict()
        m.ext[:objectives] = Dict()

        # Extract sets
        ND = m.ext[:sets][:ND]
        NT = m.ext[:sets][:NT]
        N_acd = m.ext[:sets][:N_acd]
        N_act = m.ext[:sets][:N_act]

        B_ac_fr = m.ext[:sets][:B_ac_fr]
        B_ac_frd = m.ext[:sets][:B_ac_frd] 
        B_ac_frt = m.ext[:sets][:B_ac_frt] 
        B_arcs_frd = m.ext[:sets][:B_arcs_frd]
        B_arcs_frt = m.ext[:sets][:B_arcs_frt]
        B_arcs_tod = m.ext[:sets][:B_arcs_tod]
        B_arcs_tot = m.ext[:sets][:B_arcs_tot]

        G = m.ext[:sets][:G]
        GS = m.ext[:sets][:GS]
        GC = m.ext[:sets][:GC]
        GUPD = m.ext[:sets][:GUPD] 
        GDOD = m.ext[:sets][:GDOD]       
        GUP = m.ext[:sets][:GUP]    
        GUPT = m.ext[:sets][:GUPT] 
        GDOT = m.ext[:sets][:GDOT]
        GD = m.ext[:sets][:GD] 
        GDODC = m.ext[:sets][:GDODC]
        GUPDC = m.ext[:sets][:GUPDC] 
        GDODS = m.ext[:sets][:GDODS] 
        GUPDS = m.ext[:sets][:GUPDS] 
        GUPTS = m.ext[:sets][:GUPTS]
        GUPTC = m.ext[:sets][:GUPTC]
        GDOTS = m.ext[:sets][:GDOTS]
        GDOTC = m.ext[:sets][:GDOTC]
        GUPS = m.ext[:sets][:GUPS] 
        GUPC = m.ext[:sets][:GUPC] 
        GDTS = m.ext[:sets][:GDTS] 
        GDTC = m.ext[:sets][:GDTC] 
        G_ac = m.ext[:sets][:G_ac]
        G_acupd = m.ext[:sets][:G_acupd]
        G_acdod = m.ext[:sets][:G_acdod]
        G_acupt = m.ext[:sets][:G_acupt]
        G_acupts = m.ext[:sets][:G_acupts]
        G_acuptc = m.ext[:sets][:G_acuptc]
        G_acdot = m.ext[:sets][:G_acdot]

        L = m.ext[:sets][:L]
        L_ac = m.ext[:sets][:L_ac]

        K = m.ext[:sets][:K]

        # Extract parameters
        vmmin = m.ext[:parameters][:vmmin]
        vmmax = m.ext[:parameters][:vmmax]
        pd = m.ext[:parameters][:pd]
        qd = m.ext[:parameters][:qd]
        pmax = m.ext[:parameters][:pmax]
        qbmax = m.ext[:parameters][:qbmax]
        pbmax = m.ext[:parameters][:pbmax]
        qmax = m.ext[:parameters][:qmax]
        br = m.ext[:parameters][:br]
        bx = m.ext[:parameters][:bx]
        cg = m.ext[:parameters][:cg]
        M = m.ext[:parameters][:M]

        k1 = m.ext[:parameters][:k1]
        M1 = m.ext[:parameters][:M1] 
        gamma = m.ext[:parameters][:gamma]
        GG = m.ext[:parameters][:GG] 

        ## create variables 
        # wholesale market variables
        pw = m.ext[:variables][:pw] = @variable(m, [g=G], lower_bound = 0 , base_name = "pw")
        mu_wsmup = m.ext[:variables][:mu_wsmup] = @variable(m, [g in G], lower_bound = 0, base_name = "mu_wsmup")
        lambdaw = m.ext[:variables][:lambdaw] = @variable(m, base_name = "lambdaw")
        choedw = m.ext[:variables][:choedw] = @variable(m, [gs=GS], lower_bound = 0, upper_bound = 100, base_name = "choedw")
        
        b11 = m.ext[:variables][:b11] = @variable(m, [g in G], Bin, base_name = "b11")
        b12 = m.ext[:variables][:b12] = @variable(m, [g in G], Bin, base_name = "b12")

        # flexibility market variables
        pb = m.ext[:variables][:pb] = @variable(m, [(b,i,j) in B_ac_fr], base_name = "pb") # from side active power flow (i->j)
        qbd = m.ext[:variables][:qbd] = @variable(m, [(b,i,j) in B_ac_frd], base_name = "qb") # from side reactive power flow (i->j)
        qg = m.ext[:variables][:qg] = @variable(m, [g=GD], base_name = "qg")
        wm = m.ext[:variables][:wm] = @variable(m, [i=ND], base_name = "wm") # voltage magnitude
        
        lambdaf = m.ext[:variables][:lambdaf] = @variable(m, [i=ND], base_name = "lambdaf")
        pfup = m.ext[:variables][:pfup] = @variable(m, [g=GUPC], lower_bound = 0 , base_name = "pfup")
        pfdo = m.ext[:variables][:pfdo] = @variable(m, [g=GDOD], lower_bound = 0 , base_name = "pfdo")
        choedfup = m.ext[:variables][:choedfup] = @variable(m, [gs=GUPS], lower_bound = 0, upper_bound = 100, base_name = "choedfup")
        choedfdo = m.ext[:variables][:choedfdo] = @variable(m, [gs=GDODS], lower_bound = 0, upper_bound = 100, base_name = "choedfdo")

        xfup = m.ext[:variables][:xfup] = @variable(m, [gs=GUPS,k=K], Bin, base_name = "xfup")
        pfups = m.ext[:expressions][:pfups] = @expression(m, [g=GUPS], gamma*sum((2^k)*xfup[g,k] for k in K))

        mu_fup = m.ext[:variables][:mu_fup] = @variable(m, [g in GUP], lower_bound = 0, base_name = "mu_fup")
        mu_fdo = m.ext[:variables][:mu_fdo] = @variable(m, [g in GDOD], lower_bound = 0, base_name = "mu_fdo")
        mu_qup = m.ext[:variables][:mu_qup] = @variable(m, [g in GD], lower_bound = 0, base_name = "mu_qup")
        mu_pbup = m.ext[:variables][:mu_pbup] = @variable(m, [(b,i,j) = B_ac_fr], lower_bound = 0, base_name = "mu_pbup")
        mu_qbup = m.ext[:variables][:mu_qbup] = @variable(m, [(b,i,j) = B_ac_frd], lower_bound = 0, base_name = "mu_pbup")

        mu_wlo = m.ext[:variables][:mu_wlo] = @variable(m, [i in ND], lower_bound = 0, base_name = "mu_wlo")
        mu_wup = m.ext[:variables][:mu_wup] = @variable(m, [i in ND], lower_bound = 0, base_name = "mu_wup")
        mu_f = m.ext[:variables][:mu_f] = @variable(m, [(b,i,j) = B_ac_frd], base_name = "mu_f")
        mu_wo = m.ext[:variables][:mu_wo] = @variable(m, [i in ND], base_name = "mu_wo")
        mu_qb = m.ext[:variables][:mu_qb] = @variable(m, [i in ND], base_name = "mu_qb")

        b1 = m.ext[:variables][:b1] = @variable(m, [g in GUP], Bin, base_name = "b1")
        b2 = m.ext[:variables][:b2] = @variable(m, [g in GUP], Bin, base_name = "b2")
        b3 = m.ext[:variables][:b3] = @variable(m, [g in GDOD], Bin, base_name = "b3")
        b4 = m.ext[:variables][:b4] = @variable(m, [g in GDOD], Bin, base_name = "b4")
        b5 = m.ext[:variables][:b5] = @variable(m, [g in GD], Bin, base_name = "b5")
        b6 = m.ext[:variables][:b6] = @variable(m, [g in GD], Bin, base_name = "b6")
        b8 = m.ext[:variables][:b8] = @variable(m, [i in ND], Bin, base_name = "b8")
        b9 = m.ext[:variables][:b9] = @variable(m, [i in ND], Bin, base_name = "b9")
        b14 = m.ext[:variables][:b14] = @variable(m, [(b,i,j) = B_ac_frd], Bin, base_name = "b14")
        b15 = m.ext[:variables][:b15] = @variable(m, [(b,i,j) = B_ac_frd], Bin, base_name = "b15")
        b16 = m.ext[:variables][:b16] = @variable(m, [(b,i,j) = B_ac_frd], Bin, base_name = "b16")
        b17 = m.ext[:variables][:b17] = @variable(m, [(b,i,j) = B_ac_frd], Bin, base_name = "b17")

        # redispatch market variables
        lambdar = m.ext[:variables][:lambdar] = @variable(m, [i=NT], base_name = "lambdaf")
        prup = m.ext[:variables][:prup] = @variable(m, [g=GUPT], lower_bound = 0 , base_name = "prup")
        prdo = m.ext[:variables][:prdo] = @variable(m, [g=GDOT], lower_bound = 0 , base_name = "prdo")
        choedrup = m.ext[:variables][:choedrup] = @variable(m, [gs=GUPTS], lower_bound = 0, upper_bound = 100, base_name = "choedrup")
        choedrdo = m.ext[:variables][:choedrdo] = @variable(m, [gs=GDOTS], lower_bound = 0, upper_bound = 100, base_name = "choedrdo")

        xrup = m.ext[:variables][:xrup] = @variable(m, [gs=GUPTS,k=K], Bin, base_name = "xrup")
        prups = m.ext[:expressions][:prups] = @expression(m, [g=GUPTS], gamma*sum((2^k)*xrup[g,k] for k in K))

        mu_rup = m.ext[:variables][:mu_fup] = @variable(m, [g in GUPT], lower_bound = 0, base_name = "mu_fup")
        mu_rdo = m.ext[:variables][:mu_fdo] = @variable(m, [g in GDOT], lower_bound = 0, base_name = "mu_fdo")

        br1 = m.ext[:variables][:br1] = @variable(m, [g in GUPT], Bin, base_name = "b1")
        br2 = m.ext[:variables][:br2] = @variable(m, [g in GUPT], Bin, base_name = "b2")
        br3 = m.ext[:variables][:br3] = @variable(m, [g in GDOT], Bin, base_name = "b3")
        br4 = m.ext[:variables][:br4] = @variable(m, [g in GDOT], Bin, base_name = "b4")

        b14r = m.ext[:variables][:b14r] = @variable(m, [(b,i,j) = B_ac_frt], Bin, base_name = "br14")
        b15r = m.ext[:variables][:b15r] = @variable(m, [(b,i,j) = B_ac_frt], Bin, base_name = "br15")

        # # Objective strategic agent

        m.ext[:objectives][:Obj] = @objective(m, Min, sum(cg[gs]*pw[gs] for gs in GS) + sum(cg[gc]*pw[gc] for gc in GC) - lambdaw*sum(pd[l] for l in L) + sum(pmax[gc]*mu_wsmup[gc] for gc in GC) - sum(lambdaf[5]*gamma*sum((2^k)*xfup[gs,k] for k in K) for gs in GUPS) + sum(cg[gs]*gamma*sum((2^k)*xfup[gs,k] for k in K) for gs in GUPS) - sum(lambdar[2]*gamma*sum((2^k)*xrup[gs,k] for k in K) for gs in GUPTS) + sum(cg[gs]*gamma*sum((2^k)*xrup[gs,k] for k in K) for gs in GUPTS))

        # # Wholesale market constraints 
        m.ext[:constraints][:F11] = @constraint(m, [g=GC], 0 <= cg[g] + mu_wsmup[g] - lambdaw)
        m.ext[:constraints][:F11a] = @constraint(m, [g=GC], pw[g] <= M*(1-b11[g]))
        m.ext[:constraints][:F11b] = @constraint(m, [g=GC],  cg[g] + mu_wsmup[g] - lambdaw <= M*b11[g])
        m.ext[:constraints][:F11s] = @constraint(m, [g=GS], 0 <= choedw[g] + mu_wsmup[g] - lambdaw)
        m.ext[:constraints][:F11as] = @constraint(m, [g=GS], pw[g] <= M*(1-b11[g]))
        m.ext[:constraints][:F11bs] = @constraint(m, [g=GS],  choedw[g] + mu_wsmup[g] - lambdaw <= M*b11[g])
        m.ext[:constraints][:F12] = @constraint(m, [g=G], 0 <= pmax[g] - pw[g])
        m.ext[:constraints][:F12a] = @constraint(m, [g=G], mu_wsmup[g] <= M*(1-b12[g]))
        m.ext[:constraints][:F12b] = @constraint(m, [g=G], pmax[g] - pw[g]<= M*b12[g])

        m.ext[:constraints][:F13] = @constraint(m, - sum(pw[g] for g in G) + sum(pd[l] for l in L) == 0)
        
        # # Flexibility market constraints 
        m.ext[:constraints][:F1td] = @constraint(m, [g=GUPTC], 0 <= cg[g] + mu_fup[g] - lambdaf[5])
        m.ext[:constraints][:F1atd] = @constraint(m, [g=GUPTC], pfup[g] <= M*(1-b1[g]))
        m.ext[:constraints][:F1btd] = @constraint(m, [g=GUPTC], cg[g] + mu_fup[g] - lambdaf[5] <= M*b1[g])
        m.ext[:constraints][:F1] = @constraint(m, [g=GUPDC], 0 <= cg[g] + mu_fup[g] - sum(lambdaf[i] for i in N_acd[g]))
        m.ext[:constraints][:F1a] = @constraint(m, [g=GUPDC], pfup[g] <= M*(1-b1[g]))
        m.ext[:constraints][:F1b] = @constraint(m, [g=GUPDC], cg[g] + mu_fup[g] - sum(lambdaf[i] for i in N_acd[g]) <= M*b1[g])
        m.ext[:constraints][:F2] = @constraint(m, [g=GUPC], 0 <= pmax[g]- pw[g] - pfup[g])
        m.ext[:constraints][:F2a] = @constraint(m, [g=GUPC], mu_fup[g]  <= M*(1-b2[g]))
        m.ext[:constraints][:F2b] = @constraint(m, [g=GUPC], pmax[g]- pw[g] - pfup[g] <= M*b2[g])

        m.ext[:constraints][:F1tds] = @constraint(m, [g=GUPTS], 0 <= choedfup[g] + mu_fup[g] - lambdaf[5])
        m.ext[:constraints][:F1atds] = @constraint(m, [g=GUPTS], gamma*sum((2^k)*xfup[g,k] for k in K) <= M*(1-b1[g]))
        m.ext[:constraints][:F1btds] = @constraint(m, [g=GUPTS], choedfup[g] + mu_fup[g] - lambdaf[5] <= M*b1[g])
        m.ext[:constraints][:F1s] = @constraint(m, [g=GUPDS], 0 <= choedfup[g] + mu_fup[g] - sum(lambdaf[i] for i in N_acd[g]))
        m.ext[:constraints][:F1as] = @constraint(m, [g=GUPDS], gamma*sum((2^k)*xfup[g,k] for k in K)<= M*(1-b1[g]))
        m.ext[:constraints][:F1bs] = @constraint(m, [g=GUPDS], choedfup[g] + mu_fup[g] - sum(lambdaf[i] for i in N_acd[g]) <= M*b1[g])
        m.ext[:constraints][:F2s] = @constraint(m, [g=GUPS], 0 <= pmax[g]- pw[g] - gamma*sum((2^k)*xfup[g,k] for k in K))
        m.ext[:constraints][:F2as] = @constraint(m, [g=GUPS], mu_fup[g]  <= M*(1-b2[g]))
        m.ext[:constraints][:F2bs] = @constraint(m, [g=GUPS], pmax[g]- pw[g] - gamma*sum((2^k)*xfup[g,k] for k in K) <= M*b2[g])

        m.ext[:constraints][:F3] = @constraint(m, [g=GDODC], 0 <= -cg[g] + mu_fdo[g] + sum(lambdaf[i] for i in N_acd[g]))
        m.ext[:constraints][:F3a] = @constraint(m, [g=GDODC], pfdo[g] <= M*(1-b3[g]))
        m.ext[:constraints][:F3b] = @constraint(m, [g=GDODC], -cg[g] + mu_fdo[g] + sum(lambdaf[i] for i in N_acd[g]) <= M*b3[g])
        m.ext[:constraints][:F3s] = @constraint(m, [g=GDODS], 0 <= -choedfdo[g] + mu_fdo[g] + sum(lambdaf[i] for i in N_acd[g]))
        m.ext[:constraints][:F3as] = @constraint(m, [g=GDODS], pfdo[g] <= M*(1-b3[g]))
        m.ext[:constraints][:F3bs] = @constraint(m, [g=GDODS], -choedfdo[g] + mu_fdo[g] + sum(lambdaf[i] for i in N_acd[g]) <= M*b3[g])
        m.ext[:constraints][:F4] = @constraint(m, [g=GDOD], 0 <= pw[g] - pfdo[g])
        m.ext[:constraints][:F4a] = @constraint(m, [g=GDOD], mu_fdo[g]  <= M*(1-b4[g]))
        m.ext[:constraints][:F4b] = @constraint(m, [g=GDOD], pw[g] - pfdo[g] <= M*b4[g])

        m.ext[:constraints][:F5] = @constraint(m, [g=GD], 0 <= qg[g] + qmax[g])
        m.ext[:constraints][:F5a] = @constraint(m, [g=GD], 0 <= mu_qup[g] + sum(mu_qb[i] for i in N_acd[g]))
        m.ext[:constraints][:F5b] = @constraint(m, [g=GD], qg[g] + qmax[g] <= M*(1-b5[g]))
        m.ext[:constraints][:F5c] = @constraint(m, [g=GD], mu_qup[g] + sum(mu_qb[i] for i in N_acd[g]) <= M*b5[g])
        m.ext[:constraints][:F6] = @constraint(m, [g=GD], 0 <= qmax[g]- qg[g])
        m.ext[:constraints][:F6a] = @constraint(m, [g=GD], mu_qup[g]  <= M*(1-b6[g]))
        m.ext[:constraints][:F6b] = @constraint(m, [g=GD], qmax[g]- qg[g] <= M*b6[g])

        m.ext[:constraints][:F14] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= lambdaf[i] - lambdaf[j] - 2*br[b]*mu_f[(b,i,j)] + mu_pbup[(b,i,j)])
        m.ext[:constraints][:F14a] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= pb[(b, i, j)] + pbmax[b])
        m.ext[:constraints][:F14b] = @constraint(m, [(b,i,j) = B_ac_frd], lambdaf[i] - lambdaf[j] - 2*br[b]*mu_f[(b,i,j)] + mu_pbup[(b,i,j)] <= M*(1-b14[(b,i,j)]))
        m.ext[:constraints][:F14c] = @constraint(m, [(b,i,j) = B_ac_frd], pb[(b, i, j)] + pbmax[b] <= M*b14[(b,i,j)])
        m.ext[:constraints][:F15] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= - pb[(b, i, j)] + pbmax[b])
        m.ext[:constraints][:F15a] = @constraint(m, [(b,i,j) = B_ac_frd], mu_pbup[(b,i,j)] <= M*(1-b15[(b,i,j)]))
        m.ext[:constraints][:F15b] = @constraint(m, [(b,i,j) = B_ac_frd], - pb[(b, i, j)] + pbmax[b] <= M*b15[(b,i,j)])

        m.ext[:constraints][:F16] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= - mu_qb[i] + mu_qb[j] - 2*bx[b]*mu_f[(b,i,j)] + mu_qbup[(b,i,j)])
        m.ext[:constraints][:F16a] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= qbd[(b, i, j)] + qbmax[b])
        m.ext[:constraints][:F16b] = @constraint(m, [(b,i,j) = B_ac_frd], - mu_qb[i] + mu_qb[j] - 2*bx[b]*mu_f[(b,i,j)] + mu_qbup[(b,i,j)] <= M*(1-b16[(b,i,j)]))
        m.ext[:constraints][:F16c] = @constraint(m, [(b,i,j) = B_ac_frd], qbd[(b, i, j)] + qbmax[b] <= M*b16[(b,i,j)])
        m.ext[:constraints][:F17] = @constraint(m, [(b,i,j) = B_ac_frd], 0 <= - qbd[(b, i, j)] + qbmax[b])
        m.ext[:constraints][:F17a] = @constraint(m, [(b,i,j) = B_ac_frd], mu_qbup[(b,i,j)] <= M*(1-b17[(b,i,j)]))
        m.ext[:constraints][:F17b] = @constraint(m, [(b,i,j) = B_ac_frd], - qbd[(b, i, j)] + qbmax[b] <= M*b17[(b,i,j)])

        m.ext[:constraints][:F8] = @constraint(m, [i=ND], 0 <= wm[i] - vmmin[i]^2)
        m.ext[:constraints][:F8a] = @constraint(m, [i=ND], mu_wlo[i]<= M*b8[i])
        m.ext[:constraints][:F8b] = @constraint(m, [i=ND], wm[i] - vmmin[i]^2 <= M*(1-b8[i]))
        m.ext[:constraints][:F9] = @constraint(m, [i=ND], 0 <= vmmax[i]^2- wm[i])
        m.ext[:constraints][:F9a] = @constraint(m, [i=ND], mu_wup[i]<= M*b9[i])
        m.ext[:constraints][:F9b] = @constraint(m, [i=ND], vmmax[i]^2- wm[i] <= M*(1-b9[i]))

        m.ext[:constraints][:magnij] = @constraint(m, [(b,i,j) = B_ac_frd], wm[i] - wm[j] == 2*(br[b]*pb[(b, i, j)]+bx[b]*qbd[(b, i, j)] ))
        m.ext[:constraints][:vmref] = @constraint(m, wm[5] == 1)
        m.ext[:constraints][:muwo5] = @constraint(m, mu_wo[6] == 0)
        m.ext[:constraints][:muwo6] = @constraint(m, mu_wo[7] == 0)
        m.ext[:constraints][:muwo7] = @constraint(m, mu_wo[8] == 0)
        m.ext[:constraints][:muwo8] = @constraint(m, mu_wo[9] == 0)
        m.ext[:constraints][:afgw] = @constraint(m, [i = ND], sum(mu_f[(b,i,j)] for (b,i,j) in B_arcs_frd[i]) - sum(mu_f[(b,i,j)] for (b,i,j) in B_arcs_tod[i]) - mu_wlo[i] + mu_wup[i] + mu_wo[i] == 0)

        m.ext[:constraints][:p_balanced5] = @constraint(m, -( sum(pw[g] + pfup[g] for g in GUPTC) + sum(pw[g] + gamma*sum((2^k)*xfup[g,k] for k in K) for g in GUPTS) - pd[1]) -sum(pw[g] for g in G_ac[5]) - sum(pfup[g] for g in G_acupd[5]) + sum(pfdo[g] for g in G_acdod[5]) + sum(pd[l] for l in L_ac[5]) + sum(pb[(b,i,j)] for (b,i,j) in B_arcs_frd[5]) - sum(pb[(b,i,j)] for (b,i,j) in B_arcs_tod[5]) == 0)
        m.ext[:constraints][:p_balance] = @constraint(m, [i in 6:9], -sum(pw[g] for g in G_ac[i]) - sum(pfup[g] for g in G_acupd[i]) + sum(pfdo[g] for g in G_acdod[i]) + sum(pd[l] for l in L_ac[i]) + sum(pb[(b,i,j)] for (b,i,j) in B_arcs_frd[i]) - sum(pb[(b,i,j)] for (b,i,j) in B_arcs_tod[i]) == 0)
        m.ext[:constraints][:q_balance] = @constraint(m, [i in ND], + sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) == sum(qbd[(b,i,j)] for (b,i,j) in B_arcs_frd[i]) - sum(qbd[(b,i,j)] for (b,i,j) in B_arcs_tod[i]))

        # # Redispatch market constraints 
        m.ext[:constraints][:Fr1] = @constraint(m, [g=GUPTC], 0 <= cg[g] + mu_rup[g] - sum(lambdar[i] for i in N_act[g]))
        m.ext[:constraints][:Fr1a] = @constraint(m, [g=GUPTC], prup[g] <= M*(1-br1[g]))
        m.ext[:constraints][:Fr1b] = @constraint(m, [g=GUPTC], cg[g] + mu_rup[g] - sum(lambdar[i] for i in N_act[g]) <= M*br1[g])
        m.ext[:constraints][:Fr2up] = @constraint(m, [g=GUPTC], 0 <= pmax[g]- pw[g] - pfup[g] - prup[g])
        m.ext[:constraints][:Fr2a] = @constraint(m, [g=GUPTC], mu_rup[g]  <= M*(1-br2[g]))
        m.ext[:constraints][:Fr2bup] = @constraint(m, [g=GUPTC], pmax[g]- pw[g] - pfup[g] - prup[g] <= M*br2[g])

        m.ext[:constraints][:Fr1s] = @constraint(m, [g=GUPTS], 0 <= choedrup[g] + mu_rup[g] - sum(lambdar[i] for i in N_act[g]))
        m.ext[:constraints][:Fr1as] = @constraint(m, [g=GUPTS], gamma*sum((2^k)*xrup[g,k] for k in K) <= M*(1-br1[g]))
        m.ext[:constraints][:Fr1bs] = @constraint(m, [g=GUPTS], choedrup[g] + mu_rup[g] - sum(lambdar[i] for i in N_act[g]) <= M*br1[g])
        m.ext[:constraints][:Fr2ups] = @constraint(m, [g=GUPTS], 0 <= pmax[g]- pw[g] - gamma*sum((2^k)*xfup[g,k] for k in K) - gamma*sum((2^k)*xrup[g,k] for k in K))
        m.ext[:constraints][:Fr2as] = @constraint(m, [g=GUPTS], mu_rup[g]  <= M*(1-br2[g]))
        m.ext[:constraints][:Fr2bups] = @constraint(m, [g=GUPTS], pmax[g]- pw[g] - gamma*sum((2^k)*xfup[g,k] for k in K) - gamma*sum((2^k)*xrup[g,k] for k in K) <= M*br2[g])


        m.ext[:constraints][:Fr3] = @constraint(m, [g=GDOTC], 0 <= -cg[g] + mu_rdo[g] + sum(lambdar[4]))
        m.ext[:constraints][:Fr3a] = @constraint(m, [g=GDOTC], prdo[g] <= M*(1-br3[g]))
        m.ext[:constraints][:Fr3b] = @constraint(m, [g=GDOTC], -cg[g] + mu_rdo[g] + sum(lambdar[4]) <= M*br3[g])
        m.ext[:constraints][:Fr3s] = @constraint(m, [g=GDOTS], 0 <= -choedrdo[g] + mu_rdo[g] + sum(lambdar[4]))
        m.ext[:constraints][:Fr3as] = @constraint(m, [g=GDOTS], prdo[g] <= M*(1-br3[g]))
        m.ext[:constraints][:Fr3bs] = @constraint(m, [g=GDOTS], -choedrdo[g] + mu_rdo[g] + sum(lambdar[4]) <= M*br3[g])
        m.ext[:constraints][:Fr4do] = @constraint(m, [g=GDOT], 0 <= pw[g] + pfup[g] - prdo[g])
        m.ext[:constraints][:Fr4a] = @constraint(m, [g=GDOT], mu_rdo[g]  <= M*(1-br4[g]))
        m.ext[:constraints][:Fr4bdo] = @constraint(m, [g=GDOT], pw[g] + pfup[g] - prdo[g] <= M*br4[g])

        m.ext[:constraints][:Fr14] = @constraint(m, [(b,i,j) = B_ac_frt], 0 <= lambdar[i] - lambdar[j] + mu_pbup[(b,i,j)])
        m.ext[:constraints][:Fr14a] = @constraint(m, [(b,i,j) = B_ac_frt], 0 <= pb[(b, i, j)] + pbmax[b])
        m.ext[:constraints][:Fr14b] = @constraint(m, [(b,i,j) = B_ac_frt], lambdar[i] - lambdar[j] + mu_pbup[(b,i,j)] <= M*(1-b14r[(b,i,j)]))
        m.ext[:constraints][:Fr14c] = @constraint(m, [(b,i,j) = B_ac_frt], pb[(b, i, j)] + pbmax[b] <= M*b14r[(b,i,j)])
        m.ext[:constraints][:Fr15] = @constraint(m, [(b,i,j) = B_ac_frt], 0 <= - pb[(b, i, j)] + pbmax[b])
        m.ext[:constraints][:Fr15a] = @constraint(m, [(b,i,j) = B_ac_frt], mu_pbup[(b,i,j)] <= M*(1-b15r[(b,i,j)]))
        m.ext[:constraints][:Fr15b] = @constraint(m, [(b,i,j) = B_ac_frt], - pb[(b, i, j)] + pbmax[b] <= M*b15r[(b,i,j)])

        m.ext[:constraints][:p_balancet4] = @constraint(m, -pb[(4,3,4)] + pb[(5,5,6)] + pb[(7,5,8)] - sum(pw[g] for g in G_ac[5]) - sum(pfup[g] for g in G_acupd[5]) + sum(prdo[g] for g in G_acdot[5]) == 0)
        m.ext[:constraints][:p_balancet] = @constraint(m, [i in 1:3], -sum(pw[g] for g in G_ac[i]) - sum(pfup[g] for g in G_acuptc[i]) - sum(gamma*sum((2^k)*xfup[g,k] for k in K) for g in G_acupts[i]) - sum(gamma*sum((2^k)*xrup[g,k] for k in K) for g in G_acupts[i]) - sum(prup[g] for g in G_acuptc[i]) + sum(pfdo[g] for g in G_acdot[i]) + sum(prdo[g] for g in G_acdot[i]) + sum(pd[l] for l in L_ac[i]) + sum(pb[(b,i,j)] for (b,i,j) in B_arcs_frt[i]) - sum(pb[(b,i,j)] for (b,i,j) in B_arcs_tot[i]) == 0) 

        return m 
    end

    # Build model
    build_model!(m)

    # Build model
    optimize!(m)

    # Print solution for inspection
    solution_summary(m, verbose=true)


    x = sum(length.(m.ext[:sets][:G]))
    fup = zeros(x)
    fdo = zeros(x)
    qg = zeros(x)
    rup = zeros(x)
    rdo = zeros(x)
    choedw = zeros(x)
    choedfdo = zeros(x)
    choedfup = zeros(x)
    choedrdo = zeros(x)
    choedrup = zeros(x)

    for g in m.ext[:sets][:GUPC]
        fup[g] = value.(m.ext[:variables][:pfup])[g]
    end
    for g in m.ext[:sets][:GUPS]
        fup[g] = value.(m.ext[:expressions][:pfups])[g]
    end
    for  g in m.ext[:sets][:GDOD]
        fdo[g] = value.(m.ext[:variables][:pfdo])[g]
    end
    for g in m.ext[:sets][:GD]
        qg[g] = value.(m.ext[:variables][:qg])[g]
    end
    for g in m.ext[:sets][:GUPTC]
        rup[g] = value.(m.ext[:variables][:prup])[g]
    end
    for g in m.ext[:sets][:GUPTS]
        rup[g] = value.(m.ext[:expressions][:prups])[g]
    end
    for  g in m.ext[:sets][:GDOT]
        rdo[g] = value.(m.ext[:variables][:prdo])[g]
    end
    for  g in m.ext[:sets][:GS]
        choedw[g] = value.(m.ext[:variables][:choedw])[g]
    end
    for  g in m.ext[:sets][:GUPS]
        choedfup[g] = value.(m.ext[:variables][:choedfup])[g]
    end
    for  g in m.ext[:sets][:GDODS]
        choedfdo[g] = value.(m.ext[:variables][:choedfdo])[g]
    end
    for  g in m.ext[:sets][:GUPTS]
        choedrup[g] = value.(m.ext[:variables][:choedrup])[g]
    end
    for  g in m.ext[:sets][:GDOTS]
        choedrdo[g] = value.(m.ext[:variables][:choedrdo])[g]
    end

    c = value.(m.ext[:parameters][:cg][k] for k in 1:50)
    qmax = value.(m.ext[:parameters][:qmax][k] for k in 1:50)
    pmax = value.(m.ext[:parameters][:pmax][k] for k in 1:50)
    pw = value.(m.ext[:variables][:pw])
    gen_results = DataFrame()
    gen_results.qg = qg
    gen_results.qmax = qmax
    gen_results.cost = c
    gen_results.choedw = choedw
    gen_results.pw = pw
    gen_results.choedfdo = choedfdo
    gen_results.choedfup = choedfup
    gen_results.fdo = fdo
    gen_results.fup = fup
    gen_results.choedrdo = choedrdo
    gen_results.choedrup = choedrup
    gen_results.rdo = rdo
    gen_results.rup = rup
    gen_results.pmax = pmax


    y = sum(length.(m.ext[:sets][:B]))
    qb = zeros(y)
    qbd = value.(m.ext[:variables][:qbd][k] for k in keys(m.ext[:variables][:qbd]))
    for b in 1:4
        qb[b+4] = qbd[b]
    end

    pb = value.(m.ext[:variables][:pb][k] for k in keys(m.ext[:variables][:pb]))
    qbd = value.(m.ext[:variables][:qbd][k] for k in keys(m.ext[:variables][:qbd]))
    pbmax = value.(m.ext[:parameters][:pbmax])
    qbmax = value.(m.ext[:parameters][:qbmax])
    branch_results = DataFrame()
    branch_results.pb = pb
    branch_results.pbmax = pbmax
    branch_results.qb = qb
    branch_results.qbmax = qbmax


    z = sum(length.(m.ext[:sets][:N]))
    lambw = zeros(z)
    lambf = zeros(z)
    lambr = zeros(z)
    wmm = zeros(z)
    lambw[1] = value.(m.ext[:variables][:lambdaw])
    for n in 1:5
        lambf[n+4] = value.(m.ext[:variables][:lambdaf])[n+4]
        wmm[n+4] = value.(m.ext[:variables][:wm])[n+4]
    end
    for n in 1:4
        lambr[n] = value.(m.ext[:variables][:lambdar])[n]
    end

    node_results = DataFrame()
    node_results.wmm = wmm
    node_results.lambdaw = lambw
    node_results.lambdaf = lambf
    node_results.lambdar = lambr

    objw = [value.(m.ext[:objectives][:Obj])]
    objective = DataFrame()
    objective.objw = objw

    XLSX.openxlsx("paperone.xlsx", mode="rw") do xf
        sheet = xf[t]
        XLSX.writetable!(sheet,collect(eachcol(objective)), names(objective); anchor_cell= XLSX.CellRef("U23"))
        XLSX.writetable!(sheet,collect(eachcol(gen_results)), names(gen_results); anchor_cell= XLSX.CellRef("D2"))
        XLSX.writetable!(sheet,collect(eachcol(branch_results)), names(branch_results); anchor_cell= XLSX.CellRef("U2"))
        XLSX.writetable!(sheet,collect(eachcol(node_results)), names(node_results); anchor_cell= XLSX.CellRef("U12"))
    end
end

    solution_summary(m, verbose=true)
