# # Step 0: Activate environment - ensure consistency accross computers
# using Pkg
# Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
# Pkg.instantiate()

# Pkg.add("Printf")
# Pkg.add("DataFrames")
# Pkg.add("JuMP")
# Pkg.add("Gurobi")
# Pkg.add("StatsPlots")
# Pkg.add("Plots")
# Pkg.add("Interact")
# Pkg.add("XLSX")
# # Pkg.add("CSV")

# ## Step 2: create model & pass data to model
# using Printf
# using DataFrames
# using JuMP, Gurobi
# using Plots
# using CSV
# using XLSX
# using Base.Sort

list_Fdso = []
Fdso_inv_costs = []
Fdso_gen_costs = []
Fdso_wsm_costs = []
Fdso_flex_costs = []
Fdso_flex_quant = []
Fdso_bal_costs_0 = []
Fdso_bal_costs_av = []
Fdso_bal_costs_short = []
Fdso_bal_costs_long = []
Fdso_flexenbal_costs_0 = []
Fdso_flexenbal_costs_av = []
Fdso_flexenbal_costs_long = []
Fdso_flexenbal_costs_short = []

Fdso_deltaflex = []
Fdso_flexcost_old =[]
flexcost_old = 0
Fdso_deltabal_0 = []
Fdso_deltabal_av = []
Fdso_deltabal_short = []
Fdso_deltabal_long = []
Fdso_balcost_old_0 =[]
Fdso_balcost_old_av =[]
Fdso_balcost_old_short =[]
Fdso_balcost_old_long =[]
balcost_old_0 = 0
balcost_old_av = 0
balcost_old_long = 0
balcost_old_short = 0
Fdso_deltagen = []
Fdso_gencost_old =[]
gencost_old = 0

for fcap in collect(0:1:71)

    tt_inv_costs = []
    tt_gen_costs = []
    tt_flex_costs = []
    tt_flex_quant = []
    tt_wsm_costs = []
    tt_bal_costs_0 = []
    tt_bal_costs_av = []
    tt_bal_costs_short = []
    tt_bal_costs_long = []

    dso_inv_costs = []
    flex_costs = []
    gen_costs = []
    flex_quant = []
    wsm_costs = []
    bal_costs_0 = []
    bal_costs_av = []
    bal_costs_short = []
    bal_costs_long = []

    dso_flexenbal_costs_0 = []
    dso_flexenbal_costs_av = []
    dso_flexenbal_costs_short = []
    dso_flexenbal_costs_long = []

    deltaflex = []
    deltagen = []
    deltabal_0 = []
    deltabal_av = []
    deltabal_short = []
    deltabal_long = []

    winddata = CSV.read(joinpath(@__DIR__, "Winddata3.csv"), DataFrame)
    costdata = CSV.read(joinpath(@__DIR__, "Costdata.csv"), DataFrame)


    hour = [h for h = 1:85]
    for tt in hour

        optimizer = Gurobi.Optimizer
        m = Model(optimizer_with_attributes(optimizer))

        # Step 2a: create sets
        m.ext[:sets] = Dict()
        jGen = 300 #generation at uncongested area
        J = m.ext[:sets][:J] = [j for j = 1:jGen]

        # step 2c: process input parameters
        flex = "balanced, no reservation"
        m.ext[:parameters] = Dict()
        
        normiCap = m.ext[:parameters][:normiCap] = winddata.VtoP

        # nb = m.ext[:parameters][:nb] = winddata.nb2019
        # Cbalav = m.ext[:parameters][:Cbalav] = costdata.balav19
        # Cbalshort = m.ext[:parameters][:Cbalshort] = costdata.balshort19
        # Cballong = m.ext[:parameters][:Cballong] = costdata.ballong19

        nb = m.ext[:parameters][:nb] = winddata.nb2020
        Cbalav = m.ext[:parameters][:Cbalav] = costdata.balav20
        Cbalshort = m.ext[:parameters][:Cbalshort] = costdata.balshort20
        Cballong = m.ext[:parameters][:Cballong] = costdata.ballong20

        # nb = m.ext[:parameters][:nb] = winddata.nb2021
        # Cbalav = m.ext[:parameters][:Cbalav] = costdata.balav21
        # Cbalshort = m.ext[:parameters][:Cbalshort] = costdata.balshort21
        # Cballong = m.ext[:parameters][:Cballong] = costdata.ballong21

        Cnatgas = m.ext[:parameters][:Cnatgas] = costdata.jcost
        # Cnatgas = m.ext[:parameters][:Cnatgas] = costdata.jcost2
        # Cnatgas = m.ext[:parameters][:Cnatgas] = costdata.jcost3
        Cbal0 = m.ext[:parameters][:Cbal0] = costdata.jcostbal
        Capnatgas = m.ext[:parameters][:Capnatgas] = costdata.jcap
        Capwind = m.ext[:parameters][:Capwind] = normiCap[tt]
        Cwind = m.ext[:parameters][:Cwind] = 1 # [wind; wind]#€/MWh
        Cflex = m.ext[:parameters][:Cflex] = 100 # [wind; wind]#€/MWh
        Fdso = m.ext[:parameters][:Fdso] = fcap # MW
        nwCost = m.ext[:parameters][:nwCost] = 150000 #€/kW = k€/MW 101 #€/kW Tim, 400€/kW Athir
        D = m.ext[:parameters][:D] = 500.25 #MW
        M = m.ext[:parameters][:M] = 200000 #20000000 (more zeros gives convergence errors)


        ## Step 3: construct your model
        function build_GEP_model!(m::Model)
            # Clear m.ext entries "variables", "expressions" and "constraints"
            m.ext[:variables] = Dict()
            m.ext[:expressions] = Dict()
            m.ext[:constraints] = Dict()

            # Extract sets
            J = m.ext[:sets][:J]

            # Extract parameters
            Cwind = m.ext[:parameters][:Cwind]
            normiCap = m.ext[:parameters][:normiCap]
            Cnatgas = m.ext[:parameters][:Cnatgas]
            Capnatgas = m.ext[:parameters][:Capnatgas]
            Capwind = m.ext[:parameters][:Capwind]
            Cflex = m.ext[:parameters][:Cflex]
            Cbal0 = m.ext[:parameters][:Cbal0]
            Cbalav = m.ext[:parameters][:Cbalav]
            Cbalshort = m.ext[:parameters][:Cbalshort]
            Cballong = m.ext[:parameters][:Cballong]
            Fdso = m.ext[:parameters][:Fdso]
            nwCost = m.ext[:parameters][:nwCost]
            D = m.ext[:parameters][:D]
            M = m.ext[:parameters][:M]

            # Create variables
            gammaw = m.ext[:variables][:gammaw] = @variable(m, base_name = "gammaw")
            # gammaf = m.ext[:variables][:gammaf] = @variable(m, base_name = "gammaf")
            gammaf = m.ext[:variables][:gammaf] = @variable(m,  lower_bound = 0, base_name = "gammaf")
            gammab = m.ext[:variables][:gammab] = @variable(m, base_name = "gammab")
            # gammab = m.ext[:variables][:gammab] = @variable(m,  lower_bound = 0, base_name = "gammab")

            qwind = m.ext[:variables][:qwind] = @variable(m, lower_bound = 0, base_name = "qwind")
            muw_wind = m.ext[:variables][:muw_wind] = @variable(m, lower_bound = 0, base_name = "muw_wind")
            qnatgas = m.ext[:variables][:qnatgas] = @variable(m, [j in J], lower_bound = 0, base_name = "qnatgas")
            muw_ng = m.ext[:variables][:muw_ng] = @variable(m, [j in J], lower_bound = 0, base_name = "muw_ng")

            qfwind = m.ext[:variables][:qfwind] = @variable(m, lower_bound = 0, base_name = "qfwind")
            muf_wind = m.ext[:variables][:muf_wind] = @variable(m, lower_bound = 0, base_name = "muf_wind")
           
            qbnatgas = m.ext[:variables][:qbnatgas] = @variable(m, [j in J], lower_bound = 0, base_name = "qbnatgas")
            mub_ng = m.ext[:variables][:mub_ng] = @variable(m, [j in J], lower_bound = 0, base_name = "mub_ng")

            # Create binary variables
            b1 = m.ext[:variables][:b1] = @variable(m, binary = true, base_name = "b1")
            b2 = m.ext[:variables][:b2] = @variable(m, binary = true, base_name = "b2")
            b3 = m.ext[:variables][:b3] = @variable(m, [j in J], Bin, base_name = "b3")
            b4 = m.ext[:variables][:b4] = @variable(m, [j in J], Bin, base_name = "b4")

            b6 = m.ext[:variables][:b6] = @variable(m, binary = true, base_name = "b6")
            b7 = m.ext[:variables][:b7] = @variable(m, binary = true, base_name = "b7")
            b8 = m.ext[:variables][:b8] = @variable(m, binary = true, base_name = "b8")

            b15 = m.ext[:variables][:b15] = @variable(m, [j in J], Bin, base_name = "b15")
            b16 = m.ext[:variables][:b16] = @variable(m, [j in J], Bin, base_name = "b16")
            b17 = m.ext[:variables][:b17] = @variable(m, binary = true, base_name = "b17")


            # Create affine expressions (= linear combinations of variables)
            wsm_costs = m.ext[:expressions][:wsm_costs] = @expression(m, sum(gammaw*qwind) + sum(gammaw*qnatgas[j] for j in J))
            gen_costs = m.ext[:expressions][:gen_costs] = @expression(m, sum(Cwind*qwind) + sum(Cnatgas[j]*qnatgas[j] for j in J) - sum(Cwind*qfwind) + sum(Cnatgas[j]*qbnatgas[j] for j in J))
            flex_costs = m.ext[:expressions][:flex_costs] = @expression(m, sum(Cflex*qfwind))
            bal_costs_0 = m.ext[:expressions][:bal_costs_0] = @expression(m, sum(Cbal0[j]*qbnatgas[j] for j in J))
            bal_costs_av = m.ext[:expressions][:bal_costs_av] = @expression(m, sum(Cbalav[j]*qbnatgas[j] for j in J))
            bal_costs_short = m.ext[:expressions][:bal_costs_short] = @expression(m, sum(Cbalshort[j]*qbnatgas[j] for j in J))
            bal_costs_long = m.ext[:expressions][:bal_costs_long] = @expression(m, sum(Cballong[j]*qbnatgas[j] for j in J))

            inv_costs = m.ext[:expressions][:inv_costs] = @expression(m, nwCost*Fdso)
            qnatgas_tot = m.ext[:expressions][:qnatgas_tot] = @expression(m, sum(qbnatgas[j] for j in J))
            # dso_costs = m.ext[:expressions][:dso_costs] = @expression(m,flex_costs + inv_costs)
            # tot_costs = m.ext[:expressions][:tot_costs] = @expression(m,wsm_costs + bal_costs + dso_costs)

            # wholesale market
            m.ext[:constraints][:F1] = @constraint(m, 0 <= Cwind + muw_wind - gammaw)
            m.ext[:constraints][:F1a] = @constraint(m, qwind <= M*(1-b1))
            m.ext[:constraints][:F1b] = @constraint(m, Cwind + muw_wind - gammaw <= M*b1)
            m.ext[:constraints][:F2] = @constraint(m, 0 <=  Capwind - qwind)
            m.ext[:constraints][:F2a] = @constraint(m, muw_wind <= M*(1-b2))
            m.ext[:constraints][:F2b] = @constraint(m, Capwind - qwind <= M*b2)

            m.ext[:constraints][:F3] = @constraint(m, [j=J], 0 <= Cnatgas[j] + muw_ng[j] - gammaw)
            m.ext[:constraints][:F3a] = @constraint(m, [j=J], qnatgas[j] <= M*(1-b3[j]))
            m.ext[:constraints][:F3b] = @constraint(m, [j=J], Cnatgas[j] + muw_ng[j] - gammaw <= M*b3[j])
            m.ext[:constraints][:F4] = @constraint(m, [j=J], 0 <= Capnatgas[j] - qnatgas[j])
            m.ext[:constraints][:F4a] = @constraint(m, [j=J], muw_ng[j] <= M*(1-b4[j]))
            m.ext[:constraints][:F4b] = @constraint(m, [j=J], Capnatgas[j] - qnatgas[j] <= M*b4[j])

            m.ext[:constraints][:F5] = @constraint(m, - qwind - sum(qnatgas[j] for j in J) + D == 0)

            # flexibility market
            m.ext[:constraints][:F6] = @constraint(m, 0 <= Cflex + muf_wind - gammaf)
            m.ext[:constraints][:F6a] = @constraint(m, qfwind <= M*(1-b6))
            m.ext[:constraints][:F6b] = @constraint(m, Cflex + muf_wind - gammaf <= M*b6)
            m.ext[:constraints][:F7] = @constraint(m, 0 <= qwind - qfwind)
            m.ext[:constraints][:F7a] = @constraint(m, muf_wind <= M*(1-b7))
            m.ext[:constraints][:F7b] = @constraint(m, qwind - qfwind <= M*b7)

            # m.ext[:constraints][:F8] = @constraint(m, qfwind - qwind + Fdso == 0)
               
            m.ext[:constraints][:F8] = @constraint(m, 0 <= qfwind - qwind + Fdso)
            m.ext[:constraints][:F8a] = @constraint(m, gammaf <= M*(1-b8))
            m.ext[:constraints][:F8b] = @constraint(m, qfwind - qwind + Fdso <= M*b8)

            #balancing market

            m.ext[:constraints][:F15] = @constraint(m, [j=J], 0 <= Cnatgas[j] + mub_ng[j] + gammab)
            m.ext[:constraints][:F15a] = @constraint(m, [j=J], qbnatgas[j] <= M*(1-b15[j]))
            m.ext[:constraints][:F15b] = @constraint(m, [j=J], Cnatgas[j] + mub_ng[j] + gammab <= M*b15[j])
            m.ext[:constraints][:F16] = @constraint(m, [j=J], 0 <= Capnatgas[j] - qnatgas[j] - qbnatgas[j])
            m.ext[:constraints][:F16a] = @constraint(m, [j=J], mub_ng[j] <= M*(1-b16[j]))
            m.ext[:constraints][:F16b] = @constraint(m, [j=J],  Capnatgas[j] - qnatgas[j] - qbnatgas[j] <= M*b16[j])
    
            m.ext[:constraints][:F17] = @constraint(m, - sum(qbnatgas[j] for j in J) + qfwind == 0)
            # m.ext[:constraints][:F17] = @constraint(m, 0 <= sum(qbnatgas[j] for j in J) + qwind - Fdso)
            # m.ext[:constraints][:F17a] = @constraint(m, gammab <= M*(1-b17))
            # m.ext[:constraints][:F17b] = @constraint(m, sum(qbnatgas[j] for j in J) + qwind - Fdso <= M*b17)

            return m
        end

        # Build your model
        build_GEP_model!(m)
        ## Step 4: solve

        status = optimize!(m)

        # Step 5: interpretation
        south = DataFrame()
        south.cost = value.(m.ext[:parameters][:Cnatgas])
        south.qnatgas = value.(m.ext[:variables][:qnatgas])
        south.qbnatgas= value.(m.ext[:variables][:qbnatgas])

        push!(tt_gen_costs, value.(m.ext[:expressions][:gen_costs]))
        push!(tt_flex_costs, value.(m.ext[:expressions][:flex_costs]))
        push!(tt_flex_quant, value.(m.ext[:variables][:qfwind]))
        push!(tt_inv_costs, value.(m.ext[:expressions][:inv_costs]))
        push!(tt_wsm_costs, value.(m.ext[:expressions][:wsm_costs]))
        push!(tt_bal_costs_0, value.(m.ext[:expressions][:bal_costs_0]))
        push!(tt_bal_costs_av, value.(m.ext[:expressions][:bal_costs_av]))
        push!(tt_bal_costs_short, value.(m.ext[:expressions][:bal_costs_short]))
        push!(tt_bal_costs_long, value.(m.ext[:expressions][:bal_costs_long]))


    end # of tt day loop

    gen_costs = sum(tt_gen_costs.*nb)
    flex_costs = sum(tt_flex_costs.*nb)
    flex_quant = sum(tt_flex_quant.*nb)
    wsm_costs = sum(tt_wsm_costs.*nb)
    bal_costs_0 = sum(tt_bal_costs_0.*nb)
    bal_costs_av = sum(tt_bal_costs_av.*nb)
    bal_costs_short = sum(tt_bal_costs_short.*nb)
    bal_costs_long = sum(tt_bal_costs_long.*nb)
    dso_inv_costs = first(tt_inv_costs)
    dso_flexenbal_costs_0 = (flex_costs .+ bal_costs_0)*1
    dso_flexenbal_costs_av = (flex_costs .+ bal_costs_av)*1
    dso_flexenbal_costs_short = (flex_costs .+ bal_costs_short)*1
    dso_flexenbal_costs_long = (flex_costs .+ bal_costs_long)*1

    deltaflex = flexcost_old - flex_costs
    push!(Fdso_deltaflex, deltaflex)
    flexcost_old = flex_costs
    push!(Fdso_flexcost_old, flexcost_old)

    deltabal_0 = balcost_old_0 - bal_costs_0
    push!(Fdso_deltabal_0, deltabal_0)
    balcost_old_0 = bal_costs_0
    push!(Fdso_balcost_old_0, balcost_old_0)

    deltabal_av = balcost_old_av - bal_costs_av
    push!(Fdso_deltabal_av, deltabal_av)
    balcost_old_av = bal_costs_av
    push!(Fdso_balcost_old_av, balcost_old_av)

    deltabal_short = balcost_old_short - bal_costs_short
    push!(Fdso_deltabal_short, deltabal_short)
    balcost_old_short = bal_costs_short
    push!(Fdso_balcost_old_short, balcost_old_short)

    deltabal_long = balcost_old_long - bal_costs_long
    push!(Fdso_deltabal_long, deltabal_long)
    balcost_old_long = bal_costs_long
    push!(Fdso_balcost_old_long, balcost_old_long)

    deltagen = gencost_old - gen_costs
    push!(Fdso_deltagen, deltagen)
    gencost_old = gen_costs
    push!(Fdso_gencost_old, gencost_old)

    push!(Fdso_inv_costs, dso_inv_costs)

    # push!(Fdso_inv_costs, dso_inv_costs)
    push!(Fdso_wsm_costs, wsm_costs)
    push!(Fdso_gen_costs, gen_costs)
    push!(Fdso_flex_costs, flex_costs)
    push!(Fdso_flex_quant, flex_quant)
    push!(Fdso_bal_costs_0, bal_costs_0)
    push!(Fdso_bal_costs_av, bal_costs_av)
    push!(Fdso_bal_costs_short, bal_costs_short)
    push!(Fdso_bal_costs_long, bal_costs_long)
    push!(Fdso_flexenbal_costs_0, dso_flexenbal_costs_0)
    push!(Fdso_flexenbal_costs_av, dso_flexenbal_costs_av)
    push!(Fdso_flexenbal_costs_short, dso_flexenbal_costs_short)
    push!(Fdso_flexenbal_costs_long, dso_flexenbal_costs_long)

end # of fcap loop

yearly_gen_bal = (Fdso_gen_costs)/1000000
yearly_wsm_bal = Fdso_wsm_costs/1000000
yearly_inv_bal = Fdso_inv_costs/1000000

#here 0.8
yearly_quant_bal = Fdso_flex_quant*1.0
yearly_flex_bal = Fdso_flex_costs/1000000*1.0
yearly_flexenbal_bal_0 = Fdso_flexenbal_costs_0/1000000*1.0
yearly_flexenbal_bal_av = Fdso_flexenbal_costs_av/1000000*1.0
yearly_flexenbal_bal_short = Fdso_flexenbal_costs_short/1000000*1.0
yearly_flexenbal_bal_long = Fdso_flexenbal_costs_long/1000000*1.0
yearly_bal_0 = Fdso_bal_costs_0/1000000*1.0
yearly_bal_av = Fdso_bal_costs_av/1000000*1.0
yearly_bal_short = Fdso_bal_costs_short/1000000*1.0
yearly_bal_long = Fdso_bal_costs_long/1000000*1.0

yearly_syst_bal = yearly_gen_bal + yearly_inv_bal
yearly_dso_bal_0 = yearly_inv_bal + yearly_flexenbal_bal_0 
yearly_dso_bal_no = yearly_inv_bal + yearly_flex_bal
yearly_dso_bal_av = yearly_inv_bal + yearly_flexenbal_bal_av 
yearly_dso_bal_short = yearly_inv_bal + yearly_flexenbal_bal_short
yearly_dso_bal_long = yearly_inv_bal + yearly_flexenbal_bal_long

yearly_deltaflex_bal = Fdso_deltaflex/1000000*1.0
popfirst!(yearly_deltaflex_bal)
yearly_deltabal_bal_0 = Fdso_deltabal_0/1000000*1.0
popfirst!(yearly_deltabal_bal_0)
yearly_deltabal_bal_av = Fdso_deltabal_av/1000000*1.0
popfirst!(yearly_deltabal_bal_av)
yearly_deltabal_bal_short = Fdso_deltabal_short/1000000*1.0
popfirst!(yearly_deltabal_bal_short)
yearly_deltabal_bal_long = Fdso_deltabal_long/1000000*1.0
popfirst!(yearly_deltabal_bal_long)

yearly_deltagen = Fdso_deltagen/1000000
popfirst!(yearly_deltagen)


yearly_deltadso_bal_0 = yearly_deltaflex_bal + yearly_deltabal_bal_0
yearly_deltadso_bal_no = yearly_deltaflex_bal
yearly_deltadso_bal_av = yearly_deltaflex_bal + yearly_deltabal_bal_av
yearly_deltadso_bal_short = yearly_deltaflex_bal + yearly_deltabal_bal_short
yearly_deltadso_bal_long = yearly_deltaflex_bal + yearly_deltabal_bal_long

# FIND OPTIMUM RESULTS

# Fopt_gen_bal = findmin(yearly_gen_bal)
# Fmin_gen_bal = sort(yearly_gen_bal; alg=Sort.PartialQuickSort(10))[1:10]
# Fopt_dso_bal_0 = findmin(yearly_dso_bal_0)
# Fmin_dso_bal_0 = sort(yearly_dso_bal_0; alg=Sort.PartialQuickSort(10))[1:10]

Fopt_syst_bal = findmin(yearly_syst_bal)
Fmin_syst_bal = sort(yearly_syst_bal; alg=Sort.PartialQuickSort(10))[1:10]
Fopt_dso_bal_av = findmin(yearly_dso_bal_av)
Fmin_dso_bal_av = sort(yearly_dso_bal_av; alg=Sort.PartialQuickSort(10))[1:10]
Fopt_dso_bal_short = findmin(yearly_dso_bal_short)
Fmin_dso_bal_short = sort(yearly_dso_bal_short; alg=Sort.PartialQuickSort(10))[1:10]
Fopt_dso_bal_long = findmin(yearly_dso_bal_long)
Fmin_dso_bal_long = sort(yearly_dso_bal_long; alg=Sort.PartialQuickSort(10))[1:10]
Fopt_dso_bal_no = findmin(yearly_dso_bal_no)
Fmin_dso_bal_no = sort(yearly_dso_bal_no; alg=Sort.PartialQuickSort(10))[1:10]

z_syst = Fopt_syst_bal[2]
z_dso_av = Fopt_dso_bal_av[2]
z_dso_short = Fopt_dso_bal_short[2]
z_dso_long = Fopt_dso_bal_long[2]
z_dso_no = Fopt_dso_bal_no[2]

y_syst = Fopt_syst_bal[1]
y_dso_av = Fopt_dso_bal_av[1]
y_dso_short = Fopt_dso_bal_short[1]
y_dso_long = Fopt_dso_bal_long[1]
y_dso_no = Fopt_dso_bal_no[1]

zMW_syst = z_syst - 1
zMW_dso_av = z_dso_av - 1
zMW_dso_short = z_dso_short -1
zMW_dso_long = z_dso_long -1
zMW_dso_no = z_dso_no - 1

println("""
Fopt_syst_bal = $(Fopt_syst_bal)
Fopt_dso_bal_av = $(Fopt_dso_bal_av)
Fopt_dso_bal_short = $(Fopt_dso_bal_short)
Fopt_dso_bal_long = $(Fopt_dso_bal_long)
Fopt_dso_bal_no = $(Fopt_dso_bal_no)

Fmin_syst_bal = $(Fmin_syst_bal)
Fmin_dso_bal_av = $(Fmin_dso_bal_av)
Fmin_dso_bal_short = $(Fmin_dso_bal_short)
Fmin_dso_bal_long = $(Fmin_dso_bal_long)
Fmin_dso_bal_no = $(Fmin_dso_bal_no)
""")

# CALCULATE PERCENTAGE OF CURTAILMENT
# totwindprod = 213303.995 # 2017
# totwindprod = 213319.9918 # 2019
totwindprod = 218000.6758 # 2020
# totwindprod = 180130.848 # 2021
totwindcurt_bal = yearly_quant_bal # MWh one day per season, for different values of Fdso
percwindcurt_bal = (totwindcurt_bal/totwindprod)*100

optpercwindcurt_syst = getindex(percwindcurt_bal, z_syst)
optpercwindcurt_dso_bal_av = getindex(percwindcurt_bal, z_dso_av)
optpercwindcurt_dso_bal_short = getindex(percwindcurt_bal, z_dso_short)
optpercwindcurt_dso_bal_long = getindex(percwindcurt_bal, z_dso_long)
optpercwindcurt_dso_bal_no = getindex(percwindcurt_bal, z_dso_no)

opttotwindcurt_syst = getindex(totwindcurt_bal, z_syst)
opttotwindcurt_dso_bal_av  = getindex(totwindcurt_bal, z_dso_av)
opttotwindcurt_dso_bal_short  = getindex(totwindcurt_bal, z_dso_short)
opttotwindcurt_dso_bal_long = getindex(totwindcurt_bal, z_dso_long)
opttotwindcurt_dso_bal_no = getindex(totwindcurt_bal, z_dso_no)

println("""
optpercwindcurt_syst = $(optpercwindcurt_syst)
optpercwindcurt_dso_bal_av = $(optpercwindcurt_dso_bal_av)
optpercwindcurt_dso_bal_short = $(optpercwindcurt_dso_bal_short)
optpercwindcurt_dso_bal_long = $(optpercwindcurt_dso_bal_long)
optpercwindcurt_dso_bal_no = $(optpercwindcurt_dso_bal_no)

opttotwindcurt_syst = $(opttotwindcurt_syst)
opttotwindcurt_dso_bal_av = $(opttotwindcurt_dso_bal_av)
opttotwindcurt_dso_bal_short = $(opttotwindcurt_dso_bal_short)
opttotwindcurt_dso_bal_long = $(opttotwindcurt_dso_bal_long)
opttotwindcurt_dso_bal_no = $(opttotwindcurt_dso_bal_no)
""")

# print results

optimum = DataFrame()
optimum.system = [zMW_syst, z_syst, Fopt_syst_bal[1], optpercwindcurt_syst, opttotwindcurt_syst]
optimum.bal_av = [zMW_dso_av, z_dso_av, Fopt_dso_bal_av[1], optpercwindcurt_dso_bal_av, opttotwindcurt_dso_bal_av]
optimum.bal_short = [zMW_dso_short, z_dso_short, Fopt_dso_bal_short[1], optpercwindcurt_dso_bal_short, opttotwindcurt_dso_bal_short]
optimum.bal_long = [zMW_dso_long, z_dso_long, Fopt_dso_bal_long[1], optpercwindcurt_dso_bal_long, opttotwindcurt_dso_bal_long]
optimum.bal_no = [zMW_dso_no, z_dso_no, Fopt_dso_bal_no[1], optpercwindcurt_dso_bal_no, opttotwindcurt_dso_bal_no]

costs = DataFrame()
costs.inv = yearly_inv_bal 
costs.flex = yearly_flex_bal 
costs.bal_av = yearly_bal_av 
costs.bal_short = yearly_bal_short
costs.bal_long = yearly_bal_long
costs.dso_no = yearly_dso_bal_no
costs.dso_av = yearly_dso_bal_av
costs.dso_short = yearly_dso_bal_short
costs.dso_long = yearly_dso_bal_long

delta = DataFrame()
delta.gen = yearly_deltagen
delta.dso_no = yearly_deltadso_bal_no 
delta.dso_av = yearly_deltadso_bal_av 
delta.dso_short = yearly_deltadso_bal_short 
delta.dso_long = yearly_deltadso_bal_long 

param = DataFrame()
param.wind_prod = [totwindprod]
param.cost_flex = [value.(m.ext[:parameters][:Cflex])]
param.cost_inv = [value.(m.ext[:parameters][:nwCost])]

balance = DataFrame()
balance.cost_av = [m.ext[:parameters][:Cbalav][1]]
balance.cost_short = [m.ext[:parameters][:Cbalshort][1]]
balance.cost_long = [m.ext[:parameters][:Cballong][1]]

XLSX.openxlsx("papertwo.xlsx", mode="rw") do xf
    sheet = xf[1]
    XLSX.writetable!(sheet,collect(eachcol(param)), names(param); anchor_cell= XLSX.CellRef("B2"))
    XLSX.writetable!(sheet,collect(eachcol(balance)), names(balance); anchor_cell= XLSX.CellRef("B5"))
    XLSX.writetable!(sheet,collect(eachcol(optimum)), names(optimum); anchor_cell= XLSX.CellRef("I2"))
    XLSX.writetable!(sheet,collect(eachcol(costs)), names(costs); anchor_cell= XLSX.CellRef("C9"))
    XLSX.writetable!(sheet,collect(eachcol(delta)), names(delta); anchor_cell= XLSX.CellRef("M9"))
end

# data = [yearly_dso_bal_av yearly_inv_bal yearly_flex_bal yearly_bal_av]
# labels = ["total" "investments" "flexibility" "average balancing"]
# # fig_1 = plot(data, color=[:darkgreen :orange :yellow], yrange =(0,11.5), hover = data, seriestype = :line, minorticks = true, label = labels)
# fig_1 = plot(data, color=[:darkgreen :orange :yellow :pink], yrange =(-2,50), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO costs [M€]", yguidefontsize=8,  hover = data, seriestype = :line, minorticks = true, label = labels)
# vline!([z_dso_av], color=:black, linestyle=:dash, label="") 
# savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsoav_components.png")
# @show fig_1

# data = [yearly_dso_bal_short yearly_inv_bal yearly_flex_bal yearly_bal_short]
# labels = ["total" "investments" "flexibility" "short balancing"]
# # fig_1 = plot(data, color=[:darkgreen :orange :yellow], yrange =(0,11.5), hover = data, seriestype = :line, minorticks = true, label = labels)
# fig_1 = plot(data, color=[:darkgreen :orange :yellow :pink], yrange =(-2,50), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO costs [M€]", yguidefontsize=8,  hover = data, seriestype = :line, minorticks = true, label = labels)
# vline!([z_dso_short], color=:black, linestyle=:dash, label="") 
# savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsoshort_components.png")
# @show fig_1

# data = [yearly_dso_bal_long yearly_inv_bal yearly_flex_bal yearly_bal_long]
# labels = ["total" "investments" "flexibility" "long balancing"]
# # fig_1 = plot(data, color=[:darkgreen :orange :yellow], yrange =(0,11.5), hover = data, seriestype = :line, minorticks = true, label = labels)
# fig_1 = plot(data, color=[:darkgreen :orange :yellow :pink], yrange =(-2,50), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO costs [M€]", yguidefontsize=8,  hover = data, seriestype = :line, minorticks = true, label = labels)
# vline!([z_dso_long], color=:black, linestyle=:dash, label="") 
# savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsolong_components.png")
# @show fig_1

data = [yearly_dso_bal_no yearly_inv_bal yearly_flex_bal]
labels = ["total costs" "investment costs" "flexibility costs" "DSO bal costs"]
fig_1 = plot(data, color=[:gold :black :green], yrange =(0,30), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO costs [M€]", yguidefontsize=8,  hover = data, seriestype = :line, minorticks = true, label = labels)
vline!([z_dso_no], color=:darkgrey, linestyle=:dash, label="") 
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsono_component.png")
@show fig_1


data = [yearly_dso_bal_no yearly_dso_bal_av yearly_dso_bal_short yearly_dso_bal_long]
labels = ["neglect balancing costs" "average balancing costs" "short balancing costs" "long balancing costs"]
fig_1 = plot(data, color=[:green :red :darkorange :purple :blue], yrange =(0,40), hover = data, xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO total costs [M€]", yguidefontsize=8, minorticks = true, label = labels)
vline!([z_dso_no], color=:darkgrey, linestyle=:dash, label="") 
vline!([z_dso_av], color=:darkgrey, linestyle=:dash, label="") 
vline!([z_dso_short], color=:darkgrey, linestyle=:dash, label="") 
vline!([z_dso_long], color=:darkgrey, linestyle=:dash, label="") 
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsotot_comparison.png")
@show fig_1

data = [yearly_deltadso_bal_no yearly_deltagen]
labels = ["neglect balancing costs" "system perspective"]
fig_1 = plot(data,color=[:green :blue],hover = data,  yrange =(0,1), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "Marginal costs [M€]", yguidefontsize=8, seriestype = :line, minorticks = true, label = labels)
vline!([z_dso_no], color=:darkgrey, linestyle=:dash, label="") 
vline!([z_syst-1], color=:darkgrey, linestyle=:dash, label="") 
hline!([0.25], color=:black, linestyle=:solid, label="investment costs") 
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsobal_deltacomparisonnone.png")
@show fig_1



plot_size = (1600, 1200)  # Set the desired figure size in pixels
plot!(size = plot_size)

data = [yearly_dso_bal_short yearly_dso_bal_av yearly_dso_bal_no yearly_dso_bal_long]
labels = ["worst case"  "average case" "neglect case"  "best case"]
fig_1 = plot(data, color=[:darkorange :red :green :purple :blue], yrange =(0,40), hover = data, xlabel = "Installed network capacity [MW]", xguidefontsize=8, ylabel = "Total DSO costs [M€]", yguidefontsize=8, minorticks = true, label = labels)
scatter!([z_dso_no], [y_dso_no], marker = :x, markersize = 2.8, color = :black, label = "optimum")
scatter!([z_dso_av], [y_dso_av], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_short], [y_dso_short], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_long], [y_dso_long], marker = :x, markersize = 2.8, color = :black, label = "")
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\230623totalcosts.png")
@show fig_1



data = [yearly_deltadso_bal_short yearly_deltadso_bal_av yearly_deltadso_bal_no  yearly_deltadso_bal_long yearly_deltagen]
labels = ["operational costs - worst case" "operational costs - average case" "operational costs - neglect case"  "operational costs - best case" "generation costs"]
# fig_1 = plot(data, hover = data, xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "Marginal total DSO costs [M€]", yguidefontsize=8, seriestype = :line, minorticks = true, label = labels)
fig_1 = plot(data,  color=[:darkorange :red :green :purple :blue],  yrange =(0,1.5),hover = data, xlabel = "Installed network capacity [MW]", xguidefontsize=8, ylabel = "Marginal costs [M€]", yguidefontsize=8, seriestype = :line, minorticks = true, label = labels)
hline!([0.25], color=:black, label="investment costs") 
scatter!([z_dso_no], [0.25], marker = :x, markersize = 2.8, color = :black, label = "optimum")
scatter!([z_dso_av], [0.25], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_short], [0.25], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_long], [0.25], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_syst-1], [0.25], marker = :x, markersize = 2.8, color = :black, label = "")
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\230623marginalcosts.png")
@show fig_1

# # final figures


plot_size = (1600, 1200)  # Set the desired figure size in pixels
plot!(size = plot_size)



data = [yearly_deltadso_bal_short yearly_deltadso_bal_av yearly_deltadso_bal_no  yearly_deltadso_bal_long yearly_deltagen]
labels = ["operational costs - worst case" "operational costs - average case" "operational costs - neglect case"  "operational costs - best case" "generation costs"]
# fig_1 = plot(data, hover = data, xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "Marginal total DSO costs [M€]", yguidefontsize=8, seriestype = :line, minorticks = true, label = labels)
fig_1 = plot(data,  color=[:darkorange :red :green :purple :blue], legend= false, xrange=(0,72), yrange =(0,1.5),hover = data, xlabel = "Installed network capacity [MW]", xguidefontsize=8, ylabel = "Marginal costs [M€]", yguidefontsize=8, seriestype = :line, minorticks = true, label = labels)
hline!([0.15], color=:black, label="investment costs") 
scatter!([z_dso_no], [0.15], marker = :x, markersize = 2.8, color = :black, label = "optimum")
scatter!([z_dso_av], [0.15], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_short], [0.15], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_long], [0.15], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_syst], [0.15], marker = :x, markersize = 2.8, color = :black, label = "")
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\marginalcosts.png")
@show fig_1

data = [yearly_dso_bal_short yearly_dso_bal_av yearly_dso_bal_no yearly_dso_bal_long]
labels = ["worst case"  "average case" "neglect case"  "best case"]
fig_1 = plot(data, color=[:darkorange :red :green :purple :blue], legend= false, xrange=(0,72), yrange =(0,40), hover = data, xlabel = "Installed network capacity [MW]", xguidefontsize=8, ylabel = "Total costs [M€]", yguidefontsize=8, minorticks = true, label = labels)
scatter!([z_dso_no], [y_dso_no], marker = :x, markersize = 2.8, color = :black, label = "optimum")
scatter!([z_dso_av], [y_dso_av], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_short], [y_dso_short], marker = :x, markersize = 2.8, color = :black, label = "")
scatter!([z_dso_long], [y_dso_long], marker = :x, markersize = 2.8, color = :black, label = "")
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\totalcosts.png")
@show fig_1

data = [yearly_dso_bal_av yearly_inv_bal yearly_flex_bal yearly_bal_av]
labels = ["total" "investments" "flexibility" "average balancing"]
# fig_1 = plot(data, color=[:darkgreen :orange :yellow], yrange =(0,11.5), hover = data, seriestype = :line, minorticks = true, label = labels)
fig_1 = plot(data, color=[:red :black :grey :grey], xrange=(0,72), yrange =(0,40), xlabel = "Installed network capacity Fdso [MW]", xguidefontsize=8, ylabel = "DSO costs [M€]", yguidefontsize=8,  hover = data, line = [:solid :solid :solid :dash], minorticks = true, label = labels)
scatter!([z_dso_av], [y_dso_av], marker = :x, markersize = 2.8, color = :black, label = "")
savefig("d:\\users\\ellen.beckstedde\\Desktop\\figures\\dsoav_components.png")
@show fig_1

