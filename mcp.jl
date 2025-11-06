using JuMP
using PATHSolver
using MPSGE
using DataFrames


include("models.jl")
include("mpsge.jl")

begin
data = ModelData()
E0 = data.initial_exports
M0 = data.initial_imports
X0 = data.initial_production
Q0 = data.initial_consumption
D0 = data.domestic_demand
P0 = data.initial_consumer_price_index
PDD0 = data.initial_domestic_price
PWE0 = data.initial_world_export_price
PWM0 = data.initial_world_import_price


MCP = Model(PATHSolver.Optimizer)

@variables(MCP, begin
    sigma in JuMP.Parameter(.2)
    omega in JuMP.Parameter(.2)
    BBAR in JuMP.Parameter(0.0)
    PWE in JuMP.Parameter(1)
    PWM in JuMP.Parameter(1)
    TM in JuMP.Parameter(0.0)
    TE in JuMP.Parameter(0.0)
    TD in JuMP.Parameter(0.0)

    X>=0, (start = 1)
    Q>=0, (start = 1)

    PFX>=0, (start = 1)
    PX>=0, (start = 1)
    PQ>=0, (start = 1)
    PDT>=0, (start = 1)

    Y>=0, (start = 100)

end)

# Zero Profit
export_price = PFX*PWE*(1 + TE)
revenue_X = (D0/(D0+E0)*(PDT)^(1+omega) + (E0)/(D0+E0)*export_price^(1+omega))^(1/(1+omega))


cost_X = PX

cd_X_PFX = E0*PWE * (revenue_X/export_price)^(-omega)
cd_X_PDT = D0 * (revenue_X/PDT)^(-omega)
cd_X_PX = X0
@constraint(MCP, zp_X,
    cost_X - revenue_X ⟂ X
)

import_price = PFX*PWM*(1 + TM)
domestic_price = PDT*(1+TD)
revenue_Q = PQ
cost_Q = @expression(MCP,
    ifelse(sigma == 1,
        import_price^(M0/(M0+D0)) * domestic_price^(D0/(M0+D0)),
        ((M0)/(M0+D0)*import_price^(1-sigma) + D0/(M0+D0)*domestic_price^(1-sigma))^(1/(1-sigma))
    )
)

cd_Q_PFX = M0*PWM * (cost_Q/import_price)^(sigma)
cd_Q_PDT = D0 * (cost_Q/domestic_price)^(sigma)
cd_Q_PQ = Q0
@constraint(MCP, zp_Q,
    cost_Q - revenue_Q ⟂ Q
)


# Market Clearance

# PX
@constraint(MCP, mc_PX,
    X0 - X*cd_X_PX ⟂ PX
)

# PQ
@constraint(MCP, mc_PQ,
    Q*cd_Q_PQ - Y/PQ ⟂ PQ
)

# PDT
@constraint(MCP, mc_PDT,
    X*cd_X_PDT - Q*cd_Q_PDT ⟂ PDT
)


# PFX
@constraint(MCP, mc_PFX,
    X*cd_X_PFX - Q*cd_Q_PFX + BBAR ⟂ PFX
)

# Income Balance
tax_revenue =  -X*cd_X_PFX*TE*PFX + Q*cd_Q_PFX*TM*PFX + Q*cd_Q_PDT*TD*PDT
@constraint(MCP, income_balance,
    Y - PX*X0 - BBAR*PFX - tax_revenue ⟂ Y
)

fix(PQ, 1; force=true)

set_attribute(MCP, "cumulative_iteration_limit", 0)
optimize!(MCP)

set_attribute(MCP, "cumulative_iteration_limit", 10_000)

end





begin
    params = NLPModelParameters(sigma=.2, omega=15, BBAR=10, PWE=1.1, PWM = 1.0, TD = .2, TE = .1, TM = .1)

    MP = MPSGE_model(data)
    set_parameter_values(MP, params)
    solve!(MP)

end


generate_report(MP)



begin

    set_parameter_value(sigma, params.sigma)
    set_parameter_value(omega, params.omega)
    set_parameter_value(BBAR, params.BBAR)
    set_parameter_value(PWE, params.PWE)
    set_parameter_value(PWM, params.PWM)
    set_parameter_value(TD, params.TD)
    set_parameter_value(TE, params.TE)
    set_parameter_value(TM, params.TM)
    
    optimize!(MCP)
end

zip(all_variables(MCP), value.(all_variables(MCP))) |> DataFrame
generate_report(MP)