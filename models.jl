struct ModelData
    initial_exports::Float64
    initial_imports::Float64
    initial_production::Float64
    initial_consumption::Float64
    domestic_demand::Float64
    initial_consumer_price_index::Float64
    initial_world_import_price::Float64
    initial_world_export_price::Float64
    initial_domestic_price::Float64
    ModelData(; 
        initial_exports=25.0,
        initial_imports=25.0,
        initial_production=100.0,
        initial_consumption=100.0,
        domestic_demand=75.0,
        initial_consumer_price_index=1.0,
        initial_world_import_price=1.0,
        initial_world_export_price=1.0,
        initial_domestic_price=1.0
    ) = new(
        initial_exports,
        initial_imports,
        initial_production,
        initial_consumption,
        domestic_demand,
        initial_consumer_price_index,
        initial_world_import_price,
        initial_world_export_price,
        initial_domestic_price
    )
end


struct NLPModelParameters
    sigma::Float64
    omega::Float64
    BBAR::Float64
    PWE::Float64
    PWM::Float64
    TM::Float64
    TE::Float64
    TD::Float64
    function NLPModelParameters(; 
        sigma=0.2,
        omega=0.2,
        BBAR=0.0,
        PWE=1.0,
        PWM=1.0,
        TM=0.0,
        TE=0.0,
        TD=0.0
    )
        return new(
            sigma,
            omega,
            BBAR,
            PWE,
            PWM,
            TM,
            TE,
            TD
        )
    end
end

function set_parameter_values(M::JuMP.Model, params::NLPModelParameters)
    set_parameter_value(M[:BBAR], params.BBAR)
    set_parameter_value(M[:sigma], params.sigma)
    set_parameter_value(M[:omega], params.omega)
    set_parameter_value(M[:PWE], params.PWE)
    set_parameter_value(M[:PWM], params.PWM)
    set_parameter_value(M[:TM], params.TM)
    set_parameter_value(M[:TE], params.TE)
    set_parameter_value(M[:TD], params.TD)
end






function NLP_model(data::ModelData)

    E0 = data.initial_exports
    M0 = data.initial_imports
    X0 = data.initial_production
    Q0 = data.initial_consumption
    D0 = data.domestic_demand
    P0 = data.initial_consumer_price_index
    PDD0 = data.initial_domestic_price
    PWE0 = data.initial_world_export_price
    PWM0 = data.initial_world_import_price
    

    GSS = Model(Ipopt.Optimizer)

    @variables(GSS, begin
        sigma in JuMP.Parameter(.2) #Elasticity of substitution
        omega in JuMP.Parameter(.2) #Elasticity of Transformation
        BBAR in JuMP.Parameter(0)   #Exogenous value of commercial balance

        PWE in JuMP.Parameter(PWE0) #World Export Price
        PWM in JuMP.Parameter(PWM0) #World Import Price
        TM in JuMP.Parameter(0) # Customs duty rate
        TE in JuMP.Parameter(0) # Export subsidy
        TD in JuMP.Parameter(0) # Domestic tax
    end)

    @expressions(GSS, begin
        rho, (1/sigma) - 1
        h, (1/omega) + 1
        alpha, 1/ (PDD0/PWE0 * (E0/D0)^(1/omega) +1)
        beta, ((PWM0/PDD0)*(M0/D0)^(1/sigma))/((PWM0/PDD0)*(M0/D0)^(1/sigma)+1)
        A, X0*(alpha*E0^h + (1-alpha)*D0^h)^(-1/h)
        B, Q0*(beta*M0^(-rho) + (1-beta)*D0^(-rho))^(1/rho)
    end)


    @variables(GSS, begin
        DD >= 0,     (start = D0) # Domestic Demand
        DS >= 0,     (start = D0) # Domestic Supply
        M >= 0,      (start = M0) # Imports
        E >= 0,      (start = E0) # Exports
        X >= 0,      (start = X0) # Compound good?
        Y >= 0,      (start = X0) # PNB
        PX >= 0,     (start = 1,) # Producer price index
        PQ >= 0,     (start = 1,) # Consumer price index
        PED >= 0,    (start = 1,) # Export price including subsidy
        PMD >= 0,    (start = 1,) # Import price including duties
        PDT >= 0,    (start = 1,) # Domestic Price including taxes
        PDD >= 0,    (start = 1,) # Domestic Price excluding taxes
        ER >= 0,     (start = 1,) # Nominal exchange rate
        
        Q,      (start = 100) # Consumption of the compound good
        GR,     (start = 0) # Government Revenue
    end)



    @objective(GSS, Min, Q)

    @constraints(GSS, begin
        OUTPUT, X == A * ((alpha*(E^h)) + ((1-alpha)*(DS^h)))^(1/h);
        CONS, Q == B*(((beta*(M^(-rho))) + ((1-beta)*(DD^(-rho))))^(-1/rho));
        EXPRAT, E == DS*((PED/PDD)*((1-alpha)/alpha))^omega;
        IMPRAT, M == DD*((PDT/PMD)*(beta/(1-beta)))^sigma;
        EXCH, PED == ER*PWE*(1+TE);
        PEXP, PX == (PED*E + PDD*DS)/X;
        PIMP, PMD == PWM*ER*(1+TM);
        PDOM, PDT*DD + PMD*M == PQ*Q;
        PDTEQ, PDT == (1+TD)*PDD;
        GREQ, GR == (TM*ER*PWM*M) + (TD*PDD*DD) - (TE*ER*PWE*E);
        YEQ, Y == (PX*X) + (ER*BBAR) + GR;
        G, X == X0;
        D, DD == DS;
        BOP, PWM*M - PWE*E == BBAR;
    end)

    fix(PQ, 1; force=true)

    return GSS
end

function NLP_report(M::JuMP.Model)
    return DataFrame([
    (
        sigma = value(M[:sigma]),
        omega = value(M[:omega]),
        BBAR = value(M[:BBAR]),
        PWM = value(M[:PWM]),
        PWE = value(M[:PWE]),
        Q = value(M[:Q]),
        PD = value(M[:PDD]),
        TCR = value(M[:PED]/M[:PDD]),
        TCRE = value(M[:PED]/M[:PDD]),
        TCRM = value(M[:PMD]/M[:PDD]),
        TCERQ = value(M[:ER]/M[:PQ]),
        TCERX = value(M[:ER]/M[:PX]),
        TCN = value(M[:ER]/M[:PQ]),
        E = value(M[:E]/M[:PQ]),
        DD = value(M[:DD]/M[:PQ]),
        M = value(M[:M]/M[:PQ]),
    )])
end