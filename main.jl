using Pkg
Pkg.activate(".")
Pkg.instantiate()

using JuMP
using Ipopt
using MPSGE
using DataFrames

include("models.jl")

data = ModelData()
M = NLP_model(data)

## Recreation of Table 5.4 and 5.5

shocks = [
    NLPModelParameters(BBAR = 0,    sigma = 0.2, omega = 0.2,),
    NLPModelParameters(BBAR = 10.0, sigma = 0.2, omega = 0.2,),
    NLPModelParameters(BBAR = 10.0, sigma = 0.5, omega = 0.5,),
    NLPModelParameters(BBAR = 10.0, sigma = 2,   omega = 2,  ),
    NLPModelParameters(BBAR = 10.0, sigma = 5,   omega = 5,  ),
    NLPModelParameters(BBAR = 10.0, sigma = 15,  omega = .2, ),
    NLPModelParameters(BBAR = 10.0, sigma = .2,  omega = 15, ),

    #NLPModelParameters(PWM = 1.1, sigma = 0.2, omega = 0.2,),
    #NLPModelParameters(PWM = 1.1, sigma = 0.5, omega = 0.5,),
    #NLPModelParameters(PWM = 1.1, sigma = 2,   omega = 2,  ),
    #NLPModelParameters(PWM = 1.1, sigma = 5,   omega = 5,  ),
    #NLPModelParameters(PWM = 1.1, sigma = 15,  omega = .2, ),
    #NLPModelParameters(PWM = 1.1, sigma = .2,  omega = 15, ),

]

report = DataFrame()

for shock in shocks
    set_parameter_values(M, shock)
    optimize!(M)
    report = vcat(report, NLP_report(M))
end


report