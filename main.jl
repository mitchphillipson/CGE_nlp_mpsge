using Pkg
Pkg.activate(".")
Pkg.instantiate()

using JuMP
using Ipopt
using MPSGE
using DataFrames

include("models.jl")
include("mpsge.jl")

begin
    data = ModelData()
    M = NLP_model(data)
    MP = MPSGE_model(data)

    params = NLPModelParameters(sigma=.2, omega=15, BBAR=10, PWE=1.2, PWM = 1.1, TD = .5, TE = .01, TM = .1)
    set_parameter_values(M, params)
    set_parameter_values(MP, params)

    optimize!(M)
    solve!(MP)
end


outerjoin(
    report(M) |> x -> stack(x, value_name = :NLP),
    report(MP) |> x -> stack(x, value_name = :MPSGE),
    on = [:variable]
) |>
x -> transform(x,
    [:NLP, :MPSGE] => ByRow((a,b) -> round(a-b, digits=6)) => :diff
)


value(M[:PED])

value(MP[:PFX])





value(MP[:PFX]*(1 + MP[:TE])*MP[:PWE]) # PED






value(MP[:PFX])/value(M[:PED])

value(MP[:PFX])
value(1 - MP[:TE])
value(MP[:PWE])




## Recreation of Table 5.4 and 5.5


data = ModelData()
M = NLP_model(data)

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

output = DataFrame()

for shock in shocks
    set_parameter_values(M, shock)
    optimize!(M)
    output = vcat(output, report(M))
end


output