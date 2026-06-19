module ECCO

# Write your package code here.
import Climatology, MITgcm, Optim, AirSeaFluxes
using DifferentiationInterface

_autodiff_Reverse() = (@warn "undefined")
_Reverse() = (@warn "undefined")
_Active() = (@warn "undefined")
_Duplicated() = (@warn "undefined")

"""
    dTdt_solve()
    
Placeholder for the `OrdinaryDiffEq` extension.

```
using ECCO, OrdinaryDiffEq
?ECCO.dTdt_solve
```
"""
function dTdt_solve()
    @warn "requires using dTdt_solve"
    error("Placeholder for dTdt_solve should never be used")
end

"""
    calc_adjoint(f = y->y[1], backend=AutoForwardDiff(), x=[0.0])

```
using ECCO, Mooncake
using ECCO.DifferentiationInterface

backend = AutoMooncake(; config=nothing)
f(x)=ECCO.AirSeaFluxes.bulkformulae(x[1],x[2],x[3],x[4]).hl
x=[300.,0.001,1.,10.]

(x,adx)=ECCO.calc_adjoint(f,backend,x)
```
"""
function calc_adjoint(f = y->y[1], backend=AutoForwardDiff(), x=[0.0])
    value_and_gradient(f, backend, x)
end

include("initial_examples.jl")
include("mountain_glacier.jl")
include("Lorenz_models.jl")
include("BudykoSellers_model.jl")

Zygote_ex1=Zygote_examples.Zygote_ex1

export Budyko_Sellers_models, Lorenz_models, glacier_model
export toy_problems, Zygote_examples

end
