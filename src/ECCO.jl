module ECCO

# Write your package code here.
import Climatology, MITgcm, Optim, Enzyme, AirSeaFluxes

include("initial_examples.jl")
include("mountain_glacier.jl")
include("Lorenz_models.jl")
include("BudykoSellers_model.jl")

Zygote_ex1=Zygote_examples.Zygote_ex1
DifferentiationInterface_ex1=DifferentiationInterface_example.DifferentiationInterface_ex1


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

export Budyko_Sellers_models, Lorenz_models, glacier_model
export toy_problems, Zygote_examples

end
