module ECCO

# Write your package code here.
import Climatology, MeshArrays, MITgcm, Optim, Enzyme, AirSeaFluxes

include("initial_examples.jl")
include("mountain_glacier.jl")
include("Lorenz_models.jl")
include("BudykoSellers_model.jl")

Zygote_ex1=Zygote_examples.Zygote_ex1

export Budyko_Sellers_models, Lorenz_models, glacier_model
export toy_problems, Zygote_examples

end
