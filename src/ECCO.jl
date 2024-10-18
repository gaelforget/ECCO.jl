module ECCO

# Write your package code here.
import Climatology, MeshArrays, MITgcm, Optim, Enzyme, AirSeaFluxes

include("initial_examples.jl")
include("mountain_glacier.jl")
include("Lorenz_models.jl")

end
