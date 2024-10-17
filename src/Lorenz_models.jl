
module Lorenz_models

import DynamicalSystems

struct Lorenz96{N} end # Structure for size type
function (obj::Lorenz96{N})(dx, x, p, t) where {N}
    F = p[1]
    # 3 edge cases explicitly (performance)
    @inbounds dx[1] = (x[2] - x[N - 1]) * x[N] - x[1] + F
    @inbounds dx[2] = (x[3] - x[N]) * x[1] - x[2] + F
    @inbounds dx[N] = (x[1] - x[N - 2]) * x[N - 1] - x[N] + F
    # then the general case
    for n in 3:(N - 1)
    @inbounds dx[n] = (x[n + 1] - x[n - 2]) * x[n - 1] - x[n] + F
    end
    return nothing
end

"""
x, y, z = Lorenz_models.L96()
plot3D(x, y, z)
"""
function L96(N = 5, F = 8.0, P = 0.01)
    # parameters and initial conditions
    #N = 5
    #F = 8.0
    u₀ = F * ones(N)
    u₀[1] += P # small perturbation

    # The Lorenz-96 model is predefined in DynamicalSystems.jl:
    ds = Systems.lorenz96(N; F = F)

    # Equivalently, to define a fast version explicitly, do:
    lor96 = Lorenz96{N}() # create struct
    ds = DynamicalSystems.ContinuousDynamicalSystem(lor96, u₀, [F])

    # And now evolve a trajectory
    dt = 0.01 # sampling time
    Tf = 30.0 # final time
    tr = DynamicalSystems.trajectory(ds, Tf; dt = dt)

    columns(tr)
end

end
