
module Lorenz_models

#Lorenz 1996 model

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
    L96(; N=5, F=8)

See https://en.wikipedia.org/wiki/Lorenz_96_model

```
using ECCO, CairoMakie
store=ECCO.Lorenz_models.L96()
lines(store[1,:]); lines!(store[2,:]); lines!(store[end,:])
current_figure()
```
"""
function L96(; N=5, F=8)
    L=Lorenz96{N}()
    x=F*ones(N); x[1] += 0.01 # small perturbation
    dx=0*x;  t=[0.0]; dt=0.01; nt=1000; store=zeros(N,nt)    
    for n=1:nt
        L(dx,x,[F],t)
        x.=x+dt*dx
        store[:,n].=x
    end
    store    
end

## Lorenz 1963 

# define the Lorenz63 attractor
@kwdef mutable struct Lorenz63
    dt::Float64 = 0.02
    σ::Float64 = 10
    ρ::Float64 = 28
    β::Float64 = 8/3
    x::Float64 = 2
    y::Float64 = 1
    z::Float64 = 1
end

function step!(l::Lorenz63)
    dx = l.σ * (l.y - l.x);         l.x += l.dt * dx
    dy = l.x * (l.ρ - l.z) - l.y;   l.y += l.dt * dy
    dz = l.x * l.y - l.β * l.z;     l.z += l.dt * dz
end

"""
    L63(; nt=10000)

See https://en.wikipedia.org/wiki/Lorenz_system

```
using ECCO, CairoMakie
x,y,z=ECCO.Lorenz_models.L63()
lines(x,y)
```
"""
function L63(; nt=10000)
    attractor = Lorenz63()
    store=Lorenz63[]
    for i=1:nt
        step!(attractor)
        push!(store,deepcopy(attractor))
    end
    x=[store[i].x for i in eachindex(store)]
    y=[store[i].y for i in eachindex(store)]
    z=[store[i].z for i in eachindex(store)]
    return x,y,z
end

end
