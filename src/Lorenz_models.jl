
module Lorenz_models

#import DynamicalSystems

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

end
