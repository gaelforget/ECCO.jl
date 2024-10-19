
module toy_problems

import Enzyme, Optim, ForwardDiff
import AirSeaFluxes

#export enzyme_ex1, enzyme_ex2, enzyme_ex3, enzyme_ex4

import Enzyme: autodiff, Reverse, Active, Duplicated
import AirSeaFluxes: simpleflux, bulkformulae

"""
    ECCO.toy_problems.enzyme_ex1()

```
using ECCO
(f,f_ad,x)=ECCO.toy_problems.enzyme_ex1()
f(x)
f_ad(x)
```
"""
function enzyme_ex1()
    mld=10 #mixed layer depth (m)
    tim=86400*30 #relaxation time scale (s)
    pisvel=mld/tim #piston velocity (m/s)

    Co=0.0 #ocean value (e.g. concentration of some compound)
    Ca=1.0 #atmospeheric value (e.g. equivalent compound concentration)
    #flx=simpleflux(Ca,Co,pisvel)

    flux(x)=simpleflux(x[1],Ca,pisvel)
    flux_ad(x)=autodiff(Reverse, flux, Active(x[1]))

    (flux,flux_ad,[Co])
end

## univariate example 2

"""
    ECCO.toy_problems.enzyme_ex2()

```
using ECCO
(f,f_ad,x)=ECCO.toy_problems.enzyme_ex2()
f(x...)
f_ad(x...)
```
"""
function enzyme_ex2()
    atemp=300.
    aqh=0.001
    speed=1.
    sst=10.
    #all=bulkformulae(atemp,aqh,speed,sst)

    f_evap(x,y)=bulkformulae(atemp,aqh,x,y).evap
    f_evap_ad(x,y)=autodiff(Reverse, f_evap, Active(x), Active(y))

    (f_evap,f_evap_ad,(speed, sst))
end

## multivariate examples

function f_tau(x::Array{Float64}, y::Array{Float64})
    y[1] = bulkformulae(x[1],x[2],x[3],x[4]).tau
    return y[1]
end        
function f_tau_ad(x=[300.,0.001,1.,10.], y=[0.0])
    bx = zeros(size(x)); by = ones(size(x))
    autodiff(Reverse, f_tau, Duplicated(x, bx), Duplicated(y, by));
    copy(bx)
end

function f_hl(x::Array{Float64}, y::Array{Float64})
    y[1] = bulkformulae(x[1],x[2],x[3],x[4]).hl
    return y[1]
end
function f_hl_ad(x=[300.,0.001,1.,10.], y=[0.0])
    bx = zeros(size(x)); by = ones(size(x))
    autodiff(Reverse, f_hl, Duplicated(x, bx), Duplicated(y, by));
    copy(bx)
end

"""
    ECCO.toy_problems.enzyme_ex3()

```
using ECCO
(f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex3()
f(x,y)
f_ad(x,y)
```    
"""
enzyme_ex3() = (f_tau,f_tau_ad,[300.,0.001,1.,10.],[0.0])

"""
    ECCO.toy_problems.enzyme_ex4()

```
using ECCO
(f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex4()
f(x,y)
f_ad(x,y)
```    
"""
enzyme_ex4() = (f_hl,f_hl_ad,[300.,0.001,1.,10.],[0.0])

## ForwardDiff

ForwardDiff_ex1() = begin
    f(x)=bulkformulae(x[1],x[2],x[3],x[4]).hl
    x=[300.,0.001,1.,10.]
    adx=ForwardDiff.gradient(f, x)
    (x,adx)
end

## Optim

"""
    ECCO.toy_problems.optim_ex1()

```
using ECCO
(f,x0,x1,result)=ECCO.toy_problems.optim_ex1()
```
"""
function optim_ex1()
    f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    x0 = [0.0, 0.0]
    result=Optim.optimize(f, x0)
    x1=Optim.minimizer(result)
    f,x0,x1,result
end

"""
    ECCO.toy_problems.optim_ex2()
"""
function optim_ex2()
    f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
    function g!(G,x)
        G[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
        G[2] = 200.0 * (x[2] - x[1]^2)
    end
    x0 = [0.0, 0.0]
    result=Optim.optimize(f,g!,x0)
    x1=Optim.minimizer(result)
    f,g!,x0,x1,result
end

"""
    ECCO.toy_problems.optim_ex3()
"""
function optim_ex3()
        h(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
        function h!(bx2, x) 
            bx = zeros(size(x))
            Enzyme.autodiff(Reverse, h, Duplicated(x, bx))
            bx2 .= bx
        end
        x0 = [0.0, 0.0]
        H = zeros(2)
        h!(H,x0)
    
        result=Optim.optimize(h,h!,x0)
        x1=Optim.minimizer(result)
        h,h!,x0,x1,result
end

end

## Zygote

module Zygote_examples

using Zygote, AirSeaFluxes

Zygote_ex1() = begin
    f(x)=bulkformulae(x[1],x[2],x[3],x[4]).hl
    x=[300.,0.001,1.,10.]
    adx=f'(x)
    (x,adx)
end

end
