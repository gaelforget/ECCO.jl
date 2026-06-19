
# Generic Adjoint Calculation

Using `DifferentiationInterface.jl`

```@example Generic
using ECCO, Mooncake
using ECCO.DifferentiationInterface

f(x)=ECCO.AirSeaFluxes.bulkformulae(x[1],x[2],x[3],x[4]).hl
x=[300.,0.001,1.,10.]
backend = AutoMooncake(; config=nothing)

ad=ECCO.calc_adjoint(f,backend,x)
```

# Optimization

```@example Optim
using ECCO

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
x0 = [0.0, 0.0]
op=ECCO.calc_optim(f,x0)
```

## More Optimization

```
(f,x0,x1,result)=toy_problems.optim_ex1()
(h,h!,x0,x1,result)=toy_problems.optim_ex2()
(h,h!,x0,x1,result)=toy_problems.optim_ex3()
```
