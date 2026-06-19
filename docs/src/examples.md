
# Generic Adjoint Calculation

Using `DifferentiationInterface.jl`

```@example Generic
using ECCO, Mooncake
using ECCO.DifferentiationInterface

f(x)=ECCO.AirSeaFluxes.bulkformulae(x[1],x[2],x[3],x[4]).hl
x=[300.,0.001,1.,10.]
backend = AutoMooncake(; config=nothing)

#(x,adx)=ECCO.DifferentiationInterface_example.ex1(backend)
(x,adx)=ECCO.calc_adjoint(f,backend,x)
```

# Optimization

```@example Optim
using ECCO
(f,x0,x1,result)=toy_problems.optim_ex1()
```

```@example Optim
h,h!,x0,x1,result=toy_problems.optim_ex2()
```

```@example Optim
h,h!,x0,x1,result=toy_problems.optim_ex3()
```
