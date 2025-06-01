

## Bulk formulae

## Lorenz

See https://en.wikipedia.org/wiki/Lorenz_96_model

```@example L1
using ECCO, CairoMakie
store=Lorenz_models.L96()
lines(store[1,:]); lines!(store[2,:]); lines!(store[end,:])
current_figure()
```

## Budyko-Sellers Energy Balance Model

See [this tutorial](https://brian-rose.github.io/ClimateLaboratoryBook/courseware/one-dim-ebm.html) for detailed explanations.

```@example B-S1
using ECCO, CairoMakie

(; Q, y) = Budyko_Sellers_models.params
Tsol,Tini,dTdt_ini,incr_t,incr=Budyko_Sellers_models.dTdt_demo(Q)

fig=Figure()
Axis(fig[1,1]); lines!(y,dTdt_ini)
Axis(fig[2,1]); lines!(y,Tini); lines!(y,Tsol,color=:red)
Axis(fig[3,1]); lines!(incr_t,incr)	
fig
```

### Parameter Choices

From B. Rose notebook :

R = 10^7 J/m2/K

From Walsh and Rackauckas :

- Figure 3. Equilibrium solutions (7) with albedo function (9) for five η-values. Note T∗ ηi (ηi) = Tc only for i = 2,5. 
  - Parameters: Q = 343, A = 202, B = 1.9, C = 3.04, αw = 0.32, αs = 0.62, Tc =−10.
- Figure 9. Equilibrium solutions of (2) with albedo function (34). Solid: η= 0.1. Dashed: η= 0.25. Dash-Dot: η= 0.4. 
  - Parameters: Q = 321,A = 167,B = 1.5,C = 2.25,M = 50,αw = 0.32,αi = 0.46,αs = 0.72,ρ= 0.35.

### Earlier Implementations

- <https://brian-rose.github.io/ClimateLaboratoryBook/courseware/one-dim-ebm.html>
- <https://github.com/ECCO-Summer-School/ESS25-Team_FLOW>
- <https://www.cise.ufl.edu/~luke.morris/2_4_2025/build/bsh/budyko_sellers_halfar/>

### Related References

- discrete and continuous -- on the Budyko-Bellers energy balance climate model with ice line coupling -- James Walsh, Christopher Rackauckas -- doi:10.3934/dcdsb.2015.20.2187
- Theory of Energy-Balance Climate Models -- Gerald R. North -- 1975 -- DOI: <https://doi.org/10.1175/1520-0469(1975)032<2033:TOEBCM>2.0.CO;2>
- Predictability in a Solvable Stochastic Climate Model -- Gerald R. North and Robert F. Cahalan -- DOI: <https://doi.org/10.1175/1520-0469(1981)038%3C0504:PIASSC%3E2.0.CO;2>

## Glacier


