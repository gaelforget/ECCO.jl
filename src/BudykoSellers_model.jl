
module Budyko_Sellers_models

using OrdinaryDiffEq, ForwardDiff, Statistics, Optim

##1. setup

"""# Budyko-Sellers Energy Balance Model

### Parameter Choices

From B. Rose notebook :

R = 10^7 J/m2/K

From Walsh and Rackauckas :

- Figure 3. Equilibrium solutions (7) with albedo function (9) for five η-values. Note T∗ ηi (ηi) = Tc only for i = 2,5. 
  - Parameters: Q = 343, A = 202, B = 1.9, C = 3.04, αw = 0.32, αs = 0.62, Tc =−10.
- Figure 9. Equilibrium solutions of (2) with albedo function (34). Solid: η= 0.1. Dashed: η= 0.25. Dash-Dot: η= 0.4. 
  - Parameters: Q = 321,A = 167,B = 1.5,C = 2.25,M = 50,αw = 0.32,αi = 0.46,αs = 0.72,ρ= 0.35.
"""
params=(y=-89:1:89, η=70, R=10^7,
		Q = 343.0, A = 202, B = 1.9, C = 3.04, 
		αw = 0.32, αs = 0.62, Tbar = -10.0, 
		ntime=86400*365*10)

#Insolation        
s2 = −0.482
p2(y) = 1/2 * (3* y.^2 .− 1)
s(y) = 1 .+ s2*p2(y./90)

#Albedo
(; αw, αs, η) = params 
α(y) = 
if abs(y)<η
    αw
elseif abs(y)>η
    αs
else
    0.5*(αw + αs)
end

#Equations
(; Q, A, B, C, Tbar, y) = params
Tini=fill(Tbar,size(y))
Ein(y;Q=Q)=Q[1]*s.(y).*(1 .−α.(y))
Eout(T) = A.+B*T
Etran(T)=C*(T.−Tbar)
#Etran = D dT/dx^2
#Eout = σT4(1−mtanh(19 T^6 ×10^(−16)))

dTdt(T;Q=Q)= (Ein(params.y,Q=Q).−Eout(T).−Etran.(T))./params.R

function dTdt_solve(par=[Q])
    f(u, p, t) = dTdt(u,Q=par[1])	
    prob = ODEProblem(f, Tini, (0.0, params.ntime))
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
end

"""
    dTdt_demo(par=[Q])

```
using ECCO, CairoMakie

(; Q, y) = Budyko_Sellers_models.params
Tsol,Tini,dTdt_ini,incr_t,incr=Budyko_Sellers_models.dTdt_demo(Q)

fig=Figure()
Axis(fig[1,1]); lines!(y,dTdt_ini)
Axis(fig[2,1]); lines!(y,Tini); lines!(y,Tsol,color=:red)
Axis(fig[3,1]); lines!(incr_t,incr)	
fig
```
"""
function dTdt_demo(par=[Q])
    sol=dTdt_solve([Q])
	Tsol = sol.u[end]
    
	incr=[maximum(abs.(sol.u[t+1]-sol.u[t]))/(sol.t[t+1]-sol.t[t]) for t in 1:length(sol.t)-1]
	incr_t=0.5*(sol.t[2:end]+sol.t[1:end-1])
    dTdt_ini=Budyko_Sellers_models.dTdt(Tini)
    Tsol,Tini,dTdt_ini,incr_t,incr
end

function fake_obs(x)
   dt=1
   result = zeros(eltype(x),179)
   result.=Tini
   for i in 1:10000
	   result .+= dt*dTdt(result,Q=x)
   end
   result
end

"""
    the_main_loop(x; nt=10000)

Do `nt` time steps of `dt*dTdt(result,Q=x)`.
"""
function the_main_loop(x; nt=10000)
   dt=1
   result = zeros(eltype(x),179)
   result.=Tini
   for i in 1:nt
	   result .+= dt*dTdt(result,Q=x)
   end
   return result
end


"""
    adjoint_demo(par=[Q])

- generate obs for cost function with `fake_obs`.
- call `the_main_loop_cost` and `adthe_main_loop_cost`
- return gradient check result

`grdchk=adjoint_demo([Q])`
"""
function adjoint_demo(par=[Q])
    Qobs=0.9*Q
    Tobs=fake_obs(Qobs)

    the_main_loop_cost(x)=mean((the_main_loop(x)-Tobs).^2)
    adthe_main_loop_cost(x)=ForwardDiff.gradient(the_main_loop_cost, x)

    the_main_loop_cost(Q)
    adQ=adthe_main_loop_cost([Q])

    #gradient check
    grdchk=(adQ,the_main_loop_cost(Q+0.5)-the_main_loop_cost(Q-0.5))
end

"""### Optimize Equilibrium solution"""


"""
    dTdt_solve_optim(Tobs)

Call `dTdt_solve` and then `Optim.optimize` (adjoint free).

`f,x0,x1,result=dTdt_solve_optim(Tobs)`
"""
function dTdt_solve_optim(Tobs)
	f(x)=mean((dTdt_solve(x).u[end]-Tobs).^2)
	x0 = [Q]
	result=Optim.optimize(f, x0)
	x1=Optim.minimizer(result)
	f,x0,x1,result
end

"""
    optim_demo(Q=Q; verbose=false)

Call `dTdt_solve` and then `dTdt_solve_optim`.

`optim_demo()`
"""
function optim_demo(Q=Q; verbose=false)
    Qobs=0.9*Q
    Tobs=dTdt_solve(Qobs).u[end]
    (f,x0,x1,Q)=dTdt_solve_optim(Tobs)
    verbose ? println("first guess=$(x0) optim=$(x1) truth=$(Qobs)") : nothing
end

#let
#	lines(dTdt_solve(x0_2).u[end],label="first guess")
#	lines!(dTdt_solve(x1_2).u[end],label="optim")
#	scatter!(dTdt_solve(Qobs).u[end],label="obs",marker=:x)
#	axislegend()
#	current_figure()
#end

"""
    dTdt_loop_optim()

Call `the_main_loop` and then `Optim.optimize`.

`f,x0,x1,result=dTdt_loop_optim()`
"""
function dTdt_loop_optim()
    Qobs=0.9*Q
    Tobs=fake_obs(Qobs)
    f(x)=mean((the_main_loop(x)-Tobs).^2)
	x0 = [Q]
	result=Optim.optimize(f, x0)
	x1=Optim.minimizer(result)
	f,x0,x1,result
end

"""
    optim_demo_loop()

Call `dTdt_loop_optim`.
"""
function optim_demo_loop()
	f,x0,x1,result=dTdt_loop_optim()
	(y0,y1)=(f(x0),f(x1))
	[x0 x1 0.9*Q]
end

end
