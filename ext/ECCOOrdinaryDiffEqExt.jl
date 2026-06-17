module ECCOOrdinaryDiffEqExt

using OrdinaryDiffEq
import ECCO 
import ECCO: dTdt_solve

"""
    dTdt_solve(Tini,params,par)

Solve     
"""
function dTdt_solve(Tini,params,par)
    f(u, p, t) = ECCO.Budyko_Sellers_models.dTdt(u,Q=par[1])	
    prob = ODEProblem(f, Tini, (0.0, params.ntime))
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
end

end


