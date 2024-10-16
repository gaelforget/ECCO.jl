using ECCO
using Test

@testset "ECCO.jl" begin

    (f,f_ad,x)=ECCO.toy_problems.enzyme_ex1()
    (f,f_ad,x)=ECCO.toy_problems.enzyme_ex2()
    (f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex3()
    (f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex4()
    fc=f(x,y)
    adx=f_ad(x,y);
    @test isa(adx,Vector)

    ECCO.toy_problems.f_tau([300.,0.001,1.,10.],[0.0])
    adx=ECCO.toy_problems.f_tau_ad([300.,0.001,1.,10.],[0.0])
    @test isa(adx,Vector)

    ##

    (f,x0,x1,result)=ECCO.toy_problems.optim_ex1()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    (f,g!,x0,x1,result)=ECCO.toy_problems.optim_ex2()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    (f,g!,x0,x1,result)=ECCO.toy_problems.optim_ex3()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    ##
    
    V=ECCO.glacier_model.integrate()
    @test isapprox(V,12.364982412145055)

    store=ECCO.Lorenz_models.L96()
    @test isapprox(store[1,end],-1.334792006090324)

    x,y,z=ECCO.Lorenz_models.L63()
    @test isapprox(x[end],-13.473664284864503)
end
