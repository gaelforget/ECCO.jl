using ECCO
using Test

@testset "ECCO.jl" begin

    ## Enzyme

    (f,f_ad,x)=toy_problems.Enzyme_ex1()
    (f,f_ad,x)=toy_problems.Enzyme_ex2()
    (f,f_ad,x,y)=toy_problems.Enzyme_ex3()
    (f,f_ad,x,y)=toy_problems.Enzyme_ex4()
    fc=f(x,y)
    adx=f_ad(x,y);
    @test isa(adx,Vector)

    toy_problems.f_tau([300.,0.001,1.,10.],[0.0])
    adx=toy_problems.f_tau_ad([300.,0.001,1.,10.],[0.0])
    @test isa(adx,Vector)

    if VERSION<v"1.11.0"
      x=[0.004]
      adx=glacier_model.adjoint_problem(x)
      @test isapprox(adx[1][1],2647.01879892872)
    end

    ## ForwardDiff

    (x,adx)=toy_problems.ForwardDiff_ex1()
    @test isapprox(adx[2],458.8925283180731)

    ## Zygote

    (x,adx)=ECCO.Zygote_ex1()
    @test isapprox(adx[2],458.8925283180731)

    ##

    (f,x0,x1,result)=toy_problems.optim_ex1()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    (f,g!,x0,x1,result)=toy_problems.optim_ex2()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    (f,g!,x0,x1,result)=toy_problems.optim_ex3()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    ##
    
    V=glacier_model.forward_problem()
    @test isapprox(V,12.36497768703625)

    store=Lorenz_models.L96()
    @test isapprox(store[1,end],-1.334792006090324)

    x,y,z=Lorenz_models.L63()
    @test isapprox(x[end],-13.473664284864503)

    Budyko_Sellers_models.adjoint_demo()
    Budyko_Sellers_models.dTdt_demo()
    Budyko_Sellers_models.optim_demo()
    (x0,x1,xobs)=Budyko_Sellers_models.optim_demo_loop()
    @test isapprox(x0,x1,rtol=1)
end
