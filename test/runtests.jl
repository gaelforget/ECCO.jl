using ECCO
using Zygote, OrdinaryDiffEq, Enzyme, Mooncake
using ECCO.DifferentiationInterface
using Test

@testset "ECCO.jl" begin

    ## Enzyme

    ad=toy_problems.Enzyme_ex1()
    ad=toy_problems.Enzyme_ex2()
    ad=toy_problems.Enzyme_ex3()
    ad=toy_problems.Enzyme_ex4()
    @test isa(ad.adx,Vector)

    toy_problems.f_tau([300.,0.001,1.,10.],[0.0])
    adx=toy_problems.f_tau_ad([300.,0.001,1.,10.],[0.0])
    @test isa(adx,Vector)

    if VERSION<v"1.11.0"
      x=[0.004]
      adx=glacier_model.adjoint_problem(x)
      @test isapprox(adx[1][1],2647.01879892872)
    end

    ## ForwardDiff

    ad=toy_problems.ForwardDiff_ex1()
    @test isapprox(ad.adx[2],458.8925283180731)

    ## Zygote

    ad=ECCO.Zygote_ex1()
    @test isapprox(ad.adx[2],458.8925283180731)

    ## DifferentiationInterface

    ad=ECCO.DifferentiationInterface_example.ex1()
    @test isapprox(ad.adx[2],458.8925283180731)

    backend = AutoMooncake(; config=nothing)
    ad=ECCO.DifferentiationInterface_example.ex1(backend)
    @test isapprox(ad.adx[2],458.8925283180731)

    ## Optimization

    op1=ECCO.calc_optim()
    op2=ECCO.calc_optim_ad()

    op=toy_problems.optim_ex1()
    dx=1e-4*(op.x0-op.x1)
    @test op.f(op.x1)<op.f(op.x1+dx)

    op=toy_problems.optim_ex2()
    dx=1e-4*(op.x0-op.x1)
    @test op.f(op.x1)<op.f(op.x1+dx)

    op=toy_problems.optim_ex3()
    dx=1e-4*(op.x0-op.x1)
    @test op.f(op.x1)<op.f(op.x1+dx)

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
