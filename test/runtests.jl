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

    (f,x0,x1,result)=ECCO.toy_problems.optim_ex1()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)

    (f,g!,x0,x1,result)=ECCO.toy_problems.optim_ex2()
    dx=1e-4*(x0-x1)
    @test f(x1)<f(x1+dx)
end
