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
    
end
