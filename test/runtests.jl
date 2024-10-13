using ECCO
using Test

@testset "ECCO.jl" begin

    (f,f_ad,x)=ECCO.toy_problems.enzyme_ex1()
    (f,f_ad,x)=ECCO.toy_problems.enzyme_ex2()
    (f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex3()
    (f,f_ad,x,y)=ECCO.toy_problems.enzyme_ex4()
    f(x,y); f_ad(x,y);
    @test isa(Tuple,y)
    
end
