using Test
using PhysPropsRules
using Unitful
@testset "spec" begin
    @test value(spec(t=42))==42
end

@testset "specs - material" begin
#incorrect vector constructors
    @test_throws ArgumentError specs(t=1,p=1,m=2)
    @test_throws ArgumentError specs(t=1,p=1,n=2)
    @test_throws ArgumentError specs(t=1,p=1,xn=2)
    @test_throws ArgumentError specs(t=1,p=1,xm=2)
end
