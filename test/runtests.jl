using Test
using PhysPropsRules
using Unitful
@testset "spec" begin
    @test value(spec(t=42))==42
end

