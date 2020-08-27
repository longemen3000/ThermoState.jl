using Test
using ThermoState
using Unitful
@testset "spec" begin
    @test value(spec(t=42))==42
end

