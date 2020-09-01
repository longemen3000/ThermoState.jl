using Test
using ThermoState
using Unitful
@testset "spec" begin
    @test value(spec(t=42))==42
    @test value(spec(t=42.0))==42.0
end

