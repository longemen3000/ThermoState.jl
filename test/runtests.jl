using Test
using ThermoState
using ThermoState.Types
using Unitful


@testset "spec" begin
    @test value(spec(t=42))==42
    @test value(spec(t=42.0))==42.0
end

@testset "mass consistency" begin
    #==
    test all possible combinations of: 
    -mass, mass fractions, 
    -moles, mol fractions,
    -mol numbers, mass numbers  
    ==#

    n0 = [1.1,2.2,3.3] #mol
    m0 = [0.011,0.044,0.099] #mass
    mw0 = [10.0,20.0,30.0]u"g/mol" #mw, #testing for correct unit comversion
    xn0 = n0 ./ sum(n0)
    xm0 = m0 ./ sum(m0)
    mass0 = sum(m0)
    moles0 = sum(n0)
    _t = spec(t=1)
    _p = spec(p=2)
    _n = spec(n=n0)
    _m = spec(m=m0)
    _xm = spec(xm=xm0)
    _xn = spec(xn=xn0)
    _mass = spec(mass=mass0)
    _moles = spec(moles=moles0)

    a1 = state(mass=mass0,xm = xm0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    a1 = state(_mass,_xm,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0



    a1 = state(moles=moles0,xm = xm0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0 #fail
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0 #fail
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0 #fail
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    a1 = state(_moles,_xm,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0 #fail
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0 #fail
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0 #fail
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0



    a1 = state(moles=moles0,xn = xn0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0
    
    a1 = state(_moles,_xn,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0



    a1 = state(mass=mass0,xm = xm0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    a1 = state(_mass,_xm,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0



    a1 = state(m = m0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    a1 = state(_m,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0


    a1 = state(n = n0,t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    a1 = state(_n,_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw0)≈ xn0
    @test mass_fraction(FromState(),a1,nothing,mw0)≈ xm0
    @test mass_number(FromState(),a1,u"kg",mw0)≈ m0
    @test mol_number(FromState(),a1,u"mol",mw0)≈ n0
    @test mass(FromState(),a1,u"kg",mw0)≈ mass0
    @test moles(FromState(),a1,u"mol",mw0)≈ moles0

    mw1 = 34.3
    moles1 = 2.3
    mass1 = moles1*mw1*0.001
    _mass = spec(mass=mass1)
    _moles= spec(moles=moles1)
    a1 = state(t=1,p=2)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [0.001*mw1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [1.0]
    @test mass(FromState(),a1,u"kg",mw1)≈ 0.001*mw1
    @test moles(FromState(),a1,u"mol",mw1)≈   1.0

    a1 = state(_t,_p)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [0.001*mw1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [1.0]
    @test mass(FromState(),a1,u"kg",mw1)≈ 0.001*mw1
    @test moles(FromState(),a1,u"mol",mw1)≈   1.0

    a1 = state(t=1,p=2,mass=mass1)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [mass1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [moles1]
    @test mass(FromState(),a1,u"kg",mw1)≈ mass1
    @test moles(FromState(),a1,u"mol",mw1)≈  moles1

    a1 = state(_t,_p,_mass)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [mass1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [moles1]
    @test mass(FromState(),a1,u"kg",mw1)≈ mass1
    @test moles(FromState(),a1,u"mol",mw1)≈  moles1


    
    a1 = state(t=1,p=2,moles=moles1)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [mass1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [moles1]
    @test mass(FromState(),a1,u"kg",mw1)≈ mass1
    @test moles(FromState(),a1,u"mol",mw1)≈  moles1

    a1 = state(_t,_p,_moles)
    @test mol_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_fraction(FromState(),a1,nothing,mw1)≈ [1.0]
    @test mass_number(FromState(),a1,u"kg",mw1)≈ [mass1]
    @test mol_number(FromState(),a1,u"mol",mw1)≈ [moles1]
    @test mass(FromState(),a1,u"kg",mw1)≈ mass1
    @test moles(FromState(),a1,u"mol",mw1)≈  moles1


end

@testset "energy consistency" begin
    #==
    test all possible combinations an energy spec with:: 
    -mass, mass fractions, 
    -moles, mol fractions,
    -mol numbers, mass numbers  
    ==#

    n0 = [1.1,2.2,3.3] #mol
    m0 = [0.011,0.044,0.099] #mass
    mw0 = [10.0,20.0,30.0] #mw, #testing for correct unit comversion
    xn0 = n0 ./ sum(n0)
    xm0 = m0 ./ sum(m0)
    mass0 = sum(m0)
    moles0 = sum(n0)
    mol_h0 = 101.433
    total_h0 = mol_h0*moles0
    mass_h0 = total_h0/mass0
    a1 = state(n = n0,t=1,mol_h = mol_h0)
    @test mol_enthalpy(FromState(),a1,u"J/mol",mw0) ≈ mol_h0
    @test mass_enthalpy(FromState(),a1,u"J/kg",mw0) ≈ mass_h0
    @test total_enthalpy(FromState(),a1,u"J",mw0) ≈ total_h0

    a1 = state(n = n0,t=1,mass_h = mass_h0)
    @test mol_enthalpy(FromState(),a1,u"J/mol",mw0) ≈ mol_h0
    @test mass_enthalpy(FromState(),a1,u"J/kg",mw0) ≈ mass_h0
    @test total_enthalpy(FromState(),a1,u"J",mw0) ≈ total_h0

    a1 = state(n = n0,t=1,total_h = total_h0)
    @test mol_enthalpy(FromState(),a1,u"J/mol",mw0) ≈ mol_h0
    @test mass_enthalpy(FromState(),a1,u"J/kg",mw0) ≈ mass_h0
    @test total_enthalpy(FromState(),a1,u"J",mw0) ≈ total_h0

   #test entropy 
   a1 = state(n = n0,t=1,mol_s = mol_h0)
   @test mol_entropy(FromState(),a1,u"J/(mol*K)",mw0) ≈ mol_h0
   @test mass_entropy(FromState(),a1,u"J/(kg*K)",mw0) ≈ mass_h0
   @test total_entropy(FromState(),a1,u"J/K",mw0) ≈ total_h0

   a1 = state(n = n0,t=1,mass_s = mass_h0)
   @test mol_entropy(FromState(),a1,u"J/(mol*K)",mw0) ≈ mol_h0
   @test mass_entropy(FromState(),a1,u"J/(kg*K)",mw0) ≈ mass_h0
   @test total_entropy(FromState(),a1,u"J/K",mw0) ≈ total_h0

   a1 = state(n = n0,t=1,total_s = total_h0)
   @test mol_entropy(FromState(),a1,u"J/(mol*K)",mw0) ≈ mol_h0
   @test mass_entropy(FromState(),a1,u"J/(kg*K)",mw0) ≈ mass_h0
   @test total_entropy(FromState(),a1,u"J/K",mw0) ≈ total_h0
end

@testset "volume consistency" begin

n0 = [1.1,2.2,3.3] #mol
m0 = [0.011,0.044,0.099] #mass
mw0 = [10.0,20.0,30.0]#mw, #testing for correct unit comversion
xn0 = n0 ./ sum(n0)
xm0 = m0 ./ sum(m0)
mass0 = sum(m0)
moles0 = sum(n0)
V0 = 400.0 #m3
mol_v0 = V0/moles0
mol_rho0 = 1/mol_v0
mass_v0 = V0/mass0
mass_rho0 = 1/mass_v0

a1 = state(n = n0,t=1,total_v = V0)
@test total_volume(FromState(),a1,u"m^3",mw0) ≈ V0
@test mol_volume(FromState(),a1,u"m^3/mol",mw0) ≈ mol_v0
@test mass_volume(FromState(),a1,u"m^3/kg",mw0) ≈ mass_v0
@test mol_density(FromState(),a1,u"mol/m^3",mw0) ≈ mol_rho0
@test mass_density(FromState(),a1,u"kg/m^3",mw0) ≈ mass_rho0

a1 = state(n = n0,t=1,mass_v = mass_v0)
@test total_volume(FromState(),a1,u"m^3",mw0) ≈ V0
@test mol_volume(FromState(),a1,u"m^3/mol",mw0) ≈ mol_v0 
@test mass_volume(FromState(),a1,u"m^3/kg",mw0) ≈ mass_v0
@test mol_density(FromState(),a1,u"mol/m^3",mw0) ≈ mol_rho0
@test mass_density(FromState(),a1,u"kg/m^3",mw0) ≈ mass_rho0

a1 = state(n = n0,t=1,mol_v = mol_v0)
@test total_volume(FromState(),a1,u"m^3",mw0) ≈ V0
@test mol_volume(FromState(),a1,u"m^3/mol",mw0) ≈ mol_v0
@test mass_volume(FromState(),a1,u"m^3/kg",mw0) ≈ mass_v0 
@test mol_density(FromState(),a1,u"mol/m^3",mw0) ≈ mol_rho0
@test mass_density(FromState(),a1,u"kg/m^3",mw0) ≈ mass_rho0

a1 = state(n = n0,t=1,mol_rho = mol_rho0)
@test total_volume(FromState(),a1,u"m^3",mw0) ≈ V0
@test mol_volume(FromState(),a1,u"m^3/mol",mw0) ≈ mol_v0
@test mass_volume(FromState(),a1,u"m^3/kg",mw0) ≈ mass_v0
@test mol_density(FromState(),a1,u"mol/m^3",mw0) ≈ mol_rho0
@test mass_density(FromState(),a1,u"kg/m^3",mw0) ≈ mass_rho0

a1 = state(n = n0,t=1,mass_rho = mass_rho0)
@test total_volume(FromState(),a1,u"m^3",mw0) ≈ V0
@test mol_volume(FromState(),a1,u"m^3/mol",mw0) ≈ mol_v0
@test mass_volume(FromState(),a1,u"m^3/kg",mw0) ≈ mass_v0
@test mol_density(FromState(),a1,u"mol/m^3",mw0) ≈ mol_rho0
@test mass_density(FromState(),a1,u"kg/m^3",mw0) ≈ mass_rho0

end

@testset "state points" begin
using ThermoState.StatePoints

@test pressure(FromState(),StandardConditions()) == 100000.0
@test pressure(FromState(),NormalConditions()) == 101325.0
@test pressure(FromState(),NormalBoilingPoint()) == 101325.0
@test temperature(FromState(),StandardConditions()) == 273.15
@test temperature(FromState(),NormalConditions()) == 293.15
end
