using Test
using ThermoState
using ThermoState.Types
using ThermoState.QuickStates
using ThermoState: value
using Unitful


@testset "spec" begin
    #value
    @test value(spec(t=42))==42

    #specification
    @test specification(spec(spec"t",1)) == spec"t"

    #@spec_str keyword
    st = spec(t=1)
    st2 = spec(spec"t",1)
    @test st==st2

    #variable spec
    st2 = spec(spec"t",VariableSpec())
    @test st2(1) == st

    #units
    st2 = spec(spec"t",1u"K",false)
    @test st != st2

    st2 = spec(spec"t",1u"K")
    @test st == st2

    
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

@testset "state errors" begin
    
    #@test_throws Exception state(t=1,T=1) #broken
    mass0 = spec(spec"mass",1)
    moles0 = spec(spec"moles",1)
    t0 = spec(spec"t",1)
    p0 = spec(spec"p",1)

    @test_throws Exception state(mass0,moles0,t0,p0)
end

@testset "variable spec" begin
    λ = VariableSpec()
    p0 = 1
    t0 = 1
    moles0 = 1
    st0 = state(t=t0,p=p0,moles=moles0)
    st1 = state(t=λ,p=p0,moles=moles0)
    @test st1(t0) == st0
    st1 = state(t=λ,p=λ,moles=moles0)
    @test st1(t0,p0) == st0
    st1 = state(t=λ,p=λ,moles=λ)
    @test st1(t0,p0,moles0) == st0

    sp0 = spec(t=t0)
    sp1 = spec(t=λ)
    @test sp1(t0) == sp0
end

@testset "state_type" begin
    st = state(t=1,p=2)
    @test state_type(st) isa SinglePT
    @test QuickStates.pt() isa SinglePT

    st = state(t=1,p=2,n=rand(5))
    @test state_type(st) isa MultiPT
    @test QuickStates.ptx() isa MultiPT
    @test QuickStates.ptn() isa MultiPT

    st = state(v=1,p=2)
    @test state_type(st) isa SinglePV
    @test QuickStates.pv() isa SinglePV


    st = state(v=1,p=2,n=rand(5))
    @test state_type(st) isa MultiPV
    @test QuickStates.pvx() isa MultiPV
    @test QuickStates.pvn() isa MultiPV

    st = state(t=2,v=1)
    @test state_type(st) isa SingleVT
    @test QuickStates.vt() isa SingleVT

    st = state(t=2,v=1,n=rand(5))
    @test state_type(st) isa MultiVT
    @test QuickStates.vtx() isa MultiVT
    @test QuickStates.vtn() isa MultiVT

    st = state(h=2,p=1)
    @test state_type(st) isa SinglePH
    @test QuickStates.ph() isa SinglePH

    st = state(h=2,p=1,n=rand(5))
    @test state_type(st) isa MultiPH
    @test QuickStates.phx() isa MultiPH
    @test QuickStates.phn() isa MultiPH

    st = state(sat=true,p=1)
    @test state_type(st) isa SingleSatP
    @test QuickStates.sat_p() isa SingleSatP

    st = state(sat=true,p=1,n=rand(5))
    @test state_type(st) isa MultiSatP
    @test QuickStates.sat_px() isa MultiSatP
    @test QuickStates.sat_pn() isa MultiSatP

    st = state(sat=true,t=1)
    @test state_type(st) isa SingleSatT
    @test QuickStates.sat_t() isa SingleSatT

    st = state(sat=true,t=1,n=rand(5))
    @test state_type(st) isa MultiSatT
    @test QuickStates.sat_tx() isa MultiSatT
    @test QuickStates.sat_tn() isa MultiSatT


end

@testset "@to_units" begin
    xn = rand(5)
    xn = xn ./ sum(xn)
    st = state(sat=true,t=1u"K",xn=xn)
    @test (@to_units temperature(FromState(),st)) == 1u"K"
    @test (@to_units temperature(FromState(),st,u"°C")) ==-272.15u"°C"
    @test (@to_units mol_fraction(FromState(),st)) == xn

end

