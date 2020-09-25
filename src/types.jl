module Types
    abstract type AbstractSpec end
    abstract type SpecModifier <: AbstractSpec  end

    struct SingleComponent <: AbstractSpec end
    struct OneMol <: AbstractSpec end

    abstract type MassBasisModifier <: SpecModifier end
    struct MOLAR <: MassBasisModifier end
    struct MASS <: MassBasisModifier end
    struct TOTAL <: MassBasisModifier end

    abstract type CompoundModifier <: SpecModifier end
    struct FRACTION <: CompoundModifier end
    struct TOTAL_AMOUNT <: CompoundModifier end

    abstract type VolumeModifier <: SpecModifier end
    struct DENSITY <: VolumeModifier end
    struct VOLUME <: VolumeModifier end

    abstract type AbstractIntensiveSpec{T1<:MassBasisModifier} <: AbstractSpec  end
    abstract type AbstractTotalSpec <: AbstractSpec  end
    abstract type AbstractFractionSpec <: AbstractSpec  end
    abstract type CategoricalSpec <: AbstractSpec  end




    abstract type AbstractEnergySpec{T1} <: AbstractIntensiveSpec{T1} end

    struct Enthalpy{T} <: AbstractEnergySpec{T} end
    struct InternalEnergy{T} <: AbstractEnergySpec{T} end
    struct Gibbs{T} <: AbstractEnergySpec{T} end
    struct Helmholtz{T} <: AbstractEnergySpec{T} end

    struct Entropy{T} <: AbstractIntensiveSpec{T} end

    struct VolumeAmount{T1,T2<:VolumeModifier} <: AbstractIntensiveSpec{T1} end

    struct Pressure <: AbstractSpec end
    struct Temperature <: AbstractSpec end

    # those are vectors, 

    struct MaterialCompounds{T1<:MassBasisModifier,T2<:CompoundModifier} <: AbstractTotalSpec end

    struct MaterialAmount{T1<:MassBasisModifier} <: AbstractTotalSpec end


    struct PhaseFractions <: AbstractFractionSpec end
    struct VolumeFraction <: AbstractFractionSpec end
    struct VaporFraction <: AbstractFractionSpec end

    struct PhaseTag <: CategoricalSpec end
    struct TwoPhaseEquilibrium <: CategoricalSpec end

    struct Options <: CategoricalSpec end

    #not defined for now
    struct MolecularWeight <: AbstractSpec end

    #variable spec, signal to make a variable specs
    struct VariableSpec  end

    #base type to define models
    abstract type ThermoModel end


    export AbstractSpec 
    export SpecModifier 
    export AbstractIntensiveSpec 
    export AbstractTotalSpec 
    export AbstractFractionSpec 
    export CategoricalSpec 
    export SingleComponent 
    export OneMol  
    export MOLAR  
    export MASS  
    export TOTAL  
    export FRACTION 
    export TOTAL_AMOUNT 
    export DENSITY  
    export VOLUME  
    export AbstractEnergySpec  
    export Enthalpy 
    export InternalEnergy
    export Gibbs  
    export Helmholtz  
    export Entropy  
    export VolumeAmount 
    export Pressure 
    export Temperature 
    export MaterialCompounds 
    export MaterialAmount 
    export PhaseFractions 
    export VaporFraction 
    export PhaseTag 
    export TwoPhaseEquilibrium 
    export Options 
    export MolecularWeight
    export VariableSpec
    export ThermoModel
end

using .Types

const KW_TO_SPEC = IdDict{Symbol,Any}(
:h =>  Enthalpy{MOLAR}()
,:g =>  Gibbs{MOLAR}()
,:a =>  Helmholtz{MOLAR}()
,:u =>  InternalEnergy{MOLAR}()

,:mol_h =>  Enthalpy{MOLAR}()
,:mol_g =>  Gibbs{MOLAR}()
,:mol_a =>  Helmholtz{MOLAR}()
,:mol_u =>  InternalEnergy{MOLAR}()

,:mass_h =>  Enthalpy{MASS}()
,:mass_g =>  Gibbs{MASS}()
,:mass_a =>  Helmholtz{MASS}()
,:mass_u =>  InternalEnergy{MASS}()

,:total_h =>  Enthalpy{TOTAL}()
,:total_g =>  Gibbs{TOTAL}()
,:total_a =>  Helmholtz{TOTAL}()
,:total_u =>  InternalEnergy{TOTAL}()

,:s =>  Entropy{MOLAR}()
,:mol_s =>  Entropy{MOLAR}()
,:mass_s =>  Entropy{MASS}()
,:total_s =>  Entropy{TOTAL}()

,:p =>  Pressure()
,:P =>  Pressure()
,:t =>  Temperature()
,:T => Temperature()

,:v =>  VolumeAmount{MOLAR,VOLUME}()
,:mol_v =>  VolumeAmount{MOLAR,VOLUME}()
,:mass_v =>  VolumeAmount{MASS,VOLUME}()
,:total_v =>  VolumeAmount{TOTAL,VOLUME}()
,:V =>  VolumeAmount{TOTAL,VOLUME}()


,:rho =>  VolumeAmount{MOLAR,DENSITY}() 
,:mol_rho =>  VolumeAmount{MOLAR,DENSITY}() 
,:mass_rho =>  VolumeAmount{MASS,DENSITY}() 

,:ρ =>  VolumeAmount{MOLAR,DENSITY}()
,:mol_ρ =>  VolumeAmount{MOLAR,DENSITY}()
,:mass_ρ =>  VolumeAmount{MASS,DENSITY}() 
,:mass =>  MaterialAmount{MASS}()
,:moles =>  MaterialAmount{MOLAR}()
,:xn =>  MaterialCompounds{MOLAR,FRACTION}()
,:xm =>   MaterialCompounds{MASS,FRACTION}()
,:n =>  MaterialCompounds{MOLAR,TOTAL_AMOUNT}()
,:m =>   MaterialCompounds{MASS,TOTAL_AMOUNT}()


,:mw =>  MolecularWeight()

,:vfrac =>  VaporFraction() #looking for better name
,:phase_fracs =>  PhaseFractions() #looking for better name

,:phase =>PhaseTag()

,:sat => TwoPhaseEquilibrium()
,:vle => TwoPhaseEquilibrium()
,:lle => TwoPhaseEquilibrium()

,:single_component => SingleComponent()
,:one_mol => OneMol()
,:options => Options()

)


const SPEC_TO_KW = IdDict{Any,Symbol}(
Enthalpy{MOLAR}() => :molar_h
,Gibbs{MOLAR}() => :molar_g
,Helmholtz{MOLAR}() => :molar_a
,InternalEnergy{MOLAR}() => :molar_u

,Enthalpy{MASS}() => :mass_h
,Gibbs{MASS}()  => :mass_g
,Helmholtz{MASS}() => :mass_a
,InternalEnergy{MASS}() => :mass_u

,Enthalpy{TOTAL}() => :total_h
,Gibbs{TOTAL}() => :total_g
,Helmholtz{TOTAL}() => :total_a
,InternalEnergy{TOTAL}() => :total_u

,Entropy{MOLAR}() => :molar_s
,Entropy{MASS}() => :mass_s
,Entropy{TOTAL}() => :total_s

,Pressure() => :p
,Temperature() => :t

,VolumeAmount{MOLAR,VOLUME}() => :mol_v
,VolumeAmount{MASS,VOLUME}() => :mass_v
,VolumeAmount{TOTAL,VOLUME}() => :total_v


,VolumeAmount{MOLAR,DENSITY}()  => :rho 
,VolumeAmount{MOLAR,DENSITY}()  => :mol_rho
,VolumeAmount{MASS,DENSITY}()  => :mass_rho

,MaterialAmount{MASS}() => :mass
,MaterialAmount{MOLAR}() => :moles
,MaterialCompounds{MOLAR,FRACTION}() => :xn
,MaterialCompounds{MASS,FRACTION}() => :xm
,MaterialCompounds{MOLAR,TOTAL_AMOUNT}() => :n
,MaterialCompounds{MASS,TOTAL_AMOUNT}() =>  :n


,MolecularWeight() => :mw

,VaporFraction() => :vfrac#looking for better name
,PhaseFractions() => :phase_fracs #looking for better name

,PhaseTag() => :phase

,TwoPhaseEquilibrium() => :sat


,SingleComponent() => :single_component
,OneMol() => :one_mol
,Options() => :options

)

module QuickStates
using ..Types
const SinglePT = Tuple{Pressure,Temperature,SingleComponent}
const MultiPT = Tuple{Pressure,Temperature,MaterialCompounds}

pt() = (Pressure(),Temperature(),SingleComponent())
ptx() = (Pressure(),Temperature(),MaterialCompounds{MOLAR,FRACTION}())
ptn() = (Pressure(),Temperature(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SingleVT = Tuple{VolumeAmount,Temperature,SingleComponent}
const MultiVT = Tuple{VolumeAmount,Temperature,MaterialCompounds}

vt() = (VolumeAmount{MOLAR,VOLUME}(),Temperature(),SingleComponent())
vtx() = (VolumeAmount{MOLAR,VOLUME}(),Temperature(),MaterialCompounds{MOLAR,FRACTION}())
vtn() = (VolumeAmount{MOLAR,VOLUME}(),Temperature(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

ρt() = (VolumeAmount{MOLAR,DENSITY}(),Temperature(),SingleComponent())
ρtx() = (VolumeAmount{MOLAR,DENSITY}(),Temperature(),MaterialCompounds{MOLAR,FRACTION}())
ρtn() = (VolumeAmount{MOLAR,VOLUME}(),Temperature(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SinglePS = Tuple{Pressure,Entropy,SingleComponent}
const MultiPS = Tuple{Pressure,Entropy,MaterialCompounds}

st() = (Pressure(),Entropy{MOLAR}(),SingleComponent())
stx() = (Pressure(),Entropy{MOLAR}(),MaterialCompounds{MOLAR,FRACTION}())
stn() = (Pressure(),Entropy{MOLAR}(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())


const SinglePH = Tuple{Pressure,Enthalpy,SingleComponent}
const MultiPH = Tuple{Pressure,Enthalpy,MaterialCompounds}

pht() = (Pressure(),Enthalpy{MOLAR}(),SingleComponent())
phx() = (Pressure(),Enthalpy{MOLAR}(),MaterialCompounds{MOLAR,FRACTION}())
phn() = (Pressure(),Enthalpy{MOLAR}(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SingleSatT = Tuple{TwoPhaseEquilibrium,Temperature,SingleComponent}
const MultiSatT = Tuple{TwoPhaseEquilibrium,Temperature,MaterialCompounds}

sat_t() = (TwoPhaseEquilibrium(),Temperature(),SingleComponent())
sat_tx() = (TwoPhaseEquilibrium(),Temperature(),MaterialCompounds{MOLAR,FRACTION}())
sat_tn() = (TwoPhaseEquilibrium(),Temperature(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SingleΦT = Tuple{VaporFraction,Temperature,SingleComponent}
const MultiΦT = Tuple{VaporFraction,Temperature,MaterialCompounds}

ϕt() = (VaporFraction(),Temperature(),SingleComponent())
ϕtx() = (VaporFraction(),Temperature(),MaterialCompounds{MOLAR,FRACTION}())
ϕtn() = (VaporFraction(),Temperature(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SingleSatP = Tuple{TwoPhaseEquilibrium,Pressure,SingleComponent}
const MultiSatP = Tuple{TwoPhaseEquilibrium,Pressure,MaterialCompounds}

sat_p() = (TwoPhaseEquilibrium(),Pressure(),SingleComponent())
sat_px() = (TwoPhaseEquilibrium(),Pressure(),MaterialCompounds{MOLAR,FRACTION}())
sat_pn() = (TwoPhaseEquilibrium(),Pressure(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

const SingleΦP = Tuple{VaporFraction,Pressure,SingleComponent}
const MultiΦP = Tuple{VaporFraction,Pressure,MaterialCompounds}

ϕp() = (VaporFraction(),Pressure(),SingleComponent())
ϕpx() = (VaporFraction(),Pressure(),MaterialCompounds{MOLAR,FRACTION}())
ϕpn() = (VaporFraction(),Pressure(),MaterialCompounds{MOLAR,TOTAL_AMOUNT}())

export SinglePT,MultiPT
export SingleVT,MultiVT
export SinglePS,MultiPS
export SinglePH,MultiPH
export SingleSatT,MultiSatT
export SingleΦT,MultiΦT
export SingleSatP,MultiSatP
export SingleΦP,MultiΦP


end





