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

    struct Mass <: AbstractTotalSpec end
    struct Moles <: AbstractTotalSpec end

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
    export Mass 
    export Moles 
    export MaterialCompounds 
    export MaterialAmount 
    export PhaseFractions 
    export VaporFraction 
    export PhaseTag 
    export TwoPhaseEquilibrium 
    export Options 
    export MolecularWeight
    export VariableSpec
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

