
 

_stateorder(x::PhaseFractions) = 1
_stateorder(x::TwoPhaseEquilibrium) = 2
_stateorder(x::VaporFraction) = 3
_stateorder(x::Pressure) = 4
_stateorder(x::Helmholtz) = 5
_stateorder(x::Gibbs) = 6
_stateorder(x::InternalEnergy) = 7
_stateorder(x::VolumeAmount) = 8
_stateorder(x::Temperature) = 9
_stateorder(x::Enthalpy) = 10
_stateorder(x::Entropy) = 11

_stateorder(x::MaterialCompounds) = 10002
_stateorder(x::MaterialAmount ) = 10001
_stateorder(x::PhaseTag) = 10000
_stateorder(x::Options) = 10004
_stateorder(x::MolecularWeight) = 10005
_stateorder(x::VariableSpec) = 10006


function _state_from_type(x::Type{ThermodynamicState{M,C,S}}) where {M,C,S}
    tvec = S.parameters
    state_specs =  [tvec[i].parameters[1]() for i in 1:length(tvec)]
    has_compound = false
    for spec in state_specs
        if typeof(spec) <: MaterialCompounds
            has_compound = true
            comp_spec = spec 
            break
        end
    end
    if !has_compound
        comp_spec = SingleComponent()
    end
    sorted_specs = sort(state_specs,by=_stateorder)
    res = (sorted_specs[1],sorted_specs[2],comp_spec)
end

"""
    state_type(x::ThermodynamicState)

Returns a sorted tuple of thermodynamic specifications of the input state,
 the tuple has 3 elements, corresponding to two general specifications and a last amount specification,
  indicating if the state is single or multicomponent.


"""
@generated function state_type(x::ThermodynamicState)
    res = _state_from_type(x)
    return :($res)
end

const SinglePT = Tuple{Pressure,Temperature,SingleComponent}
const MultiPT = Tuple{Pressure,Temperature,MaterialCompounds}

const SingleVT = Tuple{VolumeAmount,Temperature,SingleComponent}
const MultiVT = Tuple{VolumeAmount,Temperature,MaterialCompounds}

const SinglePS = Tuple{Pressure,Entropy,SingleComponent}
const MultiPS = Tuple{Pressure,Entropy,MaterialCompounds}

const SinglePH = Tuple{Pressure,Enthalpy,SingleComponent}
const MultiPH = Tuple{Pressure,Enthalpy,MaterialCompounds}
