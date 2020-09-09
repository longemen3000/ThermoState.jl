
 

_stateorder(::Type{Enthalpy{MOLAR}}) = 3
_stateorder(x::Entropy) = 11
_stateorder(x::InternalEnergy) = 5
_stateorder(x::Gibbs) = 2
_stateorder(x::Helmholtz) = 1
_stateorder(x::Pressure) = -10
_stateorder(x::Temperature) = 10
_stateorder(x::VolumeAmount) = 6
_stateorder(x::MaterialCompounds) = 10002
_stateorder(x::MaterialAmount ) = 10001
_stateorder(x::PhaseFractions) = -3
_stateorder(x::VaporFraction) = -1
_stateorder(x::PhaseTag) = 1000
_stateorder(x::TwoPhaseEquilibrium) = -2
_stateorder(x::Options) = 10004
_stateorder(x::MolecularWeight) = 10005
_stateorder(x::VariableSpec) = 1006





function _state_from_type(x::Type{ThermodynamicState{M,C,S}}) where {M,C,S}
    tvec = S.parameters
    state_specs =  [tvec[i].parameters[1]() for i in 1:length(tvec)]
    sorted_specs = sort(state_specs,by=_stateorder)
    res = (sorted_specs[1],sorted_specs[2])
end


@generated function state_type(x::ThermodynamicState)
    res = _state_from_type(x)
    return :($res)
end