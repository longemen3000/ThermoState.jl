print_spec(io::IO,t) = print(io,t) #fallback

print_spec(io::IO,x::Type{MOLAR}) = print(io,"Molar ")
print_spec(io::IO,x::Type{MASS}) = print(io,"Mass ")
print_spec(io::IO,x::Type{TOTAL}) = print(io,"Total ")
print_spec(io::IO,x::Type{OneMol}) = print(io,"One mol")
print_spec(io::IO,x::Type{SingleComponent}) = print(io,"Single component")


print_spec(io::IO,x::Type{HumidityDewPoint}) = print(io,"Humidity dew point")
print_spec(io::IO,x::Type{HumidityRatio}) = print(io,"Humidity ratio")
print_spec(io::IO,x::Type{WetBulbTemperature}) = print(io,"Web bulb temperature")
print_spec(io::IO,x::Type{RelativeHumidity}) = print(io,"Relative humidity")
print_spec(io::IO,x::Type{MolarHumidity}) = print(io,"Humidity wet molar fraction")
print_spec(io::IO,x::Type{MassHumidity}) = print(io,"Humidity wet molar fraction")

print_spec(io::IO,::Type{Pressure}) = print(io,"Pressure")
print_spec(io::IO,::Type{Temperature}) = print(io,"Temperature")


function print_spec(io::IO,::Type{Gibbs{T2}}) where T2
    print_spec(io,T2)
    print(io,"Gibbs energy")
end

function print_spec(io::IO,::Type{VaporFraction{T2}}) where T2
    print_spec(io,T2)
    print(io,"vapor fraction")
end

function print_spec(io::IO,::Type{HumiditySpec{T2}}) where T2
    print_spec(io,T2)
end

function print_spec(io::IO,::Type{Helmholtz{T2}}) where T2
    print_spec(io,T2)
    print(io,"Helmholtz energy")
end

function print_spec(io::IO,::Type{Enthalpy{T2}}) where T2
    print_spec(io,T2)
    print(io,"enthalpy")
end

function print_spec(io::IO,::Type{Entropy{T2}}) where T2
    print_spec(io,T2)
    print(io,"entropy")
end

function print_spec(io::IO,::Type{InternalEnergy{T2}}) where T2
    print_spec(io,T2)
    print(io,"internal energy")
end


function print_spec(io::IO,::Type{MaterialAmount{MOLAR}})
    print(io,"Moles")
end

function print_spec(io::IO,::Type{MaterialAmount{MASS}})
    print(io,"Mass")
end

function print_spec(io::IO,::Type{MaterialCompounds{MASS,FRACTION}})
    print(io,"Mass fraction")
end

function print_spec(io::IO,::Type{PhaseTag})
    print(io,"Phase specification")
end

function print_spec(io::IO,::Type{MaterialCompounds{MOLAR,FRACTION}})
    print(io,"Molar fraction")
end

function print_spec(io::IO,::Type{MaterialCompounds{MASS,TOTAL_AMOUNT}})
    print(io,"Mass amounts")
end

function print_spec(io::IO,::Type{MaterialCompounds{MOLAR,TOTAL_AMOUNT}})
    print(io,"Molar amounts")
end


function print_spec(io::IO,::Type{VolumeAmount{T,VOLUME}}) where T
    print_spec(io,T)
    print(io,"volume")
end

function print_spec(io::IO,::Type{VolumeAmount{T,DENSITY}}) where T
    print_spec(io,T)
    print(io,"density")
end

function print_spec(io::IO,::Type{TwoPhaseEquilibrium})
    print(io,"Two phase equilibrium")
end

function print_spec(io::IO,::Type{Options})
    print(io,"Options")
end



function Base.show(io::IO, sp::Spec{T}) where T 
    print(io,"spec(")
    
    a = value(sp)
    units = default_units(T)

    if a isa VariableSpec
        print(io,"x₀)(")
        a = "x₀"
    end
    print(io,string(SPEC_TO_KW[specification(sp)])," = ")    

    if is_unitful(a)
        units = unit(a)
        a = ustrip(a)
    end
    print(io,a)
    if is_real(a) 
        if !(units == Unitful.NoUnits)
            printstyled(io,'[',color = :light_black)
            print(IOContext(io,:compact => true),units)
            printstyled(io,']',color = :light_black)  
        end       
    end 
    print(io,")")
end

function print_long_spec(io::IO, sp::Spec{T}) where T 
    print_spec(io,T)
    print(io," : ")
    a = value(sp)
    units = default_units(T)

    if is_unitful(a)
        units = unit(a)
        a = ustrip(a)
    end
    print(io,a)
    if is_real(a) 
        if !(units == Unitful.NoUnits)
            printstyled(io,'[',color = :light_black)
            print(IOContext(io,:compact => true),units)
            printstyled(io,']',color = :light_black)
        end        
    end       
end


#==function Base.show(::IO, sp::T) where T <:AbstractSpec
    print(io,"::" )
    print_spec(io,T)
end
==#

function Base.show(io::IO, ::MIME"text/plain", sp::ThermodynamicState)
    @nospecialize sp
    lcall = length(sp.callables)
    lsp = length(sp.specs) 
    subscripts = '₁':'₉'

    print(io,"ThermodynamicState")
    if !iszero(lcall)
        print(io,"(")
        for (i,spec) in enumerate(sp.callables)
            (i !=1) && print(io,",")
            print(io,"x",subscripts[i])
        end
    print(io,")")
    end

    if isone(lcall+lsp)
        println(io," with 1 property:")
    else
        println(io," with ",lcall+lsp, " properties:")
    end

    if !iszero(lcall)
        for (i,spec) in enumerate(sp.callables)
            (i !=1) && println(io)
            print("  ")
            print_spec(io,typeof(spec))
            print(io," : ")
            print(io,"x",subscripts[i])
        end
        !iszero(lsp) && println(io)
    end

    if !iszero(lsp)
        for (i,spec) in enumerate(sp.specs)
            (i !=1) && println(io)
            print("  ")
            print_long_spec(io,spec)
        end
    end
end

function Base.show(io::IO, sp::ThermodynamicState)
    lcall = length(sp.callables)
    lsp = length(sp.specs) 
    subscripts = '₁':'₉'

    print(io,"state(")

    if !iszero(lcall)
        for (i,spec) in enumerate(sp.callables)
            (i !=1) && print(io,",")
            print(io,"x",subscripts[i])
        end
        print(io,")(")
    
        for (i,spec) in enumerate(sp.callables)
            (i !=1) && print(io,", ")
            print(io,string(SPEC_TO_KW[spec])," = ")
            print(io,"x",subscripts[i])
        end
        !iszero(lsp) && print(io,", ")
    end

    if !iszero(lsp)
        for (i,spec) in enumerate(sp.specs)
            (i !=1) && print(io,", ")
            print(io,string(SPEC_TO_KW[specification(spec)])," = ")
            a = value(spec)
            units = default_units(typeof(specification(spec)))
            if is_unitful(a)
                units = unit(a)
                a = ustrip(a)
            end
            if is_real(a) 
                print(IOContext(io,:compact => true),a)
                if !(units == Unitful.NoUnits)
                    printstyled(io,'[',color = :light_black)
                    print(IOContext(io,:compact => true),units)
                    printstyled(io,']',color = :light_black)
                end
            else 
                print(io,a)
            end
        end
    end

    print(io,")")
end