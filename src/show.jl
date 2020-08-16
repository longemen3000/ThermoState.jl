print_spec(t) = print(t) #fallback

print_spec(x::Type{MOLAR}) = print("Molar ")
print_spec(x::Type{MASS}) = print("Mass ")
print_spec(x::Type{TOTAL}) = print("Total ")
print_spec(x::Type{OneMol}) = print("One mol ")
print_spec(x::Type{SingleComponent}) = print("Single component ")



print_spec(::Type{Pressure}) = print("Pressure")
print_spec(::Type{Temperature}) = print("Temperature")


function print_spec(::Type{Gibbs{T2}}) where T2
    print_spec(T2)
    print("Gibbs energy")
end

function print_spec(::Type{Helmholtz{T2}}) where T2
    print_spec(T2)
    print("Helmholtz energy")
end

function print_spec(::Type{Enthalpy{T2}}) where T2
    print_spec(T2)
    print("enthalpy")
end

function print_spec(::Type{Entropy{T2}}) where T2
    print_spec(T2)
    print("entropy")
end

function print_spec(::Type{InternalEnergy{T2}}) where T2
    print_spec(T2)
    print("internal energy")
end


function print_spec(::Type{MaterialAmount{MOLAR}})
    print("Moles")
end

function print_spec(::Type{MaterialAmount{MASS}})
    print("Mass")
end

function print_spec(::Type{MaterialCompounds{MASS,FRACTION}})
    print("Mass fraction")
end

function print_spec(::Type{PhaseTag})
    print("Phase specification")
end

function print_spec(::Type{MaterialCompounds{MOLAR,FRACTION}})
    print("Molar fraction")
end

function print_spec(::Type{MaterialCompounds{MASS,TOTAL_AMOUNT}})
    print("Mass amounts")
end

function print_spec(::Type{MaterialCompounds{MOLAR,TOTAL_AMOUNT}})
    print("Molar amounts")
end


function print_spec(::Type{VolumeAmount{T,VOLUME}}) where T
    print_spec(T)
    print("volume")
end

function print_spec(::Type{VolumeAmount{T,DENSITY}}) where T
    print_spec(T)
    print("density")
end

function print_spec(::Type{TwoPhaseEquilibrium})
    print("Two phase equilibrium")
end

function print_spec(::Type{Options})
    print("Options")
end

is_real(x) = false
is_real(x::Real) = true
is_real(x::Bool) = false
is_real(x::AbstractVector{T} where T<:Real) = true 

function Base.show(::IO, sp::Spec{T}) where T 
    print_spec(T)
    print(" : ")
    a = value(sp)
    if is_real(a) 
        print(a," ",default_units(T))
    else 
        print(a)
    end
end

function Base.show(io::IO, sp::Spec{Options})
    print_spec(Options)
    print(" : ")
    show(io,value(sp))

end

function Base.show(::IO, sp::T) where T <:AbstractSpec
    print("::" )
    print_spec(T)
end






function Base.show(io::IO, sp::Specs)
    println(length(sp.specs)," Property Specifications:")
    for (i,spec) in enumerate(sp.specs)
        if i !=1
            println()
        end
        print(" ")
        show(io,spec)
    end 
    
end