print_spec(t) = print(t) #fallback

print_spec(io::IO,x::Type{MOLAR}) = print(io,"Molar ")
print_spec(io::IO,x::Type{MASS}) = print(io,"Mass ")
print_spec(io::IO,x::Type{TOTAL}) = print(io,"Total ")
print_spec(io::IO,x::Type{OneMol}) = print(io,"One mol ")
print_spec(io::IO,x::Type{SingleComponent}) = print(io,"Single component ")



print_spec(io::IO,::Type{Pressure}) = print(io,"Pressure")
print_spec(io::IO,::Type{Temperature}) = print(io,"Temperature")


function print_spec(io::IO,::Type{Gibbs{T2}}) where T2
    print_spec(io,T2)
    print(io,"Gibbs energy")
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

is_real(x) = false
is_real(x::Real) = true
is_real(x::Bool) = false
is_real(x::AbstractVector{T} where T<:Real) = true 

function Base.show(io::IO, sp::Spec{T}) where T 
    compact = get(io, :compact, false)
    if compact
        print(io,"spec(",string(SPEC_TO_KW[specification(sp)])," = ")
        
        a = value(sp)
        if is_real(a) 
            print(IOContext(io,:compact => true),a," ",default_units(T))
        else 
            print(io,a)
        end
        print(io,")")
    else
        print_spec(io,T)
        print(io," : ")
        a = value(sp)
        if is_real(a) 
            print(IOContext(io,:compact => true),a," ",default_units(T))
        else 
            print(io,a)
        end
    end
end

function Base.show(io::IO, sp::Spec{Options})
    print_spec(io,Options)
    print(io," : ")
    show(io,value(sp))
end

#==function Base.show(::IO, sp::T) where T <:AbstractSpec
    print(io,"::" )
    print_spec(io,T)
end
==#

function Base.show(io::IO, sp::ThermodynamicState)
    len1 = length(sp.specs)
    if len1 > 0
        if len1 == 1
            p1 = " Constant Property Specification:"
        else
            p1 = " Constant Property Specifications:"
        end
        println(io,len1,p1)
        
        for (i,spec) in enumerate(sp.specs)
            if i !=1
                println()
            end
            print(io," ")
            show(io,spec)
        end 
        
    end
    
    len2 = length(sp.callables)
    if len2 > 0
        if len1 > 0
        println()
        println()
        end
    if len2 == 1
        p2 = " Variable Property Specification:"
    else
        p2 = " Variable Property Specifications:"
    end
    
    println(len2,p2)
    
    for (i,spec) in enumerate(sp.callables)
        if i !=1
            println()
        end
        print(io," ")
        print_spec(typeof(spec))
    end 
end

end

