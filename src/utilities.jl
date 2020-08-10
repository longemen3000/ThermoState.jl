# mixing rule:
#= =
an operation on a vector of properties that returns a single sumber.
is a weighted operation on p, by the weights x.
# 

= =#
"""
    mixing_rule(op, x, p)

returns an efficient implementation of:
`sum(x[i]*x[j]*op(p[i],p[j]) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])`

""" 
function mixing_rule(op, x, p)
    N = length(x)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            res1 += p[i] * x[i]^2
            for j = 1:i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j])
            end
        end
    end
    return res1
end
# 
# Abstract mixing rule, generates a mixing rule,based on 
# an operation, so mij = xi*xj*op(pi,pj)*Aij
# example: mixing_rule(geometric_mean_rule,x,Tc,1.-K)
# 
"""
    mixing_rule(op, x, p,A)

returns an efficient implementation of:
`sum(A[i,j]*x[i]*x[j]*op(p[i],p[j]) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])`

""" 
function mixing_rule(op, x, p, A)
    N = length(x)
    checkbounds(A, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            res1 += p[i] * x[i]^2
            for j = 1:i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j]) * A[i, j]
            end
        end
    end
    return res1
end
# 
# Abstract asymetric mixing rule 
# Adds a Asymetric matrix A_asym, and a op_sim(xi,xj,Aij)
# the mayor example is the GERG2008 equation, where
# op_asym(xi,xj,Aij) = (xi+xj)/(Aij^2*xi + xj)
# 
"""
    mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)

returns an efficient implementation of:
` sum(A[i,j] * x[i] * x[j] * op(p[i],p[j]) * op_asym(x[i],x[j],A_asym[i,j])) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])` , op_asym doesn't follow this symmetry.

""" 
function mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)
    N = length(x)
    checkbounds(A, N, N)
    checkbounds(A_asym, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))

        for i = 1:N
            x[i] != 0 && begin
                res1 += p[i] * x[i]^2
                for j = 1:i - 1
                    res1 +=
                        2 *
                        x[i] *
                        x[j] *
                        op(p[i], p[j]) *
                        A[i, j] *
                        op_asym(x[i], x[j], A_asym[i, j])
                end
            end
        end
    end
    return res1
end

geometric_mean_rule(a, b) = sqrt(a * b)

arithmetic_mean_rule(a, b) = (a + b) / 2

harmonic_mean_rule(a, b) = 2 * a * b / (a + b)

cubic_mean_rule(a, b) = ((a^(1 / 3) + b^(1 / 3)) / 2)^3

_power_mean_rule(a, b, n) = ((a^(1 / n) + b^(1 / n)) / 2)^n

power_mean_rule(n) = (a, b) -> _power_mean_rule(a, b, n)


"""
    mixing_matrix!(A, op, p)

returns a matrix of size nxn , where A[i,j] = op(p[i],p[j]), modifying A inplace

""" 
function mixing_matrix!(A, op, p)
    N = length(p)
    @boundscheck size(A) == (N, N)
    @inbounds begin
        res1 = zero(eltype(p))
        for i = 1:N
            A[i, i] = p[i]
            for j = 1:i - 1
                A[i, j] = op(p[i], p[j])
                A[j, i] = op(p[i], p[j])
            end
        end
    end
    return A

end







"""
    mixing_matrix!(op, p)

returns a matrix of size nxn , where A[i,j] = op(p[i],p[j])
""" 
function mixing_matrix(op, p)
    N = length(p)
    A = Matrix{eltype(p)}(undef, N, N)
    return mixing_matrix!(A, op, p)
end

#a_in_b(t1, t2) = all(in(t2), t1)
#tuple_comparison(t1, t2) = ((length(t1) == length(t2)) && a_in_b(t1, t2) && a_in_b(t2, t1))
