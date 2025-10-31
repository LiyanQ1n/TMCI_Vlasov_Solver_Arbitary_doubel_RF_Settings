# Data analysis - RMS

"""
`xvec`: coordinate vector

`ρvec`: density vector of `ρ(x)`
"""
function getRMSWidth(xvec::Vector{Float64}, ρvec::Vector{Float64})
    ave = discrete_integrate(xvec, xvec, ρvec)
    sqrt(discrete_integrate(xvec, (xvec .-  ave).^2 .* ρvec))
end

"""
`ρ`: density function
"""
function getRMSWidth(ρ, lb::Float64, ub::Float64)
    xaverage = quadgk(ρ, lb, ub)[1]
    sqrt(quadgk(x->(x-xaverage)^2 * ρ(x), lb, ub)[1])
end