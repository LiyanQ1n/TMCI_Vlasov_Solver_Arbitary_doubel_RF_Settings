# The Hamiltonian is like an inverted bowl.
# Particles are concentrated on the top of the inverted bowl.

@doc raw"""

    Hamiltonian(z::Float64, δ::Float64, v1, r::Float64, ϕs::Float64, ϕ2s::Float64, η::Float64, h1::Int64, h::Int64, circum::Float64, centerenergy::Float64)

Get Hamiltonian with harmonic cavity but without wake.

Hamiltonian here has the form: ``H=H(z,δ;s)``.

Parameters
---
`z`: ``z=s-ct``. ``z>0`` is the head of the bunch.

`δ`:

`v1`:

`r`:

`ϕs`:

`ϕ2s`:

`ηp`:

`h1`: Defaut value is `756`, parmeter of HEPS.

`h`: Defaut value is `3`, parmeter of HEPS.

`centerenergy`: Defaut value is `6e9`, parmeter of HEPS.

"""
function hamiltonian(z::Float64, δ::Float64, v1, r::Float64, ϕs::Float64, ϕ2s::Float64, η::Float64, h1::Int64, h::Int64, circum::Float64, centerenergy::Float64)
    -0.5 * η * δ^2 + getPotenCavity(z, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
end
"""
    hamiltonian(z::Float64, δ::Float64, par::ParDBRF)

哈密顿量
"""
function hamiltonian(z::Float64, δ::Float64, par::ParDBRF)
    hamiltonian(z, δ, par.v1, par.r, par.ϕs, par.ϕ2s, par.η, par.h1, par.h, par.circum, par.centerenergy)
end





@doc raw"""
Get Hamiltonian with harmonic cavity but without wake.

Hamiltonian here has the form: ``H=H(z,δ;s)``.

Parameters
---
`z``: ``z=s-ct``. ``z>0`` is the head of the bunch.

`δ`:

`v1`:

`r`:

`ϕs`:

`ϕ2s`:

`funWakePoten`: Wake potential function. Type of `Interpolations.Extrapolation` or `Function`.

`ηp`:

`h1`: Defaut value is `756`, parmeter of HEPS.

`h`: Defaut value is `3`, parmeter of HEPS.

`centerenergy`: Defaut value is `6e9`, parmeter of HEPS.

"""
function hamiltonian(z::Float64, δ::Float64, v1, r::Float64, ϕs::Float64, ϕ2s::Float64, funWakePoten::Union{Interpolations.Extrapolation,Function}, η, h1, h, circum, centerenergy)
    hamiltonian(z, δ, v1, r, ϕs, ϕ2s, η, h1, h, circum, centerenergy)+funWakePoten(z)
end
function hamiltonian(z::Float64, δ::Float64, funWakePoten::Union{Interpolations.Extrapolation,Function}, par::ParDBRF)
    hamiltonian(z, δ, par.v1, par.r, par.ϕs, par.ϕ2s, funWakePoten, par.η, par.h1, par.h, par.circum, par.centerenergy)
end

