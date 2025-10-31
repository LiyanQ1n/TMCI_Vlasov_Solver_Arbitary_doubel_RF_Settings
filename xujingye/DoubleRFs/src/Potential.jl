"""
Get cavity potential term in Hamiltonian.

Parameters
---
`z`: Type `Float64`. longtitudinal position relative to synchronous particle.


Return
---
Cavity potential at `z`. Type of `Float64`.
"""
function getPotenCavity(z, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
	-v1/(2π*h1*centerenergy)*normedPotential(z, ϕs, ϕ2s, r, h1, h, circum)
end
function getPotenCavity(z, par::ParDBRF)
	getPotenCavity(z, par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
end






"""
Calculate potential term in Hamiltonian from voltage.

Infact, it should be called `getPotentialFromVoltage`.
But we don't need to do that for cavity.
"""
function getPotenWake(zvec::Vector{Float64}, wakevolvec::Vector{Float64}, Δz::Float64, centerenergy, circum)
    # il = 0
    # ir = 0
	# @inbounds for i in 1:(length(zvec)-1)
    #     if (zvec[i] <= 0) && (zvec[i+1]>=0)
    #         il = i
    #         ir = i+1
    #     end
    # end

    # Δsl = 0.5*abs(zvec[il])*wakevolvec[il]
    # Δsr = 0.5*abs(zvec[ir])*wakevolvec[ir]

    # pvec = zeros(size(zvec, 1))
    # @inbounds for i in eachindex(zvec)
    #     if i <= il
    #         s = -Δsl
    #         @inbounds @fastmath for j in i:il
    #             s -= (wakevolvec[j] + wakevolvec[j+1]) * 0.5 * Δz
    #         end
    #     elseif i >= ir
    #         s = Δsr
    #         @inbounds @fastmath for j in ir:i
    #             s += (wakevolvec[j] + wakevolvec[j+1]) * 0.5 * Δz
    #         end
    #     end
    #     pvec[i] = s
    # end
    pvec = Δz * cumsum(wakevolvec)
    -(pvec .- pvec[end]).*elementcharge/(centerenergy*circum)
end
function getPotenWake(zvec::Vector{Float64}, wakevolvec::Vector{Float64}, Δz::Float64, par::ParDBRF)
    getPotenWake(zvec, wakevolvec, Δz, par.centerenergy, par.circum)
end




@doc raw"""
Calculate wake potential by discrete method.
"""
function getPotenWake(zvec::Vector{Float64}, ρvec::Vector{Float64}, funwakez::Union{Interpolations.Extrapolation,Function}, Nb, Δz, centerenergy, circum)
    wakevol = getWakeVoltage(zvec, ρvec, funwakez, Δz, Nb)
    getPotenWake(zvec, wakevol, Δz, centerenergy, circum)
end
function getPotenWake(zvec::Vector{Float64}, ρvec::Vector{Float64}, funwakez::Union{Interpolations.Extrapolation,Function}, Nb, Δz, par::ParDBRF)
    getPotenWake(zvec, ρvec, funwakez, Nb, Δz, par.centerenergy, par.circum)
end



# get roots of `Hamiltonian(z, δ)==const`
# where hamiltonian is hamiltonian of cavity.
function guess_zero_points_of_hamiltonian(lb, ub, potentialVal, par::ParDBRF; length_of_tentative_z0=100000)
    ite_step = (ub-lb)/length_of_tentative_z0
    resvec = Float64[]
    for z in lb:ite_step:(ub-ite_step)
        if (getPotenCavity(z, par)-potentialVal) * (getPotenCavity(z+ite_step, par)-potentialVal) <= 0
            append!(resvec, z+0.5*ite_step)
        end
    end
    resvec
end

function newton_iterate_Potential_Onestep(zo::BigFloat, potentialVal, par::ParDBRF)::Float64
    zo - (getPotenCavity(zo, par)-potentialVal)/(∇Potential(zo, par)+1e-150)
end
function newton_iterate_Potential_Onestep(zo::Float64, potentialVal, par::ParDBRF)::Float64
    z0 = BigFloat(zo)
    z0 - (getPotenCavity(z0, par)-potentialVal)/(∇Potential(z0, par)+1e-150)
end


"""
获取Potential的零点。

Parameters
---
`z0`: Initial point.

`potentialVal`: Potential value
"""
function newton_iterate_Potential(z0, potentialVal, par::ParDBRF; atol=1e-12, max_try_time=5000)
    zn = z0
    zo = z0+10*atol
    try_time = 0
    while (try_time < max_try_time) && (abs(zn-zo)>atol)
        try_time += 1
        zo, zn = zn, newton_iterate_Potential_Onestep(zn, potentialVal, par)
    end
    zn
end

function get_zero_points_of_hamiltonian(potentialVal, minBucket, maxBucket, par::ParDBRF; length_of_tentative_z0=10000, atol=1e-15, max_iter_time=50000)
    zvec = guess_zero_points_of_hamiltonian(minBucket, maxBucket, potentialVal, par, length_of_tentative_z0=length_of_tentative_z0)
    for (i, z0) in enumerate(zvec)
        zvec[i] = newton_iterate_Potential(z0, potentialVal, par, atol=atol, max_try_time=max_iter_time)
    end
    sort(remove_repeated_and_outranged_data(zvec, minBucket, maxBucket, atol=atol))
end