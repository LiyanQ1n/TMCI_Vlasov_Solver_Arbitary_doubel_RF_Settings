# Simulation Tools: mapping equation

function _cache_one_turn_map(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    k1 = 2 * h1 * pi / circum
    k2 = 2 * h1 * pi / circum * h
    c1 = sin(ϕs) + r * sin(ϕ2s)
    c2 = v1 / centerenergy
    c3 = 1 / (centerenergy / restenergy)^2
    (; c1, c2, c3, k1, k2)
end
function one_turn_map_with_cache(z, δ, ϕs, ϕ2s, r, circum, αc, k1, k2, c1, c2, c3)
    δ = δ + c2 * normed_voltage_with_cache(z, ϕs, ϕ2s, r, k1, k2, c1)
    η = αc - c3 / (1 + δ)^2
    z = z - η * δ * circum
    z, δ
end
function one_turn_map_with_cache(z, δ, ϕs, ϕ2s, r, circum, αc, c)
    one_turn_map_with_cache(z, δ, ϕs, ϕ2s, r, circum, αc, c.k1, c.k2, c.c1, c.c2, c.c3)
end



"""
mapping equation without wake. Here `δ=δ(t), z=z(t)`. Not `δ=δ(s), z=z(s)`.
"""
function one_turn_map(z::Real, δ::Real, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    one_turn_map_with_cache(z, δ, par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
end
function one_turn_map!(zvec::Vector{Real}, δvec::Vector{Real}, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    n = length(zvec)
    if n != length(δvec)
        error("Error! zvec and δvec should be same length!")
    end
    for i = 1:n
        zvec[i], δvec[i] = one_turn_map_with_cache(zvec[i], δvec[i], par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
    end
    zvec, δvec
end
function one_turn_map(zvec::Vector{Real}, δvec::Vector{Real}, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    n = length(zvec)
    if n != length(δvec)
        error("Error! zvec and δvec should be same length!")
    end
    reszvec = zeros(n)
    resδvec = zeros(n)
    for i = 1:n
        reszvec[i], resδvec[i] = one_turn_map_with_cache(zvec[i], δvec[i], par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
    end
    reszvec, resδvec
end



"""
给定粒子的z,δ坐标矢量，计算one turn map的矩阵。

每行为每代的全体粒子。每列都是同一粒子演化路径。

Parameters
---
`zvec`:

`δvec`:

`max_iterate_time`:

RETURN
---
全体粒子z和δ的历史矩阵.

"""
function one_turn_map_history(zvec::Vector{Float64}, δvec::Vector{Float64}, max_iterate_time, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    zmat = zeros(max_iterate_time, size(zvec, 1))
    δmat = zeros(max_iterate_time, size(δvec, 1))
    zmat[1, :] = zvec
    δmat[1, :] = δvec
    @sync for colindex in eachindex(zvec) # each col represent a particle
        Threads.@spawn begin
            z = zmat[1, colindex]
            δ = δmat[1, colindex]
            for rowindex in 1:(max_iterate_time-1)
                z, δ = one_turn_map_with_cache(z, δ, par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
                zmat[rowindex+1, colindex] = z
                δmat[rowindex+1, colindex] = δ
            end
        end
    end
    zmat, δmat
end



"""
给定粒子的z,δ坐标矢量，迭代计算one turn map的粒子z,δ坐标矢量。

Parameters
---
`zvec`:

`δvec`:

`max_iterate_time`:

RETURN
---
全体粒子`zvec`和`δvec`的.
"""
function one_turn_map_repeat(zvec::Vector{Float64}, δvec::Vector{Float64}, max_iterate_time, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    reszvec = zvec[:]
    resδvec = δvec[:]
    @sync for idx in eachindex(zvec)
        Threads.@spawn begin
            for iter_time in 1:max_iterate_time
                reszvec[idx], resδvec[idx] = one_turn_map_with_cache(reszvec[idx], resδvec[idx], par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
            end
        end
    end
    reszvec, resδvec
end
function one_turn_map_repeat!(zvec::Vector{Float64}, δvec::Vector{Float64}, max_iterate_time, par::ParDBRF)
    c = _cache_one_turn_map(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy)
    @sync for idx in eachindex(zvec)
        Threads.@spawn begin
            for iter_time in 1:max_iterate_time
                zvec[idx], δvec[idx] = one_turn_map_with_cache(zvec[idx], δvec[idx], par.ϕs, par.ϕ2s, par.r, par.circum, par.αc, c)
            end
        end
    end
    zvec, δvec
end