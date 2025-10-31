@doc raw"""
```text
getCavityVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)
```

双高频腔压。

`z`: `Float64`. 相对于参考粒子的纵向位置。``z>0``表示头部。单位[m]。

`v1`: `Float64`. 主腔腔压。单位[V]。

`r`: `Float64`. 高频腔和主腔电压之比。

`ϕs`: `Float64`。主腔相位。单位[rad]。

`ϕ2s`: `Float64`。高频腔相位。单位[rad]。

`h1`: `Int64`。主腔谐波数，为主腔频率``\omega_\text{Primary}``与回旋频率``\omega_0``的比值``h_1 = \omega_\text{Primary}/\omega_0``。

`h`: `Int64`. 高频腔和主腔谐波数之比。

`circum`: `Float64`. 储存环周长。单位[m]。
"""
function getCavityVoltage(z::Float64, v1::Float64, r::Float64, ϕs::Float64, ϕ2s::Float64, h1::Int, h::Int, circum::Float64)
    v1*normedVoltage(z, ϕs, ϕ2s, r, h1, h, circum)
end

@doc raw"""
```text
getCavityVoltage(z, par::ParDBRF)
```

双高频腔压。

`z`: `Float64`. 相对于参考粒子的纵向位置。``z>0``表示头部。单位[m]。

`par`: `ParDBRF`. 打包的参数结构体，不过`ParDBRF`中打包的参数比这里需要的参数多，详情见`ParDBRF`。如果没有那么多参数的话，还是用另一个版本吧。
"""
function getCavityVoltage(z, par::ParDBRF)
    getCavityVoltage(z, par.v1, par.r, par.ϕs, par.ϕ2s, par.h1, par.h, par.circum)
end




@doc raw"""
Calculate wake voltage by wake field.

Infact, wake voltage here is often called wake potential. But it has the dimension of voltage `V`.

```math
- N_b e \int_z^\infty dz' ρ(z') W_1'(z-z')
```

where ``W'_1(z)`` is the longtitudinal wake function.

---

Parameters
---
`zvec`: `Vector{Float64}`. Longtitudinal position relative to synchronous particle.

`ρvec`: `Vector{Float64}`. Normalized density at `z`.

`funz2wakez`: Wake function at `z`.

`Δz`: Step of `zvec`.

`Nb`: number of electrons.

---

Return
---
Type of `Vector`. Wake voltage according to `zvec`.
"""
function getWakeVoltage(zvec::Vector{Float64}, ρvec::Vector{Float64}, funz2wakez::Union{Interpolations.Extrapolation,Function}, Δz, Nb)
	volvec = zeros(size(zvec, 1))
	@sync for i in 1:length(zvec)
		Threads.@spawn begin
			s = 0.0
			@fastmath for j in i:length(zvec)
				s += funz2wakez((i-j)*Δz) * ρvec[j] * Δz
			end
			volvec[i] = s
		end
	end
	-Nb*elementcharge*volvec
end






"""
Calculate wake voltage by impedance at positive half frequency domain.

The impedance at negative half frequency domain should be folded to the positive half.

It has problem here(see "Test"). And the speed 
"""
function getWakeVoltage(zvec::Vector{Float64}, ρvec::Vector{Float64}, ωvec::Vector{Float64}, reimpedancevec::Vector{Float64}, imimpedancevec::Vector{Float64}, Δω, Δz, Nb)
	volvec = zeros(size(zvec, 1))
	@sync for i in 1:(length(zvec)-1)
		Threads.@spawn begin
			s = 0.0
			@fastmath @inbounds for ωindex in 1:length(reimpedancevec)
				ω = ωvec[ωindex]
				tmp1 = 0.0
				tmp2 = 0.0
				@fastmath @inbounds for j in (i+1):length(zvec)
					ϕ = ω*(i-j)*Δz/clight
					tmp1 += ρvec[j]*cos(ϕ)
					tmp2 += ρvec[j]*sin(ϕ)
				end
				s += reimpedancevec[ωindex]*tmp1 - imimpedancevec[ωindex]*tmp2
			end
			volvec[i] = s
		end
	end
	-Nb*elementcharge*Δω*Δz/π/2 .* volvec	# Why /2? Because 我们拿到的正半频域的阻抗，通常是将负半频域的幅值折叠过去的。
end




function getTotalVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum, funWakeVoltage::Union{Interpolations.Extrapolation,Function})
    getCavityVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)+funWakeVoltage(z)
end
function getTotalVoltage(z, par::ParDBRF, funWakeVoltage::Union{Interpolations.Extrapolation,Function})
	getTotalVoltage(z, par.v1, par.r, par.ϕs, par.ϕ2s, par.h1, par.h, par.circum, funWakeVoltage)
end





function getTotalVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)
    getCavityVoltage(z, v1, r, ϕs, ϕ2s, h1, h, circum)
end
function getTotalVoltage(z, par::ParDBRF)
    getTotalVoltage(z, par.v1, par.r, par.ϕs, par.ϕ2s, par.h1, par.h, par.circum)
end