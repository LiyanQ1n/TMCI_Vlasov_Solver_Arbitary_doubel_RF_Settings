# """
# Similar to Bessel function in Alex. Chao, 1993, Eq.(6.75).
# Actually in single RF condition, this function equals `i^{-l} J_l(ωr/c)`
# """
# function funS(lvec::Vector{Int64}, ωvec::Vector, Jvec::Vector, zJϕ::Matrix, ϕJϕ::Matrix, par::ParDBRF)
#     llen = length(lvec)
#     ωlen = length(ωvec)
#     Jlen = length(Jvec)
#     m = zeros(ComplexF64, llen, ωlen, Jlen)
#     @inbounds for lidx=1:llen
#         @sync @inbounds for Jidx=1:Jlen
#             Threads.@spawn begin
#                 @inbounds for ωidx=1:ωlen
#                     tmpϕvec = ϕJϕ[Jidx, :]
#                     tmpres = exp.(1.0im.*(lvec[lidx].*tmpϕvec .- zJϕ[Jidx,:]./clight.*(ωvec[ωidx] .- par.ωξ)))
#                     m[lidx, ωidx, Jidx] = discrete_integrate(tmpϕvec, tmpres)/2π
#                 end
#             end
#         end
#     end
#     m
# end

"""
Return S(ω, J, l)
"""
function funS(lvec, ωvec, zθJ, ϕθJ, ωξ, nω, nJ, lm)
    resmat = zeros(ComplexF64, nω, nJ, 2*lm+1)
    for lidx = 1:2*lm+1, Jidx = 1:nJ, ωidx = 1:nω
        tmpres = exp.(1.0im .* (lvec[lidx] .* view(ϕθJ, :, Jidx) .- 1 / clight * (ωvec[ωidx] - ωξ) .* view(zθJ, :, Jidx)))
        resmat[ωidx, Jidx, lidx] = discrete_integrate(view(ϕθJ, :, Jidx), tmpres) / 2π
    end
    resmat
end


function funS1(lvec, ωvec, zθJ, ϕθJ, ωξ, nω, nJ, lm)
    resmat = zeros(ComplexF64, nω, nJ, 2 * lm + 1)
    for lidx = 1:2*lm+1, Jidx = 1:nJ, ωidx = 1:nω
        tmpres = exp.(1.0im .* (lvec[lidx] .* ϕθJ[:, Jidx] .- 1 / clight * (ωvec[ωidx] - ωξ) .* zθJ[:, Jidx]))
        resmat[ωidx, Jidx, lidx] = discrete_integrate(ϕθJ[:, Jidx], tmpres) / 2π
    end
    resmat
end

function funS2(lvec, ωvec, zθJ, ϕθJ, ωξ, nω, nJ, lm)
    resmat = zeros(ComplexF64, nω, nJ, 2 * lm + 1)
    @inbounds for lidx = 1:2*lm+1, Jidx = 1:nJ, ωidx = 1:nω
        tmpres = exp.(1.0im .* (lvec[lidx] .* view(ϕθJ, :, Jidx) .- 1 / clight * (ωvec[ωidx] - ωξ) .* view(zθJ, :, Jidx)))
        resmat[ωidx, Jidx, lidx] = discrete_integrate(view(ϕθJ, :, Jidx), tmpres) / 2π
    end
    resmat
end

function funS3(lvec, ωvec, zθJ, ϕθJ, ωξ, nω, nJ, lm)
    resmat = zeros(ComplexF64, nω, nJ, 2 * lm + 1)
    @sync @inbounds for lidx = 1:2*lm+1
        Threads.@spawn @inbounds for Jidx = 1:nJ, ωidx = 1:nω
            tmpres = exp.(1.0im .* (lvec[lidx] .* view(ϕθJ, :, Jidx) .- 1 / clight * (ωvec[ωidx] - ωξ) .* view(zθJ, :, Jidx)))
            resmat[ωidx, Jidx, lidx] = calculate_integral(view(ϕθJ, :, Jidx), tmpres) / 2π
        end
    end
    resmat
end