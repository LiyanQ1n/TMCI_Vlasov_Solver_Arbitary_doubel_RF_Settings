"""
逐圈计算angle variable
"""
function getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ::Int64, nθ::Int64, par::ParDBRF)
    ϕJθ = zeros(nJ, nθ)
    @sync @inbounds for i in 1:nJ
        Threads.@spawn begin
            ϕJθ[i,:] .= getϕJθ(zJθ1[i,:], EJ1[i], JJ1[i], zJθ2[i,:], EJ2[i], JJ2[i], nθ, par)
        end
    end
    ϕJθ
end
function getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ::Int64, nθ::Int64, func_potential, η::Float64)
    ϕJθ = zeros(nJ, nθ)
    @sync @inbounds for i in 1:nJ
        Threads.@spawn begin
            ϕJθ[i,:] .= getϕJθ(zJθ1[i,:], EJ1[i], JJ1[i], zJθ2[i,:], EJ2[i], JJ2[i], nθ, func_potential, η)
        end
    end
    ϕJθ
end



"""
计算一圈的angle variable
"""
function getϕJθ(zθ1, E1, J1, zθ2, E2, J2, nθ::Int64, par::ParDBRF)
    F2ZE1 = getF2ZEUpHalfPlane(zθ1, E1, nθ, par)
    F2ZE2 = getF2ZEUpHalfPlane(zθ2, E2, nθ, par)

    nθcri = div(nθ+1,2)
    resvec = zeros(nθ)
    
    resvec[nθcri] = π
    resvec[nθ] = 2π
    
    # Interpolations seems faster than Dierckx
    # But Spline1D more accurate when nθ is small
    interpfun1 = linear_interpolation(zθ1[nθcri:-1:1], F2ZE1[end:-1:1], extrapolation_bc=Line()) # Spline1D(zθ1[nθcri-1:-1:2], F2ZE1[end:-1:1])
    interpfun2 = linear_interpolation(zθ2[nθcri:-1:1], F2ZE2[end:-1:1], extrapolation_bc=Line()) # Spline1D(zθ2[nθcri-1:-1:2], F2ZE2[end:-1:1])
    
    @inbounds for i in 2:(nθcri-1)
        resvec[i] = (interpfun1(zθ1[i])-interpfun2(zθ1[i]))/(J1-J2) # 计算ϕ
    end
    @inbounds for i in (nθcri+1):(nθ-1)
        resvec[i] = 2π - resvec[nθ+1-i]
    end
    
    resvec
end
function getϕJθ(zθ1, E1, J1, zθ2, E2, J2, nθ::Int64, func_potential, η::Float64)
    F2ZE1 = getF2ZEUpHalfPlane(zθ1, E1, nθ, func_potential, η)
    F2ZE2 = getF2ZEUpHalfPlane(zθ2, E2, nθ, func_potential, η)

    nθcri = div(nθ+1,2)
    resvec = zeros(nθ)
    
    resvec[nθcri] = π
    resvec[nθ] = 2π
    
    # Interpolations seems faster than Dierckx
    # But Spline1D more accurate when nθ is small
    interpfun1 = Spline1D(zθ1[nθcri:-1:1], F2ZE1[end:-1:1])
    interpfun2 = Spline1D(zθ2[nθcri:-1:1], F2ZE2[end:-1:1])
    
    # linear_interpolation(zθ1[nθcri:-1:1], F2ZE1[end:-1:1], extrapolation_bc=Line())
    # linear_interpolation(zθ2[nθcri:-1:1], F2ZE2[end:-1:1], extrapolation_bc=Line())

    @inbounds for i in 2:(nθcri-1)
        resvec[i] = (interpfun1(zθ1[i])-interpfun2(zθ1[i]))/(J1-J2) # 计算ϕ
    end
    @inbounds for i in (nθcri+1):(nθ-1)
        resvec[i] = 2π - resvec[nθ+1-i]
    end
    
    resvec
end



function getF2ZEUpHalfPlane(zθ::Vector, E, nθ, par::ParDBRF)
    nθcri = div(nθ+1,2)
    F2ZE_vec = zeros(nθcri)
    for i in 1:nθcri
        F2ZE_vec[i]=integrate_of_δ(zθ[i], zθ[end], E, par)
    end
    F2ZE_vec
end
function getF2ZEUpHalfPlane(zθ::Vector, E, nθ, func_potential, η::Float64)
    nθcri = div(nθ+1,2)
    F2ZE_vec = zeros(nθcri)
    for i in 1:nθcri
        F2ZE_vec[i]=integrate_of_δ(zθ[i], zθ[end], E, func_potential, η)
    end
    F2ZE_vec
end