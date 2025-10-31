"""
对J离散化方法
"""
function calM1WithoutNb(ℱj1j2l1l2::Array{ComplexF64, 4}, ΔJJ::Vector, ψJ::Vector, nJ::Int, lm::Int, par::ParDBRF)::Matrix
    matM1 = zeros(ComplexF64, (2lm+1)*nJ, (2lm+1)*nJ)
    @inbounds for l1=-lm:lm, l2=-lm:lm
        @inbounds for j1=1:nJ
            rowidx = (lm+l1)*nJ+j1
            @inbounds for j2=1:nJ
                colidx = (lm+l2)*nJ+j2
                # 逐元素赋值
                matM1[rowidx, colidx] = ψJ[j1]*ℱj1j2l1l2[j1,j2,l1+lm+1,l2+lm+1]*ΔJJ[j2]
            end
        end
    end
    -im*clight*r0*par.ω0/(4π*par.γ*par.ωβ) .* matM1
end

"""
对J离散化方法
"""
function calM2(νsJ::Vector, nJ, lm, par::ParDBRF)::Matrix{Float64}
    matM2 = zeros((2lm+1)*nJ, (2lm+1)*nJ)
    @inbounds for l1=-lm:lm, l2=-lm:lm
        @inbounds for j1=1:nJ, j2=1:nJ
            rowidx = (lm+l1)*nJ+j1
            colidx = (lm+l2)*nJ+j2
            # 逐元素赋值
            matM2[rowidx, colidx] = ((l1 != l2)||(j1 != j2)) ? 0.0 : (l1 * νsJ[j2] * par.ω0)
        end
    end
    matM2
end

function calM2(νs1, νs2, nJ, lm, Δϕ, par::ParDBRF)::Matrix{Float64}
    matM2 = zeros((2lm+1)*nJ, (2lm+1)*nJ)
    @inbounds for l1=-lm:lm, l2=-lm:lm
        @inbounds for j1=1:nJ, j2=1:nJ
            rowidx = (lm+l1)*nJ+j1
            colidx = (lm+l2)*nJ+j2
            # 逐元素赋值
            matM2[rowidx, colidx] = ((l1 != l2)||(j1 != j2)) ? 0.0 : (l1 * νsJ[j2] * par.ω0)
        end
    end
    matM2
end

"""
`Δϕ`为
"""
function h1(Δϕ, l1, l2)
    (l1 == l2) ? Δϕ/π : (sin(l1-l2)*Δϕ)/(l1-l2)/π
end

function h2(Δϕ, l1, l2)
    (l1 == l2) ? (π-Δϕ)/π : -(sin(l1-l2)*Δϕ)/(l1-l2)/π
end