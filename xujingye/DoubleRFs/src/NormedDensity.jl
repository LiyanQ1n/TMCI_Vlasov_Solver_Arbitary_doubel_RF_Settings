"""
从倒扣的Hamiltonian的Potential Term中，得到密度分布。

Parameters
---
`zvec`: `z`的采样点，可以不均匀。

`potenvec`: 和`zvec`对应的Potential Term in Hamiltonian。

`zL`， `zR`: 密度分布的左边界和右边界。 一般来说，不需要管，但是我似乎遇到过zL,zR以外tmpexp影响特别大的情况。

---

Return
---
与`zvec`对应的密度分布。
"""
function getDistFromPotential(zvec::Vector, potenvec::Vector, zL, zR, ηp, σδ::Float64)
    potenvec = potenvec .- maximum(potenvec)
    tmpexp = exp.(potenvec/(ηp * σδ^2))
    changeOuterZero!(tmpexp, zvec, zL, zR)  # 将超出边界的部分，设置为0，防止参与求和在某些条件下产生影响；还有使最后一步自然而然将外围清零
    tmpexp./((2*sum(tmpexp)-tmpexp[1]-tmpexp[end])*(zvec[2]-zvec[1])/2)
end
function getDistFromPotential(zvec::Vector, potenvec::Vector, zL, zR, par::ParDBRF)
    getDistFromPotential(zvec, potenvec, zL, zR, ηp, par.σδ)
end



"""
Get double rf density without wake.
"""
function getNormedDensity(ϕs, ϕ2s, v1, r, σδ, αc, h1, h, circum, centerenergy)
    ηp = αc - (restenergy^2)/(centerenergy^2)
    zL, zR = getBucketBounds(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc)
    println("BucketBounds: ", zL, " ", zR)
    zvec = collect(range(zL, zR, length=100001))
    potenvec = getPotenCavity.(zvec, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    ρvec = getDistFromPotential(zvec, potenvec, zL, zR, ηp, σδ)
    zvec, ρvec
end
function getNormedDensity(par::ParDBRF)
    getNormedDensity(par.ϕs, par.ϕ2s, par.v1, par.r, par.σδ, par.αc, par.h1, par.h, par.circum, par.centerenergy)
end
function getNormedDensity(n::Int64, ϕs, ϕ2s, v1, r, σδ, αc, h1, h, circum, centerenergy)
    ηp = αc - (restenergy^2)/(centerenergy^2)
    zL, zR = getBucketBounds(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc)
    println("BucketBounds: ", zL, " ", zR)
    zvec = collect(range(zL, zR, length=n))
    potenvec = getPotenCavity.(zvec, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    ρvec = getDistFromPotential(zvec, potenvec, zL, zR, ηp, σδ)
    zvec, ρvec
end
function getNormedDensity(n::Int64, par::ParDBRF)
    getNormedDensity(n, par.ϕs, par.ϕ2s, par.v1, par.r, par.σδ, par.αc, par.h1, par.h, par.circum, par.centerenergy)
end





"""
收敛判据

Parameters
---
ρnow: 当前得到的单束团密度分布

ρpre: 先前密度分布

zvec: 上面密度采样点的z坐标矢量。因为需要对总偏差进行积分。

χ: 权重参数。

Return:
---
收敛度。其值越接近0越好。
"""
function convergecriterion(ρnow, ρpre, zvec, χ, zL, zR)
	# fun = LinearInterpolation(zvec, abs.(ρnow .- ρpre))
	# res = quadgk(fun, zL, zR)[1]
    res = discrete_integrate(zvec, abs.(ρnow .- ρpre))
	res/χ
end



"""
计算束团Bucket右侧(z>0)的边界。

由于尾场是从前往后传递，因此右侧potential不会受到wake potential的影响。只需要计算一次就可以重复利用。

Return
---
由右侧边界的坐标zR，和临界potential构成的Tuple.
"""
function getRBound(zvec::Vector{Float64}, potenVecCav::Vector{Float64},ϕs,ϕ2s,v1,r,h1,h,circum,centerenergy,σz)
    @inbounds for i in 2:(length(zvec)-1)
        if (zvec[i] > σz) && (potenVecCav[i+1] > potenVecCav[i]) && (potenVecCav[i-1]>potenVecCav[i])
            # 此时基本确定就在这zvec[i-1]和zvec[i+1]之间，感觉差不多就用zvec[i]了
            return zvec[i], getPotenCavity(zvec[i], ϕs,ϕ2s,v1,r,h1,h,circum,centerenergy)
        end
    end
    error("没有找到头部决定粒子存在的potential边界!")
end
function getRBound(zvec::Vector{Float64}, potenVecCav::Vector{Float64}, σz::Float64, par::ParDBRF)
    getRBound(zvec, potenVecCav, par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy, σz)
end




"""
计算束团Bucket左侧的边界。

由于需要考虑wake potential的影响，因此每次更新密度分布，都需要重新计算左边界（也就是尾部）。

"""
function getLBound(zvec::Vector{Float64}, potenVecTotal::Vector{Float64}, criPoten::Float64, σz::Float64)
    @inbounds for i in 2:(length(zvec))
        if (potenVecTotal[i]>=criPoten) && (potenVecTotal[i-1]<=criPoten) && (zvec[i]<σz)
            return 0.5*(zvec[i]+zvec[i-1])
        end
    end
    error("没有找到尾部决定粒子存在的potential边界!")
end


function getρVecNow(ρvec::Vector, zvec::Vector, potenVecCavity::Vector, funwakez::Union{Interpolations.Extrapolation, Function}, zR, criPotential, Nb, centerenergy, circum, η, σδ, σz, Δz)
    potenVecWake=getPotenWake(zvec, ρvec, funwakez, Nb, Δz, centerenergy, circum)
    potenVecTotal = potenVecCavity .+ potenVecWake
    zL = getLBound(zvec, potenVecTotal, criPotential, σz)
    getDistFromPotential(zvec, potenVecTotal, zL, zR, η, σδ)
end
function getρVecNow(ρvec::Vector, zvec::Vector, potenVecCavity::Vector, funwakez::Union{Interpolations.Extrapolation, Function}, zR, criPotential, Nb, Δz, σz, par::ParDBRF)
    getρVecNow(ρvec, zvec, potenVecCavity, funwakez, zR, criPotential, Nb, par.centerenergy, par.circum, par.η, par.σδ, σz, Δz)
end





"""
Get double rf density with wake.

Parameters
---
`v1`: voltage of primary cavity.

`v2`: voltage of harmonic cavity.

`ϕs`: synchronous phase of primary cavity.

`ϕ2s`: synchronous phase of harmonic cavity.

`Nb`: electron number

`σz`: 

χ: 需要根据经验确定初始值。程序提供了密度迭代的情况。可以根据收敛速度调整。默认为0.5。

Σ: 足够宽，使得即使在尾场下改变分布，此处密度也为0。默认为0.5m。 需要包含potential右侧的临界值，但是不包括左侧临界值。
"""
function getNormedDensity(v1, v2, ϕs, ϕ2s, Nb, σz, ηp, σδ, funwakez::Union{Interpolations.Extrapolation, Function}, circum, centerenergy, h1, h; χ::Float64=0.7, Σ=0.7)
    Δz=0.02σz
    r = v2/v1
    
    # 全程不变的量
    zvec = collect(-Σ:Δz:Σ)
    potenVecCavity = getPotenCavity.(zvec, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    zR, criPotential = getRBound(zvec, potenVecCavity, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, σz) # 头部粒子存在与否的边界，及该点临界potential
    
    zL = getLBound(zvec, potenVecCavity, criPotential, σz)
    #    由于右侧(头部)不会受到尾场影响，因此只需要求一次右侧z和临界的Hamiltonian即可——但是wake potential需要头部为0。
    #    左侧（尾部）则为了保险起见每次更新。


    ρvec=getDistFromPotential(zvec, potenVecCavity, zL, zR, ηp, σδ) # 计算初始密度分布

    i = 1
    # @save "zvec_ρvec_"*string(i)*".bson" zvec ρvec

    # 第一次迭代
    ρvecpre = ρvec
    ρvecnow = getρVecNow(ρvec, zvec, potenVecCavity, funwakez, zR, criPotential, Nb, centerenergy, circum, ηp, σδ, σz, Δz)


    ρvec = (1-χ).*ρvecpre .+ χ.*ρvecnow   # 这是利用默认权重χ， 直接迭代的分布， 下一行则是相应的初始判据
    cri = convergecriterion(ρvec, ρvecpre, zvec, χ, zL, zR)
    i = 2

    # 开始循环迭代
    while χ >= 0.001
        println(χ)
        ρvecpre = ρvec
        ρvecnow = getρVecNow(ρvec, zvec, potenVecCavity, funwakez, zR, criPotential, Nb, centerenergy, circum, ηp, σδ, σz, Δz)

        # 找到合适的χ，使分布收敛
        suc = false
        while suc == false
            i += 1
            ρvec = (1-χ).*ρvecpre .+ χ.*ρvecnow
            # if i%100 == 0
            #     @save "zvec_ρvec_"*string(i)*".bson" zvec ρvec
            # end
            tmpcri = convergecriterion(ρvec, ρvecpre,  zvec, χ, zL, zR)
            println("Iteration: $i . Criterion: $tmpcri")
            if tmpcri < cri     # 找到的χ能否让判据收敛
                suc = true
                cri = tmpcri
            elseif χ <= 0.001   # 如果χ过小，则不用再寻找合适的χ了
                break
            else
                χ = χ*0.9
                println(χ)
            end
        end
    end
    zvec, ρvec
end
function getNormedDensity(Nb, σz, funwakez::Union{Interpolations.Extrapolation, Function}, par::ParDBRF; χ::Float64=0.5, Σ=0.5)
    getNormedDensity(par.v1, par.v2, par.ϕs, par.ϕ2s, Nb, σz, par.η, par.σδ, funwakez, par.circum, par.centerenergy, par.h1, par.h, χ=χ, Σ=Σ)
end
