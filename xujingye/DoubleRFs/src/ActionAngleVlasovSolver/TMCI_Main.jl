# Solve Vlasov equation by action-angle variable
# 


@doc raw"""
    scanNb(ratiovec, Nb, JJ, ΔJJ, ψJ, νsJ, zJθ, ΦJθ, ωvec, imp_trans, lm, par::ParDBRF)

对单束团电荷数进行扫描。

---

Parameters
---
`ratiovec`: `Vector` 需要扫描的“流强与最大流强的比例”。

`Nb`: 束团电子数。

`JJ`: 作用量采样点。

`ΔJJ`: 采样作用量步长。之所以需要提前算，是因为过拉伸时，头部和尾部、尾部和外围作用量之间存在跳跃。

`ψJ`: 与`JJ`对应的径向密度分布.

`νsJ`: `JJ`处的纵向振荡工作点. 纵向工作点的定义为`ωs(J)/ω0`，其中`ωs(J)`通过正则变换$\dot{\phi}=\partial H / \partial J$确定。

`zJθ`,`ΦJθ`: samples of `z` and angle variable. First and second dimension/index correspond to `J` and `θ`, action and angle variable samples, respectively.

`ωvec`,`imp_trans`: Real ω `Vector` and complex impedance `Vector` of impedance data. Here `ωvec` should be all frequency, not just positive.

`lm`: 最大角向模数.

`par`: `ParDBRF`, 双高频系统参数.

---

Return
---
`ratios`: 扫描电荷量与束团电荷量的比值`ratio`的列表。 和下面两个都能形成数值对。由于一个电荷量存在多个频率，以及多个增长率。为了保证数据配对，这里将比值相应复制了多份。

`growthRates`: `Vector`. 与`ratios`对应的每秒增长率.

`tuneShifts`: `Vector`. 与`ratios`对应的频移.
"""
function scanNb(ratiovec::Vector, Nb, JJ, ΔJJ, ψJ, νsJ, zJθ, ΦJθ, ωvec, imp_trans, lm, par::ParDBRF)
    nJ = length(JJ)
    matSlωJ = funS(collect(-lm:lm), ωvec, JJ, zJθ, ΦJθ, par)
    ℱj1j2l1l2 = calℱj1j2l1l2(matSlωJ, imp_trans, ωvec, nJ, lm)
    average_ωs = get_average_νs(JJ, νsJ, ψJ) * par.ω0

    ratios, tuneShifts, growthRates = Float64[], Float64[], Float64[]
    M1o = Nb/(30clight) .* calM1WithoutNb(ℱj1j2l1l2, ΔJJ, ψJ, nJ, lm, par)
    M2 = calM2(νsJ, nJ, lm, par)
    for ratio in ratiovec
        println(ratio)
        tmpM1 = ratio .* M1o
        m = tmpM1 .+ M2
        ΔΩeigvals = eigvals(m)
        t = real.(ΔΩeigvals)
        g = imag.(ΔΩeigvals)
        append!(tuneShifts, t)
        append!(growthRates, g)
        append!(ratios, [ratio for i in t])
    end
    ratios, growthRates, tuneShifts ./ average_ωs
end



@doc raw"""
    getDatas(nJ::Int64, nθ::Int64, times_of_σz::Float64, funψzδ, par::ParDBRF)

对双高频相空间进行采样。

Parameters
---
`nJ`: 径向采样点数。

`nθ`: 角向采样点数。需要满足一定的公式。你输入的不合适程序会给出建议。根据经验，以`99`结尾都可以，`199`,`299`,...等等。

`times_of_σz`: Sample radius from extreme point.

`funψzδ`: Density distribution of `z`, `δ`.

`par`: parameters of double RFs. `struct`: `ParDBRF`

Return
---
`EJ`: 被采样的哈密顿量。

`JJ`: 与哈密顿量对应的作用量。

`ΔJJ`: 采样作用量步长。这里不是简单将`JJ`相邻两个相减，因为过拉伸时，头部和尾部、尾部和外围作用量之间存在跳跃，所以需要独自计算。

`zJθ`: 每个相空间采样点的`z`坐标。

`δJθ`: 每个相空间采样点的`δ`坐标。

`ϕJθ`: 每个相空间采样点的角变量`ϕ`坐标。

`νsJ`: `JJ`处的纵向振荡工作点. 纵向工作点的定义为`ωs(J)/ω0`，其中`ωs(J)`通过正则变换$\dot{\phi}=\partial H / \partial J$确定。

`ψJ`: 与`JJ`对应的径向密度分布.

`nJ`: `JJ`的长度。返回的类型有两种，只有一个SFP时，就是`JJ`的长度；有两个SFP时，分为头部、尾部、外围的采样点长度。
"""
function getDatas(nJ::Int64, nθ::Int64, times_of_σz::Float64, funψzδ, par::ParDBRF)
    tmpnθ = number_angle_samples(nθ)
    if tmpnθ != nθ
        error("建议nθ从$nθ 更改为$tmpnθ")
    end

    minBucketBound, maxBucketBound = getBucketBounds(par)
    sfps, ufps, _ = getFPs(par)

    if length(sfps) == 2
        # 适合过拉伸
        hamils = hamiltonian.(sfps, 0.0, par)
        ufp_between_SFPs = (sfps[1] < ufps[1] < sfps[2]) ? ufps[1] : ufps[2]
        # bound of sample points, only need bound one side
        lb_cal = max(minBucketBound+1e-5, ufp_between_SFPs-times_of_σz)
        ub_cal = min(maxBucketBound-1e-5, ufp_between_SFPs+times_of_σz)

        if hamils[1] > hamils[2]
            println("两个SFPs: $(sfps)， 左侧Hamiltonian较高")
            Δz = (sfps[1] - lb_cal)/nJ
            # Sample 1
            zs = collect(range(sfps[1]-Δz, step=-Δz, length=nJ))
            nJ, nJ_H, nJ_T, nJ_O, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, par)
            EJ1, JJ1, zJθ1 = get_E_J_z(zbounds, nJ, nθ, par)
            # Sample 2
            zs = zs .+ (sfps[1] .- zs)./1e6
            _, _, _, _, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, par)
            EJ2, JJ2, zJθ2 = get_E_J_z(zbounds, nJ, nθ, par)
        else
            println("两个SFPs: $(sfps)， 右侧Hamiltonian较高")
            Δz = (ub_cal - sfps[2])/nJ
            # Sample 1
            zs = collect(range(sfps[2]+Δz, step=Δz, length=nJ))
            nJ, nJ_H, nJ_T, nJ_O, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, par)
            EJ1, JJ1, zJθ1 = get_E_J_z(zbounds,         nJ, nθ, par)
            # Sample 2
            zs = zs .- (zs .- sfps[2]) ./ 1e6
            _, _, _, _, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, par)
            EJ2, JJ2, zJθ2 = get_E_J_z(zbounds, nJ, nθ, par)
        end
        ΔJJ1 = getΔJJ(JJ1, nJ, nJ_H, nJ_T)
        δJθ = getδJθ(zJθ1, EJ1, nJ, nθ, par)
        ψJ = getψJ(nJ, nθ, zJθ1, δJθ, funψzδ)
        νsJ = getνsJ(JJ1, JJ2, EJ1, EJ2, nJ, par.circum)
        ϕJθ = getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ, nθ, par)
        return EJ1, JJ1, ΔJJ1, zJθ1, δJθ, ϕJθ, νsJ, ψJ, (nJ_H, nJ_T, nJ_O)
    elseif length(sfps) == 1
        println("一个SFP: $(sfps)")
        # 适用单高频、欠拉伸、最优拉伸
        lb_cal = max(minBucketBound+1e-5, sfps[1]-times_of_σz)
        Δz = (sfps[1] - lb_cal)/nJ
        zs = collect(range(sfps[1]-Δz, step=-Δz, length=nJ))
        EJ1, JJ1, zJθ1 = get_E_J_z(zs, nJ, nθ, minBucketBound, maxBucketBound, par)
        zs = zs .+ (sfps[1] .- zs)./1e6
        EJ2, JJ2, zJθ2 = get_E_J_z(zs, nJ, nθ, minBucketBound, maxBucketBound, par)
        ΔJJ1 = getΔJJ(JJ1, nJ)
        δJθ = getδJθ(zJθ1, EJ1, nJ, nθ, par)
        ψJ = getψJ(nJ, nθ, zJθ1, δJθ, funψzδ)
        νsJ = getνsJ(JJ1, JJ2, EJ1, EJ2, nJ, par.circum)
        ϕJθ = getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ, nθ, par)
        return EJ1, JJ1, ΔJJ1, zJθ1, δJθ, ϕJθ, νsJ, ψJ, nJ
    end
end

@doc raw"""
    getDatas(nJ::Int64, nθ::Int64, times_of_σz::Float64, zvec::Vector, ρvec::Vector, par::ParDBRF)

获取双高频系统相空间采样点数据。

Parameters
---
`nJ`: 径向采样点。

`nθ`: 角向采样点数。因为需要考虑过拉伸一个环和两个环，角向采样点数要一致。所以这里不管欠拉伸还是过拉伸，统一要求`nθ`满足特定公式。可以先`number_angle_samples(nθ)=nθ`判断一下。根据经验，`199,299,399`都是可以的。

`times_of_σz`: `σz`的倍数。

Return
---
`EJ`: 被采样的哈密顿量。

`JJ`: 与哈密顿量对应的作用量。

`ΔJJ`: 采样作用量步长。这里不是简单将`JJ`相邻两个相减，因为过拉伸时，头部和尾部、尾部和外围作用量之间存在跳跃，所以需要独自计算。

`zJθ`: 每个相空间采样点的`z`坐标。

`δJθ`: 每个相空间采样点的`δ`坐标。

`ϕJθ`: 每个相空间采样点的角变量`ϕ`坐标。

`νsJ`: `JJ`处的纵向振荡工作点. 纵向工作点的定义为`ωs(J)/ω0`，其中`ωs(J)`通过正则变换$\dot{\phi}=\partial H / \partial J$确定。

`ψJ`: 与`JJ`对应的径向密度分布.

`nJ`: `JJ`的长度。返回的类型有两种，只有一个SFP时，就是`JJ`的长度；有两个SFP时，分为头部、尾部、外围的采样点长度。
"""
function getDatas(nJ::Int64, nθ::Int64, times_of_σz::Float64, zvec::Vector, ρvec::Vector, par::ParDBRF)
    funψzδ(z, δ)=LinearInterpolation(zvec, ρvec)(z)*1/(sqrt(2π)*par.σδ)*exp(-0.5*(δ/par.σδ)^2)
    getDatas(nJ, nθ, times_of_σz, funψzδ, par)
end


function getDatas(nJ::Int64, nθ::Int64, times_of_σz::Float64, funψzδ, func_potential, minBucketBound, maxBucketBound, sfps::Vector{Float64}, ufps::Vector{Float64}, circum::Float64, η::Float64)
    tmpnθ = number_angle_samples(nθ)
    if tmpnθ != nθ
        error("建议nθ从$nθ 更改为$tmpnθ")
    end

    if length(sfps) == 2    # 如果存在两个bucket
        # 适合过拉伸
        hamils = func_potential.(sfps)
        ufp_between_SFPs = (sfps[1] < ufps[1] < sfps[2]) ? ufps[1] : ufps[2]
        # bound of sample points, only need bound one side
        lb_cal, ub_cal = minBucketBound, maxBucketBound

        if hamils[1] > hamils[2]
            println("两个SFPs: $(sfps)， 左侧Hamiltonian较高")
            Δz = (sfps[1] - lb_cal)/nJ
            # Sample 1
            zs = collect(range(sfps[1]-Δz, step=-Δz, length=nJ))
            nJ, nJ_H, nJ_T, nJ_O, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, func_potential)
            EJ1, JJ1, zJθ1 = get_E_J_z(zbounds, nJ, nθ, func_potential, η)
            # Sample 2
            zs = zs .+ (sfps[1] .- zs)./1e6
            _, _, _, _, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, func_potential)
            EJ2, JJ2, zJθ2 = get_E_J_z(zbounds, nJ, nθ, func_potential, η)
        else
            println("两个SFPs: $(sfps)， 右侧Hamiltonian较高")
            Δz = (ub_cal - sfps[2])/nJ
            # Sample 1
            zs = collect(range(sfps[2]+Δz, step=Δz, length=nJ))
            nJ, nJ_H, nJ_T, nJ_O, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, func_potential)
            EJ1, JJ1, zJθ1 = get_E_J_z(zbounds,         nJ, nθ, func_potential, η)
            # Sample 2
            zs = zs .- (zs .- sfps[2]) ./ 1e6
            _, _, _, _, zbounds = redefine_nJ_and_zs_ChuntaoLin(zs, minBucketBound, maxBucketBound, func_potential)
            EJ2, JJ2, zJθ2 = get_E_J_z(zbounds, nJ, nθ, func_potential, η)
        end
        ΔJJ1 = getΔJJ(JJ1, nJ, nJ_H, nJ_T)
        δJθ = getδJθ(zJθ1, EJ1, nJ, nθ, func_potential, η)
        ψJ = getψJ(nJ, nθ, zJθ1, δJθ, funψzδ)
        νsJ = getνsJ(JJ1, JJ2, EJ1, EJ2, nJ, circum)
        ϕJθ = getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ, nθ, func_potential, η)
        return EJ1, JJ1, ΔJJ1, zJθ1, δJθ, ϕJθ, νsJ, ψJ, (nJ_H, nJ_T, nJ_O)
    elseif length(sfps) == 1
        println("一个SFP: $(sfps)")
        # 适用单高频、欠拉伸、最优拉伸
        lb_cal = max(minBucketBound+1e-5, sfps[1]-times_of_σz)
        Δz = (sfps[1] - lb_cal)/nJ
        zs = collect(range(sfps[1]-Δz, step=-Δz, length=nJ))
        EJ1, JJ1, zJθ1 = get_E_J_z(zs, nJ, nθ, minBucketBound, maxBucketBound, func_potential, η)
        zs = zs .+ (sfps[1] .- zs)./1e6
        EJ2, JJ2, zJθ2 = get_E_J_z(zs, nJ, nθ, minBucketBound, maxBucketBound, func_potential, η)
        ΔJJ1 = getΔJJ(JJ1, nJ)
        δJθ = getδJθ(zJθ1, EJ1, nJ, nθ, func_potential, η)
        ψJ = getψJ(nJ, nθ, zJθ1, δJθ, funψzδ)
        νsJ = getνsJ(JJ1, JJ2, EJ1, EJ2, nJ, circum)
        ϕJθ = getϕJθ(zJθ1, EJ1, JJ1, zJθ2, EJ2, JJ2, nJ, nθ, func_potential, η)
        return EJ1, JJ1, ΔJJ1, zJθ1, δJθ, ϕJθ, νsJ, ψJ, nJ
    end
end





function get_average_νs(JJ, νsJ, ψJ)
    2π * discrete_integrate(JJ, νsJ, ψJ)
end


"""
生成合适的角向采样点数.

角向采样点需要能够让一个圈的等高线和两个圈的等高线点数相等， 因此需要一定的条件。
"""
function number_angle_samples(nθ::Int64)
    m = Int64(round((nθ - 7)/4))
    n = 2+2m
    2*n+3
end


function get_E_J_z(zs, nJ, nθ, minBucket, maxBucket, par::ParDBRF)
    zJθ = zeros(nJ, nθ)
    Js = zeros(nJ)
    Es = getPotenCavity.(zs, par)
    @sync @inbounds for i=1:nJ
        Threads.@spawn begin
            zbounds = get_zero_points_of_hamiltonian(Es[i], minBucket, maxBucket, par)
            zJθ[i, :] = getZSamplesOfE(zbounds, nθ)
            Js[i] = getJOfE(zbounds, Es[i], par)
        end
    end
    Es, Js, zJθ
end
function get_E_J_z(zs, nJ, nθ, minBucket, maxBucket, func_potential, η)
    zJθ = zeros(nJ, nθ)
    Js = zeros(nJ)
    Es = func_potential.(zs)
    @sync @inbounds for i=1:nJ
        Threads.@spawn begin
            zbounds = get_solution(minBucket, maxBucket, func_potential, Es[i])
            zJθ[i, :] = getZSamplesOfE(zbounds, nθ)
            Js[i] = getJOfE(zbounds, Es[i], func_potential, η)
        end
    end
    Es, Js, zJθ
end
function get_E_J_z(zbounds, nJ, nθ, par::ParDBRF)
    zJθ = zeros(nJ, nθ)
    Js = zeros(nJ)
    Es = zeros(nJ)
    @sync @inbounds for i=1:nJ
        Threads.@spawn begin
            zJθ[i,:] = getZSamplesOfE(zbounds[i], nθ)
            Es[i] = getPotenCavity(zbounds[i][1], par)
            Js[i] = getJOfE(zbounds[i], Es[i], par)
        end
    end
    Es, Js, zJθ
end
"""
`fp`: function of potential
"""
function get_E_J_z(zbounds, nJ, nθ, func_potential, η)
    zJθ = zeros(nJ, nθ)
    Js = zeros(nJ)
    Es = zeros(nJ)
    @sync @inbounds for i=1:nJ
        Threads.@spawn begin
            zJθ[i,:] = getZSamplesOfE(zbounds[i], nθ)
            Es[i] = func_potential(zbounds[i][1])
            Js[i] = getJOfE(zbounds[i], Es[i], func_potential, η)
        end
    end
    Es, Js, zJθ
end

"""
与林椿涛师兄的采样方法相对应的zbounds
"""
function redefine_nJ_and_zs_ChuntaoLin(zs, minBucket, maxBucket, par::ParDBRF)
    Es = getPotenCavity.(zs, par)
    zero_points = get_zero_points_of_hamiltonian.(Es, minBucket, maxBucket, par)
    zbounds_H = Vector[]
    zbounds_T = Vector[]
    zbounds_O = Vector[]
    # Redefine new nJ
    nJ = length(filter(x->length(x)==2, zero_points)) + length(filter(x->length(x)==4, zero_points))
    for i in eachindex(Es)
        if length(zero_points[i]) == 4                          # 哈密顿等高线处于两个bucket内
            push!(zbounds_T, zero_points[i][1:2])
            push!(zbounds_H, zero_points[i][3:4])
        elseif (length(zero_points[i]) == 2) & (Es[i] > 0)      # 哈密顿等高线处于一个更深的bucket内
            if zero_points[i][2] < 0
                push!(zbounds_T, zero_points[i])
            else
                push!(zbounds_H, zero_points[i])
            end
        elseif length(zero_points[i]) == 2                      # 哈密顿量不在bucket中
            push!(zbounds_O, zero_points[i])
        end
    end
    nJ_H = length(zbounds_H)
    nJ_T = length(zbounds_T)
    nJ_O = length(zbounds_O)
    nJ = nJ_H + nJ_T + nJ_O
    nJ, nJ_H, nJ_T, nJ_O, vcat(zbounds_H, zbounds_T, zbounds_O)
end
function redefine_nJ_and_zs_ChuntaoLin(zs, minBucket, maxBucket, fp)
    Es = fp.(zs)
    zero_points = get_solution.(minBucket, maxBucket, fp, Es)
    zbounds_H = Vector[]
    zbounds_T = Vector[]
    zbounds_O = Vector[]
    # Redefine new nJ
    nJ = length(filter(x->length(x)==2, zero_points)) + length(filter(x->length(x)==4, zero_points))
    for i in eachindex(Es)
        if length(zero_points[i]) == 4                          # 哈密顿等高线处于两个bucket内
            push!(zbounds_T, zero_points[i][1:2])
            push!(zbounds_H, zero_points[i][3:4])
        elseif (length(zero_points[i]) == 2) & (Es[i] > 0)      # 哈密顿等高线处于一个更深的bucket内
            if zero_points[i][2] < 0
                push!(zbounds_T, zero_points[i])
            else
                push!(zbounds_H, zero_points[i])
            end
        elseif length(zero_points[i]) == 2                      # 哈密顿量不在bucket中
            push!(zbounds_O, zero_points[i])
        end
    end
    nJ_H = length(zbounds_H)
    nJ_T = length(zbounds_T)
    nJ_O = length(zbounds_O)
    nJ = nJ_H + nJ_T + nJ_O
    nJ, nJ_H, nJ_T, nJ_O, vcat(zbounds_H, zbounds_T, zbounds_O)
end
"""
徐景晔的采样方法相应的对应方法
"""
function get_zbounds_JingyeXu(zs, minBucket, maxBucket, func_potetial)
    Es = func_potetial.(zs)
    zero_points = get_solution.(minBucket, maxBucket, func_potetial, Es)
    zero_points = [zs for zs in zero_points if length(zs) % 2 == 0]
    zsL = [zs[1] for zs in zero_points]
    zsR = [zs[end] for zs in zero_points]
    n1 = n2 = 0
    for i in eachindex(zero_points)
        if length(zero_points[i]) == 4
            break
        elseif length(zero_points[i]) == 2
            n1 += 1
        end
    end
    for i = n1+1:length(zero_points)
        if length(zero_points[i]) == 2
            break
        elseif length(zero_points[i]) == 4
            n2 += 1
        end
    end
    n3 = length(zero_points) - n1 - n2
    zsL, zsR, (n1+n2+n3), n1, n2, n3, zero_points
end


function getνsJ(JJ1, JJ2, EJ1, EJ2, nJ, circum)
    [circum/(2π) * (EJ2[i] - EJ1[i])/(JJ1[i] - JJ2[i]+1e-150) for i in 1:nJ]
end


function getψJ(nJ, nθ, zJθ, δJθ, funψzδ)
    ψJ = zeros(nJ)
    @inbounds for Jidx=1:nJ
        res = 0.0
        @inbounds for θidx=1:nθ
            res += funψzδ(zJθ[Jidx, θidx], δJθ[Jidx, θidx])
        end
        ψJ[Jidx] = res/nθ
    end
    ψJ
end

"""
Calculate `δ` at `z`.
"""
function getz2δ(z, E, par::ParDBRF)
    zPoten = getPotenCavity(z, par)
    zPoten <= E ? zero(typeof(z*E)) : sqrt(2 / par.η * (zPoten - E))
end
function getz2δ(z, E, func_potential, η)
    tmp_poten = func_potential(z)
    tmp_poten <= E ? zero(typeof(z*E)) : sqrt(2 / η * (tmp_poten - E))
end



function integrate_of_δ(z, zmax, E, par::ParDBRF)
    quadgk(x -> getz2δ(x, E, par), z, zmax)[1]  # 这个地方是否要提升精度？
end
function integrate_of_δ(z, zmax, E, func_potential, η)
    quadgk(x -> getz2δ(x, E, func_potential, η), z, zmax)[1]
end




function getJOfE(zBounds, E, par::ParDBRF)
    if length(zBounds) == 2
        return 1/π*integrate_of_δ(zBounds[1], zBounds[2], E, par)
    elseif length(zBounds) == 4
        return 1/π*integrate_of_δ(zBounds[1], zBounds[2], E, par) + 1/π*integrate_of_δ(zBounds[3], zBounds[4], E, par)
    else
        error("求出的边界为：", zBounds, " 相应的E=", E, ". Potential的解应该有2或4.")
    end
end
function getJOfE(zBounds, E, func_potential, η)
    if length(zBounds) == 2
        return 1/π*integrate_of_δ(zBounds[1], zBounds[2], E, func_potential, η)
    elseif length(zBounds) == 4
        return 1/π*integrate_of_δ(zBounds[1], zBounds[2], E, func_potential, η) + 1/π*integrate_of_δ(zBounds[3], zBounds[4], E, func_potential, η)
    else
        error("求出的边界为：", zBounds, " 相应的E=", E, ". Potential的解应该有2或4.")
    end
end



"""
Similar to Bessel function in Alex. Chao, 1993, Eq.(6.75).
Actually in single RF condition, this function equals `i^{-l} J_l(ωr/c)`
"""
function funS(lvec::Vector{Int64}, ωvec::Vector, Jvec::Vector, zJϕ::Matrix, ϕJϕ::Matrix, par::ParDBRF)
    llen = length(lvec)
    ωlen = length(ωvec)
    Jlen = length(Jvec)
    m = zeros(ComplexF64, llen, ωlen, Jlen)
    @inbounds for lidx=1:llen
        @sync @inbounds for Jidx=1:Jlen
            Threads.@spawn begin
                @inbounds for ωidx=1:ωlen
                    tmpϕvec = ϕJϕ[Jidx, :]
                    tmpres = exp.(1.0im.*(lvec[lidx].*tmpϕvec .- zJϕ[Jidx,:]./clight.*(ωvec[ωidx] .- par.ωξ)))
                    m[lidx, ωidx, Jidx] = discrete_integrate(tmpϕvec, tmpres)/2π
                end
            end
        end
    end
    m
end


"""
matrix ℱ_{l,l'}(J,J')
"""
function calℱj1j2l1l2(matSlωJ::Array{ComplexF64, 3}, Zvec::Vector, ωvec::Vector, nJ::Int64, lm::Int64)
    resmat = zeros(ComplexF64, nJ, nJ, 2lm+1, 2lm+1)
    @sync @inbounds for j1=1:nJ, j2=1:nJ
        Threads.@spawn begin
            @inbounds for l1=-lm:lm, l2=-lm:lm
                idx1 = lm+l1+1
                idx2 = lm+l2+1
                resmat[j1,j2,idx1,idx2]=discrete_integrate(ωvec, Zvec, conj.(matSlωJ[idx1,:,j1]), matSlωJ[idx2,:,j2])
            end
        end
    end
    resmat
end


function getΔJJ(JJ::Vector, nJ::Int)
    ΔJJ=zeros(nJ)
    ΔJJ[1]=JJ[2]-JJ[1]
    ΔJJ[end]=JJ[end]-JJ[end-1]
    @inbounds for i=2:(nJ-1)
        ΔJJ[i] = 0.5*(JJ[i+1]-JJ[i-1])
    end
    ΔJJ
end
function getΔJJ(JJ::Vector, nJ::Int, nJH::Int, nJT::Int)
    ΔJJ=zeros(nJ)
    
    if nJH >= 2
        ΔJJ[1] = JJ[2] - JJ[1]
        ΔJJ[nJH] = JJ[nJH] - JJ[nJH-1]
        for i=2:nJH-1
            ΔJJ[i] = 0.5 * (JJ[i+1] - JJ[i-1])
        end
    end

    if nJT >= 2
        ΔJJ[nJH+1] = JJ[nJH+2] - JJ[nJH+1]
        ΔJJ[nJH+nJT] = JJ[nJH+nJT] - JJ[nJH+nJT-1]
        for i=nJH+2:nJH+nJT-1
            ΔJJ[i] = 0.5 * (JJ[i+1] - JJ[i-1])
        end
    end
    
    if nJ - nJH - nJT >= 2
        ΔJJ[nJH+nJT+1] = JJ[nJH+nJT+2] - JJ[nJH+nJT+1]
        ΔJJ[nJ] = JJ[nJ] - JJ[nJ-1]
        for i=nJH+nJT+2:nJ-1
            ΔJJ[i] = 0.5 * (JJ[i+1] - JJ[i-1])
        end
    end
    
    ΔJJ
end