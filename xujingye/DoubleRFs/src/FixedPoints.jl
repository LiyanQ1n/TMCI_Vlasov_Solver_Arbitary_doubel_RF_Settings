# Need to prevent the over-stretching condition that two sub-bucket are separated compeletly.


@doc raw"""
    classify_FP(z, ϕs, ϕ2s, r, h1, h, circum)::String

验证到底是SFP还是UFP。

输入的z要确保是FP。
"""
function classify_FP(z, ϕs, ϕ2s, r, h1, h, circum)::String
    res1 = ∇normedVoltage(z, ϕs, ϕ2s, r, h1, h, circum)
    println(res1)
    if res1 < -1e-20        # When something about fixed point is wrong, try to adjust this
        return "UFP"
    elseif res1 > 1e-20
        return "SFP"
    else
        res2 = ∇²normedVoltage(z, ϕs, ϕ2s, r, h1, h, circum)
        if res2 < 0
            return "UFP"
        else
            return "SFP"
        end
    end
end
@doc raw"""
    classify_FP(z, par::ParDBRF)::String

验证到底是SFP还是UFP。

输入的z要确保是FP。
"""
function classify_FP(z, par::ParDBRF)::String
    classify_FP(z, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end



"""
为了更高精度用BigFloat作为基础函数
"""
function newton_iterate_FP_onestep(zn::BigFloat, ϕs, ϕ2s, r, h1, h, circum)::Float64
    zn - normedVoltage(zn, ϕs, ϕ2s, r, h1, h, circum)/∇normedVoltage(zn, ϕs, ϕ2s, r, h1, h, circum)
end
function newton_iterate_FP_onestep(zn::Float64, ϕs, ϕ2s, r, h1, h, circum)
    newton_iterate_FP_onestep(BigFloat(zn), ϕs, ϕ2s, r, h1, h, circum)
end
function newton_iterate_FP_onestep(zn::Float64, par)
    newton_iterate_FP_onestep(zn, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end



"""
牛顿迭代法，计算精确的不动点
"""
function newton_iterate_FP(z0, ϕs, ϕ2s, r, h1, h, circum; tol=1e-12, max_try_time=50)
    zn = z0
    zo = z0+10*tol
    try_time = 0
    while (try_time < max_try_time) && (abs(zn-zo)>tol)
        try_time += 1
        zo, zn = zn, newton_iterate_FP_onestep(zn, ϕs, ϕ2s, r, h1, h, circum)
    end
    zn
end
function newton_iterate_FP(z0, par::ParDBRF; tol=1e-12, max_try_time=50)
    newton_iterate_FP(z0, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum, tol=tol, max_try_time=max_try_time)
end



"""
得到大概的不动点的区间
"""
function get_more_FPs(ϕs, ϕ2s, r, h1, h, circum; max_try_time)
    λ1 = circum/h1
    ite_step = λ1/10000
    resvec = [0.0]
    for z in -λ1*0.99:ite_step:(0.99*λ1-ite_step)
        if (normedVoltage(z, ϕs, ϕ2s, r, h1, h, circum)*normedVoltage(z+ite_step, ϕs, ϕ2s, r, h1, h, circum) < 0) && (z*(z+ite_step) > 0)
            tmp_z0 = z+0.5*ite_step
            tmp_z0 = newton_iterate_FP(tmp_z0, ϕs, ϕ2s, r, h1, h, circum; max_try_time=max_try_time)
            append!(resvec, tmp_z0)
        end
    end
    resvec
end
function get_more_FPs(par::ParDBRF; max_try_time)
    get_more_FPs(par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum; max_try_time=max_try_time)
end



"""
计算不动点。

RETURN
---
返回元素: 第一个为SFP的Vector，第二个为UFP的Vector，第三个为FP（前两者的综合）的Vector

"""
function getFPs(ϕs, ϕ2s, r, h1, h, circum, centerenergy, αc; tol=1.0e-12, max_try_time=50)
    η = αc - (restenergy / centerenergy)^2
    FPs = get_more_FPs(ϕs, ϕ2s, r, h1, h, circum; max_try_time=max_try_time)
    sort!(FPs)
    println(FPs)
    remove_approx_data!(FPs, atol=tol)
    println(FPs)
    filter_FPs_within_a_bucket(FPs, ϕs, ϕ2s, r, h1, h, circum, η; tol=tol)
end
function getFPs(par::ParDBRF; tol = 1.0e-12, max_try_time = 50)
    getFPs(par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum, par.centerenergy, par.αc, tol=tol, max_try_time=max_try_time)
end


"""
当得到的不动点数量很多时， 筛选出包含0的最大的bucket.

特别的， 对于过拉伸， 既会包含两个小bucket的不动点， 也会得到最大bucket的不动点。
"""
function filter_FPs_within_a_bucket(fps, ϕs, ϕ2s, r, h1, h, circum, η; tol)
    cri = classify_FP(0.0, ϕs, ϕ2s, r, h1, h, circum)
    zero_idx = findfirst(x->abs(x) < tol, fps)
    if cri=="UFP"
        return _filter_FPs_within_a_bucket_over_stretching(fps, zero_idx, η)
    else
        return _filter_FPs_within_a_bucket_under_stretching(fps, zero_idx, η)
    end
end
function _filter_FPs_within_a_bucket_over_stretching(fps, zero_idx, η)
    if η < 0
        return [fps[zero_idx-1], fps[zero_idx+1]], [fps[zero_idx-2], fps[zero_idx]], fps[zero_idx-2:zero_idx+1]
    else
        return [fps[zero_idx-1], fps[zero_idx+1]], [fps[zero_idx], fps[zero_idx+2]], fps[zero_idx-1:zero_idx+2]
    end
end
function _filter_FPs_within_a_bucket_under_stretching(fps, zero_idx, η)
    if η < 0
        return [fps[zero_idx]], [fps[zero_idx-1]], fps[zero_idx-1:zero_idx]
    else
        return [fps[zero_idx]], [fps[zero_idx+1]], fps[zero_idx:zero_idx+1]
    end
end





"""
计算bucket上下限。

只适合η>0的情形。
"""
function getBucketBoundsPositiveη(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc; tol=1.0e-12, max_try_time=50)
    SFPs, UFPs, FPs = getFPs(ϕs, ϕ2s, r, h1, h, circum, centerenergy, αc, tol=tol, max_try_time=max_try_time)
    println("SFPs: ", SFPs)
    println("UFPs: ", UFPs)
    println("FPs: ", FPs)
    ubound = maximum(UFPs)
    c = _cache_one_turn_map(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    succ = false
    while succ == false # 不断循环直到成功
        zo, δo = ubound - 1e-3, 0.0 # 开局离开SFP一点。因为之前遇到一个巧合无论怎么迭代都在SFP.
        zn, δn = one_turn_map_with_cache(zo, δo, ϕs, ϕ2s, r, circum, αc, c)
        trytime = 0
        while ((δn * δo > 0) || (trytime < 10)) && (trytime < 1e8)   # 重复足够次数后找到零点会退出；超过一定次数也会退出。
            trytime += 1
            zo, δo = zn, δn
            zn, δn = one_turn_map_with_cache(zo, δo, ϕs, ϕ2s, r, circum, αc, c)
        end
        if trytime >= 1e7    # 循环足够多次没找到，就换个值找
            ubound = ubound - 1e-3
            println(ubound)
        else
            succ = true     # 判定成功
            return zn, ubound
        end
    end
end


function getBucketBoundsNegativeη(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc; tol = 1.0e-12, max_try_time = 50)
    lbound = minimum(getFPs(ϕs, ϕ2s, r, h1, h, circum, centerenergy, αc, tol=tol, max_try_time=max_try_time)[2])
    c = _cache_one_turn_map(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    succ = false
    while succ == false # 不断循环直到成功
        zo, δo = lbound+1e-3, 0.0
        zn, δn = one_turn_map_with_cache(zo, δo, ϕs, ϕ2s, r, circum, αc, c)
        trytime = 0
        while ((δn*δo>0) || (trytime < 100)) && (trytime < 1e8)   # 重复足够次数后找到零点会退出；超过一定次数也会退出。
            trytime += 1
            zo, δo = zn, δn
            zn, δn = one_turn_map_with_cache(zo, δo, ϕs, ϕ2s, r, circum, αc, c)
        end
        if trytime >= 1e7    # 循环足够多次没找到，就换个值找
            lbound = lbound + 1e-3
        else
            succ = true     # 判定成功
            return lbound, zn
        end
    end
end





"""
    getBucketBounds(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc; tol = 1.0e-12, max_try_time = 50)

计算相稳定区边界

---

Parameters
---
`ϕ2`: 参考粒子主腔相位。单位弧度。

`ϕ2s`: 参考粒子谐波腔相位。单位弧度。

`v1`: 主腔电压。

`r`: 谐波腔与主腔电压之比。

`h1`: 主腔的谐波数.

`h`: 两个腔谐波数之比.

`circum`: 储存环周长.

`centerenergy`: 参考粒子能量.

`αc`: 动量压缩因子.

Return
---
`zL,zR`: bucket的上下限。
"""
function getBucketBounds(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc; tol = 1.0e-12, max_try_time = 50)
    γ = centerenergy/restenergy
    η = αc - 1/γ^2
    if η > 0
        return getBucketBoundsPositiveη(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc, tol=tol, max_try_time=max_try_time)
    else
        return getBucketBoundsNegativeη(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc, tol=tol, max_try_time=max_try_time)
    end
end

"""
    getBucketBounds(par::ParDBRF; tol = 1.0e-12, max_try_time = 50)

计算相稳定区边界

---

Parameters
---
`par`: `ParDBRF`。包含加速器与双高频系统参数的参数结构体。

Return
---
`zL,zR`: bucket的上下限。
"""
function getBucketBounds(par::ParDBRF; tol = 1.0e-12, max_try_time = 50)
    getBucketBounds(par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy, par.αc, tol=tol, max_try_time=max_try_time)
end