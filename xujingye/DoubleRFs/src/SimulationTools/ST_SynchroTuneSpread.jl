# Simulation Tools: Calculate tune spread
# TODO:
#   z在tracking的时候，会经历zmin和zmax，因此一次跟踪可以计算两个点。这样更快。


function get_δvec_for_tune_fft(z0, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc; max_iterate_time=10000)
    succ = false
    last2circindex = 0
    zvec = zeros(max_iterate_time)
    δvec = zeros(max_iterate_time)
    zvec[1]=z0
    δvec[1]=0.0
    time_of_zeroδ = 0
    staindex=2
    c = _cache_one_turn_map(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy)
    while succ == false
        z, δ = zvec[staindex-1], δvec[staindex-1]
        for rowindex in staindex:max_iterate_time
            # 得到演化数组
            z, δ = one_turn_map_with_cache(z, δ, ϕs, ϕ2s, r, circum, αc, c)
            zvec[rowindex]=z
            δvec[rowindex]=δ
            # 两次δ符号相反,记穿过一次z轴。穿过两次则为一圈。
            # 如果穿过z轴次数为2*n倍，则说明为n圈。
            if δvec[rowindex]*δvec[rowindex-1]<0    
                time_of_zeroδ+=1                # δ穿过零点的次数
                if time_of_zeroδ%2==0
                    last2circindex = rowindex   # 每穿过零点两次，算为一个周期。我们记录其index。
                end
            end
        end

        # 上面的结果分几种情况，一般来说δ越多越好。
        if time_of_zeroδ >= 40      # 在轨迹中循环10圈，多复制几圈，返回的结果绕了相空间100圈。
            succ = true
            return vec(repeat(δvec[1:last2circindex-1], 10))
        elseif (maximum(abs.(δvec)) > 0.1) || (maximum(abs.(zvec)) > 1)         # 如果跑出了bucket外面去了。这里z,δ的上下限人为取为2,0.2的矩形格子。
            return Nothing
        elseif size(δvec, 1) >= 80000000                                        # 如果δ的长度超过8e7，还没有10圈，频率估计极低，则别管了，也能分析出很低的频率
            return δvec
        else
            temp = Int64(max_iterate_time*ceil(20/(time_of_zeroδ+1)))
            if temp > max_iterate_time
                max_iterate_time = temp
            else
                max_iterate_time = temp*2
            end
            # println(z0, "扩列:", max_iterate_time)
            # 复制已经计算的部分到新数组
            tempzvec = copy(zvec)
            tempδvec = copy(δvec)
            zvec = zeros(max_iterate_time)
            δvec = zeros(max_iterate_time)
            zvec[1:size(tempzvec,1)] = tempzvec
            δvec[1:size(tempzvec,1)] = tempδvec
            staindex = size(tempzvec,1)+1
        end
    end
end
function get_δvec_for_tune_fft(z0, par::ParDBRF; max_iterate_time=10000)
    get_δvec_for_tune_fft(z0, par.ϕs, par.ϕ2s, par.v1, par.r, par.h1, par.h, par.circum, par.centerenergy, par.αc, max_iterate_time=max_iterate_time)
end





"""
对大于幅值`cri_ratio`的频率取平均，作为基频
"""
function get_base_frequency(fre_vec::Vector{Float64}, amp_vec::Vector{Float64}; cri_ratio=0.95)
    cri_amp = cri_ratio*maximum(amp_vec)  # 最大幅度的频率
    sum_fre = 0.0
    count = 0
    for i in 2:size(fre_vec,1)
        if amp_vec[i] > cri_amp
            sum_fre += fre_vec[i]
            count += 1
        elseif (amp_vec[i-1] > cri_amp) && (amp_vec[i] <= cri_amp)
            return sum_fre/count
        end
    end
end



"""
得到频谱。

Param
---
`yvec`: 需要分析的矢量。矢量长度。应该是2的整数次方倍。

`fr`: revolution frequency.
"""
function get_frequency_dist(yvec::Vector{Float64}, fr)
    l = length(yvec)
    Y = abs.(fft(yvec))/l
    Y[2:end-1] = Y[2:end-1]*2
    fvec = collect(0:l-1)*fr/l
    if iseven(l)
        return fvec[1:Int64(l/2)], Y[1:Int64(l/2)]
    else
        return fvec[1:Int64((l-1)/2)], Y[1:Int64((l-1)/2)]
    end
end



"""
利用fft方法，计算tune分布。

Parameters
---
`numPoints`: 采样点的数量。

`v1`: ``V_1``, voltage of main cavity.

`r`: ``r=V_2/V_1``, where ``V_2`` is the voltage of harmonic cavity.

`ϕs`: Synchronous phase of main cavity.

`ϕ2s`: Synchronous phase of harmonic cavity.

`centerenergy`: ``E=γ m_e c^2 ``.

`αc`: 

`h1`:

`h`: ``h=h_2/h_1``.

`circum`:


RETURN
---
Tuple of two `Vector{Float64}`。Vector of ``z`` and Vector of tune。

"""
function getSynchrotronTunes(numPoints, v1, r, ϕs, ϕ2s, centerenergy, αc, h1, h, circum)
    fr = clight/circum     # revolution frequency
    # 得到需要计算的z
    lbound, ubound = getBucketBounds(ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc)
    zvec = collect(range(lbound+1e-15, ubound-1e-15, length=numPoints))  # 这种划分方法显然不合适
    frevec = zeros(size(zvec, 1))
    # 得到需要计算的δ
    for i in eachindex(zvec)
        z0 = zvec[i]
        δvec = get_δvec_for_tune_fft(z0, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc)
        if typeof(δvec) <: AbstractArray    # 避免返回的nothing值
            fvec, yvec = get_frequency_dist(δvec, fr)
            frevec[i] = get_base_frequency(fvec, yvec)
        end
    end
    zvec, frevec/fr
end
function getSynchrotronTunes(numPoints, par::ParDBRF)
    getSynchrotronTunes(numPoints, par.v1, par.r, par.ϕs, par.ϕ2s, par.centerenergy, par.αc, par.h1, par.h, par.circum)
end