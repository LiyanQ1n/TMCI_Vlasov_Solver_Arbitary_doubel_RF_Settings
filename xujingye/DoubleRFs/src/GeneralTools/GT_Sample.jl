"""
Sampling according to density distribution.

Parameters
---
`n`: number of samples.

`densityFun`: function of density distribution.

`xmin`,`xmax`: low bound and up bound of `densityFun`.

`ymax`: maximum value of `densityFun`.
"""
function samplingAccordingDensityFun(n::Int64, densityFun, xmin::Float64, xmax::Float64, ymax::Float64)
    resvec = zeros(n)
    for i in 1:n
        suc = false
        while !suc
            xrand = randRange(xmin, xmax)
            yrand = randRange(0, ymax)
            if densityFun(xrand) > yrand
                resvec[i] = xrand
                suc = true
            end
        end
    end
    resvec
end


"""
get density distribution from sample vector.

Parameters
---
`samples`: a `Vector` of samples.

`n`: number of domains.

`r`: window radius.

Return
---
`xs`: `Vector` of domain centers.

`nums`: numbers in each domains.
"""
function getDensityDistFromSamples(samples::Vector, n::Int, r::Float64)
    xmin = minimum(samples) - 1e-6
    xmax = maximum(samples) + 1e-6
    xs = range(xmin, xmax, length=n)
    nums = zeros(Int64, n)
    for i in 1:n
        nums[i] = countBetween(samples, xs[i] .- r, xs[i] .+ r)
    end
    xs, nums
end
function countBetween(xvec::Vector, x1, x2)
    res = 0
    for x in xvec
        if x1 < x < x2
            res += 1
        end
    end
    res
end







"""

"""
function extrapolateAndResample(JReSampleVec::Vector{Float64}, JDataVec::Vector{Float64}, ψDataVec::Vector)
    nDataHalf = Int64(ceil(length(JDataVec)/2))
    nReSampleLength = length(JReSampleVec)
    lb = JDataVec[2]
    ub = JDataVec[end-1]
    
    lbound_linear_fun = curve_fit(LinearFit, JDataVec[1:10], ψDataVec[1:10])
    ubound_exp_fun = curve_fit(ExpFit, JDataVec[nDataHalf:end], ψDataVec[nDataHalf:end])
    interp_fun = Spline1D(JDataVec, ψDataVec)
    
    ψReSampleVec = zeros(eltype(ψDataVec), nReSampleLength)
    for i in 1:nReSampleLength
        tmpJ = JReSampleVec[i]
        if tmpJ > ub
            ψReSampleVec[i] = ubound_exp_fun(tmpJ)
        elseif tmpJ < lb
            ψReSampleVec[i] = lbound_linear_fun(tmpJ)
        else
            ψReSampleVec[i] = interp_fun(tmpJ)
        end
    end
    ψReSampleVec
end
function extrapolateAndResample(JReSampleVec::Vector{Float64}, JDataVec::Vector{Float64}, ψDataVec::Vector, outVal::Real)
    nReSampleLength = length(JReSampleVec)
    lb = JDataVec[1]
    ub = JDataVec[end]
    
    lbound_linear_fun = curve_fit(LinearFit, JDataVec[1:10], ψDataVec[1:10])
    interp_fun = Spline1D(JDataVec, ψDataVec)
    
    ψReSampleVec = zeros(eltype(ψDataVec), nReSampleLength)
    for i in 1:nReSampleLength
        tmpJ = JReSampleVec[i]
        if tmpJ > ub
            ψReSampleVec[i] = outVal
        elseif tmpJ < lb
            ψReSampleVec[i] = lbound_linear_fun(tmpJ)
        else
            ψReSampleVec[i] = interp_fun(tmpJ)
        end
    end
    ψReSampleVec
end




"""
Adaptive sample to make integrate of `f` precise.
"""
function adaptiveSample(f, initialDomainLength::Int64, a, b, relativeError; minΔx=1e-20)
    allDomains = collect(range(a, b, length=initialDomainLength))
    firstIndexToCal = 1

    while true
        lb, ub = allDomains[firstIndexToCal], allDomains[firstIndexToCal+1]
        h = ub - lb
        fl = f(lb)
        fc = f(lb+0.5h)
        fr = f(ub)
        devi = abs(fl+fr-2fc)
        if (devi/(abs(max(fl, fc, fr))+1e-60) <= relativeError) || (devi <= 1e-100)
            lengthToCal = length(allDomains) - firstIndexToCal + 1
            if lengthToCal > 2
                firstIndexToCal += 1
            elseif lengthToCal==2
                break
            end
        else
            c = 0.5*(lb+ub)
            if (c==lb) || (c==ub) || abs(ub-lb)<=minΔx
                firstIndexToCal += 1
                if length(allDomains) == firstIndexToCal
                    break
                else
                    continue
                end
            end
            insert!(allDomains, firstIndexToCal+1, c)
        end
    end
    allDomains
end
"""
Adaptive sample to make integrate of `f` precise.
"""
function adaptiveSample(f, initialDomains::Vector, relativeError; minΔx=1e-20)
    allDomains = initialDomains[:]
    firstIndexToCal = 1

    while true
        lb, ub = allDomains[firstIndexToCal], allDomains[firstIndexToCal+1]
        h = ub - lb
        fl = f(lb)
        fc = f(lb+0.5h)
        fr = f(ub)
        devi = abs(fl+fr-2fc)
        if (devi/(abs(max(fl, fc, fr))+1e-60) <= relativeError) || (devi <= 1e-100)
            lengthToCal = length(allDomains) - firstIndexToCal + 1
            if lengthToCal > 2
                firstIndexToCal += 1
            elseif lengthToCal==2
                break
            end
        else
            c = 0.5*(lb+ub)
            if (c==lb) || (c==ub) || abs(ub-lb)<=minΔx
                firstIndexToCal += 1
                if length(allDomains) == firstIndexToCal
                    break
                else
                    continue
                end
            end
            insert!(allDomains, firstIndexToCal+1, c)
        end
    end
    allDomains
end





# ================================= 手动调节的采样方法 ========================
# 变化范围从均匀采样(α=0.0), 到步长指数增长(α~1.0).
"""
步长按指数增长的采样。
"""
function samplingExpStepIncrease(n::Int64, a::Float64, b::Float64; α=0.0)
    xvec = sumOfExpStep(n, α)
    (b-a) .* xvec .+ a
end
"""
以指数增长的步长间隔的求和。
"""
function sumOfExpStep(n::Int64, α::Float64)
    if n<=1
        error("Number of samples too small")
    end

    xvec = zeros(Float64, n)
    for i in 2:n
        xvec[i] = xvec[i-1] + exp(α*20*(i-2)/n) # for different n, (i-2/n) always range from 0 to 1
    end

    xvec./maximum(xvec)     # Results always range from 0 to 1
end



"""
步长按照指数增长的采样。
"""
function samplingExpStepIncrease(n::Int64, a::Float64, b::Float64, c::Float64; α=0.0)
    nbc = getBestNbc(n, a, b, c, α=α)
    nab = n-nbc
    
    xbcvec = samplingExpStepIncrease(nbc+1, b, c, α=α)
    xabvec = xbcvec .- (b-a)
    xbcvec = (xbcvec .- xbcvec[1]) .* ((c-b)/(xbcvec[end]-xbcvec[1])) .+ b
    
    vcat(xabvec[1:nab], xbcvec[2:end])
end
function getBestNbc(n, a, b, c; α=0.0)
    nlargelist = [n]
    nsmalllist = [2]
    succ = false
    while true
        minlargelist = minimum(nlargelist)
        maxsmalllist = maximum(nsmalllist)
        if minlargelist - maxsmalllist <= 2
            return Int(round(0.5*(minlargelist+maxsmalllist)))
        end
        
        nbc = Int(round(0.5*(minlargelist+maxsmalllist)))
        n1, n2 = trynbc(nbc, a, b, c, α=α)
        if 2n1+n2>n
            append!(nlargelist, nbc)
        elseif 2n1+n2<n
            append!(nsmalllist, nbc)
        else
            return nbc
        end
    end
end
function trynbc(nbc, a, b, c; α=0.0)
    countn1n2(samplingExpStepIncrease(nbc, 0.0, c-b, α=α), 0.0, b-a, c-b)
end
function countn1n2(xbcvec::Vector{Float64}, a, b, c)
    n1 = 0
    n2 = 0
    for x in xbcvec
        if b <= x <= c
            n2 += 1
        elseif a <= x < b
            n1 += 1
        end
    end
    n1, n2
end






function samplingExpStepDecrease(n::Int64, b::Float64, a::Float64; α=0.0)
    - samplingExpStepIncrease(n, -b, -a, α=α)
end
function samplingExpStepDecrease(n::Int64, c::Float64, b::Float64, a::Float64; α=0.0)
    - samplingExpStepIncrease(n, -c, -b, -a, α=α)
end

"""
n: 采样点数
`a`,`b`: 递增或者递减的参数
α: 0~1比较合适
"""
function samplingExpStep(n::Int64, a::Float64, b::Float64; α=0.0)
    if a < b
        return samplingExpStepIncrease(n, a, b, α=abs(α))
    elseif b < a
        return samplingExpStepDecrease(n, a, b, α=abs(α))
    end
end
"""
n: 采样点数
`a`,`b`,`c`: 递增或者递减的参数
α: 0~1比较合适
"""
function samplingExpStep(n::Int64, a::Float64, b::Float64, c::Float64; α=0.0)
    if a < b < c
        return samplingExpStepIncrease(n, a, b, c, α=abs(α))
    elseif a > b > c
        return samplingExpStepDecrease(n, a, b, c, α=abs(α))
    else
        error("输入参数不是递增或递减")
    end
end