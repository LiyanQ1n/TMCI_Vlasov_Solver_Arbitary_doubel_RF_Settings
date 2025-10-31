# nθ 需要满足
# nθ = 4*[round((nθ-7)/4)]+7



"""
根据Hamilton等高线的零点数量，进行角向采样。
"""
function getZSamplesOfE(zBounds, nCirclePoint)
    if length(zBounds) == 4
        return getZSamplesOfTwoCircles(zBounds[1], zBounds[2], zBounds[3], zBounds[4], nCirclePoint)
    else
        return getZSamplesOfOneCircle(zBounds[1], zBounds[2], nCirclePoint)
    end
end


function getZSamplesOfOneCircle(zmin, zmax, nθ)
    samplesUpHalfCircle = getZSamplesOfUpHalfCircle(zmin, zmax, div(nθ+1,2))    # 上半平面， 包括左端
    samplesDoHalfCircle = samplesUpHalfCircle[end-1:-1:1]                       # 下半平面， 没有左端
    cat(samplesUpHalfCircle, samplesDoHalfCircle, dims=1)
end


function getZSamplesOfTwoCircles(zminL, zmaxL, zminR, zmaxR, nθ)
    nLCircle = div(nθ-1, 2)
    nRCircle = div(nθ+1, 2)
    samplesUpHalfCircleL = getZSamplesOfUpHalfCircle(zminL, zmaxL, div(nLCircle+1, 2))
    samplesDoHalfCircleL = samplesUpHalfCircleL[end-1:-1:1]
    samplesUpHalfCircleR = getZSamplesOfUpHalfCircle(zminR, zmaxR, div(nRCircle, 2))
    samplesDoHalfCircleR = samplesUpHalfCircleR[end:-1:1]
    cat(samplesUpHalfCircleR, samplesUpHalfCircleL, samplesDoHalfCircleL, samplesDoHalfCircleR, dims=1)
end


"""
返回上半平面， 从右往左的采样.
"""
function getZSamplesOfUpHalfCircle(zmin, zmax, sample_num)
    zcenter = 0.5 * (zmin + zmax)
    rz = 0.5 * (zmax - zmin)
    angleVector = range(0, π, length=sample_num)
    zcenter .+ rz .* cos.(angleVector) 
end