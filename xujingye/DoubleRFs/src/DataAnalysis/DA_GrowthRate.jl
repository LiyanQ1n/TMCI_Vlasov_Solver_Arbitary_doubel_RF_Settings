# Data Analysis: Growth Rate

"""
Calculate growth rate according to turns, RMS width every turn, particle number every turn.

`turnVec`: turns recording data.

`σWidthVec`: RMS width.

`numParticleVec`: particle numbers.
"""
function getGrowthRatePerTurn(turnVec::Vector{Int64}, σWidthVec::Vector{Float64}, particlenumVec::Vector{Int64}, circum::Float64)
    turnVec, σWidthVec = removeDataParticleLoss(turnVec, σWidthVec, particlenumVec)
    turnVec, σWidthVec = removeDataAfterExpGrowth(turnVec, σWidthVec)
    dot(transpose(turnVec), log.(σWidthVec)) / dot(transpose(turnVec), turnVec) / circum
end
getGrowthRatePerTurn(turnVec::Vector{Int64}, σWidthVec::Vector{Float64}, par::ParDBRF) = getGrowthRatePerTurn(turnVec,σWidthVec,par.circum)


function movingAverage(xvec, nwidth)
    ndata = length(xvec)
    resvec = zeros(ndata)
    for i in nwidth:ndata
        res = 0.0
        for j in (i-nwidth+1):i
            res += xvec[j]
        end
        resvec[i] = res/nwidth
    end
    resvec
end

function removeDataParticleLoss(turnVec, σWidthVec, particlenumVec)
    iniParticleNum = particlenumVec[1]
    for i in 1:(length(turnVec)-1)
        if particlenumVec[i+1] < iniParticleNum
            return turnVec[1:i], σWidthVec[1:i], particlenumVec[1:i]
        end
    end
    turnVec, σWidthVec
end

function removeDataAfterExpGrowth(turnVec, σWidthVec)
    σWidthVec = σWidthVec ./ σWidthVec[1]
    movingAverageVec20 = movingAverage(σWidthVec, 20)
    movingAverageVec100 = movingAverage(σWidthVec, 100)
    for i in 1000:length(turnVec)
        if (movingAverageVec20[i]>10000) && (movingAverageVec20[i-1] < movingAverageVec100[i-1]) && (movingAverageVec20[i] < movingAverageVec100[i])
            return turnVec[1:i-1], σWidthVec[1:i-1]
        end
    end
    turnVec, σWidthVec
end

function selectDataIncrease(turnVec, σWidthVec)
    suc = false
    for i in 2:length(turnVec)
        nothing
    end
end