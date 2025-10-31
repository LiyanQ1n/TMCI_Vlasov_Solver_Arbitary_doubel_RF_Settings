function getδJθ(zJθ::Matrix, Es::Vector, nJ::Int, nθ::Int, par::ParDBRF)
    δJθ = zeros(nJ, nθ)
    for i in 1:nJ
        δJθ[i,:] = getδJθ(zJθ[i,:], Es[i], nθ, par)
    end
    δJθ
end
function getδJθ(zJθ::Matrix, Es::Vector, nJ::Int, nθ::Int, func_potential, η::Float64)
    δJθ = zeros(nJ, nθ)
    for i in 1:nJ
        δJθ[i,:] = getδJθ(zJθ[i,:], Es[i], nθ, func_potential, η)
    end
    δJθ
end


"""
得到单圈的δ数据
"""
function getδJθ(zθ::Vector, E, nθ::Int, par::ParDBRF)
    nθcri = div(nθ+1,2)
    resvec = zeros(nθ)

    for i in 2:nθcri-1
        resvec[i] = getz2δ(zθ[i], E, par)
    end

    for i in (nθcri+1):(nθ-1)
        resvec[i] = -resvec[nθ+1-i]
    end
    resvec
end
function getδJθ(zθ::Vector, E, nθ::Int, func_potential, η::Float64)
    nθcri = div(nθ+1,2)
    resvec = zeros(nθ)

    for i in 2:nθcri-1
        resvec[i] = getz2δ(zθ[i], E, func_potential, η)
    end

    for i in (nθcri+1):(nθ-1)
        resvec[i] = -resvec[nθ+1-i]
    end
    resvec
end