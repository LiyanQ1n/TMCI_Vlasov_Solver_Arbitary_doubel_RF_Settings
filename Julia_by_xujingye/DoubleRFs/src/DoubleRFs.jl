module DoubleRFs

using FFTW, Interpolations, LinearAlgebra
using IntervalArithmetic: Interval, mid, radius
using IntervalRootFinding: roots
using QuadGK: quadgk
using LoopVectorization
import CurveFit: curve_fit, ExpFit, LinearFit, Polynomial

# 临时
# import BSON: @save


export ParDBRF
# Theory
export getFPs, getBucketBounds, getNormedDensity, getNormedDensityWithRandomχ, getNormedDensityWithTaperΧ, trajectoriesForPlot
export hamiltonian, getPotenCavity, getPotenWake, getWakeVoltage, getCavityVoltage

# General Tools
export adaptiveSample, extrapolateAndResample, samplingExpStep
export discrete_integrate, matrix_derivative
export curve_fit, FractionalOrderPolynomial
export randRange
export samplingAccordingDensityFun, getDensityDistFromSamples
export remove_repeated_and_outranged_data!
# Data Analysis
export getGrowthRatePerTurn, getRMSWidth
# Simulation Tools
export getSynchrotronTunes
export one_turn_map, one_turn_map!, one_turn_map_repeat
export one_turn_map_repeat!, one_turn_map_history
# Vlasov
export getDatas, scanNb # Vlasov analysis with radial discritization


const clight = 299792458.0
const elementcharge = 1.6021766208e-19
const r0 = 2.818e-15
const restenergy = 0.51099906e6

struct ParDBRF # params of double rf
    # ===== 这些参数需要提前准备 =====
    v1::Float64
    r::Float64
    ϕs::Float64
    ϕ2s::Float64
    h1::Int64
    h::Int64
    circum::Float64
    centerenergy::Float64
    αc::Float64
    ξ::Float64
    νy::Float64
    σδ::Float64
    # ===== 下面的参数会自动计算 ======
    v2::Float64
    ωβ::Float64
    ω0::Float64
    ωξ::Float64
    γ::Float64
    η::Float64
end
@doc raw"""
    ParDBRF(v1, r, ϕs, ϕ2s, h1, h, circum, centerenergy, αc, ξ, νy, σδ)::ParDBRF

为了避免反复输入很多参数，将双高频系统和环的参数塞进统一结构体。

Parameters:

`v1`: 主腔电压，单位[V]

`r`: 高频腔和主腔电压之比，单位[1]

`ϕs`: 主腔相位，单位[rad]

`ϕ2s`: 高频腔相位，单位[rad]

`h1`: 主腔周期数，单位[1]

`h`: 高频腔周期数，单位[1]

`circum`: 轨迹周期，单位[m]

`centerenergy`: 质子中心能量，单位[eV]

`αc`: 质子中心能量对应的相对论因子

`ξ`: 质子中心能量对应的相对论因子

`νy`: 横向振荡工作点，即，单位[1]

`σδ`: 束团能散

Return:

参数结构体`ParDBRF`，函数和结构体同名，是利用了多重分派特性。
"""
function ParDBRF(v1, r, ϕs, ϕ2s, h1, h, circum, centerenergy, αc, ξ, νy, σδ)::ParDBRF
    ParDBRF(v1, r, ϕs, ϕ2s, h1, h, circum, centerenergy, αc, ξ, νy, σδ, v1*r, νy * 2 * π * 299792458.0 / circum, 2 * π * 299792458.0 / circum, νy * 2 * π * 299792458.0 / circum * ξ / (αc - 1 / (centerenergy / 0.511e6)^2), centerenergy / 0.511e6, αc - 1 / (centerenergy / 0.511e6)^2)
end
Base.broadcastable(x::ParDBRF) = Ref(x)


# Theory
include("./Hamiltonian.jl")
include("./NormedDensity.jl")
include("./Potential.jl")
include("./FixedPoints.jl")
include("./Voltage.jl")
# Action Angle Vlasov
include("./ActionAngleVlasovSolver/Sample/AngleVariable.jl")
include("./ActionAngleVlasovSolver/Sample/AngularSampling.jl")
include("./ActionAngleVlasovSolver/Sample/δJθ.jl")
include("./ActionAngleVlasovSolver/TMCI_InteractionMatrix.jl")
include("./ActionAngleVlasovSolver/TMCI_Main.jl")
# (GT)General Tools
include("./GeneralTools/GT_Sum.jl")
include("./GeneralTools/GT_Pad.jl")
include("./GeneralTools/GT_Int.jl")
include("./GeneralTools/GT_MatrixDerivative.jl")
include("./GeneralTools/GT_Random.jl")
include("./GeneralTools/GT_Sample.jl")
include("./GeneralTools/GT_CurveFit.jl")
include("./GeneralTools/GT_RemoveDataList.jl")
include("./GeneralTools/GT_ZeroPoint.jl")
# (ST)Simulation Tools
# Some tools for problems without mathematic solution.
include("./SimulationTools/ST_MappingEquation.jl")
include("./SimulationTools/ST_SynchroTuneSpread.jl")
# (DA)Data Analysis
include("./DataAnalysis/DA_GrowthRate.jl")
include("./DataAnalysis/DA_RMS.jl")
# 


"""
返回两个矩阵。第一个为zmat，第二个为δmat。
"""
function trajectoriesForPlot(max_iterate_time, par::ParDBRF)
    sfps, ufps, _ = getFPs(par)
    if (length(sfps) == 2) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, sfps[1] + 1e-2, (sfps[1] + ufps[2]) / 2, ufps[2] - 1e-6, ufps[2] + 1e-6, (sfps[2] + ufps[2]) / 2, sfps[2] + 1e-2]
    elseif (length(sfps) == 2) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6, ufps[1] + 1e-6, (ufps[1] + sfps[2]) / 2, sfps[2] + 1e-2, ufps[2] - 1e-6]
    elseif (length(sfps) == 1) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, (ufps[1] + sfps[1]) / 2, sfps[1] + 1e-2]
    elseif (length(sfps) == 1) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6]
    end
    δvec = zeros(length(zvec))
    one_turn_map_history(zvec, δvec, max_iterate_time, par)
end
function trajectoriesForPlot(zvec0::Vector{Float64}, δvec0::Vector{Float64}, max_iterate_time, par::ParDBRF)
    sfps, ufps, _ = getFPs(par)
    if (length(sfps) == 2) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, sfps[1] + 1e-2, (sfps[1] + ufps[2]) / 2, ufps[2] - 1e-6, ufps[2] + 1e-6, (sfps[2] + ufps[2]) / 2, sfps[2] + 1e-2]
    elseif (length(sfps) == 2) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6, ufps[1] + 1e-6, (ufps[1] + sfps[2]) / 2, sfps[2] + 1e-2, ufps[2] - 1e-6]
    elseif (length(sfps) == 1) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, (ufps[1] + sfps[1]) / 2, sfps[1] + 1e-2]
    elseif (length(sfps) == 1) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6]
    end
    δvec = zeros(length(zvec))
    append!(zvec, zvec0)
    append!(δvec, δvec0)
    one_turn_map_history(zvec, δvec, max_iterate_time, par)
end
"""
历史遗留问题了。这是为了照顾到我的部分旧代码保留的。因为在python中定义一个ParDBRF传入并不是很方便。
"""
function trajectoriesForPlot(zvec0::Vector{Float64}, δvec0::Vector{Float64}, ϕs, ϕ2s, v1, r, h1, h, circum, centerenergy, αc, max_iterate_time)
    par = ParDBRF(v1, r, ϕs, ϕ2s, h1, h, circum, centerenergy, αc, 0.0, 100.0, 1e-4)
    sfps, ufps, _ = getFPs(par)
    if (length(sfps) == 2) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, sfps[1] + 1e-2, (sfps[1] + ufps[2]) / 2, ufps[2] - 1e-6, ufps[2] + 1e-6, (sfps[2] + ufps[2]) / 2, sfps[2] + 1e-2]
    elseif (length(sfps) == 2) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6, ufps[1] + 1e-6, (ufps[1] + sfps[2]) / 2, sfps[2] + 1e-2, ufps[2] - 1e-6]
    elseif (length(sfps) == 1) && (par.η < 0)
        zvec = [ufps[1] + 1e-6, (ufps[1] + sfps[1]) / 2, sfps[1] + 1e-2]
    elseif (length(sfps) == 1) && (par.η > 0)
        zvec = [sfps[1] + 1e-2, (sfps[1] + ufps[1]) / 2, ufps[1] - 1e-6]
    end
    δvec = zeros(length(zvec))
    append!(zvec, zvec0)
    append!(δvec, δvec0)
    one_turn_map_history(zvec, δvec, max_iterate_time, par)
end



function _cache_normedVoltage(ϕs, ϕ2s, r, h1, h, circum)
    k1 = 2 * h1 * pi / circum
    k2 = 2 * h1 * pi / circum * h
    c1 = sin(ϕs) + r * sin(ϕ2s)
    (; k1, k2, c1)
end
function normed_voltage_with_cache(z, ϕs, ϕ2s, r, c)
    sin(ϕs - c.k1 * z) + r * sin(ϕ2s - c.k2 * z) - c.c1
end
"""
`k1`: wave number of primary cavity. `2 * h1 * pi / circum`;

`k2`: wave number of harmonic cavity. `2 * h2 * pi / circum`;

`c1`: normed voltage seen by reference particle. `sin(ϕs) + r * sin(ϕ2s)`.
"""
function normed_voltage_with_cache(z, ϕs, ϕ2s, r, k1, k2, c1)
    sin(ϕs - k1 * z) + r * sin(ϕ2s - k2 * z) - c1
end

function normedVoltage(z::Number, ϕs, ϕ2s, r, h1, h, circum)
    c = _cache_normedVoltage(ϕs, ϕ2s, r, h1, h, circum)
    normed_voltage_with_cache(z, ϕs, ϕ2s, r, c)
end
function normedVoltage(z::Number, par::ParDBRF)
    normedVoltage(z, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end
function normedVoltage(zvec::Vector, ϕs, ϕ2s, r, h1, h, circum)
    n = length(zvec)
    resvec = zeros(eltype(zvec), n)
    c = _cache_normedVoltage(ϕs, ϕ2s, r, h1, h, circum)
    for i = 1:n
        resvec[i] = normed_voltage_with_cache(zvec[i], ϕs, ϕ2s, r, c)
    end
    resvec
end




@doc raw"""
    ∇normedVoltage(z, ϕs::Float64, ϕ2s::Float64, r::Float64, h1::Int64, h::Int64, circum::Float64)::Float64

归一化电压的1阶导数。
"""
function ∇normedVoltage(z, ϕs::Float64, ϕ2s::Float64, r::Float64, h1::Int64, h::Int64, circum::Float64)::Float64
    res = cos(ϕs - 2 * h1 * π * z / circum) + r * h * cos(ϕ2s - 2 * h1 * h * π * z / circum)
    -2 * h1 * π / circum * res
end
@doc raw"""
    ∇normedVoltage(z, par::ParDBRF)::Float64

归一化电压的1阶导数。
"""
function ∇normedVoltage(z, par::ParDBRF)::Float64
    ∇normedVoltage(z, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end



"""
2nd derivative of normalized voltage.
"""
function ∇²normedVoltage(z, ϕs::Float64, ϕ2s::Float64, r::Float64, h1::Int64, h::Int64, circum::Float64)
    res = sin(ϕs - 2 * h1 * π * z / circum) + r * h^2 * sin(ϕ2s - 2 * h1 * h * π * z / circum)
    -(2 * h1 * π / circum)^2 * res
end
function ∇²normedVoltage(z, par::ParDBRF)
    ∇²normedVoltage(z, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end



function _cache_normed_potential(r, ϕs, ϕ2s, h1, h, circum)
    k1 = 2π * h1 / circum
    k2 = 2π * h1 * h / circum
    c1 = (sin(ϕs) + r * sin(ϕ2s)) * k1
    c2 = (r==0.) ? 0. : r/h     # 避免不考虑谐波腔时，h=0，此时h作为分母会报错
    c3 = cos(ϕs) + c2 * cos(ϕ2s)
    (; k1, k2, c1, c2, c3)
end
function normedPotential_with_cache(z, ϕs, ϕ2s, c)
    cos(ϕs - c.k1 * z) + c.c2 * cos(ϕ2s - c.k2 * z) - c.c1 * z - c.c3
end
function normedPotential(z, ϕs::Float64, ϕ2s::Float64, r::Float64, h1::Int64, h::Int64, circum::Float64)
    c = _cache_normed_potential(r, ϕs, ϕ2s, h1, h, circum)
    normedPotential_with_cache(z, ϕs, ϕ2s, c)
end
function normedPotential(zvec::Vector, ϕs::Float64, ϕ2s::Float64, r::Float64, h1::Int64, h::Int64, circum::Float64)
    c = _cache_normed_potential(r, ϕs, ϕ2s, h1, h, circum)
    n = length(zvec)
    resvec = zeros(n)
    for i=1:n
        resvec[i] = normedPotential_with_cache(zvec[i], ϕs, ϕ2s, c)
    end
    resvec
end
function normedPotential(z::Number, par::ParDBRF)
    normedPotential(z, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end
function normedPotential(zvec::Vector, par::ParDBRF)
    normedPotential(zvec, par.ϕs, par.ϕ2s, par.r, par.h1, par.h, par.circum)
end



function ∇Potential(z, v1::Float64, centerenergy::Float64, circum::Float64, par::ParDBRF)
    -v1 / centerenergy / circum * normedVoltage(z, par)
end
function ∇Potential(z, par::ParDBRF)
    -par.v1 / par.centerenergy / par.circum * normedVoltage(z, par)
end




end # module
