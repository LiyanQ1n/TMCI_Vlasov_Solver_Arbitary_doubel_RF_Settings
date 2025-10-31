computername = ENV["COMPUTERNAME"]
if computername == "家里电脑"
    push!(LOAD_PATH, "E:/Documents/JuliaProgram/MyModules/")
elseif computername == "宿舍电脑"
    push!(LOAD_PATH, "F:/Documents/JuliaProgramme/MyModules")
elseif computername=="DESKTOP-FDIN0I7"
    push!(LOAD_PATH, "D:/Documents/Julia/MyModules")
elseif computername=="TINA"
    push!(LOAD_PATH, "D:/JuliaProgramme/MyModules")
end

using DoubleRFs, CSV, DataFrames, Interpolations, BenchmarkTools, Plots


function wake(z::Float64)::Float64
    if z > 0.0
        return 0.0
    elseif z==0.0
        return α*Rs
    else
        return 2α*Rs*exp(α*z/clight)*(cos(ωbar*z/clight)+α/ωbar*sin(ωbar*z/clight))
    end
end

function gaussianρ(z, σz)
    1/sqrt(2*π)/σz*exp.(-0.5 .* (z./σz).^2)
end

function reimpedance(ω)
    res = 1+Q^2*(ωR/ω-ω/ωR)^2
    2*Rs/res
end

function imimpedance(ω)
    -2.0*Rs*Q*(ωR/ω-ω/ωR)/(1+Q^2*(ωR/ω-ω/ωR)^2)
end

const σz = 0.00252
const Nb = 5e-9/1.6e-19
const Δz = 0.01σz
const Rs = 891898
const ωR = 2π*56e9
const Q = 10
const α = ωR/2Q
const ωbar = sqrt(ωR^2 - α^2)
const clight = 3e8
const Δω = 1e10
ωvec = collect(Δω:Δω:5e12)
zvec = collect(-0.2:Δz:0.2)
ρvec = gaussianρ(zvec, σz)
reZvec = reimpedance.(ωvec)
imZvec = imimpedance.(ωvec)

println("正在计算Method1")
wakevol1 = getWakeVoltage(zvec, ρvec, wake, Δz, Nb)
println("正在计算Method2")
wakevol2 = getWakeVoltage(zvec, ρvec, ωvec, reZvec, imZvec, Δω, Δz, Nb)

gr()
plt = plot(zvec, wakevol1, label="Method1")
plot!(plt, zvec, wakevol2, label="Method2")
display(plt)