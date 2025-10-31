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

using DoubleRFs, CSV, DataFrames, Interpolations, BenchmarkTools, Plots, QuadGK


function gaussianρ(z, σz)
    1/sqrt(2*π)/σz*exp.(-0.5 .* (z./σz).^2)
end

const σz = 0.00252
const Nb = 5e-9/1.6e-19
const Δz = 0.01σz
const v1 = 3639505.3705870000
const ϕs = 1.9203327420
const v2 = 659030.4556070120
const ϕ2s = 5.3423283818
const h1 = 756
const h = 3
circum = 1360.4
centerenergy = 6e9

zvec = collect(-0.2:Δz:0.2)
ρvec = gaussianρ(zvec, σz)

cavvolvec = getCavityVoltage.(zvec, v1, v2/v1, ϕs, ϕ2s, h1, h, circum)
cavpotenvec = getPotenCavity.(zvec, ϕs, ϕ2s, v1, v2/v1, h1, h, circum, centerenergy)
cavpotenvecbypotenwake = getPotenWake(zvec, cavvolvec, Δz, centerenergy, circum)

gr()
plt=plot(zvec, cavpotenvec .- maximum(cavpotenvec), label="Standard Cavity Potential")
plot!(plt, zvec, cavpotenvecbypotenwake, label="Cavity Potential from Voltage")
display(plt)