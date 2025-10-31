"""
对矢量index求离散导数
"""
function vector_1st_derivative(v::Vector{Float64})
    n = length(v)
    dv = zeros(n)
    dv[1] = v[2] - v[1]
    dv[end] = v[end] - v[end-1]
    for i in 2:(n-1)
       dv[i] = 0.5 * (v[i+1] - v[i-1])
    end
    dv
end
function vector_1st_derivative(v::Vector{ComplexF64})
    n = length(v)
    dv = zeros(ComplexF64, n)
    dv[1] = v[2] - v[1]
    dv[end] = v[end] - v[end-1]
    for i in 2:(n-1)
       dv[i] = 0.5 * (v[i+1] - v[i-1])
    end
    dv
end

"""
对矩阵第一个index求离散导数
"""
function matrix_1st_derivative(mat::Matrix{Float64})
    nrow, ncol = size(mat)
    dmat = zeros(Float64, nrow, ncol)
    dmat[1,:] .= mat[2,:] .- mat[1,:]
    dmat[end,:] .= mat[end,:] .- mat[end-1,:]
    for row in 2:(nrow-1)
        dmat[row,:] .= 0.5 .* (mat[row+1,:] .- mat[row-1,:])
    end
    dmat
end
function matrix_1st_derivative(mat::Matrix{ComplexF64})
    nrow, ncol = size(mat)
    dmat = zeros(ComplexF64, nrow, ncol)
    dmat[1,:] .= mat[2,:] .- mat[1,:]
    dmat[end,:] .= mat[end,:] .- mat[end-1,:]
    for row in 2:(nrow-1)
        dmat[row,:] .= 0.5 .* (mat[row+1,:] .- mat[row-1,:])
    end
    dmat
end

"""
对矩阵第二个index求离散导数
"""
function matrix_2nd_derivative(mat::Matrix{Float64})
    nrow, ncol = size(mat)
    dmat = zeros(nrow, ncol)
    dmat[:,1] .= mat[:,2] .- mat[:,1]
    dmat[:,end] .= mat[:,end] .- mat[:,end-1]
    for col in 2:(ncol-1)
        dmat[:,col] .= 0.5 .* (mat[:,col+1] .- mat[:,col-1])
    end
    dmat
end
function matrix_2nd_derivative(mat::Matrix{ComplexF64})
    nrow, ncol = size(mat)
    dmat = zeros(ComplexF64, nrow, ncol)
    dmat[:,1] .= mat[:,2] .- mat[:,1]
    dmat[:,end] .= mat[:,end] .- mat[:,end-1]
    for col in 2:(ncol-1)
        dmat[:,col] .= 0.5 .* (mat[:,col+1] .- mat[:,col-1])
    end
    dmat
end

"""
derivative of matrix/vector.

Parameters
---
mat: `Matrix`.

dims: `Int`/`Tuple` of `Int`/`Vector` of `Int`.
"""
function matrix_derivative(mat::Matrix, dims)
    if (length(dims) == 1) && (dims[1]) == 1
        return matrix_1st_derivative(mat)
    elseif (length(dims) == 1) && (dims[1]) == 2
        return matrix_2nd_derivative(mat)
    elseif (length(dims) == 2)
        return matrix_1st_derivative(mat), matrix_2nd_derivative(mat)
    end
end
"""
derivative of matrix/vector.
"""
function matrix_derivative(mat::Vector)
    vector_1st_derivative(mat)
end