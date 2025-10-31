"""
插值之前，对一维自变量进行扩充。
"""
function pad_coor_vec_before_interp(coor_vector_to_expand::Vector{Float64})
    tempvec = zeros(Float64, length(coor_vector_to_expand)+4)
    smallest_number = coor_vector_to_expand[1]
    largest_number = coor_vector_to_expand[end]
    tempvec[1] = smallest_number - 1e10
    tempvec[2] = prevfloat(smallest_number)
    tempvec[end-1] = nextfloat(largest_number)
    tempvec[end] = largest_number + 1e10
    tempvec[3:(end-2)] = coor_vector_to_expand
    tempvec
end

"""
插值之前的处理——对一维函数值补零。
"""
function pad_funval_vec_before_interp(funval_vector_to_expand::Union{Vector{Float64},Vector{ComplexF64}})
    tempvec = zeros(eltype(funval_vector_to_expand), length(funval_vector_to_expand)+4)
    tempvec[3:(end-2)] = funval_vector_to_expand
    tempvec
end


"""
将超出边界的设置为0
"""
function changeOuterZero!(fvec::Vector{Float64}, zvec::Vector{Float64}, zL, zR)
    @inbounds for i in eachindex(zvec)
        if (zvec[i] < zL) || (zvec[i] > zR)
            fvec[i] = 0.0
        end
    end
end