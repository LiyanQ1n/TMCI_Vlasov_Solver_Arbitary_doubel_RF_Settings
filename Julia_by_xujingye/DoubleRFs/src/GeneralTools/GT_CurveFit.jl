struct FractionalOrderPolynomial
    coeffs::Vector  # 系数列表
    orders::Vector  # 阶数列表
end


function (p::FractionalOrderPolynomial)(x)
    v = 0.0
    for i = 1:length(p.coeffs)
        v += p.coeffs[i] * x^p.orders[i]
    end
    return v
end


"""
对于正整数阶/分数阶多项式求多项式。

Parameters
---

`xvec`, `yvec`: `x` and `y` of one dimensional sample points.

`ordervec`: 整数阶/分数阶的矢量。

Return
---
与`ordervec`对应的阶数的系数。


# The speed is slower than `Interpolations`

# Warning. Something wrong?
"""
function curve_fit(xvec::Vector{BigFloat}, yvec::Vector{BigFloat}, ordervec::Union{Vector{Int64},Vector{Float64}})
    # Data matrix
    # [x_1^0 x_1^1 ... x_1^n]
    # [x_2^0 x_2^1 ... x_2^n]
    # ...
    # [x_N^0 x_N^1 ... x_N^n]
    X = zeros(BigFloat, size(xvec, 1), size(ordervec,1))
    for col in eachindex(ordervec)
        X[:, col] .= xvec .^ ordervec[col]
    end
    factors=inv(transpose(X) * X .+ diagm([1e-200 for i in eachindex(ordervec)])) * transpose(X) * yvec

    FractionalOrderPolynomial(factors, ordervec)
end

function curve_fit(xvec::Vector{Float64}, yvec::Vector{Float64}, ordervec::Union{Vector{Int64},Vector{Float64}})
    curve_fit(big.(xvec), big.(yvec), ordervec)
end