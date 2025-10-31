function discrete_integrate(xvec::Vector, yvec::Vector)
    s = zero(promote_type(eltype(xvec), eltype(yvec)))
    @tturbo for i in 1:length(xvec)-1
        s += (yvec[i] + yvec[i+1]) * (xvec[i+1] - xvec[i]) / 2
    end
    s
end



function discrete_integrate(xvec::Vector, y1vec::Vector, y2vec::Vector)
    s = zero(promote_type(eltype(xvec), eltype(y1vec), eltype(y2vec)))
    @tturbo for i in 1:length(xvec)-1
        s += (y1vec[i+1] * y2vec[i+1] + y1vec[i] * y2vec[i]) * (xvec[i+1] - xvec[i]) / 2
    end
    s
end

function discrete_integrate(xvec::Vector, y1vec::Vector, y2vec::Vector, y3vec::Vector)
    s = zero(promote_type(eltype(xvec), eltype(y1vec), eltype(y2vec), eltype(y3vec)))
    @tturbo for i in 1:length(xvec)-1
        s += (y1vec[i+1] * y2vec[i+1] * y3vec[i+1] + y1vec[i] * y2vec[i] * y3vec[i]) * (xvec[i+1] - xvec[i]) / 2
    end
    s
end





"""
对二维矩阵进行离散积分。

Param
---
`x1vec`: 需要积分的第一个变量的`Vector`;

`x2vec`: 需要积分的第二个变量的`Vector`;

`fmat`: 矩阵, 第一个维度为`x1vec`， 第二个维度为`x2vec`.

Return
---
二维矢量积分的结果。
"""
function discrete_integrate(x1vec::Vector, x2vec::Vector, fmat::Matrix)
    nrow = size(fmat, 1)
    intVecEachRow=zeros(typeof(x1vec[1]*fmat[1]*x2vec[1]), nrow)
    @inbounds for rowidx in 1:nrow
        intVecEachRow[rowidx] = discrete_integrate(x2vec, fmat[rowidx,:])
    end
    discrete_integrate(x1vec, intVecEachRow)
end




"""
对二维矩阵进行带权重的离散积分。

Param
---
`x1vec`: 需要积分的第一个变量的`Vector`;

`x2vec`: 需要积分的第二个变量的`Vector`;

`w1vec`: 需要积分的第一个变量的权重`Vector`;

`w2vec`: 需要积分的第二个变量的权重`Vector`;

`fmat`: 矩阵, 第一个维度为`x1vec`， 第二个维度为`x2vec`.

Return
---
二维矢量积分的结果。
"""
function discrete_integrate(x1vec::Vector, x2vec::Vector, w1vec::Vector, w2vec::Vector, fmat::Matrix)
    discrete_integrate(x1vec, x2vec, fmat .* reshape(w2vec, 1, :) .* w1vec)
end