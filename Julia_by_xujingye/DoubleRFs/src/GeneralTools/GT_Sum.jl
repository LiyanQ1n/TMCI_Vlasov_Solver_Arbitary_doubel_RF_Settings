function kahansum(xvec::Vector)
    # 更精确的求和算法
    # 参考(https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
    s = zero(eltype(xvec))
    c = zero(eltype(xvec))
    for i in eachindex(xvec)
        y=xvec[i]-c
        t=s+y
        c=t-s-y
        s=t
    end
    s
end