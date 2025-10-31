"""
Random sampling between `xmin`, `xmax`.
"""
function randRange(n, xmin, xmax)
    r = 0.5*(xmax - xmin)
    c = 0.5*(xmax + xmin)
    (rand(n).-0.5).*2r .+ c
end
function randRange(xmin, xmax)
    r = 0.5*(xmax - xmin)
    c = 0.5*(xmax + xmin)
    (rand()-0.5)*2r + c
end