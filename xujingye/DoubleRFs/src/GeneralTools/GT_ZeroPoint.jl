function newton_iterate_onestep(x::Float64, f, Δx=1e-9)
    x - 2 * Δx * f(x) / (f(x+Δx) - f(x-Δx))
end

function get_zero_point(x::Float64, f, tol=1e-6, maxiter=2000, Δx=1e-6)
    for _ = 1:maxiter
        xold = x
        x = newton_iterate_onestep(x, f, 1e-3 * Δx)
        if abs(x - xold) < tol
            return x
        end
    end
    NaN
end
function get_zero_points(xmin, xmax, f, tol=1e-6, maxiter=2000, Δx=1e-9)
    zero_vec = Float64[]
    ypre = f(xmin)
    tmpΔx = 1000 * Δx
    for x in xmin+tmpΔx:tmpΔx:xmax
        ynow = f(x)
        if ypre * ynow <= 0
            tmpx = get_zero_point(x - 0.5tmpΔx, f, tol, maxiter, Δx)
            if !isnan(tmpx)
                push!(zero_vec, tmpx)
            end
        end
        ypre = ynow
    end
    remove_approx_data!(zero_vec; atol=Δx)
    remove_outranged_data!(zero_vec, xmin, xmax)
end

"""
在`(xmin, xmax)`范围内，求解`f(x) = f_target`

Parameters
---

`xmin`, `xmax`: 解的下限和上限。

`f`: 目标函数。

`f_target`: 目标函数的值。

`tol`: 解的精度。

`maxiter`: 最大迭代次数。

`Δx`: 搜索步长。
"""
function get_solution(xmin, xmax, f, f_target, tol=1e-6, maxiter=2000, Δx=1e-9)
    get_zero_points(xmin, xmax, x -> f(x) - f_target, tol, maxiter, Δx)
end