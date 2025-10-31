"""
删除过于接近的重复数据
"""
function remove_approx_data!(lst::Vector; atol=1e-10)
    boolvec = zeros(Bool, length(lst))
    @inbounds for i in length(lst):-1:1
        @inbounds for j in i-1:-1:1
            if isapprox(lst[i], lst[j], atol=atol)
                boolvec[i] = true
                break
            end
        end
    end
    deleteat!(lst, boolvec)
end
function remove_approx_data(lst::Vector; atol=1e-10)
    tmplst = lst[:]
    remove_approx_data!(tmplst; atol=atol)
end

"""
移除超出范围的数据
"""
function remove_outranged_data!(lst::Vector, lb::Float64, ub::Float64)
    boolvec = zeros(Bool, length(lst))
    @inbounds for i in length(lst):-1:1
        if (lst[i] < lb) || (lst[i] > ub)
            boolvec[i] = true
        end
    end
    deleteat!(lst, boolvec)
end
function remove_outranged_data(lst::Vector, lb::Float64, ub::Float64)
    tmplst = lst[:]
    boolvec = zeros(Bool, length(tmplst))
    @inbounds for i in length(tmplst):-1:1
        if (lst[i] < lb) || (lst[i] > ub)
            boolvec[i] = true
        end
    end
    deleteat!(tmplst, boolvec)
end



"""
移除超出范围且重复的数据
"""
function remove_repeated_and_outranged_data!(lst::Vector, lb::Float64, ub::Float64; atol=1e-10)
    boolvec = zeros(Bool, length(lst))
    @inbounds for i in length(lst):-1:1
        if (lst[i] < lb) || (lst[i] > ub)
            boolvec[i] = true
            continue
        end

        @inbounds for j in i-1:-1:1
            if isapprox(lst[i], lst[j], atol=atol)
                boolvec[i] = true
                break
            end
        end
    end
    deleteat!(lst, boolvec)
end
function remove_repeated_and_outranged_data(lst::Vector, lb::Float64, ub::Float64; atol=1e-10)
    tmplst = lst[:]
    remove_repeated_and_outranged_data!(tmplst, lb, ub; atol=atol)
end