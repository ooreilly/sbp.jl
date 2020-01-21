module MMS
"""
Module for working with the method of manufactured solutions.

"""

using DataStructures

"""
Extract slices from an Array of ordered lengths n1, n2, n3. The sum of lengths
cannot exceed the array length. The dictionary `vars` is used to label each slice.

Example:
        ex = extract_vars([1,2,3,4,5], OrderedDict("a" => 2, "b" => 3))
        #ex["a"] = [1,2]
        #ex["b"] = [3,4,5]

"""
function extract_vars(solution::Array, vars::OrderedDict)
        extracted_vars = OrderedDict()

        offset = 0
        for (key, value) in vars
                extracted_vars[key] = solution[1+offset:offset+value]
                offset += value
        end

        return extracted_vars
end

"""
Return the l2-error for each variable stored in a dict. 
The l2-error is computed as \$\\|e\\| = \\sqrt{e^T H e}\$, 
where e is the error vector, and H is a norm matrix.
"""
function l2_errors(solution::Dict, mms::Dict, norms::Dict)
        err = Dict()
        for (var, value) in solution
                e = (solution[var] .- mms[var])   
                err[var] = sqrt(e' * norms[var] * e)
        end
        return err
end

function abs_errors(solution::Dict, mms::Dict)
        err = Dict()
        for (var, value) in solution
                err[var] = abs.(solution[var] - mms[var])
        end
        return err
end

"""
Return the functional error for each variable stored in a dict.  Consider the
linear functional \$ F(u) =  \\int u dx \$. The (absolute) error is defined as \$e = |F(u) -
\\sum_j (H * v)_j| \$, where H is a quadrature rule and v is the numerical
solution.
"""
function functional_errors(solution::Dict, mms::Dict, norms::Dict)
        err = Dict()
        for (var, value) in solution
                err[var] = abs(mms[var] - sum(norms[var] * solution[var]))
        end
        return err
end

"""
Compute the convergence rate for a list of errors that are computed with
subsequent grid refinement of a factor of two. 

"""
function log2_convergence_rates(err::Array)
        
    q = -log2.(err[2:end]./err[1:end-1])
    prepend!(q, Inf)
    return q
end

end
