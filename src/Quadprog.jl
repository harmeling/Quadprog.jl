module Quadprog
# written by Stefan Harmeling, 2014-11-12

export quadprog

using Ipopt  # a julia wrapper for Ipopt

# References for Ipopt:
# A. Wächter and L. T. Biegler, 
# ​On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm for Large-Scale Nonlinear Programming, 
# Mathematical Programming 106(1), pp. 25-57, 2006

function quadprog(Q, c;
                  A   =  Array(Float64, (0, length(c))),
                  b   =  Array(Float64, 0),
                  Aeq =  Array(Float64, (0, length(c))),
                  beq =  Array(Float64, 0),
                  lb  = -Inf*ones(length(c)),
                  ub  =  Inf*ones(length(c)),
                  x0  =  zeros(length(c)))

    # DESCRIPTION:
    # min  0.5 * x' * Q * x + c' * x
    # s.t. A   * x <= b
    #      Aeq * x == beq
    #      lb <= x <= ub
    # init x with x0

    if ~issym(Q)           error("Q must be symmetric") end
    n   = length(x0)       # number of variables
    m   = size(A,   1)     # number of inequality constraints
    meq = size(Aeq, 1)     # number of equality constraints
    if length(b)   != m    error("size missmatch") end
    if length(beq) != meq  error("size missmatch") end
    if size(A, 2) != size(Aeq, 2) != length(ub) != length(lb) != length(x0) != n
        error("size missmatch") 
    end
    
    # setup the Ipopt problem
    g_L = [-Inf*ones(m), beq]   # lower constraint bounds
    g_U = [b, beq]              # upper constraint bounds
    eval_f(x) = (0.5*x'*Q*x + c'*x)[1]
    eval_grad_f(x, grad_f) = begin grad_f[:] = Q*x + c end
    eval_g(x, g) = begin g[:] = [A*x, Aeq*x] end
    eval_jac_g(x, mode, rows, cols, values) = begin
        # grad_g = [A; Aeq]
        if mode == :Structure
            rows[:] = (1:(m+meq)) * ones(Int64, n)'
            cols[:] = ones(Int64, m+meq) * (1:n)'
        else
            values[:] = [A; Aeq]
        end
    end
    eval_h(x, mode, rows, cols, obj_factor, lambda, values) = begin
        if mode == :Structure
            rows[:] = [1:n] * ones(Int64, n)'
            cols[:] = ones(Int64, n) * [1:n]'
        else
            values[:] = -obj_factor * Q
        end
    end
    prob = createProblem(n, lb, ub, m+meq, g_L, g_U, (m+meq)*n, n*n,
                         eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)
    addOption(prob, "print_level", 0)
    addOption(prob, "tol", 1e-12)
    prob.x = x0
    status = solveProblem(prob)
    return (prob.x, prob.obj_val, status)
end

end # of module
