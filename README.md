# Quadprog.jl

A wrapper around Ipopt.jl to solve quadratic programming problems.
This package provides the function `quadprog` which calls the
`Ipopt.jl` library for Interior Point Nonlinear Optimization.  No
other solver is required.

[![Build Status](https://travis-ci.org/harmeling/Quadprog.jl.svg?branch=master)](https://travis-ci.org/harmeling/Quadprog.jl)

## Installation

Install via the package manager,

```
Pkg.add("Quadprog")
```

## Quadratic programming

The provided function `quadprog` solves the problem

```
min  0.5 * x' * Q * x + c' * x
s.t. A   * x <= b
     Aeq * x == beq
     lb <= x <= ub
start optimization at x0
```

where `Q` is a symmetric matrix.

Only the parameters `Q` and `c` are required.  The others are keyword
arguments that could be omitted and/or given in any order.


The full signature of `quadprog` with its defaults is

```
quadprog(Q, c;
         A   =  Array(Float64, (0, length(c))),
         b   =  Array(Float64, 0),
         Aeq =  Array(Float64, (0, length(c))),
         beq =  Array(Float64, 0),
         lb  = -Inf*ones(length(c)),
         ub  =  Inf*ones(length(c)),
         x0  =  zeros(length(c)))
```

## Example 1

To solve

```
min 0.5 * [x1 x2]' * [3 1; 1 4] * [x1 x2] + [5 6]' * x
```

without any constraints you call

```
(x, fx, status) = quadprog([3 1; 1 4], [5, 6])
```

The answer you should get is:

```
([-1.27273,-1.18182],-6.7272727272727275,1)
```

That means the solution is `x==[-1.27273,-1.18182]`, the function value at the solution
is `fx=-6.7272727272727275` and the status flag is `1`.
See [Ipopt.jl](https://github.com/JuliaOpt/Ipopt.jl) for a description of the status flags.


## Example 2

To solve the dual of a linear support vector machine (SVM), i.e.

```
max    sum(alpha) - 0.5 * alpha' * diagm(y) * X * X' * diagm(y) * alpha
s.t.   0 <= alpha <= C\\
       y' * alpha == 0
```

where `X` contains the training locations as rows and `y` contains the
training labels being either `1.0` or `-1.0`, you call

```
yX = diagm(y)*X
K  = yX * yX'
(alpha, val, flag) = quadprog(K, -ones(n),
                              lb=zeros(n), ub=C*ones(n),
                              Aeq=y', beq=zeros(1),
                              x0=zeros(n))
```

with `x0` being the initial value.
