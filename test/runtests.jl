using Quadprog
using Base.Test

# simple test example
(x, fx, status) = quadprog([3 1; 1 4], [5, 6])
@test norm(x  - [-1.27273,-1.18182]) < 1e-5
@test norm(fx - -6.7272727272727275) < 1e-5
@test status == 1
