# Collocation toolbox for delay differential equations (DDEs) and differential-algebraic equations (DAEs) and their adjoints


## Unified formulation
By rewriting DDEs (delay differential equations) as DAEs (differential-algebraic equations) with piecewise continuous algebraic equations, we obtain a unified framework for the formulation of DDEs/DAEs and DDAEs (delay differential-algebraic equations).


## Collocation approximation
Under the unified framework, we have semi-explicit differential equations with algebraic equations. Lagrange interpolation with different orders is applied to differential and algebraic variables. Specifically, the order of differential variables is higher than the one of algebraic variables by one. Then the differential equations and algebraic conditions are collocated at collcation nodes. The discretization is consistent with the case of ordinary differential equation presented in [1].



## Automatic adjoint construction
The adjoint of the above differential and algebraic equations is constructed as well. Such construction fits well the staged construction paradigm discussed in [2]. 


## Future work
The current implementation does not support adaptive mesh. Users have to adjust the mesh manually. 




## References
[1] Dankowicz, H., & Schilder, F. (2013). Recipes for continuation (Vol. 11). SIAM.

[2] Li, M., & Dankowicz, H. (2018). Staged Construction of Adjoints for Constrained Optimization of Integro-Differential Boundary-Value Problems. SIAM Journal on Applied Dynamical Systems, 17(2), 1117-1151.

