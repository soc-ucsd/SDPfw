# Block factor-width-two cone program
Reformulation of SDPs using block factor-width two matrices

The function factorwidth.m approximates an SDP in the standard primal vectorized form using block factor-width-two matrices

		minimize 	c'x					
    subject to	Ax = b,					
		x \in K				

where the conic constraint `x \in K` is cartesian products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

Only the PSD cones are approximated.

## Related publications
Details can be found in the following papers:
1. Zheng, Y. Sootla, A., & Papachristodoulou, A. (2019). Block factor-width-two matrices and their applications in semidefinite and sum-of-squares optimization, in final preparation.
2. Sootla, A., Zheng, Y., & Papachristodoulou, A. (2019). Block factor-width-two matrices in semidefinite programming. arXiv preprint arXiv:1903.04938.
3. Miller, J., Zheng, Y., Sznaier, M., Papachristodoulou, A. (2019). Decomposed Structured Subsets for Semidefinite and Sum-of-Squares Optimization. arXiv preprint arXiv:1911.12859
