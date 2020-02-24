# Block factor-width-two and Decomposed cone programs
Reformulation of SDPs using structured subsets of the PSD cone. Includes block factor-width two and (scaled) diagonally dominant matrices.

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

The function decomposed_subset.m approximates an SDP in the standard primal vectorized form using general structured subsets. 

The arguments to decomposed_subset includes a parameter 'cone', which is the approximating cone. Allowable values are 'dd', 'sdd', 'psd', or an integer k (block factor-width 2 with block size k). 'cone' may also be a cell indicating which approximating cone should be used for each PSD cone. 

decomposed_recover.m returns the optimum matrix solution X^* given the approximated problem from decomposed_subset.m

basis_change.m performs the Change of Basis algorithm given a model and an initial feasible point. model.basis is the Cholesky basis used for transformation

## Related publications
Details can be found in the following papers:
1. Zheng, Y. Sootla, A., & Papachristodoulou, A. (2019). Block factor-width-two matrices and their applications in semidefinite and sum-of-squares optimization, in final preparation.
2. Sootla, A., Zheng, Y., & Papachristodoulou, A. (2019). Block factor-width-two matrices in semidefinite programming. arXiv preprint arXiv:1903.04938.
3. Our work