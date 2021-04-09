# Block factor-width-two and Decomposed cone programs
This repository contains a set of MATLAB scripts that reformulate standard SDPs using structured subsets of the positive semidefinite (PSD) cone. We consider two general classes of structured subsets of the PSD cone

* Block factor-width two matrices (including diagonally dominant and scaled-diagonally dominant matrices); see our paper: [Block Factor-width-two Matrices and Their Applications to Semidefinite and Sum-of-squares Optimization](https://arxiv.org/abs/1909.11076)
* Decomposed structured subsets (a combination of chordal decomposition and factor-width decomposition); see our paper: [Decomposed Structured Subsets for Semidefinite and Sum-of-Squares Optimization](https://arxiv.org/abs/1911.12859)

## Block factor-width-two program

The function factorwidth.m approximates an SDP in the standard primal vectorized form using block factor-width-two matrices

                minimize 	c'x					
	  (1)   subject to	Ax = b,					
	                 	x \in K				

where the conic constraint `x \in K` is cartesian products of the following cones:

* R^n (free variables)
* Non-negative orthant
* Second-order cone
* Positive semidefinite cone

Only the PSD cones are approximated.    

The main function _factorwidth.m_ is called with the syntax

	>> [Anew, bnew, cnew, Knew, info] = factorwidth(A,b,c,K,opts);
	
Input data
* _A, b, c, K_ are SDP data in seudmi form
* _opts.bfw_     1 or 0,  block factor-width-two decomposition
* _opts.nop_     integer, number of blocks in the partion alpha
* _opts.size_    alternative to nop, number of entries in each block
* _opts.socp_    1 or 0,  reformualte 2 by 2 PSD cone with a second-order cone
* _opts.dual_    1 or 0, whether this should be dual or primal block
                    factorwidth two cone

The output data _Anew, bnew, cnew, Knew_ are new SDP data in sedumi form, which can be passed to SeDuMi directly. See _test_sedumi.m_ and _test_mosek.m_ for examples.


## Decomposed Structured Subsets

The function decomposed_subset.m approximates an SDP in the standard primal vectorized form using general structured subsets. 

The arguments to decomposed_subset includes a parameter 'cone' for the approximating subset of the PSD cone. Allowable values are 'dd', 'sdd', 'psd', or an integer k (block factor-width 2 with block size k). 'cone' may also be a cell indicating which approximating cone should be used for each PSD cone (in K.s). 

decomposed_recover.m returns the optimum matrix solution X^* given the approximated problem from decomposed_subset.m

basis_change.m performs the Change of Basis algorithm given a model and an initial feasible point. model.basis is the Cholesky basis used for transformation

## Related publications
Details can be found in the following papers:
1. Zheng, Y. Sootla, A., & Papachristodoulou, A. (2019). [Block Factor-width-two Matrices and Their Applications to Semidefinite and Sum-of-squares Optimization](https://arxiv.org/abs/1909.11076), under review.
2. Sootla, A., Zheng, Y., & Papachristodoulou, A. (2019). [Block factor-width-two matrices in semidefinite programming](https://arxiv.org/abs/1903.04938). In 2019 18th European Control Conference (ECC) (pp. 1981-1986). IEEE.
3. Miller, J., Zheng, Y., Sznaier, M., Papachristodoulou, A. (2019). [Decomposed Structured Subsets for Semidefinite and Sum-of-Squares Optimization](https://arxiv.org/abs/1911.12859) arXiv preprint arXiv:1911.12859
