# Block factor-width-two and Decomposed cone programs
This repository contains a set of MATLAB scripts that reformulate standard SDPs using structured subsets of the positive semidefinite (PSD) cone. We consider two general classes of structured subsets of the PSD cone

* Block factor-width two matrices (including diagonally dominant and scaled-diagonally dominant matrices); see our paper: [Block Factor-width-two Matrices and Their Applications to Semidefinite and Sum-of-squares Optimization](https://arxiv.org/abs/1909.11076)
* Decomposed structured subsets (a combination of chordal decomposition and factor-width decomposition); see our paper: [Decomposed Structured Subsets for Semidefinite and Sum-of-Squares Optimization](https://arxiv.org/abs/1911.12859)

## Block factor-width-two program

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

factorwidth.m is called with the syntax

	>> [Anew, bnew, cnew, Knew, info] = factorwidth(A,b,c,K,opts);
	
Input data
* A, b, c, K are SDP data in seudmi form
* opts.bfw     1 or 0,  block factor-width-two decomposition
* opts.nop     integer, number of blocks in the partion alpha
* opts.size    alternative to nop, number of entries in each block
* opts.socp    1 or 0,  reformualte 2 by 2 PSD cone with a second-order cone
* opts.dual    1 or 0, whether this should be dual or primal block
                    factorwidth two cone

The output data Anew, bnew, cnew, Knew are new SDP data in sedumi form, which can be passed to SeDuMi directly. See test_sedumi.m and test_mosek.m for examples.


## Decomposed Structured Subsets

The function decomposed_subset.m approximates an SDP in the standard primal vectorized form using general structured subsets. 

The arguments to decomposed_subset includes a parameter 'cone' for the approximating subset of the PSD cone. Allowable values are 'dd', 'sdd', 'psd', or an integer k (block factor-width 2 with block size k). 'cone' may also be a cell indicating which approximating cone should be used for each PSD cone (in K.s). 

decomposed_recover.m returns the optimum matrix solution X^* given the approximated problem from decomposed_subset.m

basis_change.m performs the Change of Basis algorithm given a model and an initial feasible point. model.basis is the Cholesky basis used for transformation

## Related publications
Details can be found in the following papers:
1. Zheng, Y. Sootla, A., & Papachristodoulou, A. (2019). Block factor-width-two matrices and their applications in semidefinite and sum-of-squares optimization, in final preparation.
2. Sootla, A., Zheng, Y., & Papachristodoulou, A. (2019). Block factor-width-two matrices in semidefinite programming. arXiv preprint arXiv:1903.04938.
3. Miller, J., Zheng, Y., Sznaier, M., Papachristodoulou, A. (2019). Decomposed Structured Subsets for Semidefinite and Sum-of-Squares Optimization. arXiv preprint arXiv:1911.12859
