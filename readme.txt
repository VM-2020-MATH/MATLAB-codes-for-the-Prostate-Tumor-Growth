These codes are provided to solve numerically the Prostate tumor growth 
mathematical model arising in biology in two-dimensional space. 
The RBF-FD scheme is implemented 
to discretize the spatial variables, and a semi-implicit backward differential 
formula of first-order (SBDF1) is used to deal with the time variable. Besides,
the biconjugate gradient stabilized (BiCGSTAB) method is applied to solve
the new fully discrete sacheme per time step.

Authors: Vahid Mohammadi 1, Mehdi Dehghan 1, Stefano De Marchi 2

1.Department of Applied Mathematics,
  Faculty of Mathematics and Computer Sciences, 
  Amirkabir University of Technology

2. Department of Mathematics "Tullio Levi-Civita",
   University of Padua, Italy

Contacts:  v.mohammadi@aut.ac.ir, v.mohammadi.aut@gmail.com (Vahid Mohammadi), 
           mdehghan@aut.ac.ir, mdehghan.aut@gmail.com (Mehdi Dehghan), 
           demarchi@math.unipd.it (Stefano De Marchi)

The scrip Main_PCa_RBF_FD.m provides an example of solving 
the proposed mathematical model by the developed numerical 
scheme.

The folder also contains the following subroutines:

1. weights_RBF_FD.m computes the differentiation matrix

(Hint: This subroutine was written by N. Flyer and co-authors, See
N. Flyer, G.A. Barnett, and L.J. Wicker. 
Enhancing finite differences with radial ba-sis functions: 
Experiments on the Navier–Stokes equations. Journal of Computational
Physics, 316:39–62, 2016.
)

2. points.m generate the points which are used.

3.Eigenvalues_Diff.m computes the eigenvalues of the differentiation matrix

  