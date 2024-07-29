The program in this folder is to obtain the quasi periodic solution of van der Pol-Duffing (vdPD) system with two external excitation through the ETMRM (example 1 in the paper).

This folder mainly contains the following files,
Main file: Enhanced_Tikhonov_harmonic_analysis_method.m
The parameters of the algorithm can be modified and set in this document, including the order of series solution retained, the initial value of iteration, etc.

Residual calculation and sensitivity calculation subroutine:cal_residual.m
The function of this document is to calculate the residual vector and sensitivity matrix of the system.

The subroutines related to the ERSA algorithm are as follows:csvd.m; l_corner.m; l_curve.m; lcfun.m; plot_lc.m; tikhonov.m.
These files are subroutines that iteratively solve the optimal solution of the objective function through the ERSA algorithm.