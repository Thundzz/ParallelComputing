#ifndef JACOBI_H
#define JACOBI_H

void jacobi(int maxiter, double eps, double Aii, double Cx, double Cy,
		    int Nx, int N, double *RHS, double *U, double *Uold);


#endif