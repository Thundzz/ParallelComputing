#include "jacobi.h"
#include "tools.h"

#include <math.h>
#include <stdio.h>
#include <mpi.h>


void jacobi(int maxiter, double eps, double Aii, double Cx, double Cy, 
            int Nx, int N, double *RHS, double *U, double *Uold){
  int j,M,l;
  int myrank, nb_procs;
  double invAii;
  double err, err_buf;
  MPI_Status s1, s2;
  MPI_Request r1, r2;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);

  invAii = 1.0/Aii;
  j = 0;
  M = N/Nx;
  err = 100.0;
  
  int fst = (myrank * (M / nb_procs)) * Nx + 1;
  int lst = ((myrank + 1) * (M / nb_procs)) * Nx;
  if(myrank == nb_procs-1){
    lst = N;
  }

  while( (err > eps) && (j < maxiter) ){
    // Communications pour que chacun
    // récupère les morceaux de U qui lui manquent.
    if(myrank != nb_procs - 1){
      MPI_Isend(&U[lst - Nx + 1], Nx, MPI_DOUBLE,
                myrank+1, 99, MPI_COMM_WORLD, &r1);
      MPI_Irecv(&U[lst+1], Nx, MPI_DOUBLE,
                myrank+1, 99, MPI_COMM_WORLD, &r2);
    }
    if(myrank != 0){
      MPI_Isend(&U[fst], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r2);
      MPI_Irecv(&U[fst-Nx], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r1);
    }

    if(myrank != nb_procs - 1)
      MPI_Wait(&r2, &s1);
    if(myrank != 0)
      MPI_Wait(&r1, &s2);
    
    for(l =MAX(fst-Nx, 1) ;l<= MIN(lst+Nx,N);l++){
        Uold[l]=U[l];
    }
    matvec(0.0,Cx,Cy,Nx,M,Uold,U);
    for(l=fst;l<=lst;l++){
      U[l]=(RHS[l]-U[l])*invAii;}
        
    // calcul de l'erreur
    err = 0.0;
    for(l=fst;l<=lst;l++){
      Uold[l]=U[l]-Uold[l];
      err= err + Uold[l]*Uold[l];
    }
    err = sqrt(err);

    MPI_Allreduce(&err, &err_buf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    err = err_buf;
    j++;
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)
  {
    printf("Fin jacobi, convergence en %d iterations erreur=%lf\n",j,err);    
  }
}