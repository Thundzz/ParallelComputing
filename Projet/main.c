/* Resolution du laplacien par differences finies */
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <mpi.h>

#include "cdt_bords.h"
#include "jacobi.h"
#include "conjgrad.h"
#include "tools.h"

extern int cdt_choisie;

#define UNUSED(x) (void)(x)


int main( void )
{
  FILE *Infile, *Outfile, *Timefile;
  char FileName[40], Outname[40], Timename[40];

  /* declaration des variables de discretisation du probleme */
  int Nx,Ny,N;
  double dx,dy,posx,posy;
  double Aii,Cx,Cy;
  double *U,*Uold,*RHS;
  /* declaration des variables du probleme */
  double Lx,Ly,D;
  /* variables du solveur */
  int meth, maxiter, param_cond;
  double eps;
  int i,j,k;

  int myrank, nb_procs;

  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);


  /* lecture des variables dans le fichier param.dat */
  sprintf(FileName,"param.dat");
  Infile=fopen(FileName,"r");
  fscanf(Infile,"%d%d",&Nx,&Ny);
  fscanf(Infile,"%lf%lf%lf",&Lx,&Ly,&D);
  fscanf(Infile,"%d%d%d%lf",&param_cond,&meth,&maxiter,&eps);
  fclose(Infile);
  if(param_cond <0 || param_cond >= 2)
  {
    fprintf(stderr, "Conditions aux bords: %d non supportées\n", param_cond );
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  else{
    cdt_choisie = param_cond;
  }
  // Si je suis le premier processus, j'affiche les données
  // scannée sur la sortie standard.
  if(myrank == 0)
  {
    printf("Nx,Ny= %d\t%d\n",Nx,Ny); 
    printf("Lx,Ly,D= %lf\t%lf\t%lf\n",Lx,Ly,D);
    printf("condition aux bords= %d\n", cdt_choisie);
    printf("choix de meth,maxiter,eps= %d\t%d\t%0.8f\n",meth,maxiter,eps);
  }
  /* Calcul des termes de la matrice */
  dx  = Lx/(1.0 + Nx);
  dy  = Ly/(1.0 + Ny);
  Aii = 2.0*D/(dx*dx)+ 2.0/(dy*dy);  // Terme diagonal de la matrice 
  Cx  = -1.0*D/(dx*dx);
  Cy  = -1.0*D/(dy*dy);
  N = Nx*Ny;

  /* decalration des pointeurs du probleme */
  /* +1 car je commence a 1 pour compatibilite fortran */
  U    = (double*) calloc(N+1,sizeof(double)); 
  Uold = (double*) calloc(N+1,sizeof(double));
  RHS  = (double*) calloc(N+1,sizeof(double));


  /* Remplissage du second membre de l equation */
  RightHandSide(N, Nx, Ny, dx, dy, Cx, Cy, RHS);

  double start;
  double end;

  start = MPI_Wtime();
  /* Choix du solveur pour la resolution du systeme */
  if ( meth == 1 ){     
    jacobi(maxiter,eps,Aii,Cx,Cy,Nx,N,RHS,U,Uold);}
  else if ( meth == 2 ){
    GC(maxiter,eps,Aii,Cx,Cy,Nx,N,RHS,U);} 
  else
   printf("Choix de methode non supporte");
 end = MPI_Wtime() - start;

 MPI_Allreduce(MPI_IN_PLACE, &end, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 /* ecriture de la solution dans un fichier */ 
  sprintf(Timename,"time%d", nb_procs);
  sprintf(Outname,"sol%d", myrank);
  Outfile = fopen(Outname,"w");
  Timefile = fopen(Timename,"w");

  int first_elt = (myrank * (Ny / nb_procs)) * Nx + 1;
  int last_elt = ((myrank + 1) * (Ny / nb_procs)) * Nx;
  if(myrank == nb_procs-1){
    last_elt = N;
  }

  for( i=first_elt;i<=last_elt;i++ ){
    nloc(&j,&k,i,Nx);
    posx = k*dx;
    posy = j*dy;
    fprintf(Outfile,"%lf %lf %lf\n",posx,posy,U[i]);
  }

  if(myrank == 0)
  { 
    fprintf(Timefile,"%g\n",end);
  }
  fclose(Timefile);
  fclose(Outfile);
  MPI_Finalize();
  return 0;
}