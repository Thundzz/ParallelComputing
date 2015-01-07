/* Resolution du laplacien par differences finies */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "tools.h"

int cdt_choisie= 1;

typedef struct conditions_aux_bords{
  double (*f)(double, double, double);
  double (*g)(double, double, double);
  double (*h)(double, double, double);
} cdt_aux_bords;

double f1( double posx, double posy, double t)
{
  double function;
  function = sin(posx) + cos(posy);
  return(function);
}

double f2 ( double posx, double posy, double t)
{
  return 2*(posy - posy*posy + posx - posx*posx);
}

double func_zero(double posx, double posy, double t)
{
  return 0;
}

cdt_aux_bords cdt[] = 
{
  {f2, func_zero, func_zero},
  {f1, f1, f1}
};

double func_one(double posx, double posy, double t)
{
  return 1;
}



void nloc(int *i, int *j, int n, int Nx)
{
  int q,r;

  q = n/Nx;
  r = n - q*Nx;  
  if ( r == 0 ){
    *i = q;
    *j = Nx;
  }
  else{
    *i = 1+q;
    *j = r;
  }
  return;
}






void RightHandSide(int N, int Nx, int M, double dx, double dy, double Cx, double Cy,double *RHS)
{
  int i,j,l,k;
  double posx,posy;

  M = N/Nx ; /* # de lignes */

  double (*f)(double, double, double) = cdt[cdt_choisie].f;
  double (*g)(double, double, double) = cdt[cdt_choisie].g;
  double (*h)(double, double, double) = cdt[cdt_choisie].h;

  for( i = 1; i<=N; i++ ){
    nloc(&j,&k,i,Nx);
    posx = k*dx;
    posy = j*dy;
    RHS[i] = f(posx,posy,0.0);
  }

  /* premiere ligne condition de bord du bas */
  
  for( i = 1; i<= Nx; i++ ){
	nloc(&j,&k,i,Nx);
	posx = k*dx;
	posy = j*dy;
	RHS[i] = RHS[i]-g(posx,0.0,0.0)*Cy;
      }

  /* derniere ligne condition de bord du haut */
  l = 1;
  for( i = N-Nx+1;i<=N;i++ ){ 
      nloc(&j,&k,i,Nx);
      posx = k*dx;
      posy = j*dy;
      RHS[i] =RHS[i]-g(posx,1.0,0.0)*Cy;
    }

  /* Bords droit et gauche */
    /*Ligne du bas*/
  RHS[1]  = RHS[1]  -h(0.0,dy,0.0)*Cx;
  RHS[Nx] = RHS[Nx] -h(1.0,dy,0.0)*Cx;

  /*Ligne du milieux*/
  j = 1+Nx;
  for( i = 2; i<= M-1; i++ ){
    nloc(&k,&l,j,Nx);
    RHS[j] = RHS[j] -h(0.0,k*dy,0.0)*Cx;
    RHS[j+Nx-1] = RHS[j+Nx-1] -h(1.0,k*dy,0.0)*Cx;
    j = 1 + (i)*Nx;
  }
  /*ligne du haut*/
  nloc(&k,&l,N,Nx);
  RHS[N-Nx+1] = RHS[N-Nx+1] -h(0.0,k*dy,0.0)*Cx;
  RHS[N] = RHS[N] -h(1.0,k*dy,0.0)*Cx;
}

/* Produit matrice vecteur pour une matrice tridiagonale par bloc
! Produit matrice vecteur dans le cas ou A pentadiagonale de la forme:
!
! A = B C             matrice pentadiagonale (m*Nx,m*Nx)
!     C B C
!       C B C
!         . . .
!          . . .
!            C B
! avec 
! C = Cy Id            matrice diagonale (Nx,Nx)
! 
! B = Aii Cx           matrice tridiagonale (Nx,Nx)
!     Cx  Aii Cx
!          .   .   .
!              Cx Aii
*/
void matvec(double Aii,double Cx,double Cy,int Nx,int Ny,double *Uold,double *U){
  int     i,j,k;
  int myrank;
  int nb_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  int first_elt = (myrank * (Ny/ nb_procs)) * Nx + 1;
  int last_elt = ((myrank + 1) * (Ny / nb_procs)) * Nx;
  int my_first_line = first_elt/Nx + 1;
  int my_last_line = last_elt/Nx;

  i = first_elt;

  if(my_first_line == 1){
/*Premier bloc*/
    U[1] = Aii*Uold[1] + Cx*Uold[2] + Cy*Uold[1+Nx];
    i = 1;
    for( j = 1;j<= Nx-2;j++ ){
      i = 1+j;
      U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cx*Uold[i+1] + Cy*Uold[i+Nx];
    }
    U[1+(Nx-1)] = Aii*Uold[1+(Nx-1)] + Cx*Uold[1+(Nx-1)-1] + Cy*Uold[1+(Nx-1)+Nx];
    i = 1 + (Nx-1);
  }
/*bloc general, il y a m-2 blocs generaux */
  for( k = my_first_line; k<= MIN(Ny-1, my_last_line); k++){
if(k!=1){ /* First line already done */
/*Premiere ligne*/
    i = (k-1)*Nx+1;
    U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cy*Uold[i-Nx] + Cy*Uold[i+Nx] ;
/*ligne generale*/
    for( j = 1;j<= Nx-2;j++){
      i = i+1;
      U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cx*Uold[i+1] + Cy*Uold[i+Nx] + Cy*Uold[i-Nx];
    }
/*Derniere ligne*/
    i = i + 1;
    U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cy*Uold[i-Nx] + Cy*Uold[i+Nx];
  }
}
i = i+1;
if(my_last_line == Ny){
/*Dernier bloc*/
  U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cy*Uold[i-Nx];
  for( j = 1;j<= Nx-2;j++){
    i = i+1;
    U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
  }
  i = i+1;
  U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
}
}


/* solveur de jacobi */    
void jacobi(int maxiter, double eps, double Aii, double Cx, double Cy, int Nx, int N, double *RHS, double *U, double *Uold){
  int i,j,M,l;
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
    //Communications
    if(myrank != nb_procs - 1){
      MPI_Isend(&U[lst - Nx + 1], Nx, MPI_DOUBLE, myrank+1, 99, MPI_COMM_WORLD, &r1);
      MPI_Irecv(&U[lst+1], Nx, MPI_DOUBLE, myrank+1, 99, MPI_COMM_WORLD, &r2);
    }
    if(myrank != 0){
      MPI_Isend(&U[fst], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r2);
      MPI_Irecv(&U[fst-Nx], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r1);
    }

    if(myrank != nb_procs - 1)
      MPI_Wait(&r2, &s1);
    if(myrank != 0)
      MPI_Wait(&r1, &s2);
    
    for(l = 1;l<=N;l++){
      if( l>=MAX(fst-Nx, 1)&&  l<=MIN(lst+Nx,N))
        Uold[l]=U[l];
      else
      {
        Uold[l]=0;
        U[l] = 0;
      }
    }
    matvec(0.0,Cx,Cy,Nx,M,Uold,U);
    for(l=fst;l<=lst;l++){
      U[l]=(RHS[l]-U[l])*invAii;}
        
      /*calcul de l erreur*/
    err = 0.0;
    for(l=fst;l<=lst;l++){
      Uold[l]=U[l]-Uold[l];
      err= err + Uold[l]*Uold[l];
    }
    err = sqrt(err);
    if(j %10 == 0)
      printf("rank:%d -> %g\n", myrank, err);
    MPI_Allreduce(&err, &err_buf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    err = err_buf;
    j++;
  }

  printf("fin jacobi, convergence en, %d, iterations, erreur=,%lf\n",j,err);
}

/* solveur du gradient conjugue */
void GC(int maxiter, double eps, double Aii,double Cx,double Cy,int Nx,int N,double *RHS,double *U)
{
  int myrank, nb_procs;
  int l, i, M;
  double residu, drl, dwl, alpha, beta;
  double *r, *kappa, *d, *W;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);


  l = 0;
  M = N/Nx;
  residu = 0.0;
  r     = (double*) calloc(N+1,sizeof(double)); /* +1 car je commence a 1 pour compatibilite fortran */
  kappa = (double*) calloc(N+1,sizeof(double));
  d     = (double*) calloc(N+1,sizeof(double));
  W     = (double*) calloc(N+1,sizeof(double));

  /*initialisation Gradient Conjugue*/
  for ( i=1;i<=N;i++ ){
    kappa[i] = U[i];}
  matvec(Aii,Cx,Cy,Nx,M,kappa,r);
  
  for(i=1;i<=N;i++){
    r[i]     = r[i] - RHS[i];
    residu = residu + r[i]*r[i];
    d[i]=r[i];
   }    
  
  /* boucle du Gradient conjugue */
   while( (l<=maxiter) && (residu >= eps)){    
     matvec(Aii,Cx,Cy,Nx,M,d,W);
     drl = 0.0;
     dwl = 0.0;
     for( i=1; i<=N; i++ ){
       drl = drl + d[i]*r[i];
       dwl = dwl + d[i]*W[i];
       }
      
     alpha = drl/dwl;
     for(i=1; i<=N; i++ ){
       kappa[i] = kappa[i] - alpha*d[i];
       r[i] = r[i] - alpha*W[i];}
     beta = 0.0;
     for(i=1;i<=N;i++){
       beta = beta + (r[i]*r[i]);}
     beta = beta / residu;
     residu = 0.0;
     for( i = 1; i<= N; i++){
       d[i] = r[i] + beta*d[i];   
       residu    = residu + r[i]*r[i];  
     }
     l++;
   }
   for(i=1;i<=N;i++){
     U[i] = kappa[i]; /* copie de la solution dans U */
   }

   printf("le Gradient Conjugue a converge en, %d iteration, residu= %0.12f\n",l,residu);

}/*fin du Gradient Conjugue*/



void main( void )
{
  FILE *Infile, *Outfile;
  char FileName[40], Outname[40];

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
  int i,j,k,M;

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

  printf("Nx,Ny= %d\t%d\n",Nx,Ny); 
  printf("Lx,Ly,D= %lf\t%lf\t%lf\n",Lx,Ly,D);
  printf("condition aux bords= %d\n", cdt_choisie);
  printf("choix de meth,maxiter,eps= %d\t%d\t%0.8f\n",meth,maxiter,eps);
 
  /* Calcul des termes de la matrice */
  dx  = Lx/(1.0 + Nx);
  dy  = Ly/(1.0 + Ny);
  Aii = 2.0*D/(dx*dx)+ 2.0/(dy*dy); /* Terme diagonal de la matrice */
  Cx  = -1.0*D/(dx*dx);
  Cy  = -1.0*D/(dy*dy);
  N = Nx*Ny;

  /* decalration des pointeurs du probleme */
  U    = (double*) calloc(N+1,sizeof(double)); /* +1 car je commence a 1 pour compatibilite fortran */
  Uold = (double*) calloc(N+1,sizeof(double));
  RHS  = (double*) calloc(N+1,sizeof(double));


  /* Remplissage du second membre de l equation */
    RightHandSide(N, Nx, M, dx, dy, Cx, Cy, RHS);


    /* Choix du solveur pour la resolution du systeme */
   if ( meth == 1 ){     
      jacobi(maxiter,eps,Aii,Cx,Cy,Nx,N,RHS,U,Uold);}
   else if ( meth == 2 ){
      GC(maxiter,eps,Aii,Cx,Cy,Nx,N,RHS,U);} 
   else
     printf("Choix de methode non supporte");

/* ecriture de la solution dans un fichier */ 
  sprintf(Outname,"sol%d", myrank);
  Outfile = fopen(Outname,"w");

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
  fclose(Outfile);
  MPI_Finalize();
}