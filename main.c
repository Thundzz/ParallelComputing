/* Resolution du laplacien par differences finies */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



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

double f( double posx, double posy, double t)
{

  double function;

  function = sin(posx) + cos(posy);

  return(function);
}



void RightHandSide(int N, int Nx, int M, double dx, double dy, double Cx, double Cy,double *RHS)
{
  int i,j,l,k;
  double posx,posy;

  M = N/Nx ; /* # de lignes */

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
	RHS[i] = RHS[i]-f(posx,0.0,0.0)*Cy;
      }

  /* derniere ligne condition de bord du haut */
  l = 1;
  for( i = N-Nx+1;i<=N;i++ ){ 
      nloc(&j,&k,i,Nx);
      posx = k*dx;
      posy = j*dy;
      RHS[i] =RHS[i]-f(posx,1.0,0.0)*Cy;
    }

  /* Bords droit et gauche */
    /*Ligne du bas*/
  RHS[1]  = RHS[1]  -f(0.0,dy,0.0)*Cx;
  RHS[Nx] = RHS[Nx] -f(1.0,dy,0.0)*Cx;

  /*Ligne du milieux*/
  j = 1+Nx;
  for( i = 2; i<= M-1; i++ ){
    nloc(&k,&l,j,Nx);
    RHS[j] = RHS[j] -f(0.0,k*dy,0.0)*Cx;
    RHS[j+Nx-1] = RHS[j+Nx-1] -f(1.0,k*dy,0.0)*Cx;
    j = 1 + (i)*Nx;
  }
  /*ligne du haut*/
  nloc(&k,&l,N,Nx);
  RHS[N-Nx+1] = RHS[N-Nx+1] -f(0.0,k*dy,0.0)*Cx;
  RHS[N] = RHS[N] -f(1.0,k*dy,0.0)*Cx;
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
void matvec(double Aii,double Cx,double Cy,int Nx,int m,double *Uold,double *U){
  int     i,j,k;

  /*Premier bloc*/
  U[1] = Aii*Uold[1] + Cx*Uold[2] + Cy*Uold[1+Nx];
  i = 1;
  for( j = 1;j<= Nx-2;j++ ){
    i = 1+j;
    U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cx*Uold[i+1] + Cy*Uold[i+Nx];
  } 
  U[1+(Nx-1)] = Aii*Uold[1+(Nx-1)] + Cx*Uold[1+(Nx-1)-1] + Cy*Uold[1+(Nx-1)+Nx];

  i = 1 + (Nx-1);

  /*bloc general, il y a m-2 blocs generaux */
  for( k = 1; k<=(m-2); k++){
    /*Premiere ligne*/
    i = i+1;
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
  i = i+1;
  
 /*Dernier bloc*/
  U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cy*Uold[i-Nx];
  for( j = 1;j<= Nx-2;j++){
    i = i+1;
    U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
    }
  i = i+1;
  U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
  return;
  }


/* solveur de jacobi */    
void jacobi(int maxiter, double eps, double Aii, double Cx, double Cy, int Nx, int N, double *RHS, double *U, double *Uold){
  int i,j,M,l;
  double invAii;
  double err;

  invAii = 1.0/Aii;
  j = 0;
  M = N/Nx;
  err = 100.0;

  while( (err > eps) && (j < maxiter) ){
    for(l=1;l<=N;l++){
      Uold[l]=U[l];}
    matvec(0.0,Cx,Cy,Nx,M,Uold,U);
    for(l=1;l<=N;l++){
      U[l]=(RHS[l]-U[l])*invAii;}
        
      /*calcul de l erreur*/
    err = 0.0;
    for(l=1;l<=N;l++){
      Uold[l]=U[l]-Uold[l];
      err= err + sqrt(Uold[l]*Uold[l]);
    }
    j++;
  }
printf("fin jacobi, convergence en, %d, iterations, erreur=,%lf",j,err);
}

/* solveur du gradient conjugue */
void GC(int maxiter, double eps, double Aii,double Cx,double Cy,int Nx,int N,double *RHS,double *U)
{
  int l, i, M;
  double residu, drl, dwl, alpha, beta;
  double *r, *kappa, *d, *W;
  
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
  /* variables du solveurs */
  int meth, maxiter;
  double eps;
  int i,j,k,M;

  /* lecture des variables dans le fichier param.dat */
  sprintf(FileName,"param.dat");
  Infile=fopen(FileName,"r");
  fscanf(Infile,"%d%d",&Nx,&Ny);
  fscanf(Infile,"%lf%lf%lf",&Lx,&Ly,&D);
  fscanf(Infile,"%d%d%lf",&meth,&maxiter,&eps);
  fclose(Infile);

  printf("Nx,Ny= %d\t%d\n",Nx,Ny); 
  printf("Lx,Ly,D= %lf\t%lf\t%lf\n",Lx,Ly,D);
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
  sprintf(Outname,"sol");
  Outfile = fopen(Outname,"w");

  for( i=1;i<=N;i++ ){
    nloc(&j,&k,i,Nx);
    posx = k*dx;
    posy = j*dy;
    fprintf(Outfile,"%lf %lf %lf\n",posx,posy,U[i]);
  }

  fclose(Outfile);
}

