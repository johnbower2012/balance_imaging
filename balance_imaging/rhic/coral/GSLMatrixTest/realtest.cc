using namespace std;

#include "utilities.h"

int main(){
  CGSLMatrix_Real *gslmatrix;

  double **A;
  double *x,*y,*yy;
  int i,j,k,l,dim;

  printf("Enter dimension to check : ");
  scanf("%d",&dim);

  gslmatrix=new CGSLMatrix_Real(dim);
  A=new double*[dim];
  for(i=0;i<dim;i++) A[i]=new double[dim];
  x=new double[dim];
  y=new double[dim];
  yy=new double[dim];

  printf("___ Checking Linear Eq. Solver (yy[i] should = i+1) ____\n");
  for(i=0;i<dim;i++){
    y[i]=double(1+i);
    for(j=0;j<dim;j++) A[i][j]=double(3+i)/double(4+j+2*i);
  }

  gslmatrix->SolveLinearEqs(y,A,x);

  for(i=0;i<dim;i++){
    yy[i]=0.0;
    for(j=0;j<dim;j++) yy[i]+=A[i][j]*x[j];
    printf("x[%d]=%g, yy[%d]=%g\n",i,x[i],i,yy[i]);
  }

  printf("_________ Check Eigen.. Solver and Inversion ______________\n");

  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++) A[i][j]=double(dim+i+j)/double(1+(i-j)*(i-j));

  double **eigenvec;
  double *eigenval;
  eigenval=new double[dim];
  eigenvec=new double*[dim];
  for(i=0;i<dim;i++) eigenvec[i]=new double[dim];
  
  gslmatrix->EigenFind(A,eigenvec,eigenval);
  
  for(i=0;i<dim;i++){
    printf("eigenval=%10g : ",eigenval[i]);
    printf("eigenvec=(");
    for(j=0;j<dim;j++) printf("%10g ",eigenvec[i][j]);
    printf(")\n");
  }

  printf("______ Check that Udagger*A*U is diagonalized ________\n");
  double **D;
  D=new double*[dim];
  for(i=0;i<dim;i++) D[i]=new double[dim];
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      D[i][j]=0.0;
      for(k=0;k<dim;k++){
	for(l=0;l<dim;l++){
	  D[i][j]+=eigenvec[k][i]*A[k][l]*eigenvec[l][j];
	}
      }
      printf("%8.5f ",D[i][j]);
    }
    printf("\n");
  }

  printf("______ Check Matrix Inversion ______\n");

  double **Ainv;
  Ainv=new double*[dim];
  for(i=0;i<dim;i++) Ainv[i]=new double[dim];
  gslmatrix->Invert(A,Ainv);
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++) printf("%9.2e ",Ainv[i][j]);
    printf("\n");
  }

  double **Unity;
  Unity=new double*[dim];
  for(i=0;i<dim;i++){
    Unity[i]=new double[dim];
  }

  printf("__________ Ainv*A________\n");
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      Unity[i][j]=0.0;
      for(k=0;k<dim;k++) Unity[i][j]+=Ainv[i][k]*A[k][j];
      printf("%9.5f ",Unity[i][j]);
    }
    printf("\n");
  }

  return 0;
}

