#include <cstdlib>
#include <cmath>
#include <cstdio>
using namespace std;

#include "gslmatrix.h"

int main(){
  CGSLMatrix_Complex *gslmatrix;
  complex<double> ci(0.0,1.0);
  complex<double> **A;
  complex<double> *x,*y,*yy;
  int i,j,dim;
  
  printf("Enter dimension to check : ");
  scanf("%d",&dim);
  
  gslmatrix=new CGSLMatrix_Complex(dim);
  A=new complex<double>*[dim];
  for(i=0;i<dim;i++) A[i]=new complex<double>[dim];
  x=new complex<double>[dim];
  y=new complex<double>[dim];
  yy=new complex<double>[dim];
  

  printf("___ Checking Linear Eq. Solver {yy[i] =? (i+1,dim-i)} ____\n");
  for(i=0;i<dim;i++){
    y[i]=double(1+i)+ci*double(dim-i);
    for(j=0;j<dim;j++){
      A[i][j]=double(3+i)/double(4+j+2*i);
      A[i][j]+=double(dim+2-i)/double(4+i+2*j);
    }
  }
    
  gslmatrix->SolveLinearEqs(y,A,x);

  for(i=0;i<dim;i++){
    yy[i]=0.0;
    for(j=0;j<dim;j++) yy[i]+=A[i][j]*x[j];
    printf("x[%d]=(%g,%g), yy[%d]=(%g,%g)\n",i,real(x[i]),imag(x[i]),
	   i,real(y[i]),imag(y[i]));
  }

  printf("_________ Check Eigen.. Solver and Inversion ______________\n");

  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      A[i][j]=double(dim+i+j)/double(1+(i-j)*(i-j));
      A[i][j]+=ci*double(i-j)/double(1+(i-j)*(i-j));
    }
  }

  complex<double> **eigenvec;
  double *eigenval;
  eigenval=new double[dim];
  eigenvec=new complex<double>*[dim];
  for(i=0;i<dim;i++) eigenvec[i]=new complex<double>[dim];
  
  gslmatrix->EigenFind(A,eigenvec,eigenval);
  
  for(i=0;i<dim;i++){
    printf("eigenval=%10g : ",eigenval[i]);
    printf("eigenvec=(");
    for(j=0;j<dim;j++) printf("(%g,%g) ",
			      real(eigenvec[i][j]),imag(eigenvec[i][j]));
    printf(")\n");
  }

  complex<double> **Ainv;
  Ainv=new complex<double>*[dim];
  for(i=0;i<dim;i++){
    Ainv[i]=new complex<double>[dim];
  }
  

  gslmatrix->Invert(A,Ainv);
  
  printf("__________ Ainv________\n");
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++) printf("(%g,%g)",real(Ainv[i][j]),imag(Ainv[i][j]));
    printf("\n");
  }

  complex<double> **U;
  U=new complex<double>*[dim];
  for(i=0;i<dim;i++){
    U[i]=new complex<double>[dim];
  }

  printf("__________ Ainv*A________\n");
  int k;
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      U[i][j]=0.0;
      for(k=0;k<dim;k++) U[i][j]+=(Ainv[i][k])*A[k][j];
      printf("(%9.6f,%9.6f) ",real(U[i][j]),imag(U[i][j]));
    }
    printf("\n");
  }


  return 0;
}

