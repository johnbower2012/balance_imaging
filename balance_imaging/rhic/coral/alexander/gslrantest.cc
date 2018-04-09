#include "gslrandom.h"

using namespace std;

int main(){
  CRandom random(987654321);
  int i,j,n=1000000;
  int ix;
  double x,x2bar=0.0,xbar=0.0;
  //printf("How many? : ");
  //scanf("%d",&n);

  for(i=0;i<n;i++){
    //ix=random.iran(100);
    x=random.gauss();
    xbar+=x;
    x2bar+=x*x;
    //printf("x=%10.7f\n",x);
  }
  xbar=x/double(n);
  x2bar=x2bar/double(n);
  printf("xbar=%10.7f =? 0, x2bar=%g =? 1.0\n",xbar,x2bar);
  
  return 0;
}

#include "gslrandom.cc"
