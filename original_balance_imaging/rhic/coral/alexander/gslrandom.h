#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class CRandom{
 public:
  // returns from zero to 1
  double ran(void);
  // returns integer from zero to imax-1
  unsigned long int iran(unsigned long int imax); 
  // weights with exp(-x^2/2)
  double gauss(void);
  void gauss2(double *g1,double *g2);
  CRandom(int seed);
  void reset(int seed);
  void thermal(double m,double T,double *p);
 private:
  gsl_rng *randy;
};
