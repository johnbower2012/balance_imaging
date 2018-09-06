#ifndef __INCLUDE_CRANDOM_H__
#define __INCLUDE_CRANDOM_H__

typedef double FourVector[4];

using namespace std;

class CRandom{
 public:
  //! returns from zero to 1
  double ran(void);
  //! returns integer from zero to imax-1
  unsigned long int iran(unsigned long int imax); 
  //! weights with exp(-x^2/2)
  double gauss(void);
  void gauss2(double *g1,double *g2);
  //! weights with exp(-x)
  double  ran_exp(void);
	CRandom(int seed);
	void reset(int seed);
	void generate_boltzmann(double mass,double T,FourVector &p);
	//void generate_boltzmann(double mass,double T,double *p);
	void generate_boltzmann_alt(double mass,double T,FourVector &p);
	//void generate_boltzmann_alt(double mass,double T,double *p);
	int GetNPoissonian(double eta);
private:
 	boost::random::mt19937 *generator;
 	boost::random::uniform_real_distribution<> random_uniform{0.0,1.0};
	boost::normal_distribution<double> random_normal{0.0,1.0};
	//boost::poisson_distribution<int> random_poisson;
};

#endif
