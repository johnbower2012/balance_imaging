#ifndef __MCMC_CLASS_H__
#define __MCMC_CLASS_H__

#include<cmath>
#include<random>
#include<chrono>
#include<armadillo>

class MCMC{
public:
  double max_log_likelihood;
  double likelihood;

  int observable_count;
  arma::vec target_value;
  int parameter_count;
  arma::mat range;

  arma::vec position;
  arma::vec test_position;
  arma::vec widths;

  unsigned seed;
  std::default_random_engine generator;

  MCMC();

  void set_target_value(arma::vec new_target_value);
  void set_observable_count(int new_observable_count);
  void set_parameter_count(int new_parameter_count);
  void set_range(arma::mat new_range);
  void set_position(arma::vec new_position);
  void set_widths(arma::vec new_widths);

  void set_seed(int seed);
  void set_seed_clock();
  double normal();
  double uniform();

  void step();
  double get_likelihood();
  double get_likelihood(arma::vec z);
  bool decide(arma::vec z);
};


#endif
