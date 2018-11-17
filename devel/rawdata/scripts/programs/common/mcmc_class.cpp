#include "mcmc_class.h"

/*SETUP--
***************/
MCMC::MCMC(){
  this->max_log_likelihood = -1000.0;
  set_seed(1);
}
void MCMC::set_observable_count(int new_observable_count){
  this->observable_count = new_observable_count;
  this->target_value = arma::zeros<arma::vec>(this->observable_count);
}
void MCMC::set_target_value(arma::vec new_target_value){
  this->target_value = new_target_value;
}
void MCMC::set_parameter_count(int new_parameter_count){
  this->parameter_count = new_parameter_count;
  this->range = arma::zeros<arma::mat>(this->parameter_count,2);

  this->position = arma::zeros<arma::vec>(this->parameter_count);
  this->widths = arma::zeros<arma::vec>(this->parameter_count);
  this->test_position = arma::zeros<arma::vec>(this->parameter_count);
}
void MCMC::set_range(arma::mat new_range){
  this->range = new_range;
}
void MCMC::set_position(arma::vec new_position){
  this->position = new_position;
}
void MCMC::set_widths(arma::vec new_widths){
  this->widths = new_widths;
}

/*RANDOM--
***************/
void MCMC::set_seed(int seed){
  generator.seed(seed);
}
void MCMC::set_seed_clock(){
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator.seed(seed);
}
double MCMC::normal(){
  std::normal_distribution<double> normal_dist(0.0,1.0);
  return normal_dist(generator);
}
double MCMC::uniform(){
  std::uniform_real_distribution<double> uniform_dist(0.0,1.0);
  return uniform_dist(generator);
}

/*STEP--
**************/
double MCMC::get_likelihood(arma::vec Z){
  double z=0.0,diff = 0.0;
  for(int i=0;i<observable_count;i++){
    z = Z(i) - target_value(i);
    diff += z*z;
  }
  return exp(-diff/2.0);
}

void MCMC::step(){
  double random, LH;
  for(int i=0;i<parameter_count;i++){
    test_position(i) = position(i) + normal()*widths(i);
  }
}      
bool MCMC::decide(arma::vec Z){
  double random, LH,ratio;
  for(int i=0;i<parameter_count;i++){
    if(test_position(i) < range(i,0)){
      return false;
    }
    else if(test_position(i) > range(i,1)){
      return false;
    }
  }

  LH = get_likelihood(Z);
  if(log(LH) > this->max_log_likelihood){
    this-> max_log_likelihood = log(LH);
    printf("log_LH: %f\n",this->max_log_likelihood);
  }

  ratio = LH/this->likelihood;
  if(ratio<1.0){
    random = uniform();
    if(random<ratio){
      this->position = this->test_position;
      this->likelihood = LH;
      return true;
    }else{
      return false;
    }
  } else{
    this->position = this->test_position;
    this->likelihood = LH;
    return true;
  }
}      
/*FILEWRITING--
**************/

