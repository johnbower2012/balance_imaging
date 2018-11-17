#include<iostream>
#include<iomanip>
#include<cmath>
#include<random>
#include<chrono>
#include<armadillo>
#include<string>
#include<fstream>
#include "lib/emulator_class.cpp"

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

std::ofstream ofile;

int main(int argc, char* argv[]){
  int mcmc_runs = 100;
  std::string datafile, expfile, hpfile, betafile, outfile;

  if(argc<7){
    printf("Improper usage. Please enter also 'model_file exp_file betafile hpfile mcmc_runs outfile' on sabme line.\n");
    exit(1);
  }else{
    datafile = argv[1];
    expfile = argv[2];
    betafile = argv[3];
    hpfile = argv[4];
    mcmc_runs = atoi(argv[5]);
    outfile = argv[6];
  }
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::normal_distribution<double> dist(0.0,1.0);
  std::default_random_engine generator(seed);


  double epsilon=1e-8;
  int parameters = 24;
  int observables = 12;
  int model_runs = 1000;
  int exp_runs = 1;
  int hyperparameteres = 3;

  arma::vec target_value = arma::zeros<arma::mat>(observables);
  arma::mat X = arma::zeros<arma::mat>(model_runs,parameters);
  arma::mat Y_model = arma::zeros<arma::mat>(model_runs,observables);
  arma::mat Y_exp = arma::zeros<arma::mat>(exp_runs,observables);
  arma::mat hyp = arma::zeros<arma::mat>(observables,hyperparameteres);
  arma::mat beta = arma::zeros<arma::mat>(observables,parameters+1);
  arma::mat X_s = arma::zeros<arma::mat>(1,parameters);

  load_data_file(datafile,X,Y_model);
  load_file(expfile,Y_exp);
  load_file(hpfile,hyp);
  load_beta_file(betafile,beta);
  target_value = Y_exp.row(0).t();

  arma::mat range = arma::zeros<arma::mat>(parameters,2);
  arma::vec position = arma::zeros<arma::vec>(parameters);
  arma::vec widths = arma::zeros<arma::vec>(parameters);

  for(int i=0;i<parameters;i++){
    //position(i) = fabs(dist(generator));
    if(i%6==0){
      widths(i) = (1.5-0.1)/10.0;
      range(i,0) = 0.1;
      range(i,1) = 1.5;
      position(i) = 0.8;
    }else{
      widths(i) = 1000;
      range(i,0) = -8000;
      range(i,1) = 8000;
      position(i) = 0.0;
    }      
  }


  emulator gauss(X,hyp,beta,epsilon);
  MCMC random;
  random.set_observable_count(observables);
  random.set_parameter_count(parameters);

  random.set_target_value(target_value);
  random.set_range(range);
  random.set_position(position);
  random.set_widths(widths);
  random.set_seed(1);
  //  random.set_seed_clock();

  arma::mat gauss_mat = arma::zeros<arma::mat>(1,parameters + observables*3);
  bool stepped;
  arma::vec z = arma::zeros<arma::vec>(observables);
  double fraction=0.0;
  int trace_count = (int) (0.2*mcmc_runs);
  arma::mat trace = arma::zeros<arma::mat>(trace_count,parameters);
  int j=0,k=0;
  int pts=mcmc_runs/10.0;

  arma::mat gauss_output = arma::mat(mcmc_runs,parameters + 3*observables);
  for(int i=0;i<mcmc_runs;i++){
    random.step();
    X_s = random.test_position.t();
    gauss_mat = gauss.emulate(X_s,Y_model);
    gauss_output.row(i) = gauss_mat.row(0);
    for(int l=0;l<observables;l++){
      z(l) = gauss_mat(0,parameters+observables*2+l);
    }
    stepped = random.decide(z);

    if(stepped==true){
      fraction+=1.0;
    }
    if(i%5==0){
      trace.row(j) = random.position.t();
      j++;
    }
    if(i%pts==0){
      printf("%d%% completed...\n",k*10);
      k++;
    }
  }
  printf("100%% completed...\n");
  printf("fraction: %f\n",fraction/mcmc_runs);
  int trace_final_count = 0.8*trace.n_rows;
  arma::mat trace_final = arma::zeros<arma::mat>(trace_final_count,parameters);
  for(int i=0;i<trace_final_count;i++){
    trace_final.row(trace_final_count-1-i) = trace.row(trace_count-1-i);
  }
  std::string title = "'uu_width','ud_width','us_width','ss_width'";
  write_csv_file(outfile,title,trace_final);

  write_file("gauss_output.dat",gauss_output);
  ofile.open("mcmc_stats.dat");
  ofile << "acc_ratio: " << fraction/mcmc_runs << '\n';
  ofile << "max_llh: " << random.max_log_likelihood << '\n';
  ofile.close();
  
  return 0;
}

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
