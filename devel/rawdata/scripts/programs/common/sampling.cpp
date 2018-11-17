#include "sampling.h"

arma::mat construct_latinhypercube_sampling(int samples, arma::mat range){
  int parameters = range.n_rows;
  arma::mat hypercube = arma::zeros<arma::mat>(samples,parameters);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,samples-1,samples);
  for(int i=0;i<parameters;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }
  return hypercube;
}
arma::mat construct_lhp_plus_chi_constraint(int samples, int ab, arma::mat range){
  int parameters = range.n_rows;
  int nmax = parameters/ab - 2;
  CCosh dist(nmax);
  Eigen::VectorXd G = Eigen::VectorXd::Zero(ab);

  arma::mat hypercube = arma::zeros<arma::mat>(samples,parameters);
  arma::vec hyperlist = arma::linspace<arma::vec>(0,samples-1,samples);
  arma::mat chi = {{ 7.2728e+02, -2.2654e+02, -8.3627e+01},
                   {-2.2654e+02,  7.2395e+02, -8.2073e+01},
                   {-8.3627e+01, -8.2073e+01,  3.0778e+02}};
  arma::vec Chi(4);
  Chi(0) = chi(0,0);
  Chi(1) = chi(0,1);
  Chi(2) = chi(0,2);
  Chi(3) = chi(2,2);
  for(int i=0;i<parameters;i++){
    hyperlist=shuffle(hyperlist);
    hypercube.col(i) = hyperlist;
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    hyperlist = hypercube.col(i);
    hyperlist = init + dx*hyperlist;
    hypercube.col(i) = hyperlist;
  }

  for(int i=0;i<samples;i++){
    for(int j=0;j<ab;j++){
      hypercube(i,j*(nmax+2)+1) = 1.0/hypercube(i,j*(nmax+2));
      for(int k=1;k<nmax+1;k++){
        hypercube(i,j*(nmax+2)+1) -= hypercube(i,j*(nmax+2)+k+1)*dist.Z(k);
      }
      hypercube(i,j*(nmax+2)+1) /= (dist.Z(0));
    }
  }

  /*                                                                                                       
  for(int i=0;i<samples;i++){                                                                              
    G = Eigen::VectorXd::Zero(ab);                                                                         
    for(int j=0;j<ab;j++){                                                                                 
      for(int k=0;k<nmax+1;k++){                                                                           
        G(j) += hypercube(i,j*(nmax+2)+k+1)*dist.Z(k);                                                     
      }                                                                                                    
      for(int k=0;k<nmax+1;k++){                                                                           
        hypercube(i,j*(nmax+2)+k+1) = hypercube(i,j*(nmax+2)+k+1)*Chi(j)/G(j);                             
      }                                                                                                    
    }                                                                                                      
  }                                                                                                        
  */
  return hypercube;
}
