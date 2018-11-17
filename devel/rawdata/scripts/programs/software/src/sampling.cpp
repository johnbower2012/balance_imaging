#include "sampling.h"


//////////////////////////////////////////////////////////////////////////////

void ArmaToEigen(Eigen::MatrixXd &eigen, arma::mat arm){
  int rows=arm.n_rows,
    cols=arm.n_cols;
  eigen = Eigen::MatrixXd::Zero(rows,cols);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      eigen(i,j) = arm(i,j);
    }
  }
}
void EigenToArma(arma::mat &arm, Eigen::MatrixXd eigen){
  int rows=eigen.rows(),
    cols=eigen.cols();
  arm = arma::zeros<arma::mat>(rows,cols);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      arm(i,j) = eigen(i,j);
    }
  }
}
Eigen::MatrixXd ArmaToEigen(arma::mat &arm){
  int rows=arm.n_rows,
    cols=arm.n_cols;
  Eigen::MatrixXd eigen = Eigen::MatrixXd::Zero(rows,cols);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      eigen(i,j) = arm(i,j);
    }
  }
  return eigen;
}
arma::mat EigenToArma(Eigen::MatrixXd &eigen){
  int rows=eigen.rows(),
    cols=eigen.cols();
  arma::mat arm = arma::zeros<arma::mat>(rows,cols);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      arm(i,j) = eigen(i,j);
    }
  }
  return arm;
}


//////////////////////////////////////////////////////////////////////////////

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
  /*
  arma::mat chi = {{ 7.2728e+02, -2.2654e+02, -8.3627e+01},
                   {-2.2654e+02,  7.2395e+02, -8.2073e+01},
                   {-8.3627e+01, -8.2073e+01,  3.0778e+02}};
  arma::vec Chi = arma::zeros<arma::mat>(4);
  Chi(0) = chi(0,0);
  Chi(1) = chi(0,1);
  Chi(2) = chi(0,2);
  Chi(3) = chi(2,2);
  */

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
void LHCSampling(Eigen::MatrixXd &hypercube, int samples, int ab, Eigen::MatrixXd range){
  int parameters = range.rows();
  int nmax = parameters/ab - 2;
  CCosh dist(nmax);
  Eigen::VectorXd G = Eigen::VectorXd::Zero(ab);
  std::vector<int> hyperlist;
  hypercube = Eigen::MatrixXd::Zero(samples,parameters);

  for(int i=0;i<samples;i++){
    hyperlist.push_back(i);
  }

  for(int i=0;i<parameters;i++){
    std::random_shuffle(hyperlist.begin(),hyperlist.end());
    for(int j=0;j<samples;j++){
      hypercube(j,i) = hyperlist[j];
    }
  }

  float init,final,dx;
  for(int i=0;i<parameters;i++){
    init = range(i,0);
    final = range(i,1);
    dx = (final-init)/(samples-1);
    for(int j=0;j<samples;j++){
      hypercube(j,i) = init + dx*hypercube(j,i);
    }
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
}


//////////////////////////////////////////////////////////////////////////////
