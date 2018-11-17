#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include "armadillo"
#include "pricomana.h"
#include "observables_class.h"
#include "file.h"

int main(int argc, char* argv[]){
  std::vector<std::string> model_name={"I211_J211.dat",\
				       "I2212_J2212.dat",\
				       "I321_J2212.dat",\
				       "I321_J321.dat"},

  std::vector<std::string> exp_name={"star_KK_cut.dat",\
				     "star_pipi_cut.dat",\
				     "star_ppbar_cut.dat",\
				     "star_pK_cut.dat"};

  int obs_file=3,
    files=4,
    ab=4,
    model_lines=24,
    exp_lines=16,
    model_runs=1000,
    exp_runs=1,
    param=ab*6,
    observables=obs_file*files;

  int i,j,k;

  std::string source_folder="../../calculations/data/",
    dest_folder="../../calculations",
    exp_folder="../../calculations/data/exp_run00",
    error_name,
    param_name,
    exp_name;


  /*********
	     LOAD FILE
	     CALCULATE OBS_MAT
  *********/
  std::vector<arma::mat> val_matrix;
  arma::mat obs_matrix = arma::zeros<arma::mat>(runs,observables);
  arma::vec delY_vec = arma::zeros<arma::vec>(lines),
    obs_error = arma::zeros<arma::vec>(observables);

  load_file(files, lines, runs, infilename, delY_vec, val_matrix);
  obs_matrix_moments(files,obs_file,val_matrix,delY_vec,obs_matrix);


  /*********
	     CONDUCT PCA
	     PRINT
  *********/

  arma::mat tval_matrix,zcov_matrix,eigvec_matrix,print_matrix;
  arma::vec eigval_vec,mean_vec;
  print_matrix = arma::zeros<arma::mat>(observables,observables+1);
  tval_matrix = arma::zeros<arma::mat>(observables,observables);
  zcov_matrix = arma::zeros<arma::mat>(observables,observables);
  eigvec_matrix = arma::zeros<arma::mat>(observables,observables);
  eigval_vec = arma::zeros<arma::vec>(observables);
  mean_vec = arma::zeros<arma::vec>(observables);

  arma::mat err = arma::zeros<arma::mat>(1,observables);
  arma::mat y_tilde = arma::zeros<arma::mat>(lines,observables);
  arma::mat parameters = arma::zeros<arma::mat>(runs,param);
  arma::mat exp_data = arma::zeros<arma::mat>(1,observables);

  //LOAD FILES
  error_name = exp_folder + "moments_error.dat";
  load_file(error_name,err);
  printf("%s loaded.\n",error_name.c_str());
  for(int i=0;i<observables;i++){
    obs_error(i) = err(0,i);
  }
  param_name = dest_folder + "moments_parameters.dat";
  load_file(param_name,parameters);
  printf("%s loaded.\n",param_name.c_str());

  exp_name = exp_folder + "moments_exp_data.dat";
  load_file(exp_name,exp_data);
  printf("%s loaded.\n",exp_name.c_str());
  

  //COMPUTE Z_MATRIX
  tilde_function(obs_matrix,obs_error,mean_vec,y_tilde);
  tval_matrix = calculate_tmatrix_function(obs_matrix, obs_error, eigval_vec, eigvec_matrix, mean_vec, zcov_matrix);

  //PRINT RESULTS
  print_fractional_sum(eigval_vec);
  eigval_vec.print("+++++ eigvalues +++++");
  eigvec_matrix.print("+++++ eigvectors +++++");

  //WRITE FILES
  for(i=0;i<observables;i++){
    print_matrix(i,0) = eigval_vec(i);
    for(j=0;j<observables;j++){
      print_matrix(i,j+1) = eigvec_matrix(j,i);
    }
  }

  std::string printname;
  std::string title;
  /*
  printname = dest_folder + "moments_model_data.dat";
  title = "#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,obs_matrix);
  printf("%s saved to file.\n",printname.c_str());

  printname = dest_folder + "moments_model_data_tilde.dat";
  title = "#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,y_tilde);
  printf("%s saved to file.\n",printname.c_str());
  */
  printname = dest_folder + "moments_errors.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,obs_error);
  printf("%s saved to file.\n",printname.c_str());

  printname = dest_folder + "moments_means.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,mean_vec);
  printf("%s saved to file.\n",printname.c_str());

  printname = dest_folder + "moments_model_eigvec.dat";
  title = "#colvec\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,eigvec_matrix);
  printf("%s saved to file.\n",printname.c_str());

  printname = dest_folder + "moments_model_pca.dat";
  title = "#eigval col -- eigvec rows\n#pipi ppbar pK KK\n#m0 m1 m2";
  print_file(printname,title,print_matrix);
  printf("%s saved to file.\n",printname.c_str());

  print_matrix = tval_matrix;
  printname = dest_folder + "moments_model_z.dat";
  title = "#y~*eigvec";
  print_file(printname,print_matrix);
  printf("%s saved to file.\n",printname.c_str());

  print_matrix = arma::zeros<arma::mat>(runs,param+observables);
  for(int i=0;i<param;i++){
    print_matrix.col(i) = parameters.col(i);
  }
  for(int i=0;i<observables;i++){
    print_matrix.col(i+param) = tval_matrix.col(i);
  }
  printname = dest_folder + "moments_plot.dat";
  title = "#param+model_z.s";
  print_file(printname,print_matrix);
  printf("%s saved to file.\n",printname.c_str());

  arma::mat exp_tilde = arma::zeros<arma::mat>(1,observables);
  for(int i=0;i<observables;i++){
    exp_tilde(0,i) = (exp_data(0,i) - mean_vec(i))/obs_error(i);
  }
  arma::mat exp_z = exp_tilde*eigvec_matrix;
  printname = dest_folder + "moments_exp_z.dat";
  title = "#exp_z";
  print_file(printname,exp_z);
  printf("%s saved to file.\n",printname.c_str());
  
  
  return 0;
}

