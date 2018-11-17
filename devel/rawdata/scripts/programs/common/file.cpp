#include "file.h"

void write_output(arma::mat output_mat, int param, int obs, std::string outfilename){
  //write to file as:
  //		param1, param2, ... , mean value, mean-2*std, mean+2*std, output1, output2,....
  int i, j, 
    test_points = output_mat.n_rows;
  std::ofstream ofile;
  
  ofile.open(outfilename);
  for(i=0;i<test_points;i++){
    for(j=0;j<param;j++){
      ofile << " " << output_mat(i,j);
    }
    for(j=0;j<obs;j++){
      ofile << " " << output_mat(i,param+j);
      ofile << " " << output_mat(i,param+j) + 2.0*sqrt(output_mat(i,param+obs+j));
      ofile << " " << output_mat(i,param+j) - 2.0*sqrt(output_mat(i,param+obs+j));
      ofile << " " << output_mat(i,param+2*obs+j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
void write_trainset(arma::mat X_mat, arma::mat Y_mat, std::string outfilename){
  //Write to file as:
  //		param1, param2, ... , training value
  int i, j,
    train = X_mat.n_rows,
    observables = Y_mat.n_cols,
    param = X_mat.n_cols;
  std::ofstream ofile;

  ofile.open(outfilename);
  for(i=0;i<train;i++){
    for(j=0;j<param;j++){
      ofile << " " << X_mat(i,j);
    }
    for(j=0;j<observables;j++){
      ofile << " " << Y_mat(i,j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
void load_data_file(std::string fileName, arma::mat &X, arma::mat &Y){
  int param = X.n_cols;
  int obs = Y.n_cols;
  int train = X.n_rows;

  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<train;i++){
    for(int j=0;j<param;j++){
      ifile >> X(i,j);
    }
    for(int j=0;j<obs;j++){
      ifile >> Y(i,j);
    }
  }
  ifile.close();
}
void load_beta_file(std::string betaName, arma::mat &beta){
  int load = beta.n_rows;
  int param = beta.n_cols;
  std::string temp;
  std::ifstream ifile;
  ifile.open(betaName);
  ifile >> temp;
  for(int i=0;i<load;i++){
    for(int j=0;j<param;j++){
      ifile >> beta(i,j);
    }
  }
  ifile.close();
}
void load_file(std::string fileName, arma::mat &file){
  int load = file.n_rows;
  int param = file.n_cols;
  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<load;i++){
    for(int j=0;j<param;j++){
      ifile >> file(i,j);
    }
  }
  ifile.close();
}
void write_file(std::string fileName, arma::mat &file){
  int write = file.n_rows;
  int param = file.n_cols;
  std::string temp;
  std::ofstream ofile;
  ofile.open(fileName);
  for(int i=0;i<write;i++){
    for(int j=0;j<param;j++){
      ofile << " " <<  file(i,j);
    }
    ofile << std::endl;
  }
  ofile.close();
}
void write_csv_file(std::string fileName, std::string title, arma::mat &file){
  int write = file.n_rows;
  int param = file.n_cols;
  std::string temp;
  std::ofstream ofile;
  ofile.open(fileName);
  ofile << title << '\n';
  for(int i=0;i<write;i++){
    for(int j=0;j<param-1;j++){
      ofile << file(i,j) << ", ";
    }
    ofile << file(i,param-1) << std::endl;
  }
  ofile.close();
}

