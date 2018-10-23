#include "file.h"

void print_file(std::string outfilename, Eigen::MatrixXd matrix){
  int rows = matrix.rows();
  int cols = matrix.cols();
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void load_file(std::string infilename, Eigen::MatrixXd &matrix, int parameters, int samples){
  matrix = Eigen::MatrixXd::Zero(samples,parameters);
  std::ifstream ifile;
  ifile.open(infilename);
  for(int i=0;i<samples;i++){
    for(int j=0;j<parameters;j++){
      ifile >> matrix(i,j);
    }
  }
  ifile.close();
}
