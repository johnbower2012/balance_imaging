#include<iostream>
#include<iomanip>
#include<fstream>
#include<armadillo>
void Fmean(arma::mat mat, arma::vec &mean){
  int rows=mat.n_rows;
  int cols=mat.n_cols;
  mean = arma::zeros<arma::vec>(cols);
  for(int col=0;col<cols;col++){
    for(int row=0;row<rows;row++){
      mean(col) += mat(row,col);
    }
    mean(col) /= (double) rows;
  }
}
void read_file(std::string fileName, arma::mat &file){
  std::ifstream ifile;
  ifile.open(fileName);
  for(int i=0;i<file.n_rows;i++){
    for(int j=0;j<file.n_cols;j++){
      ifile >> file(i,j);
    }
  }
  ifile.close();
}
void write_file(std::string fileName, std::vector<std::string> names, arma::mat &file){
  std::ofstream ofile;
  
  ofile.open(fileName);
  for(int i=0;i<names.size();i++){
    ofile << " " << names[i];
  }
  ofile << '\n';
  for(int i=0;i<file.n_rows;i++){
    for(int j=0;j<file.n_cols;j++){
      ofile << " " << file(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void write_file(std::string fileName, arma::vec &file){
  std::ofstream ofile;
  
  ofile.open(fileName);
  for(int i=0;i<file.n_elem;i++){
    ofile << " " << file(i);
  }
  ofile << '\n';
  ofile.close();
}
void extract_rows(arma::mat &imat, arma::mat &omat){
  int irows = imat.n_rows,orows = omat.n_rows,ratio = (irows-1)/(orows-1);
  for(int i=0;i<orows;i++){
    omat.row(i) = imat.row(i*ratio);
    printf("%d %d\n",i,i*ratio);
    omat.row(i).print();
    imat.row(i*ratio).print();
  }
}

int main(int argc, char* argv[]){
  int rows,cols;
  std::string ifile,ofile;
  if(argc<5){
    printf("Improper usage. Enter also 'rows cols ifile ofile' on same line.\n");
    exit(1);
  } else{
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
    ifile = argv[3];
    ofile = argv[4];
  }
  arma::mat mat = arma::zeros<arma::mat>(rows,cols);
  arma::vec mean = arma::zeros<arma::vec>(cols);
  read_file(ifile,mat);
  Fmean(mat,mean);
  write_file(ofile,mean);

  return 0;
}
