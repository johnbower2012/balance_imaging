#include<iostream>
#include<iomanip>
#include<fstream>
#include<armadillo>

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
  int irows,cols,orows;
  std::string ifile,ofile;
  if(argc<6){
    printf("Improper usage. Enter also 'irows cols orows ifile ofile' on same line.\n");
    exit(1);
  } else{
    irows = atoi(argv[1]);
    cols = atoi(argv[2]);
    orows = atoi(argv[3]);
    ifile = argv[4];
    ofile = argv[5];
  }
  std::vector<std::string> names(4);
  names[0] = "InvCoshPar_uu_width";
  names[1] = "InvCoshPar_ud_width";
  names[2] = "InvCoshPar_us_width";
  names[3] = "InvCoshPar_ss_width";


  arma::mat imat = arma::zeros<arma::mat>(irows,cols);
  arma::mat omat = arma::zeros<arma::mat>(orows,cols);
  read_file(ifile,imat);
  extract_rows(imat,omat);
  write_file(ofile,names, omat);

  return 0;
}
