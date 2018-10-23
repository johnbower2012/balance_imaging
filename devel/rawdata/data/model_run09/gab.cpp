#include<iostream>
#include<iomanip>
#include<cmath>
#include<string>
#include<armadillo>

void load_file(std::string fileName, arma::mat &file);
void write_file(std::string fileName, arma::mat &file);

int main(int argc, char* argv[]){
  int CHI=3,PAR=4,COUNT=20,SAM=100;

  std::string chiname="chi.dat";
  std::string parname="moments_parameters.dat";
  
  arma::mat chi=arma::zeros<arma::mat>(CHI,CHI);
  arma::mat par=arma::zeros<arma::mat>(COUNT,PAR);

  load_file(chiname,chi);
  load_file(parname,par);

  arma::mat gab=arma::zeros<arma::mat>(SAM,1+PAR*COUNT);
  double eta=0,width; int one,two;
  for(int i=0;i<SAM;i++){
    eta = (2.0-0.0)*(double) i/(double) SAM;
    gab(i,0) = eta;
    for(int a=0;a<PAR;a++){
      if(a==0){ one=two=0; }
      else if(a==1){ one=0;two=1; }
      else if(a==2){ one=0;two=2; }
      else{ one=2;two=2;}
      for(int j=0;j<COUNT;j++){
	width = par(j,a);
	gab(i,1+COUNT*a+j) = -chi(one,two)*2.0/width/acos(-1)/cosh(eta/width);
      }
    }
  }

  std::string gabname = "gab.dat";
  write_file(gabname,gab);

  return 0;
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
