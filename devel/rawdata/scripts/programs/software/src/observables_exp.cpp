#include<iostream>
#include<fstream>
#include<cmath>
#include "armadillo"
#include "../lib/pricomana.cpp"
#include "../lib/observables_class.cpp"

int main(int argc, char* argv[]){

  int obs_file=3,
    files=4,
    observables=obs_file*files,
    lines=16,
    runs=1;
  int i,j,k;
  std::string dest_folder="../../calculations";

  /*
  if(argc<4+files){
    std::cout << "Improper input. Enter also 'lines runs observables ifn*[7]' on same line." << std::endl;
    exit(1);
  }
  else{
    infilename = new std::string[files];
    dest_folder = argv[1];
    lines=atoi(argv[2]);
    runs=atoi(argv[3]);
    files=atoi(argv[4]);
    for(i=0;i<files;i++){
      infilename[i]=argv[5+i];
    }
  }
  printf("Arguments read in as:\n");
  printf("observables.x %s %d %d %d",dest_folder.c_str(),lines,runs,files);
  for(i=0;i<files;i++){
    printf(" %s",infilename[i].c_str());
  }
  printf("\n");
  */
  /*********
	     LOAD FILE
	     CALCULATE OBS_MAT
  *********/
  std::vector<arma::mat> val_matrix(files);
  arma::mat obs_matrix;
  arma::vec delY_vec, obs_error;
  obs_matrix = arma::zeros<arma::mat>(runs,observables);
  delY_vec = arma::zeros<arma::vec>(lines);
  obs_error = arma::zeros<arma::vec>(observables);
  //all observables that aren't pipi, reduce unc by factor of ten.

  load_file(files, lines, runs, infilename, delY_vec, val_matrix);
  obs_matrix_moments(files,obs_file,val_matrix,delY_vec,obs_matrix);
  //reorder to match model output
  arma::mat temp = obs_matrix;
  for(int i=0;i<observables;i++){
    temp.col(i) = obs_matrix.col((i+obs_file)%observables);
  }
  obs_matrix=temp;
  std::string printname = dest_folder + "moments_exp_data.dat";
  std::string title="#pipi/ppbar/pK/KK";
  print_file(printname,obs_matrix);
  printf("%s saved to file.\n",printname.c_str());

  temp=0.1*obs_matrix;
  printname = dest_folder + "moments_error.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,temp);
  printf("%s saved to file.\n",printname.c_str());

  temp = arma::zeros<arma::mat>(0,0);
  printname = dest_folder + "moments_info.dat";
  title="#pipi/ppbar/pK/KK";
  print_file(printname,title,temp);
  printf("%s saved to file.\n",printname.c_str());
  
  delete[] val_matrix;
  delete[] infilename;

  return 0;
}
