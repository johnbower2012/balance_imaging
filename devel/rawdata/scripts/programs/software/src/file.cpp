#include "file.h"

bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }
void RemoveSpaces(std::string& str){
  std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
  str.erase(new_end, str.end()); 
}
void print_file(std::string filename){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      printf("%s\n",line.c_str());
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
void print_formatted_file(std::string filename,std::string delimiter){
  std::string line;
  std::fstream myfile(filename,std::ios::in);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	printf("%s\n",line.c_str());
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}
void load_file(std::string filename, Eigen::MatrixXd &Matrix,std::string delimiter){
  std::string line,
    number;
  std::size_t pos,length;
  std::ifstream myfile(filename);
  int rows=0,row=0,
    cols=0,col=0,
    max=0;

  if(myfile.is_open()){
    while(getline(myfile,line)){
      cols=0;
      if( line[0]!='#' && line!="\0" ){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	while( (pos=line.find(delimiter)) != std::string::npos){
	  line.erase(0,pos+delimiter.length());
	  cols++;
	}
	if(cols>max){
	  max=cols;
	}
	rows++;
      }
      cols=max;
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }

  Matrix = Eigen::MatrixXd::Zero(rows,cols);
  myfile.open(filename);
  if(myfile.is_open()){
    while(getline(myfile,line)){
      if(line[0]!='#' && line!="\0"){
	RemoveSpaces(line);
	if(line.length()<delimiter.length()){
	}else{
	  if(line.substr(line.length()-delimiter.length(),line.length()) != delimiter){
	    line.append(delimiter);
	  }
	}
	if(line.find(delimiter) == 0){
	  line.erase(0,delimiter.length());
	}
	col=0;
	while( (pos=line.find(delimiter)) != std::string::npos){
	  number=line.substr(0,pos+delimiter.length());
	  Matrix(row,col) = std::stof( number );
	  line.erase(0,pos+delimiter.length());
	  col++;
	}
	row++;
      }
    }
    myfile.close();
  }else{
    printf("Unable to open file.\n");
  }
}



void print_file(std::string outfilename, Eigen::MatrixXd matrix){
  unsigned int rows = matrix.rows();
  unsigned int cols = matrix.cols();
  std::ofstream ofile;
  ofile.open(outfilename);
  for(unsigned int i=0;i<rows;i++){
    for(unsigned int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void load_file(std::string infilename, Eigen::MatrixXd &matrix, unsigned int parameters, unsigned int samples){
  matrix = Eigen::MatrixXd::Zero(samples,parameters);
  std::ifstream ifile;
  ifile.open(infilename);
  for(unsigned int i=0;i<samples;i++){
    for(unsigned int j=0;j<parameters;j++){
      ifile >> matrix(i,j);
    }
  }
  ifile.close();
}

/////////////////////////////

void load_file(int files, int lines, int runs, std::vector<std::string> infilename, arma::vec &delY_vec, std::vector<arma::mat> &val_matrix){
  std::ifstream ifile;
  int i, j, k;
  val_matrix = std::vector<arma::mat>(files);
  for(i=0;i<files;i++){
    val_matrix[i] = arma::zeros<arma::mat>(runs,lines);
    ifile.open(infilename[i]);
    printf("+++++ %s LOADED +++++\n", infilename[i].c_str());
    for(j=0;j<lines;j++){
      for(k=0;k<runs+1;k++){
	if(k==0){
	  ifile >> delY_vec(j);
	} else{
	  ifile >> val_matrix[i](k-1,j);
	}
      }
    }
    ifile.close();
  }
}
void print_file(std::string outfilename, std::string title, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i);
  }
  ofile << '\n';
  ofile.close();
}
void print_file(std::string outfilename, arma::vec vector){
  int elem = vector.n_elem;
  std::ofstream ofile;
  ofile.open(outfilename);
  for(int i=0;i<elem;i++){
    ofile << ' ' << vector(i);
  }
  ofile << '\n';
  ofile.close();
}
void print_file(std::string outfilename, std::string title, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
  std::ofstream ofile;
  ofile.open(outfilename);
  ofile << title.c_str() << '\n';
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ofile << ' ' << matrix(i,j);
    }
    ofile << '\n';
  }
  ofile.close();
}
void print_file(std::string outfilename, arma::mat matrix){
  int rows = matrix.n_rows;
  int cols = matrix.n_cols;
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
arma::mat load_range_file(std::string filename, int parameters, std::vector<std::string> &Names){
  std::ifstream ifile;
  arma::mat File = arma::zeros<arma::mat>(parameters,2);
  Names.resize(parameters);
  ifile.open(filename);
  for(int i=0;i<parameters;i++){
    ifile >> Names[i]; ifile >> Names[i];
    ifile >> File(i,0);
    ifile >> File(i,1);
  }
  ifile.close();
  return File;
}


void create_parameter_prior_files(std::vector<std::string> Names, arma::mat hypercube){
  std::ofstream ofile;
  int lhp_samples = hypercube.n_rows;
  int parameters = hypercube.n_cols;
  for(int i=0;i<lhp_samples;i++){
    char fn[50];
    sprintf(fn,"../../../model_output/run%04d/parameters.dat",i);
    ofile.open(fn);
    for(int j=0;j<parameters;j++){
      ofile << Names[j] << " " << hypercube(i,j) << '\n';
    }
    ofile.close();
  }
}
