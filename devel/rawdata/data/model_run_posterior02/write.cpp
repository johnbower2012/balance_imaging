#include<iostream>
#include<iomanip>
#include<fstream>
#include<armadillo>
#include<vector>

void read_header_file(std::string fileName, std::vector<std::string> &header, arma::mat &file);
void write_param_file(std::string fileName, std::vector<std::string> &header, arma::mat &file);
template <class type>
void print_vector(std::vector<type> vector);
void mkdir(std::string dirName);
std::string int_to_string(int number, int digits);
void mkdir_loop(std::string dirName,int start, int end, int digits);

///////////////////////////////////////////////
int main(int argc,char* argv[]){
  int parameters=4;
  int rows=20;
  std::string iname = "mcmcextract.dat";
  if(argc<4){
    printf("Improper usage. Please enter also 'filename rows cols' on same line.\n");
    exit(1);
  } else{
    iname=argv[1];
    rows=atoi(argv[2]);
    parameters=atoi(argv[3]);
  }
  std::vector<std::string> names(parameters);
  arma::mat posterior = arma::zeros<arma::mat>(rows,parameters);
  read_header_file(iname,names,posterior);

  int digits=4;
  std::string dirBase = "run";
  std::string dirName;
  std::string number;
  std::string fileName = "parameters.dat";
  std::string fullName;
  arma::mat row = arma::zeros<arma::mat>(1,parameters);
  std::ofstream ofile;
  for(int i=0;i<rows;i++){
    number = int_to_string(i,digits);
    dirName = dirBase + number;
    mkdir(dirName);
    fullName = dirName + "/" + fileName;
    ofile.open(fullName);
    for(int j=0;j<parameters;j++){
      ofile << names[j] << " " << posterior(i,j) << '\n';
    }
    ofile.close();
  }

  return 0;
}
/////////////////////////////////////////////////////////  
void read_header_file(std::string fileName, std::vector<std::string> &header, arma::mat &file){
  int rows = file.n_rows;
  int cols = file.n_cols;
  std::ifstream ifile;

  ifile.open(fileName);
  for(int i=0;i<cols;i++){
    ifile >> header[i];
  }
  for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
      ifile >> file(i,j);
    }
  }
  ifile.close();
}
void write_param_file(std::string fileName, std::vector<std::string> &header, arma::mat &file){
  std::ofstream ofile;
  int rows=file.n_rows;
  int cols=file.n_cols;
  ofile.open(fileName);
  for(int i=0;i<cols;i++){
    ofile << header[i] << " ";
    for(int j=0;j<rows;j++){
      ofile << file(j,i) << " ";
    }
    ofile << '\n';
  }
  ofile.close();
}
template <class type>
void print_vector(std::vector<type> vector){
  int size=(int)vector.size();
  for(int i=0;i<size;i++){
    printf("%s ",vector[i].c_str());
  }
  printf("\n");
}
void mkdir(std::string dirName){
  std::string cmd = "if [ ! -d " + dirName + " ]; then mkdir -p " + dirName + "; fi";
  const char *command = cmd.c_str();
  system(command);
}
std::string int_to_string(int number, int digits){
  int count=1,temp=number;
  std::string number_string;
  while(temp>9){
    temp /= 10;
    count++;
  }
  count = digits - count;
  while(count>0){
    number_string += "0";
    count--;
  }
  number_string += std::to_string(number);
  return number_string;
}
void mkdir_loop(std::string dirName,int start, int end, int digits){
  int count=digits,temp;
  std::string fullName,number;
  for(int i=start;i<end+1;i++){
    count=1;
    temp=i;
    while(temp>9){
      temp /= 10;
      count++;
    }
    count = digits - count;
    fullName = dirName;
    while(count>0){
      fullName += "0";
      count--;
    }
    number = std::to_string(i);
    fullName += number;
    printf("%s\n",fullName.c_str());
    std::string cmd = "if [ ! -d " + fullName + " ]; then mkdir -p " + fullName + "; fi";
    const char *command = cmd.c_str();
    system(command);
  }
}