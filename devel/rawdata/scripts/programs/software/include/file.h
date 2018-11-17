#ifndef __FILE_H__
#define __FILE_H__

#include<fstream>
#include<Eigen/Dense>
#include<string>
#include "armadillo"

bool BothAreSpaces(char lhs, char rhs);
void RemoveSpaces(std::string& str);
void print_file(std::string filename);
void print_formatted_file(std::string filename,std::string delimiter);
void load_file(std::string filename, Eigen::MatrixXd &Matrix,std::string delimiter);

void print_file(std::string outfilename, Eigen::MatrixXd matrix);
void load_file(std::string infilename, Eigen::MatrixXd &matrix, unsigned int parameters, unsigned int samples);
void load_file(int files, int lines, int runs, std::string *&infilename, arma::vec &delY_vec, arma::mat *&val_matrix);
void print_file(std::string outfilename, std::string title, arma::vec vector);
void print_file(std::string outfilename, arma::vec vector);
void print_file(std::string outfilename, std::string title, arma::mat matrix);
void print_file(std::string outfilename, arma::mat matrix);
arma::mat load_range_file(std::string filename, int parameters, std::vector<std::string> &Names);

void create_parameter_prior_files(std::vector<std::string> Names, arma::mat hypercube);

#endif
