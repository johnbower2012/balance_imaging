#ifndef FILE_H
#define FILE_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<Eigen/Dense>

void print_file(std::string outfilename, Eigen::MatrixXd matrix);
void load_file(std::string infilename, Eigen::MatrixXd &matrix, int parameters, int samples);

#endif
