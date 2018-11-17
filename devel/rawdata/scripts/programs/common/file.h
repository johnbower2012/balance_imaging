#ifndef __FILE_H__
#define __FILE_H__

#include <armadillo>
#include<Eigen/Dense>
#include <string>
#include <fstream>

void write_output(arma::mat output_mat, int param, int obs, std::string outfilename);
void write_trainset(arma::mat X_mat, arma::mat Y_mat, std::string outfilename);
void load_data_file(std::string fileName, arma::mat &X, arma::mat &Y);
void load_beta_file(std::string betaName, arma::mat &beta);
void load_file(std::string fileName, arma::mat &file);
void write_file(std::string fileName, arma::mat &file);
void write_csv_file(std::string fileName, std::string title, arma::mat &file);

#endif
