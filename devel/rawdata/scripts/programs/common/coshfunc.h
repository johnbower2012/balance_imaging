#ifndef __COSHFUNC_H__
#define __COSHFUNC_H__

#include<iostream>
#include<cmath>
#include<fstream>
#include<Eigen/Dense>
#include<random>
#include <cstdlib>
#include <cstdio>
#include <complex>

#include "parametermap.h"
#define PI 4.0*atan(1.0)

class CCosh{
public:
  int nmax;
  Eigen::MatrixXd A;
  Eigen::VectorXd I;
  Eigen::VectorXd Z;

  CCosh();
  CCosh(int nmaxset);

  void Set_NMax(int nmaxset);
  double GetG(int n,double x);
  double GetRandomX();
  void TestOverlap();
  void CalcOverlap(int n,int nprime);
  void CalcOverlapNumerical(int n,int nprime);
  void Integrate(double x);
  
  int seed;
  std::mt19937 rng;
  std::uniform_real_distribution<double> uniform_dist;
};

class CDistribution{
 public:
  CDistribution();
  string parsfilename;
  CParameterMap parmap;
  virtual void GenX(int ab,double &x,double &weight);
};

class CDistCosh : public CDistribution{
 public:
  Eigen::VectorXd width;
  Eigen::MatrixXd g;

  CCosh coshcalc;  

  //  CDistCosh(string parsfilename);
  CDistCosh(Eigen::VectorXd Width, Eigen::MatrixXd G);
  CDistCosh();

  void Set_WG(Eigen::VectorXd Width, Eigen::MatrixXd G);
  void GenX(int ab,double &x,double &weight);
  void Functions(int n, double x, Eigen::MatrixXd &M);
  void FunctionSet(int grid, double x, int samples, int ab, int nmax, Eigen::MatrixXd G, Eigen::MatrixXd &Functions);
  double Integrate(int ab);
  double IntegrateNumerical(int ab);
};


#endif