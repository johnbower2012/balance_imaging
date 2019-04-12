#ifndef __MCMC_H__
#define __MCMC_H__

#include<Eigen/Dense>
#include<random>
#include<chrono>
#include "emulator.h"

class CRandom{
 public:
  unsigned seed;
  std::default_random_engine generator;
  std::normal_distribution<double> normal_dist;
  std::uniform_real_distribution<double> uniform_dist;

  CRandom();
  CRandom(unsigned Seed);

  void setSeed(unsigned Seed);
  void setSeedClock();

  double normal();
  double uniform();
};

class CMCMC{
 public:
  double maxLogLikelihood;
  double Likelihood;
  CRandom random;

  int paramCount;
  int NSamples;

  Eigen::MatrixXd Range;
  Eigen::VectorXd Widths;

  Eigen::MatrixXd Position;
  Eigen::MatrixXd testPosition;

  CMCMC();
  CMCMC(Eigen::MatrixXd newRange, Eigen::VectorXd newWidths);
  CMCMC(std::string filename);
  void Construct(Eigen::MatrixXd newRange, Eigen::VectorXd newWidths);
  void setSeed(unsigned Seed);
  void setPosition();
  void setPosition(Eigen::MatrixXd newPosition);
  void setRange(Eigen::MatrixXd newRange);
  void setWidths(Eigen::VectorXd newWidths);

  void step();
  //void step(CCosh dist, Eigen::MatrixXd MinMax);
  double getLikelihood(Eigen::MatrixXd Z1, Eigen::MatrixXd Z2);
  bool decide(double likelihood);
  //double getLogLikelihoodGaussian(Eigen::MatrixXd Z);
  //bool decideGaussian(Eigen::MatrixXd Z);
  //double getLogLikelihoodLorentzian(Eigen::MatrixXd Z);
  //bool decideLorentzian(Eigen::MatrixXd Z);

  Eigen::MatrixXd Run(CEmulator *Emulator, Eigen::MatrixXd Target, int NSamp);
};

#endif
