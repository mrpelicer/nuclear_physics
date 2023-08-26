#pragma once

#include <iostream>
#include <vector>
#include <cmath>
//Gsl for integration 
#include <functional>
#include "constant.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include "ceres/ceres.h"
#include "glog/logging.h"
using ceres::NumericDiffCostFunction;
//using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

using namespace std;
  

struct particle{

  double mass=1., mass_eff=1.;
  double gamma=2.;            //degeneracy
  double spin=1./2., I3=0., Q=0.;         //isospin and charge
  double stg=0;                  //strangeness
  double chemPot=0., chemPot_eff=0., kf=0., kf2=0.;
  double density=0.,  condensate=0., Qdens=0.;//rhob, rhos 
  double energy=0., pressure=0., entropy=0.;   //thermodynamics
  double temperature=0.;
  string type="H";
  
  bool doB  =false;
  double Bfield=0;
  int inumax=0;
  double densityPP=0., densityP=0., densityM=0., densityMM=0.;
  double pressureParallel =0., pressureTransverse =0.;
  double magnetization=0.;

  bool doamm=false;
  double Kb=0.;
  int inumaxPP=0, inumaxP=0, inumaxM=0, inumaxMM=0;
  
  bool dosrc=false; //short range correlations
  double c_=0., phi_=0., delta_=0.;

  void calculateProperties();
  void calculateDensity();
  void calculateCondensate();
	void solveChemPotEff();
  // void setBaryonEff(double mub_, double muq_, double gphi0_, double gv0_, double gb0_);
  void setChemicalPotential(double mu_);
  void setBaryonEff(double chemPotEff_, double massEff_);

  void setQuarkEff(double mueff_);
  void setLepton(double muq_);

  void setBfield(bool dob_,double B_);
	void setAMM(bool doa_, double Kb_);

  double densityT0();
  double condensateT0();
  double energyT0();
  double pressureT0();
  
  double condensateT0_src();
  double energyT0_src();
  double pressureT0_src();

};

struct quark_particle : public particle{
  double omega0=0.;

  void calculateQProperties();

  double getOmega0(double mueff_, double mass_);
  double getDOmega0Dmass();
  double getDOmega0DmuEf();
  double getDOmega0Dt();

};


struct ChemPotFunctor{
public:
	ChemPotFunctor(particle & fermion_):fermion(fermion_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    particle &fermion;
};


double integrate(double (func)(double, void *), void *parametersPointer);
double integrate_src(double (func)(double, void *), void *parametersPointer);


double densityFunc(double x, void *p);
double density_condensateFunc(double x, void *p);
double energyFunc(double x, void *p);
double pressureFunc(double x, void *p);
double entropyFunc(double x, void *p);
double fermiDirac(double ener, double chemPotEff, double T);

double density_condensateFunc_src(double x, void *p);
double energyFunc_src(double x, void *p);
double pressureFunc_src(double x, void *p);


