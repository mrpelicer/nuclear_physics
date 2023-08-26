#pragma once

#include "particles.hpp"
#include "rmf_walecka.hpp"
#include "interpolator.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cuba.h>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>	
#include <gsl/gsl_spline.h>

//define cuba integration
// #define NDIM 1
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-8
#define EPSABS 1e-8
#define VERBOSE_AVG 0
#define VERBOSE_CLOG 1
//if using suave, must define:
#define LAST 4
#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define SEED 0
#define MINEVAL 0
#define MAXEVAL 1e6

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

class pasta_class{
public:	

	double rhoB, YpG, temperature, dim;
	string type_solver;
	bool firstRun=true;
  bool doSRC=false;
	int iType;
	double nup1_guess, nun1_guess, Mef1_guess, nup2_guess, nun2_guess, Mef2_guess;
	double f, Rd, Rw, VN, Vw;
	
	nlwm_class &cluster;
	nlwm_class &gas;

	pasta_class(nlwm_class &cluster_, nlwm_class &gas_);
	~pasta_class(void){};

	void setInitialCPA(double &nup1_, double &nun1_, double &mef1_,
									double &nup2_, double &nun2_, double &mef2_);

  void setInitialCPA_src(double &nup1_, double &nun1_, double &mef1_, double &yp1,
									double &nup2_, double &nun2_, double &mef2_, double &yp2_);

	void solveCPA(double rhoB_, double Yp_, double temp_);
  void solveCPA_src(double rhoB_, double Yp_, double temp_);
	void solveCPA_betaEq(double rhoB_, double temp_, particle &electron_);

	void solveCLD(double rhoB_, double Yp_, double temp_, double dim_, int it_);
		
};

struct cpaFunctor{
public:
	cpaFunctor(pasta_class & pasta_): pasta(pasta_){
	}

	template <typename T>
	bool operator()(const T* x, T* residuals) const;
  

private:
		pasta_class &pasta;
};


struct cpa_srcFunctor{
public:
	cpa_srcFunctor(pasta_class & pasta_): pasta(pasta_){
	}

	template <typename T>
	bool operator()(const T* x, T* residuals) const;
  

private:
		pasta_class &pasta;
};


struct cpaFunctor_betaEq{
public:
	cpaFunctor_betaEq(pasta_class & pasta_, particle &electron_): pasta(pasta_), electron(electron_){
	}

	template <typename T>
	bool operator()(const T* x, T* residuals) const;
  

private:
		pasta_class &pasta;
		particle &electron;
};

struct cldFunctor{
public:
	cldFunctor(pasta_class & pasta_): pasta(pasta_){
	}

	template <typename T>
	bool operator()(const T* x, T* residuals) const;
  

private:
		pasta_class &pasta;
};


double getSurfaceTension(nlwm_class &cluster_, double Yp_, double temp_);

double getSurfaceTensionDerivative(nlwm_class &cluster_, double Yp_, double temp_);

void setSurfaceParameters(nlwm_class &cluster_, double &sigma0, double &sigma1, 
			double &sa1, double &sa2, double &sa3, double &sa4, double &sa5, double &sa6,
  		double &aa0, double &aa1, double &aa2, double &aa3, double &aa4, double &aa5,
  		double &ba0, double &ba1, double &ba2, double &ba3, double &ba4, double &ba5,
  		double &ca0, double &ca1, double &ca2, double &ca3, double &ca4, double &ca5);

double getPhiFunc(double dim_, double u_);	

double getPhiFuncDerivative(double dim_, double u_);

double getRadiusD(double dim_, double beta_, double Yp_, 
									nlwm_class &cluster_, nlwm_class &gas_);
	



// transport part:
class pasta_transport_class{
public:
  explicit pasta_transport_class(particle electron_);
  particle electron;
  
  int i_dim;
  double q0=0.; // only in case we extend this for solids!
  double xr=0., enerf=0., betar=0., gammar=0.;
  double ae=0., kTF=0.;
  double q;
  std::vector<double> qv, F2v;

  void setMomentum(double q_);
  void setMomentumVec(std::vector<double> qv_);
  void setFormFactor2Vec(std::vector<double> qv_);

  double q2dielectric_function(double q_);

  virtual double getVolume(){return 0;};
  virtual double getArea(){return 0;};
  virtual double getStructureFunction(double q_){return 0;};  
};


class cDroplet : public pasta_transport_class{
public:
  using pasta_transport_class::pasta_transport_class;
  double radius, a;
  vector<double> F2v;

  void setRadius(double r_);
  double getVolume();
  double getArea();
  double getStructureFunction(double q_);

};

//================ rods ==============
class cRod : public pasta_transport_class{
public:
  
  using pasta_transport_class::pasta_transport_class;

  double radius, length;
  double min_cost, max_cost;
  vector<double> Fa2v, Fp2v;

  void setLengths(double r_, double h_);

  void setLim_CosTheta(double min_cost_, double max_cost_);

  double getVolume();
  double getArea();

  double getStructureFunction(double q_, double cost_);
  double getStructureFunction2_Axial(double q_);
  double getStructureFunction2_Trans(double q_);
};


 extern int IF2_rod_axial(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata);

extern int IF2_rod_trans(const int *ndim, const cubareal xx[],const int *ncomp, 
                       cubareal ff[], void *userdata);


//================ slabs ==============

class cSlab : public pasta_transport_class{
public:
  using pasta_transport_class::pasta_transport_class;
  
  double lx, ly, lz;
  double min_cost, max_cost, min_phi, max_phi;
  vector<double> Fa2v, Fp2v;

  void setLengths(double lx_, double ly_, double lz_); 

  void setLim_CosTheta(double min_cost_, double max_cost_);
  void setLim_Phi(double min_phi_, double max_phi_);

  double getVolume();
  double getArea();
 
  double getStructureFunction(double q_, double cost_, double phi_);
  double getStructureFunction2_Axial(double q_);
  double getStructureFunction2_Trans(double q_);

};

extern int IF2_slab_axial(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata);

extern int IF2_slab_trans(const int *ndim, const cubareal xx[],const int *ncomp, 
                       cubareal ff[], void *userdata);


double getUnintegratedLambda_ei(double q_, pasta_transport_class pasta_);
double getCoulombIntegral(pasta_transport_class shape_);

extern int Integral_Coulomb(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata);

double integrate_coulomb(double (func)(double, void *), void *parametersPointer);

double coulomb_gsl(double x, void *p);
