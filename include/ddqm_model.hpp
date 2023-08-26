#pragma once

#include "particles.hpp"

class quarks_class{
	public:
	double rhoB, temperature=0., PressureTot;
	double Yu, Yd, Ys;
	double Ye, Ym;
	int iFlavor;
	double Bfield=0.;

//Particles
	quark_particle qu, qd, qs;
	double C=0., D=0., tcrit=0.;
	double lambda=1.60581;
	double muB, muQ;
	//Do hyperons and Deltas?  
	bool firstRun		=true; //so that initializing parameters are only used in the first run, 
												// afterwards, the previous solution is used!
// Construct (initialize) the object:
	quarks_class();
  ~quarks_class(void){};
//Define the parameterization
	void setParameters(double C_, double D_, double tcrit_);
	
	void setEOS_betaEq(double rhob_, double temp_, particle &electron_, particle &muon_);
	void setEOS_symmetric(double rhob_, double temp_);
	void 	setEoSFlavorFixed(double rhob_, double temperature_, 
																					double Yu_,  double Yd_,  double Ys_);
	void setEoSFlavor_PressFixed(double press_, double temp_, 
                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_);
	void setEoSFlavor_muBFixed(double press_, double temp_, 
	                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_);
	void setEoSFlavor_muBFixed2(double mub_, double temp_, 
                                      particle &electron_, particle &muon_,
                                      double Yu_,  double Yd_,  double Ys_, 
																			double Ye_, double Ym_);

	void setEffectiveMasses(double rhob_, double temp_);

	void setFlavorNumber(int if_);

	void setDensities(double mu_u, double mu_d, double mu_s);
	double sigma(double t_);
	double getDmassDtemp();
	double getDmassDdens();

	double getBaryonDens();
	double getEnergy();
	double getPressure();
	double getEntropy();
	double getFenergy();
	double getOmega();
	
	void setTemperature(double temp_);
	void setBfield(bool dob_, double B_);
	void setAMM();

};

struct QuarkBetaEqFunctor{
public:
	QuarkBetaEqFunctor(quarks_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};


struct SymmetricFunctor{
public:
	SymmetricFunctor(quarks_class & quarks_): quarks(quarks_)
	{}

template <typename T>
	bool operator()(const T* arg, T* residuals) const;

private:
	quarks_class &quarks;
};


struct QuarkFlavor_PressFixed{
public:
	QuarkFlavor_PressFixed(quarks_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};

struct QuarkFlavor_muBFixed{
public:
	QuarkFlavor_muBFixed(quarks_class & quarks_, particle & electron_, particle & muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};


struct QuarkFlavor_muBFixed2{
public:
	QuarkFlavor_muBFixed2(quarks_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};



