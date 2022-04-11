#ifndef quark_model_h
#define quark_model_h

#include "particles.h"

class quarks_class{
	public:
	double rhoB, temperature=0., PressureTot;
	double Yu, Yd, Ys;
	int iFlavour;
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
	void setEOS_betaEq_2F(double rhob_, double temp_, particle &electron_, particle &muon_);
	void setEOS_symmetric(double rhob_, double temp_);
	void setEOS_2flavour(double rhob_, double temp_);
	void 	setEoSFlavorFixed(double rhob_, double temperature_, 
																					double Yu_,  double Yd_,  double Ys_);
	void setEoSFlavor_PressFixed(double press_, double temp_, 
                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_);

	void setEffectiveMasses(double rhob_, double temp_);


	void setDensities(double mu_u, double mu_d, double mu_s);
	double sigma(double t_);
	double getDmassDtemp();
	double getDmassDdens();

	double getEnergy();
	double getPressure();
	double getEntropy();
	double getFenergy();
	double getOmega();
	
	void setTemperature(double temp_);
	void setBfield(double B_);
	void setAMM();

};

struct ThreeFlavourBetaEqFunctor{
public:
	ThreeFlavourBetaEqFunctor(quarks_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};

struct TwoFlavourBetaEqFunctor{
public:
	TwoFlavourBetaEqFunctor(quarks_class & quarks_, particle & electron_, particle &muon_):
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


struct TwoFlavourSymFunctor{
public:
	TwoFlavourSymFunctor(quarks_class & quarks_): quarks(quarks_)
	{}

template <typename T>
	bool operator()(const T* arg, T* residuals) const;

private:
	quarks_class &quarks;
};

struct QuarkFlavour_PressFixed{
public:
	QuarkFlavour_PressFixed(quarks_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    quarks_class 	&quarks;
		particle 			&electron;
		particle 			&muon;
};


#endif
