#ifndef bag_model_h
#define bag_model_h

#include "particles.h"

class bag_model_class{
public:

	std::string parametrization, parhyp, pardelta;
//Nucleon parameters
  double Mv=0., gv=0., xsi=0., Bag=0., tcrit=0.;
  double  V0;
	double Bfield=0.;

//Hyperons and deltas parameters:
	double xvu=0., xvd=0., xvs=0.;

	int iFlavor=3;
//Thermodynamic variables:
  double rhoB, rhoS, rho3, rhoQ, Yp, temperature;
	double rhoB_integrated, rhoS_integrated, rho3_integrated, rhoq_integrated; 
	double rhoB_eff, rhoS_eff, rho3_eff;// w/ hyperon couplings
	double muB, muQ, PressureTot;

	double yN, yH, yD;
//Particles
	particle qu, qd, qs;
	bool firstRun		=true; //so that initializing parameters are only used in the first run, 
												// afterwards, the previous solution is used!
// Construct (initialize) the object:
	bag_model_class();
  ~bag_model_class(void){};
	void setParameters(double bag_, double Gv_, double Xv_, double tcrit_);
	void setParametrization(std::string parametrization_);
	void printParameters(void);

	void setThermodynamics();
	void setTemperature(double temp_);
	void setBfield(bool dob_, double B_);
	void setAMM(bool doa_);
	void setFlavorNumber(int if_);
	void setDensities(double muef_u, double muef_d, double muef_s);

//Set EOS fixed proton fraction: input density, proton fraction and temperature
	void setEOS_symmetric(double rhob_, double temp_);
	//Set EOS beta-equilibrium: input density and temperature
	void setEOS_betaEq(double rhoB_, double temp_,	particle &electron_, particle &muon_);
	void setEOS_betaEq_2F(double rhob_, double temp_, particle &electron_, particle &muon_);
	
	void setEOS_betaEq(double rhoB_, double temp_,	particle &electron_);
	void setEOS_betaEq_PressureFixed(double press_, double temp_,	
																		particle &electron_, particle &muon_);
//Set EOS for the pasta solver: input effective chemical potentials and mass.
	void setEOS_coexistence(double nup_, double nun_, double mef_);	

//Set the solver for meson equations of motion:
	void setVectorMeanFields(); 

//Input chemical potentials and meson fields and calculate density:
	void setDensities(double mub_, double muq_, double phi0_, double v0_, double b0_);
	void setDensities(double mub_, double muq_, double phi0_, double v0_, double b0_, double theta0_);

	double getBaryonDens();
	double getIsoDens();
	double getChargeDens();

	double getSigmaEffDens();
	double getOmegaEffDens();
	double getIsoEffDens();
	double getThetaEffDens();

	std::vector<double> getNucleonPotential();
	std::vector<double> getHyperonPotential();
	std::vector<double> getDeltaPotential();


//Set initial value for solvers:
	//hyperons+ deltas,b=0
	void setInitial_hd(double &mub_, double &mue_, double  &phi0_, double &v0_, double &b0_);
	//hyperons+ deltas,b>0
	void setInitial_hdb(double &mub_, double &mue_, double  &phi0_, double &v0_, double &b0_);

//Calculate reside of meson equations:
  double sigmaMeson_eom_residue(double rhoS_);
	double omegaMeson_eom_residue(double rhoB_);
	double rhoMeson_eom_residue(	double rho3_);
	double thetaMeson_eom_residue(double rhoT_);

//get the thermodynamic quantities for quarks:
  double getEnergy(void);
  double getPressure(void);
  double getEntropy(void);
	double getFenergy(void);


  double getPressureParallel(void);
  double getPressureTransverse(void);
	double getMagnetization(void);
};

// Functor for Ceres:

//Solve for the vector meson
struct VQFunctor{
public:
	VQFunctor(bag_model_class & quarks_):quarks(quarks_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		bag_model_class & quarks;
};

struct SymmetricFunctor{
public:
	SymmetricFunctor(bag_model_class & quarks_): quarks(quarks_)
	{}

template <typename T>
	bool operator()(const T* arg, T* residuals) const;

private:
	bag_model_class &quarks;
};

struct QuarkBetaEqFunctor{
public:
	QuarkBetaEqFunctor(bag_model_class & quarks_, particle & electron_, particle &muon_):
									quarks(quarks_), electron(electron_), muon(muon_)
  {}

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    bag_model_class 	&quarks;
		particle 					&electron;
		particle 					&muon;
};

#endif
