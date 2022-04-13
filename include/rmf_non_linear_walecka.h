#ifndef quantum_hadrodynamics_h
#define quantum_hadrodynamics_h

#include "particles.h"

class nlwm_class{
public:

	std::string parametrization, parhyp, pardelta;
//Nucleon parameters
  double Mn=0., Ms=0., Mv=0., Mb=0., Mt=0.; //Masses
  double gs=0., gv=0., gb=0., gt=0.; //Adimensional!
  double gs3=0., gs4=0.; //gs3 in fm^(-1) and gs4 is adimensional.
  double xsi=0., Lv=0.;
  double rho0=0., Mstar=0.;
  double phi0, V0, b0, theta0, Mef;
	double Bfield=0.;

//Hyperons and deltas parameters:
	double xsl=0., xss=0., xsx=0., xsd=0.;
	double xvl=0., xvs=0., xvx=0., xvd=0.;
	double xbl=0., xbs=0., xbx=0., xbd=0.;
	double xtl=0., xts=0., xtx=0.;

//Thermodynamic variables:
  double rhoB, rhoS, rho3, rhoQ, Yp, temperature;
	double rhoB_integrated, rhoS_integrated, rho3_integrated, rhoq_integrated; 
	double rhoB_eff, rhoS_eff, rho3_eff;// w/ hyperon couplings
	double muB, muQ, PressureTot;

	double yN, yH, yD;
//Particles
	particle proton, neutron;
	particle lambda0, sigmap, sigma0, sigmam, xi0, xim;
	particle deltapp, deltap, delta0, deltam;

	//Do hyperons and Deltas?  
	bool doHyperons	=false;
	bool doDeltas		=false;
	bool firstRun		=true; //so that initializing parameters are only used in the first run, 
												// afterwards, the previous solution is used!
// Construct (initialize) the object:
	nlwm_class(std::string parameters);
  ~nlwm_class(void){};
	void includeHyperons(bool do_, std::string parameters_);
	void includeDeltas(	 bool do_, std::string parameters_);

	void setParametrization(std::string parametrization_);
	void printParameters(void);

	void setThermodynamics();
	void setTemperature(double temp_);
	void setBfield(bool dob_, double B_);
	void setAMM(bool doa_);
	
//Set EOS fixed proton fraction: input density, proton fraction and temperature
  void setEOS_nucleons(double rhoB_, double Yp_, double temp_); //npe
//Set EOS beta-equilibrium: input density and temperature
	void setEOS_neutrons(double rhoB_, double temp_, particle &electron_, particle &muon_);
	void setEOS_betaEq(double rhoB_, double temp_,	particle &electron_, particle &muon_);
	void setEOS_betaEq(double rhoB_, double temp_,	particle &electron_);
	void setEOS_betaEq_PressureFixed(double press_, double temp_,	
																		particle &electron_, particle &muon_);
	void setEOS_betaEq_muB(double mub_, double temp_,	
																		particle &electron_, particle &muon_);
//Set EOS for the pasta solver: input effective chemical potentials and mass.
	void setEOS_coexistence(double nup_, double nun_, double mef_);	

//Set the solver for meson equations of motion:
	void setScalarMeanFields();
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

//get the thermodynamic quantities for baryons:
  double getEnergy(void);
  double getPressure(void);
  double getEntropy(void);

  double getPressureParallel(void);
  double getPressureTransverse(void);
	double getMagnetization(void);
};

// Functor for Ceres:

//Solve for the vector meson
struct VFunctor{
public:
	VFunctor(nlwm_class & baryons_):baryons(baryons_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class & baryons;
};

//Solve for scalar meson given the densities
struct SFunctor{
public:
	SFunctor(nlwm_class & baryons_):baryons(baryons_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class & baryons;
};

struct NeutronFunctor{
	public:
	NeutronFunctor(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};


struct BetaEqFunctor{
	public:
	BetaEqFunctor(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};

struct BetaEqFunctor2{
	public:
	BetaEqFunctor2(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};

struct BetaEqFunctor3{
public:
	BetaEqFunctor3(nlwm_class &baryons_, particle &electron_):baryons(baryons_), 
																														 electron(electron_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class &baryons;
		particle 	 &electron;
};

struct BetaEqFunctor_PressFixed{
	public:
	BetaEqFunctor_PressFixed(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};


struct BetaEqFunctor2_PressFixed{
	public:
	BetaEqFunctor2_PressFixed(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};

struct BetaEqFunctor_muBFixed{
	public:
	BetaEqFunctor_muBFixed(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};

struct BetaEqFunctor2_muBFixe{
	public:
	BetaEqFunctor2_muBFixe(nlwm_class &baryons_, particle &electron_, particle &muon_):
													baryons(baryons_), electron(electron_), muon(muon_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
};

#endif
