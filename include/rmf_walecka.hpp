#pragma once
#include "particles.hpp"

class nlwm_class{
public:

	std::string parametrization, parhyp, pardelta;

//Nucleon parameters
  	double Mn=0., Ms=0., Md=0., Mv=0., Mr=0., Mp=0.; //Masses
  	double gs=0., gv=0., gr=0., gp=0., gd=0.; //Adimensional!
  	double gs3=0., gs4=0.; //gs3 in fm^(-1) and gs4 is adimensional.
  	double gv4=0., Lvr=0.;
  	double rho0=0., Mstar=0.;
  	double sigma_meson=0., omega_meson=0., rho_meson=0., phi_meson=0., delta_meson=0.;	
  
	double Mef=0.;
	double Bfield=0.;

//Hyperons and deltas parameters:
	//x_{Meson Baryon} = g_{Meson Baryon}/g_{Meson Nucleon}
	double xsl=0., xss=0., xsx=0., xsd=0.;	// x_sigma_lambda, x_sigma_sigma, x_sigma_xi, x_sigma_delta
	double xdl=0., xds=0., xdx=0., xdd=0.;	// x_delta_lambda, x_delta_delta, x_delta_xi, x_delta_delta
	double xvl=0., xvs=0., xvx=0., xvd=0.;	// x_omega_lambda, x_omega_omega, x_omega_xi, x_omega_delta
	double xrl=0., xrs=0., xrx=0., xrd=0.;	// x_rho_lambda, x_rho_rho, x_rho_xi, x_rho_delta
	double xpl=0., xps=0., xpx=0.;			// x_phi_lambda, x_phi_sigma, x_phi_xi

//density dependent coupling parameters:
	double as=0., bs=0., cs=0., ds=0.;
	double ad=0., bd=0., cd=0., dd=0.;
	double av=0., bv=0., cv=0., dv=0.;
	double ar=0.;

//Thermodynamic variables:
  	double rhoB, rho3, rhoS, rhoS3, rhoQ, Yp, temperature;
	double rhoB_integrated, rhoS_integrated, rho3_integrated, rhoS3_integrated, rhoq_integrated; 
	double rhoB_eff, rhoS_eff, rho3_eff, rhoS3_eff;// w/ hyperon couplings
	double muB, muQ, muS, PressureTot;

	double yN, yH, yD;
//Particles
	particle proton, neutron;
	particle lambda0, sigmap, sigma0, sigmam, xi0, xim;
	particle deltapp, deltap, delta0, deltam;

	//Do hyperons and Deltas?  
	bool useHyperons	=false;
	bool useDeltas		=false;

	bool useDensityDependentCoupling	=false;
	
	//Do electrons and muons?
	bool useElectron=true;
	bool useMuon=true;

	bool useSigmaMeson=true;
	bool useOmegaMeson=true;
	bool useRhoMeson=true;
	bool usePhiMeson=true;
	bool useDeltaMeson=true;

	bool firstRun		=true; //so that initializing parameters are only used in the first run, 
												// afterwards, the previous solution is used!
// Construct (initialize) the object:
	nlwm_class(std::string parameters);
	nlwm_class();
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
	void setEOS_nucleons(double rhoB_, double Yp_, double temp_);
	
	//Set EOS fixed proton fraction with short range correlations
	void setEOS_src_nucleons(double rhoB_, double Yp_, double temp_); //npe
	
	//Set EOS beta-equilibrium: input density and temperature
	void setEOS_neutrons(double rhoB_, double temp_);
	void setEOS_betaEq(double rhoB_, double temp_,	particle &electron_, particle &muon_);
	void setEOS_fixedYl(double rhoB_, double temp_, double Yle_, double Ylm_, 
												particle &electron_, particle &muon_, particle &ne_, particle &nm_);

	
	//Set EOS for the pasta: input effective chemical potentials and mass
	void setEOS_coexistence(double nup_, double nun_, double mef_);	
	void setEOS_coexistence_src(double nup_, double nun_, double mef_, double yp_);	

	//Set the solver for scalar meson equations of motion given the densities
	void setSigmaMeanField();
	//Set the solver for vector meson equations of motion given the densities and effective_mass
	void setVectorMeanFields(); 
	//Set the solver for all meson equations of motion given proton and neutron densities
	void setMesonFields();


//Input chemical potentials and meson fields and calculate density:
	// void setDensities(double mub_, double muq_, double sigma_, double omega_, double rho_);
	// void setDensities(double mub_, double muq_, double sigma_, double omega_, double rho_, double phi_);
	void setDensities(double mub_, double muq_, double sigma_,  double delta_, double omega_, double rho_, double phi_, double rear_);
	// void setDensitiesDD(double mub_, double muq_, double sigma_, double omega_, double rho_, double phi_, double rear_);

	double getBaryonDens();
	double getIsospinDensity();
	double getChargeDens();

	double getSigmaSource();
	double getDeltaSource();
	double getOmegaSource();
	double getRhoSource();
	double getPhiSource();

	std::vector<double> getNucleonPotential();
	std::vector<double> getHyperonPotential();
	std::vector<double> getDeltaPotential();

	double getCoupling_sigma(double rhob_);
	double getCoupling_delta(double rhob_);
	double getCoupling_phi(double rhob_);
	double getCoupling_omega(double rhob_);
	double getCoupling_rho(double rhob_);
	
	double getDerivativeCoupling_sigma(double rhob_);
	double getDerivativeCoupling_delta(double rhob_);
	double getDerivativeCoupling_omega(double rhob_);
	double getDerivativeCoupling_rho(double rhob_);
	double getDerivativeCoupling_phi(double rhob_);

	double getRearrangementEnergy(void);

//Set initial value for solvers:
	//hyperons+ deltas,b=0
	void setInitial_hd(double &mub_, double &mue_, double  &sigma_, double &omega_, double &rho_);
	//hyperons+ deltas,b>0
	void setInitial_hdb(double &mub_, double &mue_, double  &sigma_, double &omega_, double &rho_);

//Calculate reside of meson equations:
  	double sigmaMeson_eom_residue(double rhoS_);
	double deltaMeson_eom_residue(double rhoS3_);
	double omegaMeson_eom_residue(double rhoB_);
	double rhoMeson_eom_residue(	double rho3_);
	double phiMeson_eom_residue(double rhoT_);

//get the thermodynamic quantities for baryons:
  double getEnergy(void);
  double getPressure(void);
  double getEntropy(void);

  double getPressureParallel(void);
  double getPressureTransverse(void);
	double getMagnetization(void);
};



struct MesonFieldsFunctor{
public:
	MesonFieldsFunctor(nlwm_class & baryons_):baryons(baryons_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class & baryons;
};

struct VFunctor{
public:
	VFunctor(nlwm_class & baryons_):baryons(baryons_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class & baryons;
};


//Find the scalar meson given the baryon densities/effective chemical potentials
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
	NeutronFunctor(nlwm_class &baryons_):baryons(baryons_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 		&baryons;
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

struct YlFunctor{
	public:
	YlFunctor(nlwm_class &baryons_, particle &electron_, particle &muon_,
																		particle &ne_, 			 particle &nm_,
																		double Yle_, double Ylm_):
													baryons(baryons_), electron(electron_), muon(muon_),
																						ne(ne_), nm(nm_), Yle(Yle_), Ylm(Ylm_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
		particle 		&ne;
		particle 		&nm;
		double Yle, Ylm;
};

struct YlFunctor2{
	public:
	YlFunctor2(nlwm_class &baryons_, particle &electron_, particle &muon_,
																		particle &ne_, 			 particle &nm_,
																		double Yle_, double Ylm_):
													baryons(baryons_), electron(electron_), muon(muon_),
																						ne(ne_), nm(nm_), Yle(Yle_), Ylm(Ylm_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
		particle 		&ne;
		particle 		&nm;
		double Yle, Ylm;
};


struct YlFunctorDD{
	public:
	YlFunctorDD(nlwm_class &baryons_, particle &electron_, particle &muon_,
																		particle &ne_, 			 particle &nm_,
																		double Yle_, double Ylm_):
													baryons(baryons_), electron(electron_), muon(muon_),
																						ne(ne_), nm(nm_), Yle(Yle_), Ylm(Ylm_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		nlwm_class 	&baryons;
		particle 		&electron;
		particle 		&muon;
		particle 		&ne;
		particle 		&nm;
		double Yle, Ylm;
};



