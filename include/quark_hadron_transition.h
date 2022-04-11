#ifndef quark_hadron_transition
#define quark_hadron_transition

#include "particles.h"
#include "rmf_non_linear_walecka.h"
#include "quark_model.h"
#include "interpolator.h"
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


class phasetransition_class{
public:	

	bool firstRun=true;

  double temperature;
  double muB, muQ;
	// nlwm_class   &hadrons;
	// quarks_class &quarks;

	// phasetransition_class(nlwm_class   &hadrons_, 	quarks_class &quarks_);
	// ~phasetransition_class(void){};

	// void setInitialQH(double &nup1_, double &nun1_, double &mef1_,
	// 								double &nup2_, double &nun2_, double &mef2_);

	void solveQHTransition(double mub_, double pressure_, double temperature_, 
											nlwm_class &hadrons_, quarks_class &quarks_,
											particle &electron_, particle &muon_);
	void solveEqualFlavorFraction(double rhoB_, double temp_, particle &electron_);	

		
};

struct QH_TransitionFunctor{
public:
	QH_TransitionFunctor(nlwm_class & hadrons_, quarks_class &quarks_, 
										particle &electron_, particle &muon_): 
										hadrons(hadrons_),quarks(quarks_), electron(electron_), muon(muon_)
  {}

	template <typename T>
	bool operator()(const T* x, T* residuals) const;
  

private:
	nlwm_class & hadrons;
	quarks_class &quarks;
	particle &electron;
  particle &muon;
};


#endif