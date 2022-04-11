#include "quark_hadron_transition.h"

// phasetransition_class::phasetransition_class(nlwm_class &hadrons_, quarks_class &quarks_):
//                               hadrons(hadrons_), quarks(quarks_){

// }

//=============== Functor beta-equilibrium for Ceres w/ 3 mesons: ===============

void phasetransition_class::solveQHTransition(double mub_, double pressure_,double temperature_, 
																						nlwm_class &hadrons, quarks_class &quarks,
																						particle &electron_, particle &muon_){
  temperature=temperature_;
  hadrons.setTemperature(temperature_);
  quarks.setTemperature(temperature_);
	electron_.temperature=temperature;
	muon_.temperature=temperature;
	
  double mue_, phi0_, v0_, b0_, theta0_=0.;
	
  if(firstRun){
		mub_=1.67397;
		mue_=0.0795405;
		phi0_=0.110646;
		v0_=0.120147;
		b0_=-0.00466873;
		if(hadrons.parametrization=="fsu2h" || hadrons.parametrization=="l3wr") theta0_ = -0.0118828;
	}else{
    mub_    = muB;
    mue_    = electron_.chemPot;
    phi0_   = hadrons.phi0;
		v0_     = hadrons.V0;
		b0_     = hadrons.b0;
		if(hadrons.parametrization=="fsu2h" || hadrons.parametrization=="l3wr") theta0_=hadrons.theta0;
  }

	if(hadrons.parametrization=="fsu2h" || hadrons.parametrization=="l3wr"){ 

		double x[]={mub_, mue_, phi0_, v0_, b0_, theta0_};
	
		Problem pQHT;
		CostFunction* costQHT =	new NumericDiffCostFunction<QH_TransitionFunctor, ceres::CENTRAL, 6, 6>
																									(new QH_TransitionFunctor(hadrons, quarks, electron_, muon_));
		// 
		pQHT.AddResidualBlock(costQHT, NULL, x);

		Solver::Options optionQHT; 		
		optionQHT.dense_linear_algebra_library_type=ceres::LAPACK;
		// 
		optionQHT.parameter_tolerance = 1e-7;	 		//1e-8					
		optionQHT.function_tolerance = 1e-8;			//1e-6				
		optionQHT.gradient_tolerance=1e-10;			//1e-10			
		optionQHT.max_num_iterations=1e3;	
		// 
		optionQHT.use_nonmonotonic_steps=true;	
		optionQHT.linear_solver_type= ceres::DENSE_QR;
		//optionQHT.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //default
		//optionQHT.dogleg_type = ceres::TRADITIONAL_DOGLEG;// default	
		// 
		//optionQHT.minimizer_type= ceres::LINE_SEARCH;
		//optionQHT.line_search_direction_type= 	ceres::NONLINEAR_CONJUGATE_GRADIENT;
		//optionQHT.line_search_type= ceres::ARMIJO;
		//optionQHT.trust_region_strategy_type = ceres::DOGLEG;
		//optionQHT.dogleg_type = ceres::SUBSPACE_DOGLEG;

		optionQHT.minimizer_progress_to_stdout = true;	
		Solver::Summary summaryCPA;
		Solve(optionQHT, &pQHT, &summaryCPA);

		std::cout << summaryCPA.BriefReport() << "\n";

		std::cout << mub_ << " " << mue_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ <<  " " << theta0_ << 
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]  << " " << x[5]
		<< std::endl << std::endl;

		mub_	 =x[0];
		mue_	 =x[1];
		phi0_	 =x[2];
		v0_		 =x[3];
		b0_		 =x[4];
		theta0_=x[5];
  hadrons.setDensities(x[0], -x[1],  x[2],  x[3], x[4], x[5]);
	electron_.setLepton(x[1]);
  muon_.setLepton(x[1]);

	double Yu= (2.*hadrons.proton.density + hadrons.neutron.density
							+hadrons.lambda0.density  + 2.*hadrons.sigmap.density+ hadrons.sigma0.density
							+hadrons.xi0.density)/(3.*hadrons.getBaryonDens()); 

	double Yd= (hadrons.proton.density + 2.*hadrons.neutron.density
							+hadrons.lambda0.density  +  hadrons.sigma0.density + 2.*hadrons.sigmam.density
							+hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							
	
	double Ys= (hadrons.lambda0.density  +  hadrons.sigmap.density + hadrons.sigma0.density
							+hadrons.sigmam.density +2.*hadrons.xi0.density
							+2.*hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							
							
	quarks.setEoSFlavorFixed(hadrons.getBaryonDens(), hadrons.temperature,
															Yu, Yd, Ys);


	  firstRun=false;

	}else{

	}

}
//  //=============== Functor for beta-eq transition ===============
template <typename T>
bool QH_TransitionFunctor::operator()(const T* x, T* residuals) const{
    
  hadrons.setDensities(x[0], -x[1],  x[2],  x[3], x[4], x[5]);
	electron.setLepton(x[1]);
  muon.setLepton(x[1]);

	double Yu= (2.*hadrons.proton.density + hadrons.neutron.density
							+hadrons.lambda0.density  + 2.*hadrons.sigmap.density+ hadrons.sigma0.density
							+hadrons.xi0.density)/(3.*hadrons.getBaryonDens()); 

	double Yd= (hadrons.proton.density + 2.*hadrons.neutron.density
							+hadrons.lambda0.density  +  hadrons.sigma0.density + 2.*hadrons.sigmam.density
							+hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							
	
	double Ys= (hadrons.lambda0.density  +  hadrons.sigmap.density + hadrons.sigma0.density
							+hadrons.sigmam.density +2.*hadrons.xi0.density
							+2.*hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							
							
	quarks.setEoSFlavorFixed(hadrons.getBaryonDens(), hadrons.temperature,
															Yu, Yd, Ys);

	residuals[0] = hadrons.getChargeDens() +	electron.Qdens + muon.Qdens;
	residuals[1] = hadrons.sigmaMeson_eom_residue(	hadrons.getSigmaEffDens());
	residuals[2] = hadrons.omegaMeson_eom_residue(	hadrons.getOmegaEffDens());
	residuals[3] = hadrons.rhoMeson_eom_residue(		hadrons.getIsoEffDens());
	residuals[4] = hadrons.thetaMeson_eom_residue(	hadrons.getThetaEffDens());
	residuals[5] = hadrons.getPressure() - quarks.getPressure();
	// residuals[5] = hadrons.muB - quarks.muB;

		return true;
}
