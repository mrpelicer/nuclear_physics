#include "pasta.h"

pasta_class::pasta_class(nlwm_class &cluster_, nlwm_class &gas_):cluster(cluster_), gas(gas_){

}

//=============== Set the initial guess for solving the CPA equilibrium ===============
void pasta_class::setInitialCPA(double &nup1_, double &nun1_, double &mef1_,
									double &nup2_, double &nun2_, double &mef2_){

	if(cluster.parametrization=="iufsu"){
    if(YpG>0.3){
      nup1_ = 0.669306;
      nun1_ = 0.669306;
      mef1_ = 0.609329;
      nup2_ = 0.982546; 
      nun2_ = 0.982546;
      mef2_ = 0.999;
    }else if(YpG>0.1){
      nup1_ =0.677603;     
      nun1_ =0.705576;
      mef1_ =0.638736;
      nup2_ =0.946885;
      nun2_ =0.994252;
      mef2_ =0.989749;

    }else{
      //0.70, 0.73, 0.65,0.92, 0.97, 0.95
      // nup1_ = 0.73049;
      // nun1_ = 0.778542;
      // mef1_ = 0.719429 ;
      // nup2_ = 0.819118;
      // nun2_ = 0.866568;
      // mef2_ = 0.829265;
       nup1_ =0.721149;     
       nun1_ =0.767649;
       mef1_ =0.706676;
       nup2_ =0.843373;
       nun2_ =0.892196;
       mef2_ =0.861294;

      // nup1_ = 0.7;
      // nun1_ = 0.73;
      // mef1_ = 0.65 ;
      // nup2_ = 0.92;
      // nun2_ = 0.97;
      // mef2_ = 0.95;
    }
	}else{
    nup1_ = 0.669306;
    nun1_ = 0.669306;
    mef1_ = 0.609329;
    nup2_ = 0.982546; 
    nun2_ = 0.982546;
    mef2_ = 0.999;
  }
}

//=============== Solve CPA equilibrium given dens, yp, temperature ===============
void pasta_class::solveCPA(double rhoB_, double Yp_, double temp_){
								// 
	rhoB=rhoB_;
	YpG=Yp_;
	temperature=temp_;
	cluster.setTemperature(temperature);
	gas.setTemperature(temperature);

  double nup1, nun1, mef1, nup2, nun2, mef2;
  if(firstRun) setInitialCPA(nup1, nun1, mef1, nup2, nun2, mef2);
  else{ nup1 = cluster.proton.chemPot_eff;
        nun1 = cluster.neutron.chemPot_eff;
        mef1 = cluster.proton.mass_eff;
        nup2 = gas.proton.chemPot_eff;
        nun2 = gas.neutron.chemPot_eff;
        mef2 = gas.proton.mass_eff;
  }

	double x[]={nup1, nun1, mef1, nup2, nun2, mef2};
  
	Problem pCPA;
	CostFunction* costCPA =	new NumericDiffCostFunction<cpaFunctor, ceres::CENTRAL, 6, 6>
																														(new cpaFunctor(*this));
	// 
	pCPA.AddResidualBlock(costCPA, NULL, x);

	Solver::Options optionsCPA; 		
	optionsCPA.dense_linear_algebra_library_type=ceres::LAPACK;
	// 
	optionsCPA.parameter_tolerance = 1e-8;	 		//1e-8					
	optionsCPA.function_tolerance = 1e-8;			//1e-6				
	optionsCPA.gradient_tolerance=1e-10;			//1e-10			
	optionsCPA.max_num_iterations=1e3;	
	// 
	optionsCPA.use_nonmonotonic_steps=true;	
	optionsCPA.linear_solver_type= ceres::DENSE_QR;
	//optionsCPA.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //default
	//optionsCPA.dogleg_type = ceres::TRADITIONAL_DOGLEG;// default	
// 
	//optionsCPA.minimizer_type= ceres::LINE_SEARCH;
	//optionsCPA.line_search_direction_type= 	ceres::NONLINEAR_CONJUGATE_GRADIENT;
	//optionsCPA.line_search_type= ceres::ARMIJO;
	//optionsCPA.trust_region_strategy_type = ceres::DOGLEG;
	//optionsCPA.dogleg_type = ceres::SUBSPACE_DOGLEG;

	optionsCPA.minimizer_progress_to_stdout = true;	
	Solver::Summary summaryCPA;
	Solve(optionsCPA, &pCPA, &summaryCPA);

	std::cout << summaryCPA.BriefReport() << "\n";

	 std::cout << "Solution-> "  << nup1 <<  " " << nun1 << " " << mef1 <<  " "
	 					<< nup2 <<  " " <<  nun2 << " " << mef2 << std::endl;
	 std::cout << "----> " << x[0] <<  " " << x[1] << " " << x[2] <<  " "
	 					<< x[3] <<  " " << x[4] << " " << x[5] << std::endl;


	cluster.setEOS_coexistence(x[0], x[1], x[2]);
	gas.setEOS_coexistence(x[3], x[4], x[5]);
	f= (rhoB- gas.rhoB)/(cluster.rhoB-gas.rhoB);
	 
  firstRun=false;
}
 
 //=============== Functor for CPA equilibrium w/ fixed yp ===============
template <typename T>
bool cpaFunctor::operator()(const T* x, T* residuals) const{
		pasta.cluster.setEOS_coexistence(x[0], x[1], x[2]);
		pasta.gas.setEOS_coexistence(x[3], x[4], x[5]);
		double u= (pasta.rhoB- pasta.gas.rhoB)/(pasta.cluster.rhoB - pasta.gas.rhoB);
		// 
		residuals[0] = pasta.cluster.proton.chemPot -  pasta.gas.proton.chemPot;
		residuals[1] = pasta.cluster.neutron.chemPot -  pasta.gas.neutron.chemPot;
		residuals[2] = pasta.cluster.getPressure() - pasta.gas.getPressure();
		residuals[3] = u*pasta.cluster.proton.density + (1.-u)*pasta.gas.proton.density
																																-pasta.YpG*pasta.rhoB;
		residuals[4] = pasta.cluster.sigmaMeson_eom_residue(pasta.cluster.rhoS);
		residuals[5] = pasta.gas.sigmaMeson_eom_residue(pasta.gas.rhoS);
		// 
		return true;
}
 
//=============== Solve CPA in beta equilibrium given dens, temperature ===============
void pasta_class::solveCPA_betaEq(double rhoB_, double temp_, particle &electron_){
							// 
	rhoB=rhoB_;
	temperature=temp_;
	cluster.setTemperature(temperature);
	gas.setTemperature(temperature);
	// 

  double nup1_, nun1_, mef1_, nup2_, nun2_, mef2_, mue_;
  if(firstRun){
    if(temperature<3./Mnucleon){
      nup1_= 0.741031;
      nun1_= 0.790173;
      mef1_= 0.733351;
      nup2_= 0.795927;
      nun2_= 0.843354;
      mef2_= 0.799942;
      mue_ = 0.106509;
    }else{
      nup1_ =0.724612;
      nun1_ =0.768577;
      mef1_ =0.709865;
      nup2_ =0.863269;
      nun2_ =0.913784;
      mef2_ =0.889335;
      mue_  =0.127652;
    }
  }else{
    nup1_=cluster.proton.chemPot_eff;
    nun1_=cluster.neutron.chemPot_eff;
    mef1_=cluster.neutron.mass_eff;
    nup2_=gas.proton.chemPot_eff;
    nun2_=gas.neutron.chemPot_eff;
    mef2_=gas.neutron.mass_eff;
    mue_ =electron_.chemPot;
  }

  double x[]= {nup1_, nun1_, mef1_, nup2_, nun2_, mef2_, mue_};
	Problem pCPA;
	CostFunction* costCPA =	new NumericDiffCostFunction<cpaFunctor_betaEq, ceres::CENTRAL, 7, 7>
																														(new cpaFunctor_betaEq(*this, electron_));
	// 
	pCPA.AddResidualBlock(costCPA, NULL, x);

	Solver::Options optionsCPA; 	
  						//default:		
	optionsCPA.parameter_tolerance = 1e-10;	 		//1e-8					
	optionsCPA.function_tolerance = 1e-8;			//1e-6				
	optionsCPA.gradient_tolerance=1e-10;			//1e-10			
	optionsCPA.max_num_iterations=5e3;	
	// 
	optionsCPA.linear_solver_type= ceres::DENSE_QR;

	optionsCPA.dense_linear_algebra_library_type=ceres::LAPACK;
	// 
	//optionsCPA.trust_region_strategy_type = ceres::DOGLEG;
	//optionsCPA.dogleg_type = ceres::SUBSPACE_DOGLEG;
		// 
	// optionsCPA.use_nonmonotonic_steps= true;
	//optionsCPA.update_state_every_iteration = true;
	//optionsCPA.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //default
	//optionsCPA.dogleg_type = ceres::TRADITIONAL_DOGLEG;// default	
// 
	//optionsCPA.minimizer_type= ceres::LINE_SEARCH;
	//optionsCPA.line_search_direction_type= 	ceres::NONLINEAR_CONJUGATE_GRADIENT;
	//optionsCPA.line_search_type= ceres::ARMIJO;
	//optionsCPA.trust_region_strategy_type = ceres::DOGLEG;
	//optionsCPA.dogleg_type = ceres::SUBSPACE_DOGLEG;
	// 
	// 
	optionsCPA.minimizer_progress_to_stdout = true;	
	// 
	Solver::Summary summaryCPA;
	Solve(optionsCPA, &pCPA, &summaryCPA);
	//
	std::cout << summaryCPA.BriefReport() << "\n";
// 
	std::cout << "Solution-> "  << nup1_ <<  " " << nun1_ << " "  << mef1_ << " " 
            << nup2_ <<  " " << nun2_ << " "   << mef2_ << " " << mue_ 
	          <<  std::endl;
	// 
	std::cout << "----> " << x[0] <<  " " << x[1] << " " << x[2] <<  " "
	 					<< x[3] <<  " " << x[4] << " " << x[5] << " " << x[6] << std::endl;
// 
	cluster.setEOS_coexistence(x[0], x[1], x[2]);
	gas.setEOS_coexistence(x[3], x[4], x[5]);
	f= (rhoB- gas.rhoB)/(cluster.rhoB-gas.rhoB);
	electron_.setLepton(x[6]);
	electron_.calculateProperties();
  YpG= (f*(cluster.proton.density - gas.proton.density) + gas.proton.density)
        /(f*(cluster.rhoB - gas.rhoB) + gas.rhoB);
      
  firstRun=false;

// 
}

 //=============== Functor for CPA beta-equilibrium  ===============
template <typename T>
bool cpaFunctor_betaEq::operator()(const T* x, T* residuals) const{
		pasta.cluster.setEOS_coexistence(x[0], x[1], x[2]);
		pasta.gas.setEOS_coexistence(x[3], x[4], x[5]);
		electron.setLepton(x[6]);
		electron.calculateProperties();
		//electron
		double u= (pasta.rhoB- pasta.gas.rhoB)/(pasta.cluster.rhoB - pasta.gas.rhoB);
		// 
		residuals[0] = pasta.cluster.proton.chemPot -  pasta.gas.proton.chemPot;
		residuals[1] = pasta.cluster.neutron.chemPot -  pasta.gas.neutron.chemPot;
		residuals[2] = pasta.cluster.getPressure() - pasta.gas.getPressure();
		residuals[3] = u*pasta.cluster.proton.density + (1.-u)*pasta.gas.proton.density-electron.density;
		residuals[4] = pasta.cluster.sigmaMeson_eom_residue(pasta.cluster.rhoS);
		residuals[5] = pasta.gas.sigmaMeson_eom_residue(pasta.gas.rhoS);
		residuals[6]	=electron.chemPot + pasta.cluster.proton.chemPot - pasta.cluster.neutron.chemPot;
		return true;
}

//=============== Solve LCD in equilibrium given dens, yp, temperature ===============
void pasta_class::solveCLD(double rhoB_, double Yp_, double temp_,double dim_, int itype_){

	rhoB=rhoB_;
	YpG=Yp_;
	temperature=temp_;
	dim=dim_;
	iType=itype_;
	cluster.setTemperature(temperature);
	gas.setTemperature(temperature);

  double nup1, nun1, mef1, nup2, nun2, mef2;
  if(firstRun) setInitialCPA(nup1, nun1, mef1, nup2, nun2, mef2);
  else{ nup1 = cluster.proton.chemPot_eff;
        nun1 = cluster.neutron.chemPot_eff;
        mef1 = cluster.proton.mass_eff;
        nup2 = gas.proton.chemPot_eff;
        nun2 = gas.neutron.chemPot_eff;
        mef2 = gas.proton.mass_eff;
  }

	double x[]={nup1, nun1, mef1, nup2, nun2, mef2};
  
	Problem pCLD;
	CostFunction* costCLD =	new NumericDiffCostFunction<cldFunctor, ceres::CENTRAL, 6, 6>
																														(new cldFunctor(*this));
	pCLD.AddResidualBlock(costCLD, NULL, x);
	//pCLD.SetParameterLowerBound(x, 0, 0.);
	//pCLD.SetParameterLowerBound(x, 1, 0.);
	//pCLD.SetParameterLowerBound(x, 2, 0.);
	//pCLD.SetParameterUpperBound(x, 2, 1.);
//
	//pCLD.SetParameterLowerBound(x, 3, 0.);
	//pCLD.SetParameterLowerBound(x, 4, 0.);
	//pCLD.SetParameterLowerBound(x, 5, 0.);
	//pCLD.SetParameterUpperBound(x, 5, 1.);
//
	// 
	Solver::Options optionsCLD; 							//default:	0.5 fluc	
	optionsCLD.parameter_tolerance = 1e-8; 		//1e-8			9 			
	optionsCLD.function_tolerance = 1e-6;			//1e-6			8	
	optionsCLD.gradient_tolerance=1e-10;			//1e-10			12
// 
	optionsCLD.max_num_iterations=1e3;
	// optionsCLD.use_nonmonotonic_steps=true;	
	optionsCLD.linear_solver_type= ceres::DENSE_QR;
	//optionsCLD.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
	//optionsCLD.dogleg_type = ceres::SUBSPACE_DOGLEG;
	optionsCLD.minimizer_progress_to_stdout = true;
	//optionsCLD.trust_region_strategy_type = ceres::DOGLEG;
	//optionsCLD.dogleg_type = ceres::TRADITIONAL_DOGLEG;
	// 
	Solver::Summary summaryCLD;
	Solve(optionsCLD, &pCLD, &summaryCLD);
	//
	std::cout << summaryCLD.BriefReport() << "\n";
// 
	std::cout << "rhoB= " << rhoB*pow(Mnucleon/hc, 3.) << std::endl;
	 std::cout << "Solution (" << dim << " , " << iType << ") : "  << nup1 <<  " " << nun1 << " " << mef1 <<  " "
	 					<< nup2 <<  " " << nun2 << " " << mef2 << std::endl;
	 std::cout << "----> " << x[0] <<  " " << x[1] << " " << x[2] <<  " "
	 					<< x[3] <<  " " << x[4] << " " << x[5] << std::endl;
// 
	cluster.setEOS_coexistence(x[0], x[1], x[2]);
	gas.setEOS_coexistence(x[3], x[4], x[5]);
	f= (rhoB- gas.rhoB)/(cluster.rhoB-gas.rhoB);
	// 
// 
	firstRun=false;
}

 //=============== Functor for CDL equilibrium  w/ fixed yp===============
template <typename T>
bool cldFunctor::operator()(const T* x, T* residuals) const{
		pasta.cluster.setEOS_coexistence(x[0], x[1], x[2]);
		pasta.gas.setEOS_coexistence(x[3], x[4], x[5]);
		// 
		double u= (pasta.rhoB- pasta.gas.rhoB)/(pasta.cluster.rhoB - pasta.gas.rhoB);
		double beta=u;
		double sign_=1.;
	// 
		if(pasta.iType==1){beta=(1.-u); sign_=-1.;}
		// 
		double sigma=getSurfaceTension(pasta.cluster, pasta.YpG, pasta.temperature);
		double Phi=getPhiFunc(pasta.dim, beta);
		double PhiD=getPhiFuncDerivative(pasta.dim, beta);
		// 
		double Rd=getRadiusD(pasta.dim, beta, pasta.YpG, pasta.cluster, pasta.gas);
		double Fs= sigma*pasta.dim/Rd;
		double Fc= Fs/2.;
		double Fsc= Fc+Fs;
		// 
		double interfaceMuP=0.;
		double interfacePrs=0.;
// 
		//std::cout << Rd << " " << Fc << " " << u << " " << Phi << " " << PhiD << std::endl;
		if(u>0. || u<1.){
			interfaceMuP= 2.*beta*Fc/(u*(1.-u)*(pasta.cluster.proton.density-pasta.gas.proton.density));
// 
			interfacePrs= sign_*(Fsc+beta*Fc*PhiD/Phi)
										-2.*beta*Fc*((1.-u)*pasta.cluster.proton.density + u*pasta.gas.proton.density)/
											(u*(1.-u)*(pasta.cluster.proton.density - pasta.gas.proton.density));	
	/*		
			interfacePrs= Fc*(1. +beta*PhiD/Phi- 2.*beta*pasta.gas.proton.density/
												(u*(1.-u)*(pasta.cluster.proton.density - pasta.gas.proton.density)));	
		*/	
		}
				// 
		residuals[0] = pasta.cluster.proton.chemPot - pasta.gas.proton.chemPot+ interfaceMuP;						
		residuals[1] = pasta.cluster.neutron.chemPot -  pasta.gas.neutron.chemPot ;								
		residuals[2] =-pasta.cluster.getPressure() + pasta.gas.getPressure() + interfacePrs;
		residuals[3] = u*pasta.cluster.proton.density+(1.-u)*pasta.gas.proton.density-pasta.YpG*pasta.rhoB ;
		residuals[4] = pasta.cluster.sigmaMeson_eom_residue(pasta.cluster.rhoS);
		residuals[5] = pasta.gas.sigmaMeson_eom_residue(pasta.gas.rhoS);
		// 
		return true;
}
// 
// 
double getSurfaceTension(nlwm_class &cluster_, double Yp_, double T){
// 
  double sigma0=0., sigma1, sa1, sa2, sa3, sa4, sa5, sa6;
  double aa0, aa1, aa2, aa3, aa4, aa5;
  double ba0, ba1, ba2, ba3, ba4, ba5;
  double ca0, ca1, ca2, ca3, ca4, ca5;
// 
	setSurfaceParameters( cluster_, sigma0, sigma1, sa1, sa2, sa3, sa4, sa5, sa6,
																			  aa0, aa1, aa2, aa3, aa4, aa5,
 																			  ba0, ba1, ba2, ba3, ba4, ba5,
 																			  ca0, ca1, ca2, ca3, ca4, ca5);
// 
	// 
	double x=pow(1.-2.*Yp_, 2.);
	double sigma_= sigma0*exp(-sigma1*pow(x, 1.5))
                      *(1.+ sa1*x          + sa2*pow(x, 2.)
                          + sa3*pow(x, 3.) + sa4*pow(x, 4.)
                          + sa5*pow(x, 5.) + sa6*pow(x, 6.)
						);
	if(T!=0.){
		T*=Mnucleon;
		double aT= aa0 + aa1*T +aa2*pow(T, 2.) + aa3*pow(T, 3.)
                           + aa4*pow(T, 4.) + aa5*pow(T, 5.);
// 
		double bT= ba0 + ba1*T +ba2*pow(T, 2.) + ba3*pow(T, 3.)
                           + ba4*pow(T, 4.) + ba5*pow(T, 5.);
// 
		double cT= ca0 + ca1*T +ca2*pow(T, 2.) + ca3*pow(T, 3.)
                           + ca4*pow(T, 4.) + ca5*pow(T, 5.);
// 
		sigma_*=(1.- aT*x*T - bT*pow(T, 2.)- cT*T*pow(x, 2.) );
	}
// 
   return sigma_; // abs only for very small yp
}
// 
// 
double getSurfaceTensionDerivative(nlwm_class &cluster_, double Yp_, double T){
	// 
	double x=pow(1.-2.*Yp_, 2.);
// 
  double sigma0=0., sigma1, sa1, sa2, sa3, sa4, sa5, sa6;
  double aa0, aa1, aa2, aa3, aa4, aa5;
  double ba0, ba1, ba2, ba3, ba4, ba5;
  double ca0, ca1, ca2, ca3, ca4, ca5;
	setSurfaceParameters(cluster_, sigma0, sigma1, sa1, sa2, sa3, sa4, sa5, sa6,
																			  aa0, aa1, aa2, aa3, aa4, aa5,
 																			  ba0, ba1, ba2, ba3, ba4, ba5,
 																			  ca0, ca1, ca2, ca3, ca4, ca5);
// 
	double sigma_= sigma0*exp(-sigma1*pow(x, 1.5));
	double SigmaXT=	1.;
	double DSigmaXT=0.;
	// 
	double Px= (1.+ sa1*x + sa2*pow(x, 2.)
								+ sa3*pow(x, 3.) + sa4*pow(x, 4.)
                + sa5*pow(x, 5.) + sa6*pow(x, 6.)
						);
	double DPx= (sa1 + 2.*sa2*x
									 + 3.*sa3*pow(x, 2.) + 4.*sa4*pow(x, 3.)
                	 + 5.*sa5*pow(x, 4.) + 6.*sa6*pow(x, 5.)
							);
							// 
	double Dsigma_= -4.*(1.-2.*Yp_)*sigma_*(DPx*SigmaXT
																-3.*Px*SigmaXT*sigma1*sqrt(x)/2.);
																// 
	if(T!=0.){
		T*=Mnucleon;
		double aT= aa0 + aa1*T +aa2*pow(T, 2.) + aa3*pow(T, 3.)
                           + aa4*pow(T, 4.) + aa5*pow(T, 5.);
// 
		double bT= ba0 + ba1*T +ba2*pow(T, 2.) + ba3*pow(T, 3.)
                           + ba4*pow(T, 4.) + ba5*pow(T, 5.);
// 
		double cT= ca0 + ca1*T +ca2*pow(T, 2.) + ca3*pow(T, 3.)
                           + ca4*pow(T, 4.) + ca5*pow(T, 5.);
// 
		SigmaXT= (1.- aT*x*T-bT*pow(T, 2.)- cT*T*pow(x, 2.) );
		DSigmaXT=-(aT*T+2.*cT*T*x );
		Dsigma_= -4.*(1.-2.*Yp_)*sigma_*(DPx*SigmaXT + Px*DSigmaXT
																-3.*Px*SigmaXT*sigma1*sqrt(x)/2.);
	}
// 
   return Dsigma_;
}
// 
// 
double getPhiFunc(double dim_, double u_){
// 
  double phi=0.;
  double Lu=log(u_);
// 
  if(dim_!=2.){phi=( (2.-dim_*pow(u_, 1.-2./dim_))/(dim_-2.) +u_ )/(dim_+2.);}
  else {phi=(u_-1.-Lu)/(dim_+2.);}
// 
  return phi;
}
// 
double getPhiFuncDerivative(double dim_, double u_){
	double phiD=0.;
// 
  if(dim_!=2.){phiD=( (dim_*(2./dim_-1.)*pow(u_, -2./dim_) )/(dim_-2.) +1.)/(dim_+2.);}
  else {phiD=(1.- 1/u_)/(dim_+2.);}
	return phiD;
}
// 
double getRadiusD(double dim_, double u_, double Yp_, nlwm_class &cluster_, 
																											nlwm_class &gas_){
// 
	double sigma_= getSurfaceTension(cluster_, Yp_, cluster_.temperature);
	double PhiFunc_=getPhiFunc(dim_, u_);
	double rd_=pow( (sigma_*dim_)/(4.*M_PI*pow(eGS, 2.)*pow((cluster_.proton.density
												 -gas_.proton.density), 2.)*PhiFunc_), 1./3.);

	return  rd_;
}
// 
// 
void setSurfaceParameters(nlwm_class &cluster_, double &sigma0, double &sigma1, 
			double &sa1, double &sa2, double &sa3, double &sa4, double &sa5, double &sa6,
  		double &aa0, double &aa1, double &aa2, double &aa3, double &aa4, double &aa5,
  		double &ba0, double &ba1, double &ba2, double &ba3, double &ba4, double &ba5,
  		double &ca0, double &ca1, double &ca2, double &ca3, double &ca4, double &ca5){

// parameters provided by:
  if(cluster_.parametrization=="nl3")
	{
    sigma0=1.12307/( pow(Mnucleon, 3.)/pow(hc, 2.) );
    sigma1= 20.7779;
    sa1= -5.84915;
    sa2= 138.839;
    sa3= -1631.42;
    sa4= 8900.34;
    sa5=-21592.3;
    sa6=20858.6;
// 
    aa0= 0.0121222;
    aa1= 0.01664;
    aa2= -0.00137266;
    aa3= 4.0257e-5;
    aa4=0.;
    aa5=0.;
// 
    ba0= 0.00792168;
    ba1=-8.2504e-5;
    ba2=-4.59336e-6;
    ba3=-2.81679e-7;
    ba4= 0.;
    ba5=0.;
// 
    ca0=0.;
    ca1=0.;
    ca2=0.;
    ca3=0.;
    ca4=0.;
    ca5=0.;
  }
	// ==== NL3wr ====
	if(cluster_.parametrization=="nl3wr")
	{
	  sigma0=1.12013/( pow(Mnucleon, 3.)/pow(hc, 2.) );
    sigma1= 14.0774;
    sa1=-2.15376;
    sa2= 57.8455;
    sa3=-431.365;
    sa4= 1854.81 ;
    sa5=-3653.96;
    sa6= 3214.82;
// 
    aa0=-5.80451e-5;
    aa1= 0.0233833;
    aa2=-0.00507732;
    aa3= 0.000490863;
    aa4=-1.59473e-5;
    aa5=-7.55062e-8;
// 
    ba0= 0.00725961;
    ba1= 0.000318409;
    ba2=-0.000104941;
    ba3= 1.19645e-5;
    ba4=-7.19099e-7;
    ba5= 1.62087e-8;
// 
    ca0=-0.00259094;
    ca1=-0.053756;
    ca2= 0.0114598;
    ca3=-0.000354375;
    ca4=-4.76451e-5;
    ca5= 2.28389e-6;
	}
// ===FSUGold ====
	if(cluster_.parametrization=="fsu")
  {
    sigma0= 1.1223/( pow(Mnucleon, 3.)/pow(hc, 2.) );
    sigma1=-1.45717;
    sa1=-3.17729;
    sa2=-9.5121;
    sa3= 70.5609;
    sa4=-155.641;
    sa5= 154.691;
    sa6=-58.9476;
// 
    aa0=-0.0133789;
    aa1= 0.0330912;
    aa2=-0.00786564;
    aa3= 0.000902286;
    aa4=-4.84828e-5;
    aa5= 9.56728e-7;
// 
    ba0= 0.00773356;
    ba1=-0.000240406;
    ba2= 4.52523e-5;
    ba3=-7.64893e-6;
    ba4= 5.33346e-7;
    ba5=-1.45394e-8;
// 
    ca0= 0.0408077;
    ca1=-0.0971609;
    ca2= 0.0195288;
    ca3=-0.00140166;
    ca4= 4.97386e-5;
    ca5=-1.20803e-6;
  }
	// 
// ===IU-FSU ====
  if(cluster_.parametrization=="iufsu")
  {
    sigma0 = 1.16473/( pow(Mnucleon, 3.)/pow(hc, 2.) );
    sigma1 =-0.659167;
    sa1    =-2.25482;
    sa2    =-5.64237;
    sa3    = 37.8471;
    sa4    =-81.6617;
    sa5    = 81.2696;
    sa6    =-31.0227;
// 
    aa0= 0.00404325;
    aa1= 0.00828207;
    aa2=-0.00153301;
    aa3= 7.26763e-5;
    aa4=0.;
    aa5=0.;
// 
    ba0= 0.00767923;
    ba1=-8.58068e-5;
    ba2= 4.43918e-7;
    ba3=-5.44453e-7;
    ba4= 0.;
    ba5= 0.;
// 
    ca0= 0.0066774;
    ca1=-0.0514285;
    ca2= 0.00949505;
    ca3=-0.000427613;
    ca4= 0.;
    ca5= 0.;
  }
// 
}

pasta_transport_class::pasta_transport_class(particle electron_){
  electron=electron_;
  xr    = electron.kf/electron.mass;
  enerf = sqrt(pow(electron.mass, 2.)+ pow(electron.kf, 2.));
  gammar=sqrt(1.+pow(xr, 2.));
  betar = xr/gammar;          //=kf/enerf;
  ae=pow(3./(4.*M_PI*electron.density), 1./3.);
  kTF=0.185/(ae*sqrt(betar));

}

void pasta_transport_class::setMomentum(double q_){
  q=q_;
};


void pasta_transport_class::setMomentumVec(std::vector<double> qv_){
  qv=qv_;
}

void pasta_transport_class::setFormFactor2Vec(std::vector<double> F2_){
  F2v= F2_;
  assert(qv.size()==F2v.size());
}

double pasta_transport_class::q2dielectric_function(double q_){
  double y_=q_/(2.*electron.kf);
  double t1= - 2.*pow(y_, 2.)*xr*log(xr+gammar)/(3.*gammar)  ;
  double t2= (pow(xr, 2.) +1. - 3.*pow(xr*y_, 2.))*log(fabs( (1.+y_)/(1.-y_)))/(6.*y_*xr*xr) ;
  double t3= ((2.*pow(y_*xr, 2.)-1.)*sqrt(1.+pow(xr*y_, 2.) )*
              log(fabs( (y_*gammar+sqrt(1.+pow(xr*y_, 2.) ) )/(y_*gammar-sqrt(1.+pow(xr*y_, 2.))) ))
                                                                                      /(6.*y_*xr*xr*gammar));

  double t= (y_<1e-7) ? 1. : (2./3. +t1+t2+t3) ;
  double eps_= pow(q_, 2.)+ pow(kTF, 2.)*t;
  return eps_;
}

void cDroplet::setRadius(double r_){
  radius=r_; 
}

double cDroplet::getVolume(){
  return 4.*M_PI*pow(radius, 3.)/3.;
}

double cDroplet::getArea(){
  return 2.*M_PI*pow(radius, 2.);
}

double cDroplet::getStructureFunction(double q_){
  double a_=q_*radius;
  a=a_;
  double F_= a_>1e-5 ? 3.*(sin(a_) - a_*cos(a_))/pow(a_, 3) : 1.;
  return F_;
}

void cRod::setLengths(double r_, double h_){
  radius=r_; length=h_;
}

double cRod::getVolume(){
  return M_PI*pow(radius, 2.)*length;
}

double cRod::getArea(){
  return 2.*M_PI*pow(radius, 2.) + 2*M_PI*radius*length;
}

void cRod::setLim_CosTheta(double min_cost_, double max_cost_){
  min_cost=min_cost_;
  max_cost=max_cost_;
}

double cRod::getStructureFunction(double q_, double cost_){
  double x=cost_;
  double y=sqrt(1.-x*x);
  double az= q_*x*length/2.;
  double ar= q_*y*radius;
  double Fz= (az>1e-7) ? sin(az)/az : 1.;
  double Fp= (ar>1e-7) ? 2.*gsl_sf_bessel_J1(ar)/(ar)  : 1.;
  return Fz*Fp;
}


double cRod::getStructureFunction2_Axial(double q_){
  setMomentum(q_);
  cubareal integral_rod[NCOMP],  error_rod[NCOMP], prob_rod[NCOMP];
  int neval, fail;
  // int comp;
  Vegas(1, NCOMP, IF2_rod_axial, this, NVEC,
    EPSREL, EPSABS, VERBOSE_AVG, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral_rod, error_rod, prob_rod);

  // printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
  // for(int comp = 0; comp < NCOMP; ++comp )
  //   printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
  //     (double)integral_rod[comp], (double)error_rod[comp], (double)prob_rod[comp]);

  return integral_rod[0];
}

double cRod::getStructureFunction2_Trans(double q_){
  setMomentum(q_);
  cubareal integral_rod[NCOMP],  error_rod[NCOMP], prob_rod[NCOMP];
  int neval, fail;
  // int comp;
  Vegas(1, NCOMP, IF2_rod_trans, this, NVEC,
    EPSREL, EPSABS, VERBOSE_AVG, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral_rod, error_rod, prob_rod);

  // printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
  // for(int comp = 0; comp < NCOMP; ++comp )
  //   printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
  //     (double)integral_rod[comp], (double)error_rod[comp], (double)prob_rod[comp]);

  return integral_rod[0];
}


int IF2_rod_axial(const int *ndim, const cubareal xx[],const int *ncomp, 
                     cubareal ff[], void *userdata) {
  cRod &rod= *reinterpret_cast<cRod*>(userdata);
  double range=rod.max_cost - rod.min_cost;
  double x = range*xx[0] +rod.min_cost;
  //double y = xx[1];
  //double z = xx[2];
  ff[0] = range*pow(x, 2)*pow(rod.getStructureFunction(rod.q, x), 2.)/2.;
  return 0;
}

int IF2_rod_trans(const int *ndim, const cubareal xx[],const int *ncomp, 
                     cubareal ff[], void *userdata) {
  cRod &rod= *reinterpret_cast<cRod*>(userdata);
  double range=rod.max_cost - rod.min_cost;
  double x = range*xx[0] +rod.min_cost;
  //double y = xx[1];
  //double z = xx[2];
  ff[0] = range*(1.-pow(x, 2) )*pow(rod.getStructureFunction(rod.q, x), 2.)/4.;
  return 0;
}



void cSlab::setLengths(double lx_, double ly_, double lz_){
  lx=lx_; ly=ly_; lz=lz_;
}


double cSlab::getVolume(){
  return lx*ly*lz;
}

double cSlab::getArea(){
  return 2.*(lx*ly+lx*lz+ly*lz);
}

void cSlab::setLim_CosTheta(double min_cost_, double max_cost_){
  min_cost=min_cost_;
  max_cost=max_cost_;
}

void cSlab::setLim_Phi(double min_phi_, double max_phi_){
  min_phi=min_phi_;
  max_phi=max_phi_;
}

double cSlab::getStructureFunction(double q_, double cost_, double phi_){
  double theta_= acos(cost_);
  double qx= q_*cos(phi_)*sin(theta_);
  double qy= q_*sin(phi_)*sin(theta_);
  double qz= q_*cost_;
  
  double ax= qx*lx/2.;
  double ay= qy*ly/2.;
  double az= qz*lz/2.;
  double Fx= ax>1e-7 ? sin(ax)/(ax) : 1.;
  double Fy= ay>1e-7 ? sin(ay)/(ay) : 1.;
  double Fz= az>1e-7 ? sin(az)/(az) : 1.;
  
  return Fx*Fy*Fz;
}

double cSlab::getStructureFunction2_Axial(double q_){
  setMomentum(q_);
  cubareal integral_slab[NCOMP],  error_slab[NCOMP], prob_slab[NCOMP]; 
  int  neval, fail;
  Vegas(2, NCOMP, IF2_slab_axial, this, NVEC,
    EPSREL, EPSABS, VERBOSE_AVG, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral_slab, error_slab, prob_slab);

  
  // for(int comp = 0; comp < NCOMP; ++comp )
  //    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
  //      (double)integral_slab[comp], (double)error_slab[comp], (double)prob_slab[comp]);
  return integral_slab[0];
}

double cSlab::getStructureFunction2_Trans(double q_){
  setMomentum(q_);
  cubareal integral_slab[NCOMP],  error_slab[NCOMP], prob_slab[NCOMP]; 
  int neval, fail;
  Vegas(2, NCOMP, IF2_slab_trans, this, NVEC,
    EPSREL, EPSABS, VERBOSE_AVG, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral_slab, error_slab, prob_slab);

  // for(int comp = 0; comp < NCOMP; ++comp )
  //    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
  //      (double)integral_slab[comp], (double)error_slab[comp], (double)prob_slab[comp]);
  return integral_slab[0];
}

int IF2_slab_axial(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata) {
  cSlab &slab= *reinterpret_cast<cSlab*>(userdata);
  double range_cost =slab.max_cost  - slab.min_cost;
  double range_phi  =slab.max_phi   - slab.min_phi;
  double x0 = range_cost*xx[0]  +slab.min_cost;
  double x1= range_phi*xx[1]    +slab.min_phi;
  ff[0] = range_cost*range_phi*pow(x0, 2)*pow(slab.getStructureFunction(slab.q, x0, x1), 2.)/(4.*M_PI);
  return 0;
}
int IF2_slab_trans(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata) {
  cSlab &slab= *reinterpret_cast<cSlab*>(userdata);
  double range_cost =slab.max_cost  - slab.min_cost;
  double range_phi  =slab.max_phi   - slab.min_phi;
  double x0 = range_cost*xx[0]  +slab.min_cost;
  double x1= range_phi*xx[1]    +slab.min_phi;
  ff[0] = range_cost*range_phi*(1.-pow(x0, 2))*pow(slab.getStructureFunction(slab.q, x0, x1), 2.)/(8.*M_PI);
  return 0;
}


double getUnintegratedLambda_ei(double q_, pasta_transport_class pasta_){

  double q2eps= pasta_.q2dielectric_function(q_);
  double f2_= interpolation_func(q_, pasta_.F2v, pasta_.qv);
  double ef_= hypot(pasta_.electron.kf, pasta_.electron.mass);
  return pow(q_, 3.)*(1.-q_*q_/(4.*ef_*ef_))*f2_/(q2eps*q2eps);


}

//cuba integration

double getCoulombIntegral(pasta_transport_class input_){
  cubareal integral_tau[NCOMP],  error_tau[NCOMP], prob_tau[NCOMP];
  int comp, neval, fail;
  pasta_transport_class *par_ = &input_;

  Vegas(1, 1, Integral_Coulomb, par_, NVEC,
        EPSREL, EPSABS, VERBOSE_CLOG, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
       GRIDNO, STATEFILE, SPIN,
       &neval, &fail, integral_tau, error_tau, prob_tau);
  
  printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral_tau[comp], (double)error_tau[comp], (double)prob_tau[comp]);
  
  return integral_tau[0];
}

int Integral_Coulomb(const int *ndim, const cubareal xx[],const int *ncomp, 
                      cubareal ff[], void *userdata) {
  pasta_transport_class &par= *reinterpret_cast<pasta_transport_class*>(userdata);
  double range =2.*par.electron.kf  - par.q0;
  double q_ = range*xx[0]  +par.q0;

  ff[0] = range*getUnintegratedLambda_ei(q_,par);
  return 0;
}



double integrate_coulomb(double (func)(double, void *), void *parametersPointer)
{
  double result = 0e0;
  double er_res = 0e0;
  double err_abs = 1e-5; //1e-13
  double err_rel = 1e-8; //1e-10;
  size_t max_steps = 1e8;
  pasta_transport_class &par= *reinterpret_cast<pasta_transport_class*>(parametersPointer);
  double min=par.q0;
  double max=2.*par.electron.kf;
  
  //gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
	
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);

	gsl_function My_function;
  My_function.function = func;
  My_function.params = parametersPointer;
	
	// gsl_integration_qag(&My_function, 0., max, err_abs, err_rel, max_steps, 1,w, &result, &er_res);
  double alpha=0.0;
  double beta=0.0;
  int mu=0;
  int nu=0;
  gsl_integration_qaws_table * table=gsl_integration_qaws_table_alloc(alpha,beta,mu,nu);
	gsl_integration_qaws(&My_function, min, max, table, err_abs, err_rel, max_steps, w, &result, &er_res);
  gsl_integration_qaws_table_free(table);
  
  //gsl_set_error_handler(old_handler);
	gsl_integration_workspace_free(w);
  return result;
}

double coulomb_gsl(double q_, void *p){

  pasta_transport_class &par= *reinterpret_cast<pasta_transport_class*>(p);
  // double eps= par.q2dielectric_function(q_);

  // double f2_= interpolation_func(q_, par.F2v, par.qv);
  double I1 = getUnintegratedLambda_ei(q_, par);
  return I1;
}