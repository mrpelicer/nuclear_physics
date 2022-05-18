#include "quark_model.h"

quarks_class::quarks_class(){
  qu.mass=5./Mnucleon;       qu.gamma = 6.;    qu.Q=  2./3.;
  qd.mass=10./Mnucleon;      qd.gamma = 6.;    qd.Q= -1./3.;
  qs.mass=100./Mnucleon;     qs.gamma = 6.;    qs.Q= -1./3.;  
}

void quarks_class::setParameters(double C_, double D_, double tcrit_){
  C=C_;
  D=D_;
  tcrit=tcrit_;
}


void quarks_class::setDensities(double muef_u, double muef_d, double muef_s){
	qu.chemPot_eff =muef_u;
  qd.chemPot_eff =muef_d;
  
  qu.kf2=pow(qu.chemPot_eff, 2.)- pow(qu.mass_eff, 2.);	
  qu.kf2<=0. ? qu.kf=0. : qu.kf=sqrt(qu.kf2);
  
  qd.kf2=pow(qd.chemPot_eff, 2.)- pow(qd.mass_eff, 2.);	
  qd.kf2<=0. ? qd.kf=0. : qd.kf=sqrt(qd.kf2);

  qu.calculateDensity();
  qd.calculateDensity();

  if(iFlavor==3){
    qs.chemPot_eff =muef_s;
    qs.kf2=pow(qs.chemPot_eff, 2.)- pow(qs.mass_eff, 2.);	
    qs.kf2<=0. ? qs.kf=0. : qs.kf=sqrt(qs.kf2);
    qs.calculateDensity();
  }

}

void quarks_class::setFlavorNumber(int if_){
  iFlavor=if_;
}

void quarks_class::setTemperature(double temp_){
  temperature=temp_;
  qu.temperature=temp_;
  qd.temperature=temp_;
  if(iFlavor==3) qs.temperature=temp_;
}

void quarks_class::setEffectiveMasses(double rhob_, double temp_){
  qu.mass_eff= qu.mass + (D/cbrt(rhob_) + C*cbrt(rhob_))*sigma(temp_);
  qd.mass_eff= qd.mass + (D/cbrt(rhob_) + C*cbrt(rhob_))*sigma(temp_);
  if(iFlavor==3) qs.mass_eff= qs.mass + (D/cbrt(rhob_) + C*cbrt(rhob_))*sigma(temp_);
}

double quarks_class::getBaryonDens(){
  double rho_= qu.density +   qd.density;
   if(iFlavor==3) rho_+=qs.density;
   return rho_/3.;
}
double quarks_class::getEnergy(){
  double en0= qu.energy + qd.energy;
  if(iFlavor==3) en0+= qs.energy ;

  if(temperature>Tmin_integration){
    double edm= iFlavor==3 ? qu.getDOmega0Dmass() +qd.getDOmega0Dmass() +qs.getDOmega0Dmass()
                      : qu.getDOmega0Dmass() +qd.getDOmega0Dmass();

    en0+= edm*getDmassDtemp();
  }

  return en0;
}

double quarks_class::getPressure(){
  return -getOmega();
}

double quarks_class::getEntropy(){
    double s0=0.;
  if(temperature>Tmin_integration){
    s0= qu.entropy + qd.entropy;
    if(iFlavor==3) s0+= qs.entropy;

    double sdm=0;
    sdm= iFlavor==3 ? qu.getDOmega0Dmass() +qd.getDOmega0Dmass() +qs.getDOmega0Dmass()
                      : qu.getDOmega0Dmass() +qd.getDOmega0Dmass();

    s0-= sdm*getDmassDtemp();
  }
  return s0;
}

double quarks_class::getFenergy(){
  double fe_= qu.omega0 + qd.omega0 + qu.chemPot_eff*qu.density + qd.chemPot_eff*qd.density;
  if(iFlavor==3) fe_+= qs.omega0 + qs.chemPot_eff*qs.density;

  return fe_;
}

double quarks_class::getOmega(){
  double omega_= qu.omega0 + qd.omega0;
  if(iFlavor==3) omega_+= qs.omega0 ;

  double dOmegadMefu= qu.getDOmega0Dmass();
  double dOmegadMefd= qd.getDOmega0Dmass();
  double dOmegadMefs= qs.getDOmega0Dmass();

  omega_-= (dOmegadMefu+dOmegadMefd)*getDmassDdens()*rhoB;
  if(iFlavor==3) omega_-=dOmegadMefs*getDmassDdens()*rhoB;

  return omega_;
}

double quarks_class::sigma(double temp_){
  double s_= 1.-  8*temp_*exp(-lambda*tcrit/temp_)/(lambda*tcrit);

  if(temp_>tcrit) s_=0.;
  return temp_>Tmin_integration ? s_ : 1.;
}

double quarks_class::getDmassDdens(){
  double d0_= (-D/pow(rhoB, 4./3.) + C/pow(rhoB, 2./3.))/3.;
  if(temperature>Tmin_integration) d0_*=sigma(temperature);
  return d0_;
}

double quarks_class::getDmassDtemp(){

  double dsigma_= - 8*(temperature + lambda*tcrit)
                     *exp(-lambda*tcrit/temperature)/(lambda*tcrit*temperature);
  dsigma_*=(D/cbrt(rhoB) + C*cbrt(rhoB));
  
  if(temperature>tcrit) dsigma_=0.;
  return temperature>Tmin_integration ? dsigma_ : 0.;
}


void quarks_class::setEOS_betaEq(double rhob_, double temp_, particle &electron_, particle &muon_){
  rhoB=rhob_;
  setTemperature(temp_);
  setEffectiveMasses(rhoB, temperature);
  double mueff_u, mueff_d;
  if(firstRun){
    mueff_u=500./Mnucleon;
    mueff_d=600./Mnucleon;
  }else{
    mueff_u=qu.chemPot_eff;
    mueff_d=qd.chemPot_eff;
  }
	double x[]={mueff_u, mueff_d};

	Problem pBetaEq;
	CostFunction* costBetaEq= 
							new NumericDiffCostFunction<QuarkBetaEqFunctor,ceres::CENTRAL, 2, 2>
							(new  QuarkBetaEqFunctor(*this, electron_, muon_));

	pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

	Solver::Options optionsBetaEq;
	optionsBetaEq.parameter_tolerance = 1e-10;
	optionsBetaEq.function_tolerance = 1e-10;
	optionsBetaEq.gradient_tolerance=1e-12;
	// optionsBetaEq.line_search_direction_type= ceres::STEEPEST_DESCENT;
	optionsBetaEq.use_nonmonotonic_steps= true;
	optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
	optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
  optionsBetaEq.update_state_every_iteration = true;
	//optionsBetaEq.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
	//optionsBetaEq.linear_solver_type= ceres::DENSE_QR;

	//optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
	//optionsBetaEq.dogleg_type = ceres::TRADITIONAL_DOGLEG;

	
	optionsBetaEq.minimizer_progress_to_stdout = false;
	Solver::Summary summaryBetaEq;
	optionsBetaEq.max_num_iterations=1e4;	

	//Run
	Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

	//Print if convergence was achieved.
	 std::cout << summaryBetaEq.BriefReport() << "\n";
	 std::cout  << mueff_u << " " << mueff_d << " ---> "
              << x[0] << " " << x[1] 
	            << std::endl;

  setDensities(x[0], x[1], x[1]);
  electron_.setLepton(x[1]-x[0]);
  muon_.setLepton(x[1]-x[0]);
    
  qu.calculateQProperties();
  qd.calculateQProperties();
  if(iFlavor==3) qs.calculateQProperties();
  electron_.calculateProperties();
  muon_.calculateProperties();

  double dOmegadM_= qu.getDOmega0Dmass() +qd.getDOmega0Dmass();
  if(iFlavor==3) dOmegadM_+= qs.getDOmega0Dmass();
  qu.chemPot=    qu.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
  qd.chemPot=    qd.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
  if(iFlavor==3) qs.chemPot=    qs.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;

  // cout << muB << " " << muQ << " -> ";
  muB= qu.chemPot + 2.* qd.chemPot;
  muQ= - electron_.chemPot;
  // cout << muB << " " << muQ << endl;
  firstRun=false;
}



template <typename T>
bool QuarkBetaEqFunctor::operator()(const T* x, T* residuals) const{

	quarks.qu.chemPot_eff =x[0];
  quarks.qd.chemPot_eff =x[1];
  if(quarks.iFlavor==3) quarks.qs.chemPot_eff =x[1];
  
  quarks.qu.kf2=pow(quarks.qu.chemPot_eff, 2.)- pow(quarks.qu.mass_eff, 2.);	
  quarks.qu.kf2<=0. ? quarks.qu.kf=0. : quarks.qu.kf=sqrt(quarks.qu.kf2);
  
  quarks.qd.kf2=pow(quarks.qd.chemPot_eff, 2.)- pow(quarks.qd.mass_eff, 2.);	
  quarks.qd.kf2<=0. ? quarks.qd.kf=0. : quarks.qd.kf=sqrt(quarks.qd.kf2);
  if(quarks.iFlavor==3){
    quarks.qs.kf2=pow(quarks.qs.chemPot_eff, 2.)- pow(quarks.qs.mass_eff, 2.);	
    quarks.qs.kf2<=0. ? quarks.qs.kf=0. : quarks.qs.kf=sqrt(quarks.qs.kf2);
  }

  quarks.qu.calculateDensity();
  quarks.qd.calculateDensity();
  if(quarks.iFlavor==3) quarks.qs.calculateDensity();

  electron.setLepton(x[1]-x[0]);
  muon.setLepton(x[1]-x[0]);

  double rhob_=(quarks.qu.density + quarks.qd.density +quarks.qs.density)/3.;

	residuals[0] = quarks.rhoB 		-	rhob_;
	residuals[1] = electron.Qdens+ muon.Qdens +quarks.qu.Qdens +quarks.qd.Qdens +quarks.qs.Qdens;
  
	return true;
}



void quarks_class::setEOS_symmetric(double rhob_, double temp_){
  rhoB=rhob_;
  setTemperature(temp_);
  setEffectiveMasses(rhoB, temperature);
  double mueff_=523./Mnucleon;
  
  if(qu.mass_eff>0. && qd.mass_eff>0. && qs.mass_eff>0.){
   if(firstRun){
     mueff_=523./Mnucleon; 
   }else{
     mueff_=qu.chemPot_eff;
   }
  	double x[]={mueff_};

  	Problem pSym;
  	CostFunction* costSym= 
  							new NumericDiffCostFunction<SymmetricFunctor,ceres::CENTRAL, 1, 1>
  							(new  SymmetricFunctor(*this));

  	pSym.AddResidualBlock(costSym, NULL, x);

  	Solver::Options optionsSym;
  	optionsSym.parameter_tolerance = 1e-10;
  	optionsSym.function_tolerance = 1e-10;
  	optionsSym.gradient_tolerance=1e-12;
  	// optionsSym.line_search_direction_type= ceres::STEEPEST_DESCENT;
  	optionsSym.use_nonmonotonic_steps= true;
  	optionsSym.dense_linear_algebra_library_type=ceres::LAPACK;
  	optionsSym.linear_solver_type= ceres::DENSE_QR;
    optionsSym.update_state_every_iteration = true;
  	//optionsSym.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
  	//optionsSym.linear_solver_type= ceres::DENSE_QR;

  	//optionsSym.trust_region_strategy_type = ceres::DOGLEG;
  	//optionsSym.dogleg_type = ceres::TRADITIONAL_DOGLEG;

  
  	optionsSym.minimizer_progress_to_stdout = false;
  	Solver::Summary summarySym;
  	optionsSym.max_num_iterations=1e4;	

  	//Run
  	Solve(optionsSym, &pSym, &summarySym);

  	//Print if convergence was achieved.
  	 std::cout << summarySym.BriefReport() << "\n";
  	 std::cout  << mueff_ << " ---> " << x[0] 
  	            << std::endl;


  	mueff_	=x[0];
   setDensities(x[0], x[0], x[0]);
   //qu.chemPot= qu.chemPot_eff - 

   qu.calculateQProperties();
   qd.calculateQProperties();
   qs.calculateQProperties();
   double dOmegadM_= qu.getDOmega0Dmass() +qd.getDOmega0Dmass() + qs.getDOmega0Dmass();
   qu.chemPot=    qu.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
   qd.chemPot=    qd.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
   qs.chemPot=    qs.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
   muB= qu.chemPot + 2.* qd.chemPot;
   muQ= -qu.chemPot + qd.chemPot;
  firstRun=false;
  }


} 

template <typename T>
bool SymmetricFunctor::operator()(const T* x, T* residuals) const{

	quarks.qu.chemPot_eff =x[0];
  quarks.qd.chemPot_eff =x[0];
  if(quarks.iFlavor==3) quarks.qs.chemPot_eff =x[0];
  
  quarks.qu.kf2=pow(quarks.qu.chemPot_eff, 2.)- pow(quarks.qu.mass_eff, 2.);	
  quarks.qu.kf2<=0. ? quarks.qu.kf=0. : quarks.qu.kf=sqrt(quarks.qu.kf2);
  
  quarks.qd.kf2=pow(quarks.qd.chemPot_eff, 2.)- pow(quarks.qd.mass_eff, 2.);	
  quarks.qd.kf2<=0. ? quarks.qd.kf=0. : quarks.qd.kf=sqrt(quarks.qd.kf2);
  
  if(quarks.iFlavor==3){
    quarks.qs.kf2=pow(quarks.qs.chemPot_eff, 2.)- pow(quarks.qs.mass_eff, 2.);	
    quarks.qs.kf2<=0. ? quarks.qs.kf=0. : quarks.qs.kf=sqrt(quarks.qs.kf2);
  }

  quarks.qu.calculateDensity();
  quarks.qd.calculateDensity();
  if(quarks.iFlavor==3) quarks.qs.calculateDensity();

  double rhob_=(quarks.qu.density + quarks.qd.density +quarks.qs.density)/3.;

	residuals[0] = quarks.rhoB 		-	rhob_;
  
	return true;
}


	void 	quarks_class::setEoSFlavorFixed(double rhob_, double temp_, 
																	double Yu_,  double Yd_,  double Ys_){

  rhoB= rhob_;
  setTemperature(temp_);
  setEffectiveMasses(rhoB, temperature);
  
  Yu=Yu_; Yd=Yd_; Ys=Ys_;
  qu.density= 3.*Yu_*rhoB;
  qd.density= 3.*Yd_*rhoB;
  qs.density= 3.*Ys_*rhoB;
  qu.kf= cbrt(6.*pi2*qu.density/qu.gamma);
  qu.kf2=qu.kf*qu.kf;
  qd.kf= cbrt(6.*pi2*qd.density/qd.gamma);
  qd.kf2=qd.kf*qd.kf;
  if(iFlavor==3){
    qs.kf= cbrt(6.*pi2*qs.density/qs.gamma);
    qs.kf2=qs.kf*qs.kf;
  }
  qu.solveChemPotEff();
  qd.solveChemPotEff();
  if(iFlavor==3) qs.solveChemPotEff();

  qu.calculateQProperties();
  qd.calculateQProperties();
  if(iFlavor==3) qs.calculateQProperties();

  double dOmegadM_= qu.getDOmega0Dmass() +qd.getDOmega0Dmass() + qs.getDOmega0Dmass();
  qu.chemPot=    qu.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
  qd.chemPot=    qd.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
  if(iFlavor==3) qs.chemPot=    qs.chemPot_eff +    dOmegadM_*getDmassDdens()/3.;
//  muB= (qu.energy + qd.energy + qs.energy  -temperature*(qu.entropy + qd.entropy + qs.entropy )
//            +qu.pressure + qd.pressure + qs.pressure )/rhoB;
}

void quarks_class::setEoSFlavor_PressFixed(double press_, double temp_, 
                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_){

  PressureTot=press_;
  setTemperature(temp_);
  Yu=Yu_; Yd=Yd_; Ys=Ys_;
    double rhob_;
    if(firstRun){
      rhob_=0.5*pow(hc/Mnucleon, 3);
    }else{
      rhob_=rhoB;
    }
	  double x[]={rhob_};

	  Problem p;
	  CostFunction* cost= 
	  						new NumericDiffCostFunction<QuarkFlavor_PressFixed,ceres::CENTRAL, 1, 1>
	  						(new QuarkFlavor_PressFixed(*this, electron_, muon_));

	  p.AddResidualBlock(cost, NULL, x);
    p.SetParameterLowerBound(x, 0, 0.);
    // p.SetParameterUpperBound(x, 0, 1.5*pow(hc/Mnucleon, 3));
	  Solver::Options options;
	  options.parameter_tolerance = 1e-8;
	  options.function_tolerance = 1e-10;
	  options.gradient_tolerance=1e-10;
	  options.max_num_iterations=1e4;	

    // if(PressureTot*Mnucleon*pow(Mnucleon/hc, 3.) < 50. ){
		//   // optionsBetaEq.parameter_tolerance = 1e-22;
		// 	// optionsBetaEq.function_tolerance = 1e-22;
		// 	// optionsBetaEq.gradient_tolerance=1e-25;
		// 	options.parameter_tolerance = 1e-15;
		// 	options.function_tolerance = 1e-15;
		// 	options.gradient_tolerance=1e-16;
		// 	options.max_num_iterations=1e6;	

		// }	
	  // options.line_search_direction_type= ceres::STEEPEST_DESCENT;
	  options.use_nonmonotonic_steps= true;
	  options.dense_linear_algebra_library_type=ceres::LAPACK;
	  options.linear_solver_type= ceres::DENSE_QR;
    options.update_state_every_iteration = true;
	  //options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
	  //options.linear_solver_type= ceres::DENSE_QR;

	  //options.trust_region_strategy_type = ceres::DOGLEG;
	  //options.dogleg_type = ceres::TRADITIONAL_DOGLEG;

  
	  options.minimizer_progress_to_stdout = false;
	  Solver::Summary summary;

	  //Run
	  Solve(options, &p, &summary);

	  //Print if convergence was achieved.
	   std::cout << summary.BriefReport() << "\n";
	   std::cout  << rhob_ << " ---> " << x[0] 
	              << std::endl;


	  rhoB=x[0];
    setEoSFlavorFixed(rhoB, temperature, Yu_, Yd_, Ys_);
     
    muB= (getEnergy() +electron_.energy + muon_.energy
      -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
      + getPressure()+electron_.pressure + muon_.pressure)/rhoB;

    firstRun=false;
}



template <typename T>
bool QuarkFlavor_PressFixed::operator()(const T* x, T* residuals) const{

  quarks.setEoSFlavorFixed(x[0], quarks.temperature, quarks.Yu, quarks.Yd, quarks.Ys);
  
  residuals[0] = quarks.getPressure() + electron.pressure+ muon.pressure - quarks.PressureTot;
  return true;
}


void quarks_class::setEoSFlavor_muBFixed(double mub_, double temp_, 
                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_){

  muB=mub_;
  setTemperature(temp_);
  Yu=Yu_; Yd=Yd_; Ys=Ys_;
  double rhob_;
  if(firstRun){
    rhob_=0.5*pow(hc/Mnucleon, 3);
  }else{
    rhob_=rhoB;
  }
	double x[]={rhob_};
	
  Problem p;
	CostFunction* cost= 
						new NumericDiffCostFunction<QuarkFlavor_muBFixed,ceres::CENTRAL, 1, 1>
						(new QuarkFlavor_muBFixed(*this, electron_, muon_));

	p.AddResidualBlock(cost, NULL, x);
  p.SetParameterLowerBound(x, 0, 0.);
  // p.SetParameterUpperBound(x, 0, 1.5*pow(hc/Mnucleon, 3));
	Solver::Options options;
	options.parameter_tolerance = 1e-8;
	options.function_tolerance = 1e-10;
	options.gradient_tolerance=1e-10;
	options.max_num_iterations=1e5;	

   // if(PressureTot*Mnucleon*pow(Mnucleon/hc, 3.) < 50. ){
	//   // optionsBetaEq.parameter_tolerance = 1e-22;
	// 	// optionsBetaEq.function_tolerance = 1e-22;
	// 	// optionsBetaEq.gradient_tolerance=1e-25;
	// 	options.parameter_tolerance = 1e-15;
	// 	options.function_tolerance = 1e-15;
	// 	options.gradient_tolerance=1e-16;
	// 	options.max_num_iterations=1e6;	
	// }	
	 // options.line_search_direction_type= ceres::STEEPEST_DESCENT;
	 options.use_nonmonotonic_steps= true;
	 options.dense_linear_algebra_library_type=ceres::LAPACK;
	 options.linear_solver_type= ceres::DENSE_QR;
   options.update_state_every_iteration = true;
	 //options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
	 //options.linear_solver_type= ceres::DENSE_QR;
	 //options.trust_region_strategy_type = ceres::DOGLEG;
	 //options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
  
	options.minimizer_progress_to_stdout = false;
	Solver::Summary summary;

	//Run
	Solve(options, &p, &summary);

	//Print if convergence was achieved.
	 std::cout << summary.BriefReport() << "\n";
	 std::cout  << rhob_ << " ---> " << x[0] 
	            << std::endl;
  rhoB=x[0];
  setEoSFlavorFixed(rhoB, temperature, Yu_, Yd_, Ys_);
   
  muB= (getEnergy() +electron_.energy + muon_.energy
    -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
    + getPressure()+electron_.pressure + muon_.pressure)/rhoB;
  
  firstRun=false;
}



template <typename T>
bool QuarkFlavor_muBFixed::operator()(const T* x, T* residuals) const{

  quarks.setEoSFlavorFixed(x[0], quarks.temperature, quarks.Yu, quarks.Yd, quarks.Ys);

  double mub_ = (quarks.getEnergy() +electron.energy + muon.energy
      -quarks.temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + quarks.getPressure()+electron.pressure + muon.pressure)/quarks.getBaryonDens();
  residuals[0]=  quarks.muB - mub_;

  return true;
}

void quarks_class::setEoSFlavor_muBFixed2(double mub_, double temp_, 
                                            particle &electron_, particle &muon_,
                                            double Yu_,  double Yd_,  double Ys_,
                                            double Ye_, double Ym_){
  muB=mub_;
  setTemperature(temp_);
  Yu=Yu_; Yd=Yd_; Ys=Ys_;
  Ye=Ye_; Ym=Ym_;
  double rhob_;
  if(firstRun){
    rhob_=0.5*pow(hc/Mnucleon, 3);
  }else{
    rhob_=rhoB;
  }
  double x[]={rhob_};

	Problem p;
	CostFunction* cost= 
							new NumericDiffCostFunction<QuarkFlavor_muBFixed2,ceres::CENTRAL, 1, 1>
							(new QuarkFlavor_muBFixed2(*this, electron_, muon_));
  
  p.AddResidualBlock(cost, NULL, x);
  p.SetParameterLowerBound(x, 0, 0.);
  // p.SetParameterUpperBound(x, 0, 1.5*pow(hc/Mnucleon, 3));
	Solver::Options options;
	options.parameter_tolerance = 1e-8;
	options.function_tolerance = 1e-10;
	options.gradient_tolerance=1e-10;
	options.max_num_iterations=1e5;	

   // if(PressureTot*Mnucleon*pow(Mnucleon/hc, 3.) < 50. ){
	//   // optionsBetaEq.parameter_tolerance = 1e-22;
	// 	// optionsBetaEq.function_tolerance = 1e-22;
	// 	// optionsBetaEq.gradient_tolerance=1e-25;
	// 	options.parameter_tolerance = 1e-15;
	// 	options.function_tolerance = 1e-15;
	// 	options.gradient_tolerance=1e-16;
	// 	options.max_num_iterations=1e6;	
	// }	
	// options.line_search_direction_type= ceres::STEEPEST_DESCENT;
	options.use_nonmonotonic_steps= true;
	options.dense_linear_algebra_library_type=ceres::LAPACK;
	options.linear_solver_type= ceres::DENSE_QR;
  options.update_state_every_iteration = true;
	//options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
	//options.linear_solver_type= ceres::DENSE_QR;

	//options.trust_region_strategy_type = ceres::DOGLEG;
	//options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
  options.minimizer_progress_to_stdout = false;
	Solver::Summary summary;

//Run
  Solve(options, &p, &summary);

//Print if convergence was achieved.
  std::cout << summary.BriefReport() << "\n";
  std::cout  << rhob_ << " ---> " << x[0] 
              << std::endl;
  rhoB=x[0];
  setEoSFlavorFixed(rhoB, temperature, Yu_, Yd_, Ys_);

  electron_.density= Ye*x[0];
  electron_.kf= cbrt(6.*pi2*electron_.density/electron_.gamma);
  electron_.kf2=electron_.kf*electron_.kf;
  electron_.solveChemPotEff();
  electron_.chemPot= electron_.chemPot_eff;
  electron_.calculateProperties();

  muon_.density= Ym*x[0];
  muon_.kf= cbrt(6.*pi2*muon_.density/muon_.gamma);
  muon_.kf2=muon_.kf*muon_.kf;
  muon_.solveChemPotEff();
  muon_.chemPot= muon_.chemPot_eff;
  muon_.calculateProperties();
  
  muB= (getEnergy() +electron_.energy + muon_.energy
    -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
    + getPressure()+electron_.pressure + muon_.pressure)/rhoB;
  
  firstRun=false;
}



template <typename T>
bool QuarkFlavor_muBFixed2::operator()(const T* x, T* residuals) const{

  quarks.setEoSFlavorFixed(x[0], quarks.temperature, quarks.Yu, quarks.Yd, quarks.Ys);

  electron.density= quarks.Ye*x[0];
  electron.kf= cbrt(6.*pi2*electron.density/electron.gamma);
  electron.kf2=electron.kf*electron.kf;
  electron.solveChemPotEff();
  electron.chemPot= electron.chemPot_eff;
  electron.calculateProperties();

  muon.density= quarks.Ym*x[0];
  muon.kf= cbrt(6.*pi2*muon.density/muon.gamma);
  muon.kf2=muon.kf*muon.kf;
  muon.solveChemPotEff();
  muon.chemPot= muon.chemPot_eff;
  muon.calculateProperties();

  double mub_ = (quarks.getEnergy() +electron.energy + muon.energy
      -quarks.temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + quarks.getPressure()+electron.pressure + muon.pressure)/quarks.getBaryonDens();
 
  residuals[0]=  quarks.muB - mub_;

  return true;
}