#include "bag_model.h"

bag_model_class::bag_model_class(){
	qu.mass= 4./Mnucleon; 	qu.spin	= 1./2.;	qu.gamma = 6.;    qu.Q=  2./3.; qu.type="Q";   
	qd.mass= 4./Mnucleon; 	qd.spin	= 1./2.;	qd.gamma = 6.;    qd.Q= -1./3.; qd.type="Q"; 
	qs.mass= 95./Mnucleon;	qs.spin	= 1./2.;	qs.gamma = 6.;    qs.Q= -1./3.; qs.type="Q";     
	
	qu.mass_eff= qu.mass;
	qd.mass_eff= qd.mass;
	qs.mass_eff= qs.mass;
  
  Mv=780.0/Mnucleon;
	gv=0.;
	xsi=0;//6.*0.0256;
	xvu= 1.;
	xvd= 1.;
	xvs= 1.;//1.
	Bag= 148./Mnucleon;
  tcrit=0./Mnucleon;

}

//=============== Set RMF parameters: nucleons/meson couplings and masses  ===============
void bag_model_class::setParametrization(std::string parametrization_){
		parametrization=parametrization_;

}

void bag_model_class::setParameters(double bag_, double Gv_, double xsi_, double Xv_, double tcrit_){
  Bag= bag_;
  gv=sqrt(Gv_)*Mv;
  xsi=xsi_;
  xvs=Xv_;
  tcrit=tcrit_;
}


//=============== Calculate thermodynamic properties for all quarks: ===============
void bag_model_class::setThermodynamics(){
	qu.calculateProperties();
	qd.calculateProperties();
	if(iFlavor==3) qs.calculateProperties();
}


void bag_model_class::setDensities(double muef_u, double muef_d, double muef_s){
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

void bag_model_class::setTemperature(double temp_){
  temperature=temp_;
  qu.temperature=temp_;
  qd.temperature=temp_;
  if(iFlavor==3) qs.temperature=temp_;
  if(tcrit>0) Bag*=(1. + pow(temperature/tcrit, 4.));
}

//=============== Set the magnetic field for all quarks: ===============
void bag_model_class::setBfield(bool dob_, double B_){
	Bfield=B_;
  qu.setBfield(dob_, B_);
  qd.setBfield(dob_, B_);
  if(iFlavor==3) qs.setBfield(dob_, B_);
}


//=============== Set the magnetic moment for all quarks: ===============          
void bag_model_class::setAMM(bool doa_){							//MUB			Kb
  qu.setAMM(doa_, 0.);
  qd.setAMM(doa_, 0.);
  if(iFlavor==3) qs.setAMM(doa_, 0.);
}

void bag_model_class::setFlavorNumber(int if_){
  iFlavor=if_;
}

//=============== Set nucleon EoS with fixed proton fraction and temperature  ===============
void bag_model_class::setEOS_symmetric(double rhob_, double temp_){
  rhoB=rhob_;
  setTemperature(temp_);
  double muup_, v0_;
  
  if(firstRun){
    muup_=1500./Mnucleon; 
    v0_= 0.12;
  }else{
    muup_=qu.chemPot;
    v0_=V0;
  }
 	double x[]={muup_, v0_};

  Problem pSym;
  CostFunction* costSym= 
  						new NumericDiffCostFunction<SymmetricFunctor,ceres::CENTRAL, 2, 2>
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
   std::cout  << muup_ << " ---> " << x[0] 
              << std::endl;
  muup_	=x[0];
	V0= x[1];
	qu.chemPot=x[0];
	qd.chemPot=x[0];
	if(iFlavor==3) qs.chemPot=x[0];
  setDensities(x[0] - gv*xvu*x[1], 
               x[0] - gv*xvd*x[1], 
               x[0] - gv*xvs*x[1]);

  setThermodynamics();
  muB= qu.chemPot + 2.* qd.chemPot;
  muQ= -qu.chemPot + qd.chemPot;
  firstRun=false;

}

template <typename T>
bool SymmetricFunctor::operator()(const T* x, T* residuals) const{

  quarks.V0=x[1];
	quarks.setDensities(x[0]- quarks.gv*quarks.xvu*x[1],
                      x[0]- quarks.gv*quarks.xvd*x[1],
                      x[0]- quarks.gv*quarks.xvs*x[1]);

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
  residuals[1] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());
	return true;
}


void bag_model_class::setEOS_betaEq(double rhob_, double temp_, 
                                    particle & electron_, particle &muon_){
  rhoB=rhob_;
  setTemperature(temp_);
  double muup_, mue_, v0_;
  
  if(firstRun){
    muup_=1500./Mnucleon; 
    mue_= 0.12;
    v0_= 0.12;
  }else{
    muup_=qu.chemPot;
    mue_=electron_.chemPot;
    v0_=V0;
  }
 	double x[]={muup_, mue_, v0_};

  Problem problem;
  CostFunction* cost= 
  						new NumericDiffCostFunction<QuarkBetaEqFunctor,ceres::CENTRAL, 3, 3>
  						(new  QuarkBetaEqFunctor(*this, electron_, muon_));

  problem.AddResidualBlock(cost, NULL, x);
  Solver::Options options;
  options.parameter_tolerance = 1e-10;
  options.function_tolerance = 1e-10;
  options.gradient_tolerance=1e-12;
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
  options.max_num_iterations=1e4;	
  //Run
  Solve(options, &problem, &summary);
  //Print if convergence was achieved.
   std::cout << summary.BriefReport() << "\n";
   std::cout  << muup_ << " " << mue_ << " " << v0_ << " ---> " 
              << x[0] << " " << x[1] << " " << x[2] 
              << std::endl;
	qu.chemPot =x[0];
  qd.chemPot =x[0]+x[1];
  qs.chemPot =x[0]+x[1];
  electron_.setLepton(x[1]);
  muon_.setLepton(x[1]);
  V0=x[2];
  setDensities(qu.chemPot- gv*xvu*V0,
               qd.chemPot- gv*xvd*V0,
               qs.chemPot- gv*xvs*V0);
  setThermodynamics();
  // muB= qu.chemPot + 2.* qd.chemPot;
  muB= (getEnergy() - temperature*getEntropy() + getPressure()
        +electron_.energy - temperature*electron_.entropy + electron_.pressure
        +muon_.energy - temperature*muon_.entropy + muon_.pressure)/rhoB;
  // muQ= -qu.chemPot + qd.chemPot;
  firstRun=false;

} 

template <typename T>
bool QuarkBetaEqFunctor::operator()(const T* x, T* residuals) const{


	quarks.qu.chemPot =x[0];
  quarks.qd.chemPot =x[0]+x[1];
  quarks.qs.chemPot =x[0]+x[1];
  electron.setLepton(x[1]);
  muon.setLepton(x[1]);
  quarks.V0=x[2];
  quarks.setDensities(quarks.qu.chemPot- quarks.gv*quarks.xvu*x[2],
                      quarks.qd.chemPot- quarks.gv*quarks.xvd*x[2],
                      quarks.qs.chemPot- quarks.gv*quarks.xvs*x[2]);


  double rhob_=(quarks.qu.density + quarks.qd.density +quarks.qs.density)/3.;

	residuals[0] = quarks.rhoB 		-	rhob_;
	residuals[1] = electron.Qdens+ muon.Qdens +quarks.getChargeDens();
  residuals[2] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());
	return true;
}



//=============== Solve for the vector meson fields with Ceres library: ===============
void bag_model_class::setVectorMeanFields(){
	
	double v0_= getOmegaEffDens()*gv/pow(Mv , 2) ;       
	if(xsi!=0){
		Problem pV;
    CostFunction* costV =	new NumericDiffCostFunction<VQFunctor, ceres::CENTRAL, 1, 1>
														(new VQFunctor(*this) );


    pV.AddResidualBlock(costV, NULL, &v0_);

   Solver::Options optionsV;
    // optionsV.parameter_tolerance = 1e-10;
    // optionsV.function_tolerance = 1e-8;
    // optionsV.gradient_tolerance=1e-12;
    //optionsV.trust_region_strategy_type = ceres::DOGLEG;
    //optionsV.dogleg_type = ceres::TRADITIONAL_DOGLEG;
		optionsV.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsV.linear_solver_type= ceres::DENSE_QR;

    optionsV.minimizer_progress_to_stdout = false;

    Solver::Summary summaryV;
    Solve(optionsV, &pV, &summaryV);
 //	std::cout << "v: " << summaryV.BriefReport() << "\n";
	
	}
	V0=v0_;
}
//=============== Functor Vector meson solver for Ceres: ===============
template <typename T>
bool VQFunctor::operator()(const T* x, T* residuals) const{

	quarks.V0=x[0];
  residuals[0] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());
	return true;
}

//=============== Calculate total baryon density: ===============
double bag_model_class::getBaryonDens(){
	double densb_=	qu.density + qd.density;
	if(iFlavor==3) densb_+=qs.density;
	return densb_/3.;
}

//=============== Calculate total isospin density: ===============

//=============== Calculate total charge density: ===============
double bag_model_class::getChargeDens(){
	double charge_=	qu.Qdens + qd.Qdens;
	if(iFlavor==3) charge_+=qs.Qdens;

	return charge_;
}

//=============== Calculate baryon density weighted by x_{vb} ratios: ===============
double bag_model_class::getOmegaEffDens(){
	double densv_=	xvu*qu.density + xvd*qd.density;
	if(iFlavor==3) densv_+= xvs*qs.density;
	return densv_;
}



//=============== Calculate residue of isoscalar-vector meson eom: ===============
double bag_model_class::omegaMeson_eom_residue(double rhob_){
  return gv*rhob_ - V0*pow(Mv, 2.) 
								-xsi*pow(gv, 4.)*pow(V0, 3.)/6.;
}


//=============== Calculate total baryonic energy: ===============
double bag_model_class::getEnergy(void){
	double enerMeson= pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/8.;
  double en0= qu.energy + qd.energy;
  if(iFlavor==3) en0+= qs.energy ;

  return en0+enerMeson + Bag;
}

double bag_model_class::getPressure(){
	// double press_meson=pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/24.;
  // double press_= qu.pressure + qd.pressure;
  // if(iFlavor==3) press_+= qs.pressure ;
  // press_+= press_meson - Bag;
  double press_=qu.chemPot*qu.density + qd.chemPot*qd.density
                - getEnergy() ;
  if(iFlavor==3) press_+= qs.chemPot*qs.density;
  if(temperature>Tmin_integration) press_+=temperature*getEntropy();
  return press_ ;
}

double bag_model_class::getEntropy(){
    double s0=0.;

  if(temperature>Tmin_integration){
    s0= qu.entropy + qd.entropy;
    if(iFlavor==3) s0+= qs.entropy;
	}
  return s0;
}

double bag_model_class::getFenergy(void){
  return getEnergy()- temperature*getEntropy();
}

	void 	bag_model_class::setEoSFlavorFixed(double rhob_, double temp_, double v0_, 
																	double Yu_,  double Yd_,  double Ys_){

  // iFlavor=3; 
  rhoB= rhob_;
  setTemperature(temp_);
  V0=v0_;
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
  
  qu.chemPot=    qu.chemPot_eff +   gv*xvu*V0;
  qd.chemPot=    qd.chemPot_eff +   gv*xvd*V0;
  if(iFlavor==3){
    qs.chemPot=    qs.chemPot_eff +   gv*xvs*V0;
  }

  setThermodynamics();
}

void bag_model_class::setEoSFlavor_muBFixed(double mub_, double temp_, 
                                            particle electron_, particle muon_,
                                            double Yu_,  double Yd_,  double Ys_){

  muB=mub_;
  setTemperature(temp_);
  Yu=Yu_; Yd=Yd_; Ys=Ys_;

  if(Bfield==0){

    double rhob_, v0_;

    if(firstRun){
      rhob_=1.5*pow(hc/Mnucleon, 3);
      v0_=0.03;
    }else{
      rhob_=rhoB;
      v0_=V0;
    }

	  double x[]={rhob_, v0_};
	  Problem p;
	  CostFunction* cost= 
	    						new NumericDiffCostFunction<QuarkFlavor_muBFixed,ceres::CENTRAL, 2, 2>
	    						(new QuarkFlavor_muBFixed(*this, electron_, muon_));

	    p.AddResidualBlock(cost, NULL, x);
      p.SetParameterLowerBound(x, 0, 0.);
      // p.SetParameterLowerBound(x, 1, 0.);
      // p.SetParameterUpperBound(x, 0, 2.*pow(hc/Mnucleon, 3));
	    Solver::Options options;
	    options.parameter_tolerance = 1e-8;
	    options.function_tolerance = 1e-10;
	    options.gradient_tolerance=1e-12;
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
	     std::cout  <<"quarks: rhob, v0= " <<  rhob_*pow(Mnucleon/hc, 3) << " " << v0_ 
                    << "---> " << x[0]*pow(Mnucleon/hc, 3) << " " << x[1] 
	                << std::endl;
      

	    rhoB=x[0];
      V0  =x[1];
      setEoSFlavorFixed(x[0], temperature, x[1], Yu_, Yd_, Ys_);
      cout << "mus: " << qu.chemPot << " " << qd.chemPot <<  " " << qs.chemPot << endl;
      muB = (getEnergy() +electron_.energy + muon_.energy
       -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
       + getPressure()+electron_.pressure + muon_.pressure)/getBaryonDens();
  }else{

    double rhob_, v0_, muu_, mud_, mus_;

    if(firstRun){
      rhob_= 2.0202*pow(hc/Mnucleon, 3);
      v0_=0.0368514;
      muu_=0.71131;
      mud_=0.78430;
      mus_=0.49338;
    }else{
      rhob_=rhoB;
      v0_=V0;
      muu_=qu.chemPot;
      mud_=qd.chemPot;
      mus_=qs.chemPot;
    }

	  double x[]={rhob_, v0_, muu_, mud_, mus_};
	  Problem p;
	  CostFunction* cost= 
	    						new NumericDiffCostFunction<QuarkFlavor_muBFixed_B,ceres::CENTRAL, 5, 5>
	    						(new QuarkFlavor_muBFixed_B(*this, electron_, muon_));

	    p.AddResidualBlock(cost, NULL, x);
      p.SetParameterLowerBound(x, 0, 0.);
      // p.SetParameterLowerBound(x, 2, 0.);
      // p.SetParameterLowerBound(x, 3, 0.);
      // p.SetParameterLowerBound(x, 4, 0.);
    
      // p.SetParameterLowerBound(x, 1, 0.);
      // p.SetParameterUpperBound(x, 0, 2.*pow(hc/Mnucleon, 3));
	    Solver::Options options;
	    options.parameter_tolerance = 1e-8;
	    options.function_tolerance = 1e-10;
	    options.gradient_tolerance=1e-12;
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
	     std::cout  <<"quarks: rhob, v0, muu, mud, mus= "  <<  rhob_*pow(Mnucleon/hc, 3) << " " 
          << v0_ << " " <<    muu_ << " "  << mud_ << " " << mus_
                    << "---> " << x[0]*pow(Mnucleon/hc, 3) << " " << x[1] << " " 
                      << x[2] << x[3] << " " << x[4] 
	                << std::endl;      

	  rhoB=x[0];
    V0  =x[1];
  
    qu.chemPot=x[2];
    qd.chemPot=x[3];
    qu.setQuarkEff(x[2] -xvu*gv*x[1]);
    qd.setQuarkEff(x[3] -xvd*gv*x[1]);
    if(iFlavor==3){
      qs.chemPot=x[4];
      qs.setQuarkEff(x[4] - xvs*gv*x[1]); 
    }
  
    setThermodynamics();

    muB = (getEnergy() +electron_.energy + muon_.energy
      -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
      + getPressure()+electron_.pressure + muon_.pressure)/getBaryonDens();
  }

    firstRun=false;
}



template <typename T>
bool QuarkFlavor_muBFixed::operator()(const T* x, T* residuals) const{

  quarks.setEoSFlavorFixed(x[0], quarks.temperature, x[1], quarks.Yu, quarks.Yd, quarks.Ys);

  double mub_ = (quarks.getEnergy() +electron.energy + muon.energy
      -quarks.temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + quarks.getPressure()+electron.pressure + muon.pressure)/quarks.getBaryonDens();
  // double mub_=quarks.Yu* quarks.qu.chemPot + quarks.Yd*quarks.qd.chemPot
  //                 +  quarks.Ys*quarks.qs.chemPot
	// 								+ (electron.chemPot*electron.density +muon.chemPot*muon.density)/x[0];
  residuals[0]=  quarks.muB - mub_;
  residuals[1] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());

  // residuals[0] = quarks.getPressure() + electron.pressure+ muon.pressure - quarks.PressureTot;
  return true;
}

template <typename T>
bool QuarkFlavor_muBFixed_B::operator()(const T* x, T* residuals) const{
  quarks.V0 = x[1];
  quarks.qu.chemPot=x[2];
  quarks.qd.chemPot=x[3];

  quarks.qu.setQuarkEff(x[2] - quarks.xvu*quarks.gv*x[1]);
  quarks.qd.setQuarkEff(x[3] - quarks.xvd*quarks.gv*x[1]);

  if(quarks.iFlavor==3){
    quarks.qs.chemPot=x[4];
    quarks.qs.setQuarkEff(x[4] - quarks.xvs*quarks.gv*x[1]); 
  }
  
  quarks.setThermodynamics();

  double mub_ = (quarks.getEnergy() +electron.energy + muon.energy
      -quarks.temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + quarks.getPressure()+electron.pressure + muon.pressure)/x[0];
  // double mub_=quarks.Yu* quarks.qu.chemPot + quarks.Yd*quarks.qd.chemPot
  //                 +  quarks.Ys*quarks.qs.chemPot
	// 								+ (electron.chemPot*electron.density +muon.chemPot*muon.density)/x[0];

  residuals[0]=  quarks.muB - mub_;
  residuals[1] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());
  residuals[2]= quarks.qu.density - 3.*quarks.Yu*x[0];
  residuals[3]= quarks.qd.density - 3.*quarks.Yd*x[0];
  residuals[4]= quarks.qs.density - 3.*quarks.Ys*x[0];
  
  return true;
}

void bag_model_class::setEoSFlavor_muBFixed2(double mub_, double temp_, 
                                            particle & electron_, particle & muon_,
                                            double Yu_,  double Yd_,  double Ys_, 
                                            double Ye_, double Ym_){

  muB=mub_;
  setTemperature(temp_);
  Yu=Yu_; Yd=Yd_; Ys=Ys_;
  Ye=Ye_;  Ym=Ym_;

  if(Bfield==0){

    double rhob_, v0_;

    if(firstRun){
      rhob_=2.*pow(hc/Mnucleon, 3);
      v0_=0.03;
    }else{
      rhob_=rhoB;
      v0_=V0;
    }

	  double x[]={rhob_, v0_};
	  Problem p;
	  CostFunction* cost= 
	    						new NumericDiffCostFunction<QuarkFlavor_muBFixed2,ceres::CENTRAL, 2, 2>
	    						(new QuarkFlavor_muBFixed2(*this, electron_, muon_));

	   p.AddResidualBlock(cost, NULL, x);
     p.SetParameterLowerBound(x, 0, 0.);
     // p.SetParameterLowerBound(x, 1, 0.);
     // p.SetParameterUpperBound(x, 0, 2.*pow(hc/Mnucleon, 3));
	   Solver::Options options;
	   options.parameter_tolerance = 1e-8;
	   options.function_tolerance = 1e-10;
	   options.gradient_tolerance=1e-12;
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
	  std::cout  <<"quarks: rhob, v0= " <<  rhob_*pow(Mnucleon/hc, 3) << " " << v0_ 
                 << "---> " << x[0]*pow(Mnucleon/hc, 3) << " " << x[1] 
	             << std::endl;
	  rhoB=x[0];
    V0  =x[1];
    setEoSFlavorFixed(x[0], temperature, x[1], Yu_, Yd_, Ys_);

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

   cout << "mus: " << qu.chemPot << " " << qd.chemPot <<  " " << qs.chemPot << " " 
        << electron_.chemPot << " " << muon_.chemPot << endl;

    muB = (getEnergy() +electron_.energy + muon_.energy
        -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
        + getPressure()+electron_.pressure + muon_.pressure)/getBaryonDens();
  }else{
    
    double rhob_, v0_, muu_, mud_, mus_, mue_, mum_;

    if(firstRun){
      rhob_= 2.0202*pow(hc/Mnucleon, 3);
      v0_=0.0368514;
      muu_=0.71131;
      mud_=0.78430;
      mus_=0.49338;
      mue_= 0.12;
      mum_= 0.12;
    }else{
      rhob_=rhoB;
      v0_=V0;
      muu_=qu.chemPot;
      mud_=qd.chemPot;
      mus_=qs.chemPot;
      mue_= electron_.chemPot;
      mum_= muon_.chemPot;
    }

	  double x[]={rhob_, v0_, muu_, mud_, mus_, mue_, mum_};
	  Problem p;
	  CostFunction* cost= 
	    						new NumericDiffCostFunction<QuarkFlavor_muBFixed2_B,ceres::CENTRAL, 7, 7>
	    						(new QuarkFlavor_muBFixed2_B(*this, electron_, muon_));

	    p.AddResidualBlock(cost, NULL, x);
      p.SetParameterLowerBound(x, 0, 0.);
      // p.SetParameterLowerBound(x, 2, 0.);
      // p.SetParameterLowerBound(x, 3, 0.);
      // p.SetParameterLowerBound(x, 4, 0.);
    
      // p.SetParameterLowerBound(x, 1, 0.);
      // p.SetParameterUpperBound(x, 0, 2.*pow(hc/Mnucleon, 3));
	    Solver::Options options;
	    options.parameter_tolerance = 1e-8;
	    options.function_tolerance = 1e-10;
	    options.gradient_tolerance=1e-12;
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
	     std::cout  <<"quarks: rhob, v0, muu, mud, mus= "  <<  rhob_*pow(Mnucleon/hc, 3) << " " 
        << v0_ << " " <<    muu_ << " "  << mud_ << " " << mus_ << " " << mue_ << " " << mum_
                    << "---> " << x[0]*pow(Mnucleon/hc, 3) << " " << x[1] << " " 
                      << x[2] << x[3] << " " << x[4] << " " << x[5] << " " << x[6]
	                << std::endl;      

	  rhoB=x[0];
    V0  =x[1];
  
    qu.chemPot=x[2];
    qd.chemPot=x[3];
    qu.setQuarkEff(x[2] -xvu*gv*x[1]);
    qd.setQuarkEff(x[3] -xvd*gv*x[1]);
    if(iFlavor==3){
      qs.chemPot=x[4];
      qs.setQuarkEff(x[4] - xvs*gv*x[1]); 
    }
  
    setThermodynamics();

    electron_.setLepton(x[5]);
    electron_.calculateProperties();
    muon_.setLepton(x[6]);
    muon_.calculateProperties();

    muB = (getEnergy() +electron_.energy + muon_.energy
      -temperature*(getEntropy()+electron_.entropy + muon_.entropy)
      + getPressure()+electron_.pressure + muon_.pressure)/getBaryonDens();

  }

  firstRun=false;
}

template <typename T>
bool QuarkFlavor_muBFixed2::operator()(const T* x, T* residuals) const{

  quarks.setEoSFlavorFixed(x[0], quarks.temperature, x[1], quarks.Yu, quarks.Yd, quarks.Ys);

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
  // double mub_=quarks.Yu* quarks.qu.chemPot + quarks.Yd*quarks.qd.chemPot
  //                 +  quarks.Ys*quarks.qs.chemPot
	// 								+ (electron.chemPot*electron.density +muon.chemPot*muon.density)/x[0];
  residuals[0]=  quarks.muB - mub_;
  residuals[1] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());

  // residuals[0] = quarks.getPressure() + electron.pressure+ muon.pressure - quarks.PressureTot;
  return true;
}


template <typename T>
bool QuarkFlavor_muBFixed2_B::operator()(const T* x, T* residuals) const{

  quarks.V0 = x[1];
  quarks.qu.chemPot=x[2];
  quarks.qd.chemPot=x[3];

  quarks.qu.setQuarkEff(x[2] - quarks.xvu*quarks.gv*x[1]);
  quarks.qd.setQuarkEff(x[3] - quarks.xvd*quarks.gv*x[1]);

  if(quarks.iFlavor==3){
    quarks.qs.chemPot=x[4];
    quarks.qs.setQuarkEff(x[4] - quarks.xvs*quarks.gv*x[1]); 
  }
  
  quarks.setThermodynamics();
  
  electron.setLepton(x[5]);
  electron.calculateProperties();
  
  muon.setLepton(x[6]);
  muon.calculateProperties();

  double mub_ = (quarks.getEnergy() +electron.energy + muon.energy
      -quarks.temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + quarks.getPressure()+electron.pressure + muon.pressure)/x[0];
  // double mub_=quarks.Yu* quarks.qu.chemPot + quarks.Yd*quarks.qd.chemPot
  //                 +  quarks.Ys*quarks.qs.chemPot
	// 								+ (electron.chemPot*electron.density +muon.chemPot*muon.density)/x[0];

  residuals[0]=  quarks.muB - mub_;
  residuals[1] = quarks.omegaMeson_eom_residue(quarks.getOmegaEffDens());
  residuals[2]= quarks.qu.density - 3.*quarks.Yu*x[0];
  residuals[3]= quarks.qd.density - 3.*quarks.Yd*x[0];
  residuals[4]= quarks.qs.density - 3.*quarks.Ys*x[0];
  residuals[5]= electron.density  - quarks.Ye*x[0];
  residuals[6]= muon.density      - quarks.Ym*x[0];
  
  return true;
}
