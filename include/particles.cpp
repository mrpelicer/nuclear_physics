#include "particles.h"

//Call functions to calculate particle properties
void particle::calculateProperties(){
  energy           = 0.;
  pressure         = 0.;
  entropy          = 0.;  
  
  if(density>0.){
    if(!doamm){// if do amm==false
      if(!doB || Q==0){ //if magnetic field = 0 OR uncharged particle w/ B!=0.
        if(temperature>Tmin_integration){
          energy             = integrate(energyFunc,    this);
          pressure           = integrate(pressureFunc,  this);
          entropy            = integrate(entropyFunc,   this);
        }else{
          energy             = energyT0();
          pressure           = pressureT0();//chemPot*density - energy; 
        }
        // std::cout << "chargeless?: " <<  Q << " " << kf2 << " " 
        //           << mass_eff << endl;

      
      }else{ // dob, q!=0
        //size_t in=0;
        double nu_;
        double ener_=0., presst_=0., pressp_=0., mag_=0.;
        double eB=fabs(Bfield*Q*eHL);
        // kf2=pow(chemPot_eff, 2.);
        // std::cout << "test_AA: " << Q  << " " << mass*Mnucleon <<  endl;

        //do{
        for(int in=0; in<=inumax; in++){
          nu_=(double) in;
         double Mn2= pow(mass_eff, 2.)  + 2.*nu_*eB;
         double Mn= sqrt(Mn2);

          kf2= pow(chemPot_eff, 2.) - Mn2;
          kf= (kf2<=0.) ? 0. : sqrt(kf2);
          if(spin>1.){ //if delta
            if(in==0)       gamma=2.;
            else if(in==1)  gamma=3.;
            else            gamma=4.;
            // in==0 ? gamma=1. : gamma=2.;
          }else{   //spin 1/2 particles
            in==0 ? gamma=1. : gamma=2.;
          }
          if(kf>=0.){ //&& nu_<inumax
            // if(spin<1. || in<inumax){
            ener_  +=gamma*(kf*chemPot_eff + Mn2*log((kf+chemPot_eff)/Mn) ); 
            pressp_+=gamma*(kf*chemPot_eff - Mn2*log((kf+chemPot_eff)/Mn) );
            presst_+=gamma*nu_*log((kf+chemPot_eff)/Mn );
            mag_+=gamma*nu_*log((kf+chemPot_eff)/Mn );
            // }
            //if(spin>1){
                // cout << " testp: " <<  Q << " " << nu_ << " " << inumax << " "  << endl
                //      << kf2 << " " << kf << " " << gamma << " "  << endl
                //      << pressp_ << " " << presst_ << " " << mag_ << " " 
                //      << pressp_/2. -eB*presst_ << endl;
                    //}

          } //else{cout<< "no sum: " << nu_ << " " << kf2 << " " << kf << endl;}
         // in++;
          //std::cout << nu_ << " " << in << " " << kf2 << std::endl;
      }//}while(kf2>0.);
        energy=  eB*ener_/(4.*pi2);
        pressureParallel  = eB*pressp_/(4.*pi2); //chemPot_eff*density- energy; 
        pressureTransverse= pow(eB, 2.)*presst_/(2.*pi2);
        //  magnetization= -((energy - chemPot_eff*density)
        //               +pow(eB, 2.)*mag_/(2.*pi2))/Bfield ;
        magnetization= (pressureParallel - pressureTransverse)/Bfield;
      // if(spin>1){ cout << " testF: " << inumax << " " << in << " " <<  kf << " " << Q << " "
      //               <<  energy << " " << pressureParallel << " " 
      //               << pressureTransverse << " " << magnetization
      //               << endl;}
      }

    }else{// if doamm ==true:

      double Delta= -Bfield*Kb;
      double ener_=0.,presst_=0., pressp_=0., mag_=0.;
      double sp= 2.*spin;

      if(Q==0){ 
        for(double s=-sp; s<=sp; s+=2.){
          double Ms= mass_eff + s*Delta;
          kf2= pow(chemPot_eff, 2.) -pow(Ms, 2.);
          kf2<=0. ? kf=0. : kf=sqrt(kf2);
          if(kf >0.){
            ener_+= pow(chemPot_eff, 3.)*(kf/2. +2.*s*Delta*(asin(Ms/chemPot_eff) - M_PI/2.)/3. )
                    +(s*Delta/3.- Ms/4.)*(Ms*chemPot_eff*kf + pow(Ms, 3.)*log( (chemPot_eff + kf)/Ms));
            pressp_+= kf*chemPot_eff*(2.*pow(chemPot_eff, 2.) - 5.*Ms*Ms + 8.*s*Delta*Ms)
                    + 4.*s*Delta*pow(chemPot_eff, 3.)*(atan(Ms/kf) - M_PI/2.) 
                    + pow(Ms, 3.)*(3.*Ms - 4.*s*Delta)*log( (chemPot_eff + kf)/Ms);
            presst_+=kf*chemPot_eff*( 2.*pow(chemPot_eff, 2.) 
                    - 5.*Ms*Ms + 12.*s*Delta*Ms - 12.*pow(s*Delta, 2.))
                    + 3.*pow(Ms, 2.)*pow(Ms - 2.*s*Delta, 2.)*log( (chemPot_eff + kf)/Ms);
            mag_+= s*(  chemPot_eff*kf*(Ms - 3.*s*Delta) 
                      - pow(chemPot_eff, 3.)*(atan(Ms/kf) - M_PI/2.)  
                      - pow(Ms, 2.)*(2.*Ms - 3.*s*Delta)*log( (chemPot_eff + kf)/Ms) );
          }
        }
        energy= ener_/(4.*pi2);
        pressureParallel  =pressp_/(48.*pi2);//chemPot_eff*density- energy; //
        pressureTransverse=presst_/(48.*pi2);
        //magnetization=Kb*mag_/(12.*pi2);
        magnetization= (pressureParallel - pressureTransverse)/Bfield;

      }else{
        //size_t in=0;
        double nu_;
        double eB=std::fabs(Bfield*Q*eHL);
        double spmin, spmax;
        for(int in=0; in<=inumax; in++){
          nu_=(double) in;

          if(sp<=1.){//spin 1/2 spin sum limits
            if(in==0 && Q>0){
              spmin= sp;
              spmax= sp;
            }else if(in==0 && Q<0){
              spmin=-sp;
              spmax=-sp;
            }else{
              spmin=-sp;
              spmax= sp;
            }
          }else{//spin 3/2 spin sum limits
            if(Q>0){
              if(in==0){
               spmin= 1.;
               spmax= sp;
              }else if(in==1){
                spmin=-1.;
                spmax=sp;
              }else{
                spmin=-sp;
                spmax= sp;
              }
            }else{
              if(in==0){
               spmin= -sp;
               spmax= -1.;
              }else if(in==1){
                spmin=-sp;
                spmax=1.;
              }else{
                spmin=-sp;
                spmax= sp;
              }
            }
            }

          for(double s=spmin; s<=spmax; s+=2.){
            double Mn2= pow(mass_eff, 2.)  + 2.*nu_*eB;
            double Mn= sqrt(Mn2);
            double Ms= Mn+ s*Delta;
            kf2= pow(chemPot_eff, 2.) -pow(Ms, 2.);
            kf2<=0. ? kf=0. : kf=sqrt(kf2);
            if(kf>0.){
              ener_  +=kf*chemPot_eff+ pow(Ms, 2.)*log((chemPot_eff + kf)/Ms);
              pressp_+=kf*chemPot_eff- pow(Ms, 2.)*log((chemPot_eff + kf)/Ms);
              presst_+=(eB*nu_*Ms/Mn + s*Delta*Ms)*log((chemPot_eff + kf)/Ms);
            }
          }
         // in++;
       //}while(kf2>0.);
      }
      energy = eB*ener_/(4.*pi2);
      pressureParallel  =eB*pressp_/(4.*pi2); // chemPot_eff*density- energy; //
      pressureTransverse=eB*presst_/(2.*pi2);
      magnetization= (pressureParallel - pressureTransverse)/Bfield;

      } //end loop of charged particles
    
    }//end doamm
  
  }//end loop if(dens>0)

}


void particle::calculateDensity(){

  if(temperature>Tmin_integration){
		density            = integrate(densityFunc,    this);
	}else{
		density            = densityT0();
	}

 Qdens=Q*density;
}


void particle::calculateCondensate(){
  
  if(temperature>Tmin_integration){
		condensate         = integrate(density_condensateFunc,  this);
	}else{
    condensate = condensateT0();
	}
}

void particle::setBaryonEff(double mub_, double muq_, double gbphi, double gbv0, double gbb0){
//std::cout << doB << std::endl;
  chemPot= mub_+Q*muq_;
	mass_eff	=mass - gbphi;
	chemPot_eff= chemPot - gbv0 - gbb0*I3;

  if(!doamm){

    if(!doB || Q==0.){ 
      //std::cout << "neutral: " << Q << std::endl;
      kf2= pow(chemPot_eff, 2.) -pow(mass_eff, 2.);
      kf2<=0. ? kf=0. : kf=sqrt(kf2);
      calculateDensity();
      calculateCondensate();

      double sp= 2.*spin;
      if(sp<=1.){//spin 1/2 spin sum limits
        densityP=density/2.;
        densityM=density/2.;
      }else{//spin 3/2 spin sum limits
        densityPP=density/4.;
        densityP =density/4.;
        densityM =density/4.;
        densityMM=density/4.;
      }

    }else{ //if do B and Q!= 0
      // int in=0;
      double nu_;
      double rho_=0.;
      double rhos_=0.;
      double eB=std::fabs(Bfield*Q*eHL);
      double sp= 2.*spin;
      densityPP=0.;
      densityP =0.;
      densityM =0.;
      densityMM=0.;
      inumax= floor((pow(chemPot_eff, 2.)  - pow(mass_eff, 2.))/(2.*eB));

      // do{
      for(int in=0; in<=inumax; in++){
        nu_=(double) in;
        kf2= pow(chemPot_eff, 2.) -pow(mass_eff, 2.) - 2.*eB*nu_;
        kf= (kf2<=0.) ? 0. : sqrt(kf2);
        if(spin>1.){ //if delta
          if(in==0)       gamma=2.;
          else if(in==1)  gamma=3.;
          else            gamma=4.;
        }else{   //spin 1/2 particles
          in==0 ? gamma=1. : gamma=2.;
        }

        if(kf>0.){//&& mass_eff >0
          rho_+= gamma*kf ;
          rhos_+= gamma*log((kf+chemPot_eff)/sqrt(pow(mass_eff, 2.) + 2.*eB*nu_) );
          if(sp<=1.){//spin 1/2 spin sum limits
            if(in==0 && Q>0){
             densityP+=kf;
            }else if(in==0 && Q<0){
              densityM+=kf;
            }else{
              densityP+=kf;
              densityM+=kf;
            }
          }else{//spin 3/2 spin sum limits
            if(Q>0){
              if(in==0){
                densityPP+=kf;
                densityP+=kf;
              }else if(in==1){
                densityPP+=kf;
                densityP+=kf;
                densityM+=kf;
              }else{
                densityPP+=kf;
                densityP+=kf;
                densityM+=kf;
                densityMM+=kf;
              }
            }else{
              if(in==0){
                densityM+=kf;
                densityMM+=kf;
              }else if(in==1){
                densityP+=kf;
                densityM+=kf;
                densityMM+=kf;
              }else{
                densityPP+=kf;
                densityP+=kf;
                densityM+=kf;
                densityMM+=kf;
              }
            }
          }
        }
      //  in++;      
        
      //}while(kf2>0.);
      }
      density=  (rho_>0. && rhos_>0.) ?   eB*rho_/(2.*pi2)  : 0.;
      condensate= (rho_>0.&& rhos_>0.) ? eB*rhos_*mass_eff/(2.*pi2) : 0.;
      densityPP*=eB/(2.*pi2);
      densityP*=eB/(2.*pi2);
      densityM*=eB/(2.*pi2);
      densityMM*=eB/(2.*pi2);
      if(Q!=0.) Qdens= Q*density;
      //inumax= in-2;
    }//end doB, Q!=0

  }else{ //if do amm:
    double Delta= -Bfield*Kb;
    double rho_=0.;
    double rhos_=0.;
    double sp= 2.*spin;
    densityPP=0.;
    densityP=0.;
    densityM=0.;
    densityMM=0.;

    if(Q==0){ //do amm for q=0
      for(double s=-sp; s<=sp; s+=2.){
        // std::cout << sp << " " << s << std::endl;
        double Ms= mass_eff + s*Delta;
        kf2= pow(chemPot_eff, 2.) -pow(Ms, 2.);
        kf2<=0. ? kf=0. : kf=sqrt(kf2);
        if(kf >0. && mass_eff>0.){
          double rho_tmp=(pow(kf, 3.)/3. +s*Delta*( kf*Ms 
                                    + pow(chemPot_eff, 2.)*(asin(Ms/chemPot_eff) - M_PI/2.) )/2.); 
          rho_+= rho_tmp;
          rhos_+=chemPot_eff*kf -pow(Ms, 2.)*log( (chemPot_eff + kf)/Ms);

          if(s==-3.){
              densityMM=rho_tmp;
            }else if(s==-1.){
              densityM =rho_tmp;
            }else if(s==+1.){
              densityP =rho_tmp;
            }else if(s==+3.){
              densityPP=rho_tmp;
            }
        }
      }
      density= (rho_>0.&& rhos_>0.)    ? rho_/(2.*pi2) : 0.;      
      condensate= (rho_>0.&& rhos_>0.) ?  mass_eff*rhos_/(4.*pi2) : 0.;
      densityMM*=1./(2.*pi2);
      densityM *=1./(2.*pi2);
      densityP *=1./(2.*pi2);
      densityPP*=1./(2.*pi2);
    }else{ //doamm for Q!=0
      size_t in=0;
      double nu_;
      double eB=std::fabs(Bfield*Q*eHL);
      double spmin, spmax;
      double kf2P=0., kf2M=0.;
      double kf2PP=0., kf2MM=0.;
      do{
        nu_=(double) in;
        
        if(sp<=1.){//spin 1/2 spin sum limits
          if(in==0 && Q>0){
            spmin= sp;
            spmax= sp;
          }else if(in==0 && Q<0){
            spmin=-sp;
            spmax=-sp;
          }else{
            spmin=-sp;
            spmax= sp;
          }
        }else{//spin 3/2 spin sum limits
          if(Q>0){
            if(in==0){
             spmin= 1.;
             spmax= sp;
            }else if(in==1){
              spmin=-1.;
              spmax=sp;
            }else{
              spmin=-sp;
              spmax= sp;
            }
          }else{
            if(in==0){
             spmin= -sp;
             spmax= -1.;
            }else if(in==1){
              spmin=-sp;
              spmax=1.;
            }else{
              spmin=-sp;
              spmax= sp;
            }
          }

        }
        double Mn, Ms, Mn2;
        for(double s=spmin; s<=spmax; s+=2.){
          Mn2= pow(mass_eff, 2.) + 2.*nu_*eB;
          Mn=sqrt(Mn2);
          Ms= Mn+ s*Delta;
          kf2= pow(chemPot_eff, 2.) -pow(Ms, 2.);
          kf2<=0. ? kf=0. : kf=sqrt(kf2);
          if(s==-3.){//spin 1/2 spin sum limits
            kf2MM=kf2;
          }else if(s==-1.){
            kf2M=kf2;
          }else if(s==+1.){
            kf2P=kf2;
          }else if(s==+3.){
            kf2PP=kf2;
          }

          if(kf>0. && mass_eff>0.){
            rho_+= kf;
            rhos_+= (1.+ s*Delta/Mn)*log( (chemPot_eff + kf)/Ms);
            if(s==-3.){//spin 1/2 spin sum limits
              densityMM+=kf;
            }else if(s==-1.){
              densityM+=kf;
            }else if(s==+1.){
              densityP+=kf;
            }else if(s==+3.){
              densityPP+=kf;
            }
          }
        }
            in++;
          //std::cout << nu_ << " " << in << " " << kf2 << std::endl;
      }while(kf2P>0. || kf2M>0. || kf2PP>0. || kf2MM>0.);
      density   = (rho_>0. && rhos_>0.)   ?  eB*rho_/(2.*pi2)           : 0.;      
      condensate= (rho_>0. && rhos_>0.)  ?  eB*mass_eff*rhos_/(2.*pi2) : 0.;
      densityPP*=eB/(2.*pi2);
      densityP*=eB/(2.*pi2);
      densityM*=eB/(2.*pi2);
      densityMM*=eB/(2.*pi2);
      Qdens= Q*density;
      inumax= in-2;

    }
  }

}

void particle::setLepton(double mul_){
  chemPot= mul_;
	chemPot_eff= chemPot;

  if(!doB){  
   kf2= pow(chemPot, 2.) -pow(mass, 2.);
   kf2<=0. ? kf=0. : kf=sqrt(kf2);
   calculateDensity();
  }else if(doB && !doamm){
    // size_t in=0;
    double nu_;
    double rho_=0.;
    double eB=std::fabs(Bfield*Q*eHL);
     inumax= floor((pow(chemPot, 2.)  - pow(mass, 2.))/(2.*eB));
   
     // do{
     for(int in=0; in<=inumax; in++){
       nu_= (double) in;
      kf2= pow(chemPot, 2.) -pow(mass, 2.) - 2.*eB*nu_;
      kf2<=0. ? kf=0. : kf=sqrt(kf2);
      in==0 ? gamma=1. : gamma=2.;
      if(kf>0.){
        rho_+=  gamma*kf;
        if(in==0){
          densityM+=kf;
        }else{
          densityP+=kf;
          densityM+=kf;
        }
      }
     //  in++;
     }//}while(kf2>0.);

    rho_>0.  ? density   =eB*rho_/(2.*pi2) : density=0.;
    densityP*=eB/(2.*pi2);
    densityM*=eB/(2.*pi2);
    Qdens=Q*density;
    // inumax= in-2;
    //if(inumax<0) inumax=0;

  }else{
    size_t in=0;
    double nu_;
    double eB=std::fabs(Bfield*Q*eHL);
    double spmin, spmax;
    double sp=2.*spin;
    double Delta= -Bfield*Kb;
    double rho_=0.;
    double kF2P=0., kF2M=0.;
    do{
      nu_=(double) in;    
      
      if(in==0){
        spmin=-sp;
        spmax=-sp;
      }else{
        spmin=-sp;
        spmax= sp;
      }
    
      double Mn, Ms, Mn2;
      for(double s=spmin; s<=spmax; s+=2.){
        Mn2= pow(mass, 2.) + 2.*nu_*eB;
        Mn=sqrt(Mn2);
        Ms= Mn+ s*Delta;
        kf2= pow(chemPot, 2.) -pow(Ms, 2.);
        kf2<=0. ? kf=0. : kf=sqrt(kf2);
        if(s==-1.){
          kF2M=kf2;
        }else if(s==+1.){
          kF2P=kf2;
        }

        if(kf>0.){
          rho_+= kf;
          if(s==-1.){
            densityM+=kf;
          }else if(s==+1.){
            densityP+=kf;
          }
        } 
      }
    in++;
    }while(  kF2P>0.|| kF2M>0.); //
        // cout << kF2P << " " << kF2M << " " << kf2 << endl;

    density   = (rho_>0.)   ?  eB*rho_/(2.*pi2)           : 0.;      
    densityP*=eB/(2.*pi2);
    densityM*=eB/(2.*pi2);
    Qdens= Q*density;
    inumax= in-2;


  }

}

void particle::setBfield(bool dob_, double B_){
  doB=dob_;
  Bfield=B_;
}

void particle::setAMM(bool doa_,double kb_){
  doamm=doa_;
  Kb= kb_*eHL/2.;
}

//Get the (effective) chemical potential for a particle of given density.
void particle::solveChemPotEff(){
	
  // kf= cbrt(6.*pi2*density/gamma);
	double nu= hypot(kf, mass_eff);
 
 if(temperature>Tmin_integration){
    
		Problem problemNu;
  	CostFunction* nu_function =
  								new NumericDiffCostFunction<ChemPotFunctor, ceres::CENTRAL, 1, 1>
																												(new ChemPotFunctor(*this));
  	problemNu.AddResidualBlock(nu_function, NULL, &nu);
  	problemNu.SetParameterLowerBound(&nu, 0, 0.);
  	
		Solver::Options optionsNu;
  	Solver::Summary summaryNu;
		optionsNu.gradient_tolerance=1e-12;			//1e-10			12
		optionsNu.use_nonmonotonic_steps=true;	
		optionsNu.linear_solver_type= ceres::DENSE_QR;
  	optionsNu.minimizer_progress_to_stdout = false;
		Solve(optionsNu, &problemNu, &summaryNu);
	}
	chemPot_eff= nu;
}
//Functor:
template <typename T>
bool ChemPotFunctor::operator()(const T* arg, T* residuals) const {
	fermion.chemPot_eff= arg[0];
	T rhs=  integrate(densityFunc, &fermion);
	residuals[0] = fermion.density - rhs;
	return true;
}

double integrate(double (func)(double, void *), void *parametersPointer)
{
  particle &part= *reinterpret_cast<particle*>(parametersPointer);
  double result = 0e0;
  double er_res = 0e0;
  double err_abs = 1e-13; //1e-13
  double err_rel = 1e-10; //1e-10;
  size_t max_steps = 1e8;
	
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);

	gsl_function My_function;
  My_function.function = func;
  My_function.params = parametersPointer;
	
	if(part.temperature>Tmin_integration){
		gsl_integration_qagiu(&My_function, 0., err_abs, err_rel, max_steps, w, &result, &er_res);
  }else{
		if(part.kf!=part.kf){result=std::nan("1");}
		else{
		gsl_integration_qag(&My_function, 0., part.kf, err_abs, err_rel, max_steps, 1,
                                                           w, &result, &er_res);
		}
	}
		
	gsl_integration_workspace_free(w);
  return result;
}
// ================= Thermodynamic functions - T=0 =================


double particle::densityT0(){
  double dens=0.;
  if(kf>0.) dens=gamma*pow(kf, 3.)/(6.*pi2);
  return dens;
}

double particle::condensateT0(){
  double ener= hypot(kf, mass_eff);
  double cond=0.;
  if(kf>0. && mass_eff >0)
     cond= gamma*mass_eff*(kf*ener - pow(mass_eff, 2.)*log( (kf+ener)/mass_eff )) /(4.*pi2);
  
  return cond;
}

double particle::energyT0(){
  double ener= hypot(kf, mass_eff);
  //  return  gamma*( ( 2.*pow(kf, 2.) + pow(mass_eff, 2.) )*kf*ener 
  //                 - pow(mass_eff, 4.)*log( (kf+ener)/mass_eff )) /(16.*pi2);

   return  gamma*(kf*pow(ener, 3.) + pow(kf, 3.)*ener 
                  - pow(mass_eff, 4.)*log( (kf+ener)/mass_eff )) /(16.*pi2);

}

double particle::pressureT0(){
  double ener= hypot(kf, mass_eff);
  return  gamma*( (2.*pow(kf, 3.) - 3.*pow(mass_eff, 2.)*kf)*ener 
                +  3.*pow(mass_eff, 4.)*log( (kf+ener)/mass_eff )) /(48.*pi2);
}



// ============================ DDQM particles =============================
void quark_particle::calculateQProperties(){
  energy           = 0.;
  pressure         = 0.;
  entropy          = 0.;  
  if(density>0.){

    if(temperature>Tmin_integration){
     omega0      = getOmega0(chemPot_eff, mass_eff);
     entropy     = -getDOmega0Dt();
     energy      = omega0+ chemPot_eff*density+temperature*entropy;
     pressure    = -omega0;
    }else{
      omega0     = getOmega0(chemPot_eff, mass_eff);
      energy     = omega0+ chemPot_eff*density; //+temperature*entropy
      pressure   = -omega0;
    }
  }
}

double quark_particle::getOmega0(double mueff_, double masseff_){
 double omega=0.;
 if(temperature<Tmin_integration){
  double kf2_=pow(mueff_, 2.) - pow(masseff_, 2.);
  double kf_=(kf2_>0) ? sqrt(kf2_ ) : 0.;
   if(kf_>0.)
      omega= -gamma*(mueff_*kf_*(2.*pow(mueff_, 2.) - 5.*pow(masseff_, 2.))
                     + 3.*pow(masseff_, 4.)*log((mueff_ +kf_)/masseff_ ) )/(48*pi2);
 }else{
  double chemPot_eff_Save = chemPot_eff;
  double mass_eff_Save    = mass_eff;
  chemPot_eff= mueff_;
  mass_eff= masseff_;
  omega= -integrate(pressureFunc,this);
  chemPot_eff = chemPot_eff_Save;
  mass_eff    = mass_eff_Save;
 }
  return omega;
}

double quark_particle::getDOmega0Dmass(){
  double dOmega0dM = ( getOmega0(chemPot_eff, mass_eff-2.*hdif)
                  - 8.*getOmega0(chemPot_eff, mass_eff-hdif)
                  + 8.*getOmega0(chemPot_eff, mass_eff+hdif)
                  - getOmega0(chemPot_eff, mass_eff+2.*hdif))/(12.*hdif);
  return dOmega0dM;

}
double quark_particle::getDOmega0DmuEf(){
  double dOmega0dmuef = ( getOmega0(chemPot_eff-2.*hdif, mass_eff)
                - 8.*getOmega0(chemPot_eff-hdif, mass_eff)
                + 8.*getOmega0(chemPot_eff+hdif, mass_eff)
                - getOmega0(chemPot_eff+2.*hdif, mass_eff))/(12.*hdif);
  return dOmega0dmuef;

}

double quark_particle::getDOmega0Dt(){
  double temp_save=temperature;
  temperature=temp_save+hdif;
  double omg_p1=  getOmega0(chemPot_eff, mass_eff);
  temperature=temp_save+2.*hdif;
  double omg_p2=  getOmega0(chemPot_eff, mass_eff);
  temperature=temp_save-hdif;
  double omg_m1=  getOmega0(chemPot_eff, mass_eff);
  temperature=temp_save-2.*hdif;
  double omg_m2=  getOmega0(chemPot_eff, mass_eff);

  double dOmega0dT = ( omg_m2 
                  - 8.*omg_m1
                  + 8.*omg_p1
                  - omg_p2)/(12.*hdif);

  temperature=temp_save;  

  return dOmega0dT;
}


// ================= GSL functions to be integrated - T >0 =================

double densityFunc(double x, void *p){
  particle &part_= *reinterpret_cast<particle*>(p);
  double ener= sqrt ( pow(x, 2.) + pow(part_.mass_eff, 2.) );
  double F;

  if(part_.temperature!=0){
    double fdp, fdm;
    fdp= fermiDirac(ener, part_.chemPot_eff, part_.temperature);
    fdm= fermiDirac(ener, -part_.chemPot_eff, part_.temperature);
    F=fdp-fdm;
  }else{
    F=1.;
  }

  return part_.gamma*pow(x, 2.)*F/(2.*pi2);
}

double density_condensateFunc(double x, void *p){
  particle &part_= *reinterpret_cast<particle*>(p);
  double ener= sqrt ( pow(x, 2.) + pow(part_.mass_eff, 2.) );
  double F;

  if(part_.temperature!=0){
    double fdp, fdm;
    fdp= fermiDirac(ener, part_.chemPot_eff, part_.temperature);
    fdm= fermiDirac(ener, -part_.chemPot_eff, part_.temperature);
    F=fdp+fdm;
  }else{
    F=1.;
  }

  return part_.gamma*pow(x, 2.)*part_.mass_eff*F/(2.*ener*pi2);
}

double energyFunc(double x, void *p){
  particle &part_= *reinterpret_cast<particle*>(p);
  double ener= sqrt ( pow(x, 2.) + pow(part_.mass_eff, 2.) );
  double F;

  if(part_.temperature!=0){
    double fdp, fdm;
    fdp= fermiDirac(ener, part_.chemPot_eff, part_.temperature);
    fdm= fermiDirac(ener, -part_.chemPot_eff, part_.temperature);
    F=fdp+fdm;
  }else{
    F=1.;
  }

  return part_.gamma*pow(x, 2.)*ener*F/(2.*pi2);
}

double pressureFunc(double x, void *p){
  particle &part_= *reinterpret_cast<particle*>(p);
  double ener= sqrt ( pow(x, 2.) + pow(part_.mass_eff, 2.) );
  double F;

  if(part_.temperature!=0){
    double fdp, fdm;
    fdp= fermiDirac(ener,  part_.chemPot_eff, part_.temperature);
    fdm= fermiDirac(ener, -part_.chemPot_eff, part_.temperature);
    F=fdp+fdm;
  }else{
    F=1.;
  }

  return part_.gamma*pow(x, 4.)*F/(6.*pi2*ener);
}

double entropyFunc(double x, void *p){
  particle &part_= *reinterpret_cast<particle*>(p);
  double ener= sqrt ( pow(x, 2.) + pow(part_.mass_eff, 2.) );

  double fdp= fermiDirac(ener, +part_.chemPot_eff, part_.temperature);
  double fdm= fermiDirac(ener, -part_.chemPot_eff, part_.temperature);
  double tp, tm;

  if(fdp==1. || fdp==0.){
    tp= 0.;
  }else{
    tp= fdp*log(fdp)+ (1.-fdp)*log(1.-fdp);
  }

  if(fdm==1. || fdm==0.){
    tm=0.;
  }else{
    tm= fdm*log(fdm)+ (1.-fdm)*log(1.-fdm);
  }

  double F = tp + tm;

  return -part_.gamma*pow(x, 2.)*F/(2.*pi2);
}

double fermiDirac(double ener, double chemPotEff, double T){
  double x1=(ener-chemPotEff)/T;
  return 1./(exp(x1)+1.);
}

