//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
#include "../../include/interpolator.h"
#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_interp.h>	
#include <gsl/gsl_spline.h>

int main(){
//Choose pararametrization
	std::string parametrization= "ddme2";
	nlwm_class hmg_matter(parametrization);

	bool doHyperons	=	true;
	bool doDeltas		=	true;
	bool doBfield		= false;
	bool doamm 			= false;

	std::string hyperon_params ="ddme2-a"; 
	std::string delta_params 	 ="su6";  //su6(1.), mplA_1(beta=1.1), mplA_2, prd89_1, prd89_1
	double rhoB, temperature, Bfield;
	temperature=50.;
	double Bg=3.e18; // Gauss
	
	
	double rhoBMin=5e-3/pow(Mnucleon/hc, 3);
  double rhoBMax=1./pow(Mnucleon/hc, 3);///7.5*hmg_matter.rho0;
	//0.62/pow(Mnucleon/hc, 3); fsu2h c amm ou b
  int iR=240;
  double dRho=  (rhoBMax-rhoBMin)/iR;

	//Construct a gas of p,n,e matter
	// nlwm_class derivative(parametrization);

	particle electron, muon;
	electron.mass	= Me/Mnucleon;	
  muon.mass			=	Mm/Mnucleon; 
	electron.mass_eff= electron.mass;
	muon.mass_eff= muon.mass;	
	electron.Q=-1.;
	muon.Q=-1.;
	electron.gamma=2.;
	muon.gamma=2.;
	double Bc=pow(electron.mass_eff, 2.)/eHL;  
	Bfield=Bg*Bc/4.41e13;
	//Bfield *= 1.95e-14/pow(Mnucleon, 2.);  // Conversion factor Gauss to MeV^2 to adim
	//Define thermodynamic variables
	double Energy, FreeEn, Pressure, PressureP, PressureT, Entropy, Magnetization;

	int iRLog= 79;

	if(doHyperons){ 
		//double as=1.;
		//double av=1.;
	 	hmg_matter.includeHyperons(doHyperons, hyperon_params);
	}

	if(doDeltas){
		hmg_matter.includeDeltas(		doDeltas, 	delta_params);
	}
	hmg_matter.printParameters();
	if(doBfield){
		hmg_matter.setBfield(	doBfield, Bfield);
		electron.setBfield(		doBfield, Bfield);
	 	muon.setBfield(				doBfield, Bfield);

		if(doamm){
			hmg_matter.setAMM(doamm);
			electron.setAMM(doamm, (1.15965e-3)/electron.mass);
			muon.setAMM(		doamm, (1.16592e-3)/muon.mass);
			// electron.setAMM(doamm, 0.);
			// muon.setAMM(doamm,		 0.);
		}
	}

	hmg_matter.printParameters();
	std::vector<double> enerv, enerdensv, pressv, pressPv, pressTv, rhobv;

	double yN=0., yH=0., yD=0.;

		//guess mue:
	//electron.chemPot=0.1; /* as long as it is smaller than ~0.5, it should converge 
	//																							 if not, try changing parametrization. */

	//std::vector<double> rhobv, enerv, pressv; 
	// for(size_t it= 0; it<temp_vec.size(); it++){
		 
		//Set output files
		std::cout << "T= " << temperature << std::endl;
		// std::ofstream outUbar("../data/data_hmg/potentials_" +parametrization+".txt");
		
		//outFile << "#rho_b P freeEn mun mup mue M* V0 b0 Y_p" << std::endl;
		
		std::string filename1, filename2, filename3, filename4;
		if(!doBfield){
			filename1="data/dens_beta_" +parametrization+"_noB.txt";
			// filename2="data/hmg_beta_" 	+parametrization+"_noB.txt";
			filename2="data/press_" 	+parametrization+"_T"+to_string(temperature)+".txt";
			filename3="data/spin_dens_"	+parametrization+"_noB.txt";
			filename4="data/eos_"				+parametrization+"_noB.txt";
		}else if(doBfield && !doamm){
			filename1="data/dens_beta_" 	+parametrization+"_wtB.txt";
			filename2="data/hmg_beta_"		+parametrization+"_wtB.txt";
			filename3="data/spin_dens_" 	+parametrization+"_wtB.txt";
			filename4="data/eos_"					+parametrization+"_wtB.txt";
		}else{
			filename1="data/dens_beta_" 	+parametrization+"_wtA.txt";
			filename2="data/hmg_beta_"		+parametrization+"_wtA.txt";
			filename3="data/spin_dens_"		+parametrization+"_wtA.txt";
			filename4="data/eos_"					+parametrization+"_wtA.txt";
		}
//../data/data_hmg/
		
		std::ofstream outDens(filename1);
		std::ofstream outFile(filename2);
		std::ofstream outSpin(filename3);
		std::ofstream outEos(filename4);

		//Adimensional temperature;
		temperature*=1./Mnucleon;
		electron.temperature=temperature;
		muon.temperature=temperature;

		//=== Loop over barionic density
		for(int irho=0; irho<iR; irho++){
			rhoB=rhoBMax- (double)irho*dRho;
			//  rhoB=rhoBMax+ (double)irho*dRho;
			
			//Solve self-consistently:
			hmg_matter.setEOS_betaEq(rhoB, temperature, electron, muon);
			// cout << "test: " << hmg_matter.proton.chemPot_eff + hmg_matter.getRearrangementEnergy()  << " " << hmg_matter.neutron.chemPot_eff + hmg_matter.getRearrangementEnergy() << " " 
			// 									<< hmg_matter.getRearrangementEnergy() << " " << hmg_matter.getDerivativeCoupling_sigma(rhoB) << " "
			// 									 << hmg_matter.getDerivativeCoupling_omega(rhoB) << " " << hmg_matter.getDerivativeCoupling_rho(rhoB) 
			// 									 << endl;
			// std::cout << hmg_matter.proton.mass_eff << " " << hmg_matter.delta0.mass_eff << " " << hmg_matter.deltap.mass_eff << " " <<  hmg_matter.delta0.chemPot << std::endl;
			electron.pressure	=electron.chemPot*electron.density 	- electron.energy + temperature*electron.entropy;
			muon.pressure			=muon.chemPot*muon.density 					- muon.energy			+ temperature*muon.entropy;
			//Set variables:
			Energy		= hmg_matter.getEnergy() 		+ electron.energy	 + muon.energy;
			Entropy		= hmg_matter.getEntropy() 	+ electron.entropy + muon.entropy;
			FreeEn		= Energy -temperature*Entropy;
			Pressure	= -Energy + temperature*Entropy + hmg_matter.neutron.chemPot*rhoB;
			//hmg_matter.getPressure()	+ electron.pressure+ muon.pressure;
			PressureP	= hmg_matter.getPressureParallel()+ electron.pressureParallel+ muon.pressureParallel;
			PressureT	= hmg_matter.getPressureTransverse()+ electron.pressureTransverse+ muon.pressureTransverse;
			Magnetization= hmg_matter.getMagnetization() + electron.magnetization+muon.magnetization;

			rhobv.push_back(rhoB);
			enerv.push_back(Energy/rhoB);
			enerdensv.push_back(Energy);
			pressv.push_back(Pressure);
			pressPv.push_back(PressureP);
			pressTv.push_back(PressureT);

			yN= (hmg_matter.proton.density + hmg_matter.neutron.density)/rhoB;
			yH= (hmg_matter.lambda0.density + hmg_matter.sigmam.density+ hmg_matter.sigma0.density+ hmg_matter.sigmap.density
						+ hmg_matter.xi0.density+ hmg_matter.xim.density)/rhoB;
			yD= (hmg_matter.deltapp.density + hmg_matter.deltap.density+ hmg_matter.delta0.density+ hmg_matter.deltam.density)/rhoB;
		  
			outFile << rhoB*pow(Mnucleon/hc, 3) << " "
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3) << " " 
					<< (FreeEn/rhoB - 1.)*Mnucleon  << " " 
					<< Energy*Mnucleon*pow(Mnucleon/hc, 3)  << " " <<  hmg_matter.neutron.mass_eff*Mnucleon << " "
					<< hmg_matter.neutron.chemPot*Mnucleon 			<< " " << electron.chemPot*Mnucleon << " " 
					<< yN << " " << yH << " " << yD  << " "
					<< hmg_matter.proton.chemPot_eff + hmg_matter.getRearrangementEnergy()  << " " << hmg_matter.neutron.chemPot_eff + hmg_matter.getRearrangementEnergy() << " " 
					<< hmg_matter.getRearrangementEnergy() << " " << hmg_matter.getDerivativeCoupling_sigma(rhoB) << " "
					 << hmg_matter.getDerivativeCoupling_omega(rhoB) << " " << hmg_matter.getDerivativeCoupling_rho(rhoB) 
					<<std::endl;

			outSpin << rhoB*pow(Mnucleon/hc, 3) << " " // hmg_matter.rho0 or *pow(Mnucleon/hc, 3) 
					<<	hmg_matter.proton.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.proton.densityM*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.neutron.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.neutron.densityM*pow(Mnucleon/hc, 3) << " "
					<<	electron.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	electron.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	muon.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	muon.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.sigmap.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.sigmap.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.sigma0.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.sigma0.densityM*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.sigmam.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.sigmam.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.xi0.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.xi0.densityM*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.xim.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.xim.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityPP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityMM*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.deltap.densityPP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltap.densityP*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.deltap.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltap.densityMM*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.delta0.densityPP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.delta0.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.delta0.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.delta0.densityMM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltam.densityPP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltam.densityP*pow(Mnucleon/hc, 3) << " "
					<<	hmg_matter.deltam.densityM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.deltam.densityMM*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.lambda0.densityP*pow(Mnucleon/hc, 3) << " " 
					<<	hmg_matter.lambda0.densityM*pow(Mnucleon/hc, 3)
					<<  std::endl;

			outDens << rhoB*pow(Mnucleon/hc, 3) 									<< " " // /hmg_matter.rho0 			 << " " // *pow(Mnucleon/hc, 3) << " " 
				<< hmg_matter.proton.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.neutron.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< electron.density						*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< muon.density								*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.lambda0.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.sigmap.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.sigma0.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.sigmam.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB	 << " "
				<< hmg_matter.xi0.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.xim.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.deltapp.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.deltap.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.delta0.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< hmg_matter.deltam.density	*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
				<< std::endl;

			// outEos << rhoB*pow(Mnucleon/hc, 3) << " "
					// << Energy*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_gdivcm3  << " "  
					// << Pressure*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_dyndivcm2
					// << std::endl;

				outEos << rhoB*pow(Mnucleon/hc, 3) << " "
					<< (Energy+Bfield*Bfield/2.)*Mnucleon*pow(Mnucleon/hc, 3) << " "  
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3) << " "
					<< (PressureP-Bfield*Bfield/2.)*Mnucleon*pow(Mnucleon/hc, 3) << " "
					<< (PressureT+Bfield*Bfield/2.)*Mnucleon*pow(Mnucleon/hc, 3)
				<< std::endl;
				
		}
		
	outFile.close();
	outSpin.close();
	outDens.close();

	//  double rhoBMaxLog= log10(rhoB);
	//  std::cout << rhoB*pow(Mnucleon/hc, 3) << std::endl;
	//  double rhoBMinLog= log10( (1e-14)/pow(Mnucleon/hc, 3) );
	//  double dlogR= (rhoBMaxLog-rhoBMinLog)/iRLog;
	//  for(int irl=1; irl<iRLog; irl++){
	// 		//rhoB=rhoBMax- (double)irho*dRho;
	 		// rhoB= pow(10., rhoBMaxLog- irl*dlogR);
	 			// hmg_matter.setAMM(false);
	 			// electron.setAMM(false, 0.);
	 			// muon.setAMM(false, 0.);
	 		//  if(rhoB*pow(Mnucleon/hc, 3) < 1e-4){
	  			// hmg_matter.setBfield(	false, Bfield);
	 			// electron.setBfield(		false, Bfield);
	  			// muon.setBfield(				false, Bfield);
	// 
	 		//  }
	// 		Solve self-consistently:
	 		// hmg_matter.setEOS_betaEq(rhoB, temperature, electron, muon);
 		//std::cout << hmg_matter.proton.mass_eff << " " << hmg_matter.delta0.mass_eff << " " << hmg_matter.deltap.mass_eff << " " <<  hmg_matter.delta0.chemPot << std::endl;
	 	//	electron.pressure	=electron.chemPot*electron.density 	- electron.energy;
	 	//	muon.pressure			=muon.chemPot*muon.density 					- muon.energy;
	 	//	Set variables:
	 		// Energy= hmg_matter.getEnergy() 		+ electron.energy+ muon.energy;
	 		// Pressure= hmg_matter.getPressure()+ electron.pressure+ muon.pressure;
	 		// Entropy= hmg_matter.getEntropy() 	+ electron.entropy + muon.entropy;
	 		// FreeEn= Energy -temperature*Entropy;
		 		// outEos << iRLog-irl << " " << rhoB*pow(Mnucleon/hc, 3) << " "
	 				// << Energy*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_gdivcm3  << " "  
	 				// << Pressure*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_dyndivcm2
	 				// << std::endl;
	//  }
	
	outEos.close();

	// }

  std::reverse(rhobv.begin(), rhobv.end());
	std::reverse(enerv.begin(), enerv.end());
	std::reverse(enerdensv.begin(), enerdensv.end());
	std::reverse(pressv.begin(), pressv.end());
	std::reverse(pressPv.begin(), pressPv.end());
	std::reverse(pressTv.begin(), pressTv.end());
	

	// std::vector<double> enerd1v, enerd2v, enerd3v; 
	// for(int irho=0; irho<iR; irho++){
	// 	rhoB=	rhobv[0]+irho*dRho;
	// 	enerd1v.push_back(deriv_func(rhoB, enerv, rhobv));
	// }
	// for(int irho=0; irho<iR; irho++){
	// 	rhoB=	rhobv[0]+irho*dRho;
	// 	enerd2v.push_back(deriv_func(rhoB, enerd1v, rhobv));
	// }
	// for(int irho=0; irho<iR; irho++){
	// 	rhoB=	rhobv[0]+irho*dRho;
	// 	enerd3v.push_back(deriv_func(rhoB, enerd2v, rhobv));
	// }

	// std::vector<double> k0v, Kv, q0v, k1v, q1v, cs2v, cs2Pv, cs2Tv;
	// double cs2=0., cs2p, cs2t, K_=0.;

	// for(int irho=0; irho<iR; irho++){
	// 	rhoB=	rhobv[0]+irho*dRho;
	// 	k1v.push_back(pow(3.*rhoB, 2.)*enerd2v[irho]);
	// 	q1v.push_back(pow(3.*rhoB, 3.)*enerd3v[irho]);
	// 	cs2=deriv_func(rhoB, pressv, rhobv)/deriv_func(rhoB, enerdensv, rhobv);
	// 	cs2p=deriv_func(rhoB, pressPv, rhobv)/deriv_func(rhoB, enerdensv, rhobv);
	// 	cs2t=deriv_func(rhoB, pressTv, rhobv)/deriv_func(rhoB, enerdensv, rhobv);
	// 	K_= deriv_func(rhoB, pressv, rhobv);
	// 	double h=1e-4;
	// 	double k0_=0., q0_=0., dpdr=0., dedr=0., cs2_=0.;
	// 	if( ((rhoB-2.*h) > rhobv[0]) && ((rhoB+2.*h) < rhobv[rhobv.size()-1]) ){ //http://web.media.mit.edu/~crtaylor/calculator.html
	// 		k0_= (-interpolation_func(rhoB-2.*h, enerv, rhobv)+ 16.*interpolation_func(rhoB-h, enerv, rhobv) 
	// 								- 30.*interpolation_func(rhoB, enerv, rhobv)
	// 								+ 16.*interpolation_func(rhoB+h, enerv, rhobv) -interpolation_func(rhoB+2.*h, enerv, rhobv))/(12.*h*h);
	// 		q0_= (-interpolation_func(rhoB-2.*h, enerv, rhobv)+ 2.*interpolation_func(rhoB-h, enerv, rhobv) 
	// 								- 2.*interpolation_func(rhoB+h, enerv, rhobv) +interpolation_func(rhoB+2.*h, enerv, rhobv))/(2.*h*h*h);
			
	// 		dpdr= (interpolation_func(rhoB-2.*h, pressv, rhobv) - 8.*interpolation_func(rhoB-h, pressv, rhobv) 
	// 								+8.*interpolation_func(rhoB+h, pressv, rhobv) - interpolation_func(rhoB+2.*h, pressv, rhobv))/(12.*h);
	// 		dedr= (interpolation_func(rhoB-2.*h, enerdensv, rhobv) - 8.*interpolation_func(rhoB-h, enerdensv, rhobv) 
	// 								+8.*interpolation_func(rhoB+h, enerdensv, rhobv) - interpolation_func(rhoB+2.*h, enerdensv, rhobv))/(12.*h);
	// 		cs2_= dpdr/dedr;
	// 	}else{
	// 		k0_=enerd2v[irho];
	// 		q0_=enerd3v[irho];
	// 		cs2_=cs2;
	// 	}
	// 	Kv.push_back(9.*K_);
	// 	k0v.push_back(pow(3.*rhoB, 2.)*k0_);
	// 	q0v.push_back(pow(3.*rhoB, 3.)*q0_);
	// 	cs2v.push_back(cs2_);
	// 	cs2Pv.push_back(cs2p);
	// 	cs2Tv.push_back(cs2t);
		
	// }

	// std::string outbulkfile;
	// if(!doBfield){
	// 	outbulkfile="data/bulk_"+parametrization+"_noB.txt";
	// }else if(doBfield && !doamm){
	// 	outbulkfile="data/bulk_"+parametrization+"_wtB.txt";
	// }else{
	// 	outbulkfile="data/bulk_"+parametrization+"_wtA.txt";
	// }

	// std::ofstream outBulk(outbulkfile);

	// for(int irho=0; irho<iR; irho++){
	// 	rhoB=	rhobv[0]+irho*dRho;

	// 	outBulk	<<  rhoB*pow(Mnucleon/hc, 3) << " " 
	// 						<< interpolation_func(rhoB, pressv, rhobv)*Mnucleon*pow(Mnucleon/hc, 3) << " "
	// 						<< interpolation_func(rhoB, enerdensv, rhobv)*Mnucleon*pow(Mnucleon/hc, 3) << " "
	// 						<< (interpolation_func(rhoB, enerv, rhobv)-1.)*Mnucleon << " "
	// 					 	<< interpolation_func(rhoB, Kv, rhobv)*Mnucleon << " "
	// 						<< interpolation_func(rhoB, k0v, rhobv)*Mnucleon << " "
	// 						<< interpolation_func(rhoB, cs2v, rhobv) << " " 
	// 						<< interpolation_func(rhoB, cs2Pv, rhobv) << " " 
	// 						<< interpolation_func(rhoB, cs2Tv, rhobv) << " " 
	// 						<< ( interpolation_func(rhoB, cs2Pv, rhobv) + 
	// 									2.*interpolation_func(rhoB, cs2Tv, rhobv) )/3. 
	// 						<< std::endl;
	// 													// << interpolation_func(rhoB, k0v, rhobv)*Mnucleon << " "
	// 						// << interpolation_func(rhoB, cs2v, rhobv) << " " 

	// }
	// outBulk.close();

  return 0;
}
