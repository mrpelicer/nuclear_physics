//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_walecka.h"
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

	bool doHyperons	=	false;
	bool doDeltas	=	false;
	// bool doBfield		= false;
	// bool doamm 			= false;
	// double Bg=3.e18; // Ga;uss

	std::string hyperon_params ="ddme2-a"; 
	std::string delta_params 	 ="su6";  //su6(1.), mplA_1(yl=1.1), mplA_2, prd89_1, prd89_1
	double rhoB, temperature; //, Bfield;
	temperature=0.;
	
	
	double rhoBMin=1e-2/pow(Mnucleon/hc, 3);
  double rhoBMax=.3/pow(Mnucleon/hc, 3);///7.5*hmg_matter.rho0;
	//0.62/pow(Mnucleon/hc, 3); fsu2h c amm ou b
  int iR=240;
  double dRho=  (rhoBMax-rhoBMin)/iR;

	//Construct a gas of p,n,e matter
	// nlwm_class derivative(parametrization);


	double s_nb=0.;
	cout << "Set entropy per density (S/nb):" << endl;
	cin >> s_nb ;

	double Yle=0.;
	cout << "Set fraction of e-leptons (Y_le):" << endl;
	cin >> Yle ;
	assert(Yle<=0.5);

	particle electron, neut_e;

	neut_e.gamma=1.;

	electron.mass	= Me/Mnucleon;	
	neut_e.mass=0.;
	
	electron.mass_eff= electron.mass;
	neut_e.mass_eff=0.;

	electron.Q=-1.;
	neut_e.Q=0.;
	electron.gamma=2.;
	neut_e.gamma=1.;

	// double Bc=pow(electron.mass_eff, 2.)/eHL;  

	//Define thermodynamic variables
	double Energy, FreeEn, Pressure, PressureP, PressureT, Entropy, Magnetization;
	std::vector<double> enerv, enerdensv, pressv, pressPv, pressTv, rhobv;
	double yN=0., yH=0., yD=0.;


	if(doHyperons){ 
	 	hmg_matter.includeHyperons(doHyperons, hyperon_params);
	}

	if(doDeltas){
		hmg_matter.includeDeltas(		doDeltas, 	delta_params);
	}
	hmg_matter.printParameters();
		
	std::string filename1, filename2, filename3, filename4;
	
	filename1="data/dens_yl_" +parametrization+"_noB.txt";
	filename2="data/press_" 	+parametrization+".txt";
	filename3="data/eos_"				+parametrization+"_noB.txt";
	
	std::ofstream outDens(filename1);
	std::ofstream outFile(filename2);
	std::ofstream outEos(filename3);

	//=== Loop over barionic density
	for(int irho=0; irho<iR; irho++){
		rhoB=rhoBMax- (double)irho*dRho;
		
				// hmg_matter.setEOS_betaEq(rhoB, temperature, electron, muon);
			hmg_matter.setEOS_fixedYl(rhoB, s_nb, Yle, 	electron, neut_e);

			// cout << "test: " << hmg_matter.proton.chemPot_eff + hmg_matter.getRearrangementEnergy()  << " " << hmg_matter.neutron.chemPot_eff + hmg_matter.getRearrangementEnergy() << " " 
			// 									<< hmg_matter.getRearrangementEnergy() << " " << hmg_matter.getDerivativeCoupling_sigma(rhoB) << " "
			// 									 << hmg_matter.getDerivativeCoupling_omega(rhoB) << " " << hmg_matter.getDerivativeCoupling_rho(rhoB) 
			// 									 << endl;
			// std::cout << hmg_matter.proton.mass_eff << " " << hmg_matter.delta0.mass_eff << " " << hmg_matter.deltap.mass_eff << " " <<  hmg_matter.delta0.chemPot << std::endl;
			electron.pressure	=electron.chemPot*electron.density 	- electron.energy + temperature*electron.entropy;
			// muon.pressure			=muon.chemPot*muon.density 					- muon.energy			+ temperature*muon.entropy;
			//Set variables:
			Energy		= hmg_matter.getEnergy() 		+ electron.energy	  ;// + muon.energy;
			Entropy		= hmg_matter.getEntropy() 	+ electron.entropy  ;// + muon.entropy;
			FreeEn		= Energy -temperature*Entropy;
			Pressure	= -Energy + temperature*Entropy + hmg_matter.neutron.chemPot*rhoB;
			//hmg_matter.getPressure()	+ electron.pressure+ muon.pressure;
			PressureP	= hmg_matter.getPressureParallel()+ electron.pressureParallel ;// + muon.pressureParallel;
			PressureT	= hmg_matter.getPressureTransverse()+ electron.pressureTransverse ;// + muon.pressureTransverse;
			Magnetization= hmg_matter.getMagnetization() + electron.magnetization ;// +muon.magnetization;

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
		  
			outFile << rhoB*pow(Mnucleon/hc, 3) << " " << hmg_matter.temperature*Mnucleon << " " 
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3) << " " 
					<< (FreeEn/rhoB - 1.)*Mnucleon  << " " 
					<< Energy*Mnucleon*pow(Mnucleon/hc, 3)  << " " <<  hmg_matter.neutron.mass_eff*Mnucleon << " "
					<< hmg_matter.neutron.chemPot*Mnucleon 			<< " " << electron.chemPot*Mnucleon << " " 
					<< yN << " " << yH << " " << yD  << " "
					<< hmg_matter.proton.chemPot_eff + hmg_matter.getRearrangementEnergy()  << " " << hmg_matter.neutron.chemPot_eff + hmg_matter.getRearrangementEnergy() << " " 
					<< hmg_matter.getRearrangementEnergy() << " " << hmg_matter.getDerivativeCoupling_sigma(rhoB) << " "
					 << hmg_matter.getDerivativeCoupling_omega(rhoB) << " " << hmg_matter.getDerivativeCoupling_rho(rhoB) 
					<<std::endl;

			outDens << rhoB*pow(Mnucleon/hc, 3) 									<< " " // /hmg_matter.rho0 			 << " " // *pow(Mnucleon/hc, 3) << " " 
				<< hmg_matter.proton.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.neutron.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< electron.density				 /rhoB  << " " // 		*pow(Mnucleon/hc, 3) << " "
				<< 0.							 /rhoB  << " " // 	*pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.lambda0.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.sigmap.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.sigma0.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.sigmam.density	 /rhoB	 << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.xi0.density		 /rhoB  << " " // 	*pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.xim.density		 /rhoB  << " " // 	*pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.deltapp.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.deltap.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.delta0.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< hmg_matter.deltam.density	 /rhoB  << " " // *pow(Mnucleon/hc, 3) << " "
				<< neut_e.density				 /rhoB  << " " // 		*pow(Mnucleon/hc, 3) << " "
				<< std::endl;

			// outEos << rhoB*pow(Mnucleon/hc, 3) << " "
					// << Energy*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_gdivcm3  << " "  
					// << Pressure*Mnucleon*pow(Mnucleon/hc, 3)*MeVdivfm3_to_dyndivcm2
					// << std::endl;

				outEos << rhoB*pow(Mnucleon/hc, 3) << " "
					<< Energy*Mnucleon*pow(Mnucleon/hc, 3) << " "  
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3)
				<< std::endl;
				
		}
		
	outFile.close();
	outDens.close();
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