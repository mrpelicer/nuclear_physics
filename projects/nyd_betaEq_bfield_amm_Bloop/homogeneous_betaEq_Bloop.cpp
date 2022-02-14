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


//===========================main codea =============================

int main(){
//Choose pararametrization
	std::string parametrization= "l3wr";

//Construct a gas of p,n,e matter
	nlwm_class hmg_matter(parametrization);
	// nlwm_class derivatives(parametrization);
	particle electron, muon;
	electron.mass	= Me/hmg_matter.Mn;	
  muon.mass			=	Mm/hmg_matter.Mn; 
	electron.mass_eff= electron.mass;
	muon.mass_eff= muon.mass;	
	electron.Q=-1.;
	muon.Q=-1.;
	
	double rhoB=0.7/pow(hmg_matter.Mn/hc, 3);
	double temperature, Bfield;

//Set Bfield variables
	double BgMin=1.e17;
	double lBgMin=log10(BgMin);
	double BgMax=3.e18; // Gauss
	double lBgMax=log10(BgMax);
	double Bc=pow(electron.mass_eff, 2.)/eHL; 
	double Bg; 
	//Bfield *= 1.95e-14/pow(hmg_matter.Mn, 2.);  // Conversion factor Gauss to MeV^2 to adim
	
	int ibMax=1000;
	std::vector<double> enerv, enerdensv, pressv, bfieldv, cs2v;

	double dlogB= (lBgMax-lBgMin)/ibMax;

	std::vector<double> temp_vec={0.};
	//Define thermodynamic variables
	double Energy, FreeEn, Pressure, Entropy;
	//double BindEn

	bool doHyperons	=	true;
	bool doDeltas		=	false;
	bool doBfield		= false; 
	bool doamm 			= false;
	if(doHyperons){ 
		//double as=1.;
		//double av=1.;
		std::string hyperon_params ="l3wr3";  //gm (Glendenning), su3 (código do Kauan) 
																					// ou fsu2h, fsu2h_1, fsu2h_2

		//hmg_matter.includeHyperons(doHyperons,	as,	av);
	 	hmg_matter.includeHyperons(doHyperons, hyperon_params);
		// derivatives.includeHyperons(doHyperons, hyperon_params);
	}

	if(doDeltas){
		std::string delta_params 	 ="su6";  //su6(1.), mplA_1(beta=1.1), mplA_2, prd89_1, prd89_1
		hmg_matter.includeDeltas(		doDeltas, 	delta_params);
		// derivatives.includeDeltas(		doDeltas, 	delta_params);
	}


	hmg_matter.printParameters();
	double yN=0., yH=0., yD=0.;

	temperature=temp_vec[0];
		 		
		std::string filename1, filename2, filename3, filename4;
		if(!doBfield){
			filename1="data/dens_beta_" 	+parametrization+"_noB.txt";
			filename2="data/hmg_beta_" 	+parametrization+"_noB.txt";
			filename3="data/spin_dens_"	+parametrization+"_noB.txt";
			filename4="data/fraction_"		+parametrization+"_noB.txt";
		}else if(doBfield && !doamm){
			filename1="data/dens_beta_" 	+parametrization+"_wtB.txt";
			filename2="data/hmg_beta_"		+parametrization+"_wtB.txt";
			filename3="data/spin_dens_" 	+parametrization+"_wtB.txt";
			filename4="data/fraction_"		+parametrization+"_wtB.txt";
		}else{
			filename1="data/dens_beta_" 	+parametrization+"_wtA.txt";
			filename2="data/hmg_beta_"		+parametrization+"_wtA.txt";
			filename3="data/spin_dens_"	+parametrization+"_wtA.txt";
			filename4="data/fraction_"		+parametrization+"_wtA.txt";
		}
//../data/data_hmg/
		
		std::ofstream outDens(filename1);
		std::ofstream outFile(filename2);
		std::ofstream outSpin(filename3);
		std::ofstream outFrac(filename4);
		//Adimensional temperature;
		temperature*=1./hmg_matter.Mn;
		electron.temperature=temperature;
		muon.temperature=temperature;

		//=== Loop over magnetic field
		for(int ib=0; ib<ibMax; ib++){
			Bg= pow(10., lBgMax- ib*dlogB);
			Bfield=Bg*Bc/4.41e13;				


			if(doBfield){
				hmg_matter.setBfield(doBfield, Bfield);
				// derivatives.setBfield(doBfield, Bfield);
				electron.setBfield(doBfield, Bfield);
	 			muon.setBfield(doBfield, Bfield);
				if(doamm){
					hmg_matter.setAMM(doamm);
					// derivatives.setAMM(doamm);
					// electron.setAMM(0.);
					// muon.setAMM(0.);
					electron.setAMM(doamm, (1.15965e-3)/electron.mass);
					muon.setAMM(doamm, (1.16592e-3)/muon.mass);
				}
		}

		//Solve self-consistently:
			std::cout << "rhob= " << rhoB*pow(hmg_matter.Mn/hc, 3) << " , B= " << Bfield/(Bc/4.41e13) << " " << lBgMax- ib*dlogB << std::endl;

			hmg_matter.setEOS_betaEq(rhoB, temperature, electron, muon);

			electron.pressure	=electron.chemPot*electron.density 	- electron.energy;
			muon.pressure			=muon.chemPot*muon.density 					- muon.energy;
			//Set variables:
			Energy= hmg_matter.getEnergy() 		+ electron.energy+ muon.energy;
			Pressure= hmg_matter.getPressure()+ electron.pressure+ muon.pressure;
			Entropy= hmg_matter.getEntropy() 	+ electron.entropy + muon.entropy;
			FreeEn= Energy -temperature*Entropy;

			// particle electron_tmp= electron;
			// particle muon_tmp= muon;
			// double hdev=1e-4;
			//std::cout << "TESTAO: " << rhoB << " " << rhoB+2.*hdev << std::endl;
			// derivatives.setEOS_betaEq(rhoB+2.*hdev, temperature, electron_tmp, muon_tmp);
			// electron_tmp.pressure	=electron_tmp.chemPot*electron_tmp.density 	- electron_tmp.energy;
			// muon_tmp.pressure			=muon_tmp.chemPot*muon_tmp.density 					- muon_tmp.energy;
			// double enerPP=	derivatives.getEnergy() 		+ electron_tmp.energy+ muon_tmp.energy;
			// double PressPP= derivatives.getPressure() 		+ electron_tmp.pressure+ muon_tmp.pressure;
// 
			// derivatives.setEOS_betaEq(rhoB+hdev, temperature, electron_tmp, muon_tmp);
			// electron_tmp.pressure	=electron_tmp.chemPot*electron_tmp.density 	- electron_tmp.energy;
			// muon_tmp.pressure			=muon_tmp.chemPot*muon_tmp.density 					- muon_tmp.energy;
			// double enerP=	derivatives.getEnergy() 		+ electron_tmp.energy+ muon_tmp.energy;
			// double PressP= derivatives.getPressure() 		+ electron_tmp.pressure+ muon_tmp.pressure;
// 
			// derivatives.setEOS_betaEq(rhoB-hdev, temperature, electron_tmp, muon_tmp);
			// electron_tmp.pressure	=electron_tmp.chemPot*electron_tmp.density 	- electron_tmp.energy;
			// muon_tmp.pressure			=muon_tmp.chemPot*muon_tmp.density 					- muon_tmp.energy;
			// double enerM=	derivatives.getEnergy() 		+ electron_tmp.energy+ muon_tmp.energy;
			// double PressM= derivatives.getPressure() 		+ electron_tmp.pressure+ muon_tmp.pressure;
// 
			// derivatives.setEOS_betaEq(rhoB-2.*hdev, temperature, electron_tmp, muon_tmp);
			// electron_tmp.pressure	=electron_tmp.chemPot*electron_tmp.density 	- electron_tmp.energy;
			// muon_tmp.pressure			=muon_tmp.chemPot*muon_tmp.density 					- muon_tmp.energy;
			// double enerMM=	derivatives.getEnergy() 		+ electron_tmp.energy+ muon_tmp.energy;
			// double PressMM= derivatives.getPressure() 		+ electron_tmp.pressure+ muon_tmp.pressure;
// 
			// double dpdr=(PressMM - 8.*PressM +8.*PressP -PressPP)/(12.*hdev);
			// double dedr=(enerMM - 8.*enerM +8.*enerP -enerPP)/(12.*hdev);
// 
			// double cs2= dpdr/dedr;
			bfieldv.push_back(Bfield);
			pressv.push_back(Pressure);
			enerdensv.push_back(Energy);
			enerv.push_back(Energy/rhoB);
			// cs2v.push_back(dpdr/dedr);

			yN= (hmg_matter.proton.density + hmg_matter.neutron.density)/rhoB;
			yH= (hmg_matter.lambda0.density + hmg_matter.sigmam.density+ hmg_matter.sigma0.density+ hmg_matter.sigmap.density
						+ hmg_matter.xi0.density+ hmg_matter.xim.density)/rhoB;
			yD= (hmg_matter.deltapp.density + hmg_matter.deltap.density+ hmg_matter.delta0.density+ hmg_matter.deltam.density)/rhoB;
		  
			outFile << Bg << " " // hmg_matter.rho0 or *pow(hmg_matter.Mn/hc, 3) 
					<< Pressure*hmg_matter.Mn*pow(hmg_matter.Mn/hc, 3) << " " 
					<< (FreeEn/rhoB - 1.)*hmg_matter.Mn  << " " 
					<< Energy*hmg_matter.Mn*pow(hmg_matter.Mn/hc, 3)  << " " <<  hmg_matter.neutron.mass_eff*hmg_matter.Mn << " "
					<< hmg_matter.neutron.chemPot*hmg_matter.Mn 			<< " " << electron.chemPot*hmg_matter.Mn << " " 
					<< yN << " " << yH << " " << yD
					<<std::endl;

			outSpin << Bg << " " // hmg_matter.rho0 or *pow(hmg_matter.Mn/hc, 3) 
					<<	hmg_matter.proton.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.proton.densityM*pow(hmg_matter.Mn/hc, 3) << " "
					<<	electron.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	electron.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	muon.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	muon.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.sigmap.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.sigmap.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.sigmam.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.sigmam.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.xim.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.xim.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityPP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltapp.densityMM*pow(hmg_matter.Mn/hc, 3) << " "
					<<	hmg_matter.deltap.densityPP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltap.densityP*pow(hmg_matter.Mn/hc, 3) << " "
					<<	hmg_matter.deltap.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltap.densityMM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltam.densityPP*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltam.densityP*pow(hmg_matter.Mn/hc, 3) << " "
					<<	hmg_matter.deltam.densityM*pow(hmg_matter.Mn/hc, 3) << " " 
					<<	hmg_matter.deltam.densityMM*pow(hmg_matter.Mn/hc, 3)
					<<  std::endl;

			outDens << Bg									<< " " // /hmg_matter.rho0 			 << " " // *pow(hmg_matter.Mn/hc, 3) << " " 
				<< hmg_matter.proton.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.neutron.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< electron.density						*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< muon.density								*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.lambda0.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.sigmap.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.sigma0.density	*pow(hmg_matter.Mn/hc, 3)	 << " "		//	/rhoB
				<< hmg_matter.sigmam.density	*pow(hmg_matter.Mn/hc, 3)	 << " "		//	/rhoB
				<< hmg_matter.xi0.density			*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.xim.density			*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.deltapp.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.deltap.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.delta0.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< hmg_matter.deltam.density	*pow(hmg_matter.Mn/hc, 3)  << " "		//	/rhoB
				<< std::endl;
		
			outFrac << Bg		<< " "
				<< yN << " " << yH << " " << yD 
				<< std::endl;


		}

		std::reverse(bfieldv.begin(), bfieldv.end());
		std::reverse(pressv.begin(), pressv.end());
		std::reverse(enerdensv.begin(), enerdensv.end());
		std::reverse(enerv.begin(), enerv.end());

		
		
	outFile.close();
	outSpin.close();
	outDens.close();


	std::string outMagnetization;
	if(!doBfield){
		outMagnetization="magnetization_"+parametrization+"_noB.txt";
	}else if(doBfield && !doamm){
		outMagnetization="magnetization_"+parametrization+"_wtB.txt";
	}else{
		outMagnetization="magnetization_"+parametrization+"_wtA.txt";
	}
	std::ofstream outMag(outMagnetization);

	//=== Loop over magnetic field
		for(int ib=0; ib<ibMax; ib++){
			Bg= pow(10., lBgMax- ib*dlogB);
			Bfield=Bg*Bc/4.41e13;				
		
		outMag	<< Bg  << " " 
							<< interpolation_func(Bfield, pressv, bfieldv)*hmg_matter.Mn*pow(hmg_matter.Mn/hc, 3) << " " 
							<< deriv_func(Bfield, pressv, bfieldv)*hmg_matter.Mn*pow(hmg_matter.Mn/hc, 3)*4.41e13/Bc
							<< std::endl; 

	}
	outMag.close();

  return 0;
}