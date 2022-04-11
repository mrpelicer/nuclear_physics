//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
#include "../../include/quark_model.h"
#include "../../include/quark_hadron_transition.h"
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
	std::string parametrization= "gm1";
	std::string hyperon_params ="gm"; 

	nlwm_class hadrons(parametrization);

	bool doHyperons	=	true;
	bool doDeltas		=	false;

	if(doHyperons) hadrons.includeHyperons(doHyperons, hyperon_params);

	quarks_class quarks;

	double tcrit=170./Mnucleon;
	double C= 0.;
	double D= pow(165., 2.)/pow(Mnucleon, 2);
	quarks.setParameters(C, D, tcrit);

	// phasetransition_class qht;
	//Construct a gas of p,n,e matter

	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff=Me/Mnucleon;
	electron.Q=-1.;
	electron.gamma=2.;

	particle muon;
	muon.mass= Mm/Mnucleon;
	muon.mass_eff=Mm/Mnucleon;
	muon.Q=-1.;
	muon.gamma=2.;


//Set system variables
	double muB, temperature;
	// double muBMin=0.;
  // double muBMax=1700./Mnucleon;
  // int imuMax=100;
  // double dMu=  (muBMax-muBMin)/imuMax;
	double tempMin=0./Mnucleon;
  double tempMax=150./Mnucleon;
  int itMax=1;
  double dt=  (tempMax-tempMin)/itMax;
	
	double rhoB;
	double rhoBMax=1.0*pow(hc/Mnucleon, 3);
	double rhoBMin=0.02*pow(hc/Mnucleon, 3);
	int iR=200;
	double dRho= (rhoBMax - rhoBMin)/iR;
	//Define thermodynamic variables
	double Energy, FreeEn, Pressure;

	std::vector<double> enerv, enerdensv, pressv, rhobv;


		//=== Loop over barionic density
	// for(int imu=0; imu<imuMax; imu++){
		// muB= muBMax- imu*dMu;
	for(int it=0; it<itMax; it++){
 		temperature=tempMin+ it*dt;
		std::string filename1;
		filename1="data/diagram_H_" +parametrization+
													"_Q_"+to_string(C)+"_"+to_string(sqrt(D)*Mnucleon)+
			 											 "_T"+to_string(temperature*Mnucleon)+".txt";
	
		std::ofstream outFile(filename1);

		electron.temperature=temperature;
		muon.temperature=temperature;
		cout << temperature*Mnucleon << endl;

		for(int irho=0; irho<iR; irho++){
			// qht.solveQHTransition(temperature, hadrons, quarks, electron, muon);
			rhoB=(rhoBMax- (double)irho*dRho);

			hadrons.setEOS_betaEq(rhoB, temperature, electron, muon);

			double Yu= (2.*hadrons.proton.density + hadrons.neutron.density
								+hadrons.lambda0.density  + 2.*hadrons.sigmap.density+ hadrons.sigma0.density
								+hadrons.xi0.density)/(3.*rhoB); 

			double Yd= (hadrons.proton.density + 2.*hadrons.neutron.density
									+hadrons.lambda0.density  +  hadrons.sigma0.density + 2.*hadrons.sigmam.density
									+hadrons.xim.density)/(3.*rhoB); 							

			double Ys= (hadrons.lambda0.density  
									+  hadrons.sigmap.density + hadrons.sigma0.density
									+hadrons.sigmam.density +2.*hadrons.xi0.density
									+2.*hadrons.xim.density)/(3.*rhoB); 							


			quarks.setEoSFlavorFixed(rhoB, hadrons.temperature,
																Yu, Yd, Ys);

				cout << rhoB*pow(Mnucleon/hc, 3) << " " 
				<< temperature*Mnucleon << " "
				<< hadrons.muB*Mnucleon << " "<< quarks.muB*Mnucleon << " " 
				<< hadrons.getPressure()*Mnucleon*pow(Mnucleon/hc, 3) << " " << quarks.getPressure()*Mnucleon*pow(Mnucleon/hc, 3)
				<<std::endl;

				outFile << rhoB*pow(Mnucleon/hc, 3) << " " 
				<< temperature*Mnucleon << " "
				<< hadrons.muB*Mnucleon << " "<< quarks.muB*Mnucleon << " " 
				<< hadrons.getPressure()*Mnucleon*pow(Mnucleon/hc, 3) << " " << quarks.getPressure()*Mnucleon*pow(Mnucleon/hc, 3)
				<<std::endl;


		}
			outFile.close();

	}
		


  return 0;
}