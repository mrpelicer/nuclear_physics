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
#include <vector>

using namespace std;

//===========================main codea =============================
int main(int argc, char** argv)
{
 
 string parametrization;

//Set system variables
	double rhoB, Yp, temperature;

	cout << "Specify the parametrization, proton fraction and temperature: " << endl;
	cin >> parametrization >> Yp >> temperature;
	cout << "You chose: " 
				<< parametrization  << " parametrization " << endl 
				<< "Yp= "<< Yp << endl  
				<< "T = " << temperature << endl;

	ofstream outFile("data/"+parametrization+"_Yp"+to_string(Yp)+".txt");		

// 	//Construct a gas of p,n,e matter
	nlwm_class qhd(parametrization);
 	//qhd.printParameters();
  particle electron;
	electron.mass_eff= Me/Mnucleon;
	electron.mass= electron.mass_eff;
	electron.Q=-1.;

	double rhoBMax=1./pow(Mnucleon/hc, 3);
	double rhoBMin=(1e-3)/pow(Mnucleon/hc, 3);
	int iR=1000;
	double dRho= (rhoBMax-rhoBMin)/iR;

// 	//Define thermodynamic variables
	double Energy, FreeEn, Pressure, Entropy;
	
	//Adimensional temperature;
	if(temperature>0.)	temperature*=1./Mnucleon;
		
		//=== Loop over barionic density
		for(int irho=0; irho<=iR; irho++){
			rhoB=(rhoBMax- (double)irho*dRho);

			qhd.setEOS_nucleons(rhoB, Yp, temperature);

			electron.density=Yp*rhoB;
			electron.temperature=temperature;
			electron.kf=pow(3.*pi2*electron.density, 1./3.);
			electron.solveChemPotEff();
			electron.chemPot=electron.chemPot_eff;
			electron.calculateProperties();
			 cout << qhd.proton.chemPot << " " << qhd.neutron.chemPot << endl;
			Energy= qhd.getEnergy()		 	;//+electron.energy; 		
			Pressure= qhd.getPressure()	;//+electron.pressure;
			Entropy= qhd.getEntropy()	 	;//+electron.entropy;
			FreeEn= Energy -temperature*Entropy;
			
			outFile << rhoB*pow(Mnucleon/hc, 3) << " " ///qhd.rho0 << " " // 
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3) << " " 
					<< (FreeEn/rhoB - 1.)*Mnucleon  << " " 
					<< Mnucleon*qhd.proton.chemPot  << " "
					<< Mnucleon*qhd.neutron.chemPot << " " 
					<< qhd.proton.mass_eff << " " 
					<< qhd.neutron.mass_eff << " " 
					<< qhd.phi0 << " " 
					<< qhd.V0 << " " 
					<< qhd.b0 << " " 
					<< qhd.muB << " " 
					<< qhd.muQ << " "
					<< qhd.proton.chemPot_eff  << " "
					<< qhd.neutron.chemPot_eff << " " 
					<< endl;

		}
		
	outFile.close();

  return 0;
 }

