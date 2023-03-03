//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/ddqm_model.h"
#include "../../include/interpolator.h"


#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <functional> 
#include <fstream>
// #include <gsl/gsl_interp.h>	
// #include <gsl/gsl_spline.h>

using namespace std;
int main(){

	double rhoB, temperature;
  double rhoBMax=1.5*pow(hc/Mnucleon, 3);
	double rhoBMin=0.05*pow(hc/Mnucleon, 3);
  int iR=200;
  double dRho= (rhoBMax - rhoBMin)/iR;

	cout << "Specify the temperature (MeV): " << endl;
	cin >> temperature ;
	cout << "You chose: " << endl
				<< "T = " << temperature << " MeV" <<  endl;

	temperature*=1./Mnucleon;
	quarks_class quarks;
	quarks.setFlavorNumber(3);
	quarks_class twoFlv;
	twoFlv.setFlavorNumber(2);

	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff=Me/Mnucleon;
	electron.Q=-1.;
	electron.temperature=temperature;
	electron.gamma=2.;

	particle muon;
	muon.mass= Mm/Mnucleon;
	muon.mass_eff=Mm/Mnucleon;
	muon.Q=-1.;
	muon.temperature=temperature;
	muon.gamma=2.;

// 	//Define thermodynamic variables
	double Energy, Pressure, FreeEn;
	string outStrS= "data/stability_ms"+to_string(quarks.qs.mass*Mnucleon)+
																	"_T"+std::to_string(temperature*Mnucleon)+".txt";
	ofstream  outStability(outStrS);
	//Set output file

//to do a single point, i*M=1 and set *Max
	double tcrit=170./Mnucleon;

	double Cmin=-1;
	// double Cmax=0.2;
	double Cmax=0.55;
	int iCM= 1;
	double dC= (Cmax-Cmin)/iCM;
	double Dmin=pow(120., 2.)/pow(Mnucleon, 2);
	// double Dmax=pow(195., 2.)/pow(Mnucleon, 2);
	double Dmax=pow(135.	, 2.)/pow(Mnucleon, 2);
	int iDM= 1;
	double dD= (Dmax-Dmin)/iDM;

	for(int ic=0; ic< iCM; ic++){
	for(int id=0; id< iDM; id++){
	
		double C= Cmax-ic*dC;
		double D= Dmax-id*dD;
		quarks.setParameters(C, D, tcrit);
		twoFlv.setParameters(C, D, tcrit);
		//Set output file		
		string outBeta= "data/eos_"+to_string(C)+"_"+to_string(sqrt(D)*Mnucleon)+
			 											 "_T"+to_string(temperature*Mnucleon)+".txt";
		ofstream outFile(outBeta);

		Eigen::VectorXd rhobv(iR), Pressv(iR), pressAbsv(iR), Ener3Fv(iR), Ener2Fv(iR), MupEfv(iR);
			//=== Loop over barionic density		
		for(int irho=0; irho<iR; irho++){
			rhoB=(rhoBMax- (double)irho*dRho);

			quarks.setEOS_betaEq(rhoB, temperature, electron, muon);

			particle electron_2f= electron;
			particle muon_2f= 		muon;
			twoFlv.setEOS_betaEq(rhoB, temperature, electron_2f, muon_2f);

			Pressure= quarks.getPressure() 			+ electron.pressure + muon.pressure;
			Energy 	= quarks.getEnergy()		 		+ electron.energy + muon.energy;
			FreeEn 	= quarks.getFenergy()				+ electron.energy - temperature*electron.entropy 
																					+ muon.energy - temperature*muon.entropy;

			rhobv(irho) = rhoB;
			Pressv(irho)= Pressure;
			pressAbsv(irho) = fabs(Pressure);
			if(temperature<Tmin_integration){
				Ener3Fv(irho)= Energy/rhoB;
				Ener2Fv(irho)	= twoFlv.getEnergy()/rhoB;
			}else{
				Ener3Fv(irho)= FreeEn/rhoB;
				Ener2Fv(irho)	= twoFlv.getFenergy()/rhoB;
			}
			MupEfv(irho)= quarks.qu.mass_eff;

			if(quarks.qu.mass_eff>0. && quarks.qd.mass_eff>0. && quarks.qs.mass_eff>0. ){
				outFile << rhoB*pow(Mnucleon/hc, 3.) << " "
									<< quarks.qu.density *pow(Mnucleon/hc, 3.) << " "
									<< quarks.qd.density *pow(Mnucleon/hc, 3.) << " "
									<< quarks.qs.density *pow(Mnucleon/hc, 3.) << " "
									<< electron.density	 *pow(Mnucleon/hc, 3.) << " " 
									<< muon.density			 *pow(Mnucleon/hc, 3.)  << " " 
									<< Energy*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
									<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
									<< quarks.muB*Mnucleon << " " 
									<< Ener3Fv(irho)*Mnucleon << " " 
									<< Ener2Fv(irho)*Mnucleon << " " 
									<< endl;
			}
		}
			
	Eigen::MatrixXd::Index minIE3F, minIE2F, minIMup, minIP;
				
	double minE3F= Ener3Fv.minCoeff(&minIE3F);
	double minE2F= Ener2Fv.minCoeff(&minIE2F);
	double minMup= MupEfv.minCoeff(&minIMup);
	double minPre= pressAbsv.minCoeff(&minIP);
	(void) minPre;
	double ePmin= Ener3Fv(minIP);
	cout << minE3F*Mnucleon << " " << minE2F*Mnucleon << endl;


	int iStable=0;
	if(minMup<0.){
		iStable=5;
	}else{
		if(minE2F*Mnucleon< 930. ){
			//stable 2F matter -- inexistent!
			iStable=1;
		}else{
			if(minE3F*Mnucleon<= 930.){
				iStable=2; //stable SQM
			}else if(minE3F*Mnucleon<= 939.){
				iStable=3; //metastable SQM
			}else{
				iStable=4; //unstable SQM
			}		
		}
	}

	cout << 	"iS=1 : Stable 2f,    2: Stable SQM,   3: metastable SQM,    4: unstable SQM" << endl
				<< "iS= " << to_string(iStable) <<  " for C= " << to_string(C) 
																				<< ", D= " << to_string(sqrt(D)*Mnucleon) <<endl;

	outStability << C << " " << sqrt(D)*Mnucleon << " " << iStable << " " 
							 <<  minE2F*Mnucleon << " " << minE3F*Mnucleon  
							 << " " << ePmin*Mnucleon*pow(Mnucleon/hc, 3.) << " " << rhobv(minIP)*pow(Mnucleon/hc, 3.)
							 << endl;
	outFile.close();
	}
	}

	outStability.close();

	return 0;
}