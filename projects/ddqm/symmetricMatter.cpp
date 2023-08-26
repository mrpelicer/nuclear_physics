//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.hpp"
#include "../../include/particles.hpp"
#include "../../include/ddqm_model.hpp"
#include "../../include/interpolator.hpp"

#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <functional> 
#include <fstream>
// #include <gsl/gsl_interp.h>	
// #include <gsl/gsl_spline.h>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

using namespace std;
int main(){

// //Set system variables
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

	
	//Define thermodynamic variables
	double Energy, Pressure, FreeEn;
	std::string outStrS= "data/stability_ms"+std::to_string(quarks.qs.mass*Mnucleon)+
											"_T"+std::to_string(temperature*Mnucleon)+".txt";
	std::ofstream  outStability(outStrS);


//to do a single point, iC/iD=1 and set CMax/DMax to the desired value.
	double tcrit=170./Mnucleon;
	double Cmin=-0.8;
	double Cmax=0.2;
	int iCM= 1;
	double dC= (Cmax-Cmin)/iCM;
	double Dmin=pow(150., 2.)/pow(Mnucleon, 2);
	double Dmax=pow(195., 2.)/pow(Mnucleon, 2);
	int iDM= 1;
	double dD= (Dmax-Dmin)/iDM;

	for(int ic=0; ic< iCM; ic++){
	for(int id=0; id< iDM; id++){
		
		double C= Cmax-ic*dC;
		double D= Dmax-id*dD;
		//Set output file		
		string outSym= "data/eos_"+std::to_string(C)+"_"+std::to_string(sqrt(D)*Mnucleon)+
			 											 "_T"+std::to_string(temperature*Mnucleon)+".txt";
		//string outSym= "symEq.txt";
		ofstream outFile(outSym);

		Eigen::VectorXd rhobv(iR), Pressv(iR), pressAbsv(iR), Ener3Fv(iR), Ener2Fv(iR), MupEfv(iR);
				
		quarks.setParameters(C, D, tcrit);
		twoFlv.setParameters(C, D, tcrit);

		for(int irho=0; irho<iR; irho++){
			rhoB=(rhoBMax- (double)irho*dRho);

			quarks.setEOS_symmetric(rhoB, temperature);
			twoFlv.setEOS_symmetric(rhoB, temperature);
			
			Pressure= quarks.getPressure();
			Energy 	= quarks.getEnergy();
			FreeEn 	= quarks.getFenergy();

			rhobv(irho) = rhoB;
			Pressv(irho)= Pressure;
			pressAbsv(irho) = std::fabs(Pressure);
			if(temperature<Tmin_integration){
				Ener3Fv(irho)= Energy/rhoB;
				Ener2Fv(irho)	= twoFlv.getEnergy()/rhoB;
			}else{
				Ener3Fv(irho)= FreeEn/rhoB;
				Ener2Fv(irho)	= twoFlv.getFenergy()/rhoB;
			}
			MupEfv(irho)= quarks.qu.mass_eff;
			std::cout << rhoB*pow(Mnucleon/hc, 3.) << " " 
								<< Ener3Fv(irho)*Mnucleon << " " << Ener2Fv(irho)*Mnucleon << " " 
								<< twoFlv.getEnergy()*Mnucleon*pow(Mnucleon/hc, 3) << " " 
								<< Energy*Mnucleon*pow(Mnucleon/hc, 3) << " " 
								<< twoFlv.qu.mass_eff*Mnucleon<< " " << twoFlv.qd.mass_eff*Mnucleon << " "  
								<< twoFlv.qs.mass_eff*Mnucleon
								<< std::endl;

		if(quarks.qu.mass_eff>0. && quarks.qd.mass_eff>0. && quarks.qs.mass_eff>0. ){
			outFile << rhoB*pow(Mnucleon/hc, 3.) << " "
								<< quarks.qu.density*pow(Mnucleon/hc, 3.) << " "
								<< quarks.qd.density*pow(Mnucleon/hc, 3.) << " "
								<< quarks.qs.density*pow(Mnucleon/hc, 3.) << " "
								<< 0. << " " << 0. << " "
								<< Energy*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
								<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
								<< quarks.muB*Mnucleon << " " 
								<< Ener3Fv(irho)*Mnucleon << " " 
								<< Ener2Fv(irho)*Mnucleon << " " 
								<< std::endl;
		}	
		}
		
	Eigen::MatrixXd::Index minIE3F, minIE2F, minIMup, minIP;
		
	double minE3F= Ener3Fv.minCoeff(&minIE3F);
	double minE2F= Ener2Fv.minCoeff(&minIE2F);
	double minMup= MupEfv.minCoeff(&minIMup);
	double minPre= pressAbsv.minCoeff(&minIP);
	(void) minPre;
	double ePmin= Ener3Fv(minIP);
	std::cout << minE3F*Mnucleon << " " << minE2F*Mnucleon << std::endl;

	int iStable=0;
	if(minMup<0.){
		iStable=5;
	}else{
		if(minE2F*Mnucleon< 930. ){
			//stable 2F mattter
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
				<< "iS= " << to_string(iStable) << endl;

	outStability << C << " " << sqrt(D)*Mnucleon << " " << iStable << " " 
							 <<  minE2F*Mnucleon << " " << minE3F*Mnucleon  
							 << " " << ePmin*Mnucleon*pow(Mnucleon/hc, 3.) << " " << rhobv(minIP)*pow(Mnucleon/hc, 3.)
							 << std::endl;
	outFile.close();
	}
	}

	outStability.close();
	return 0;
}