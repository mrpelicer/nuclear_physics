//Nuclear matter properties using Mean Feild Theory (NLWM) at T=0.
#include "../../include/constant.hpp"
#include "../../include/particles.hpp"
#include "../../include/bag_model.hpp"
#include "../../include/interpolator.hpp"

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
	double rhoBMin=0.*pow(hc/Mnucleon, 3);
  int iR=150;
  double dRho= (rhoBMax - rhoBMin)/iR;

	cout << "Specify the temperature (MeV): " << endl;
	cin >> temperature ;
	cout << "You chose: " << endl
				<< "T = " << temperature << " MeV" <<  endl;

	temperature*=1./Mnucleon;
	bag_model_class quarks;
	bag_model_class twoFlv;

// 	//Define thermodynamic variables
	double Energy, Pressure, FreeEn;


	double Bag=pow(148./Mnucleon, 4);	//MeV
	double Xv=0; 											//adim
	double Gv=0.*pow(Mnucleon/hc, 2); 	// fm^2
	double xsi=0.;
	double tcrit=0/Mnucleon;					//MeV
	
	//double Bmin=pow(145./Mnucleon, 4);
	//double Bmax=pow(154./Mnucleon, 4);
	//int iBM= 200;
	//double dB= (Bmax-Bmin)/iBM;

	string outStrS= "data/stability_bag"+to_string(Bag*pow(Mnucleon, 4))+
																	"_T"+std::to_string(temperature*Mnucleon)+".txt";
	ofstream  outStability(outStrS);
	//Set output file

//to do a single point, i*M=1 and set *Max

	string outSym= "data/eos_sym_T"+to_string(temperature*Mnucleon)+".txt";
	ofstream outFile(outSym);

	//for(int ib=0; ib< iBM; ib++){
	
		//Bag= Bmax-iB*dB;

		quarks.setParameters(Bag, Gv, xsi, Xv, tcrit);
		quarks.setFlavorNumber(3);
		twoFlv.setParameters(Bag, Gv, xsi, Xv, tcrit);
		twoFlv.setFlavorNumber(2);
		//Set output file		
		Eigen::VectorXd rhobv(iR), Pressv(iR), pressAbsv(iR), Ener3Fv(iR), Ener2Fv(iR), MupEfv(iR);
			
			//=== Loop over barionic density		
		for(int irho=0; irho<iR; irho++){
			rhoB=(rhoBMax- (double)irho*dRho);

			quarks.setEOS_symmetric(rhoB, temperature);
			twoFlv.setEOS_symmetric(rhoB, temperature);

			Pressure= quarks.getPressure();
			Energy 	= quarks.getEnergy();
			FreeEn 	= quarks.getFenergy();
																		

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

				outFile << rhoB*pow(Mnucleon/hc, 3.) << " "
									<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
									<< Energy*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
									<< quarks.muB*Mnucleon << " " 
									<< quarks.qu.density *pow(Mnucleon/hc, 3.) << " "
									<< quarks.qd.density *pow(Mnucleon/hc, 3.) << " "
									<< quarks.qs.density *pow(Mnucleon/hc, 3.) << " "
									<< 0 << " " 
									<< 0  << " " 
									<< Ener3Fv(irho)*Mnucleon << " " 
									<< Ener2Fv(irho)*Mnucleon << " " 
									<< endl;
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
				<< "iS= " << to_string(iStable) <<  " for B^1/4= " << to_string(pow(Bag, 1/4.)*Mnucleon) 
																				<< ", Gv= " << to_string(Gv)
																				<< ", Xv= " << to_string(Xv) <<endl;

	outStability << iStable << " " 
							 <<  minE2F*Mnucleon << " " << minE3F*Mnucleon  
							 << " " << ePmin*Mnucleon*pow(Mnucleon/hc, 3.) << " " << rhobv(minIP)*pow(Mnucleon/hc, 3.)
							 << endl;
	outFile.close();

	outStability.close();

	return 0;
}