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
#include <vector>

using namespace std;

//===========================main codea =============================
int main(){
// //Choose pararametrization
	string parametrization;
	
	cout << "Choose parametrization" << endl;
	cin >> parametrization;
// 	//Construct a gas of p,n,e matter
	nlwm_class qhd(parametrization);
// 	qhd.printParameters();
	qhd.printParameters();
//Set system variables
	double rhoB,Yle=0.5, temperature=0.;
  	double rhoBMax=1./pow(Mnucleon/hc, 3);
  	int iR=200;
  	double dRho= rhoBMax/iR;

// 	//Define thermodynamic variables
	double Energy, FreeEn, Pressure, Entropy;
	double Bind0En, Esym, Lsym, K0, cs2;
	// VectorXd ener(iR), esym(iR), enerDens(iR), press(iR);
	vector<double> enerv, enerdensv, esymv, pressv, rhobv;


	bool doHyperons	=	false;
	bool doDeltas		=	false;

	string hyperon_params;  //gm (Glendenning), su3 (código do Kauan) 
	cout << "Choose hyperon parametrization" << endl;
	cin >> hyperon_params;

 	qhd.includeHyperons(doHyperons, hyperon_params);
	string delta_params;  //su6(1.), mplA_1(beta=1.1), mplA_2, prd89_1, prd89_1
	cout << "Choose delta parametrization" << endl;
	cin >> delta_params;
	qhd.includeDeltas(		doDeltas, 	delta_params);

	vector<double> Unv, Upv, Ul0v, Usmv, Us0v, Uspv, Uxmv, Ux0v, Udmv, Ud0v, Udpv, Udppv;

	cout << "Yp, T= " << Yle << " " << temperature << endl;
		
	string outStr= "data/"+parametrization+".txt";
	ofstream outFile(outStr);
			
	outFile << "#rho_b P freeEn mun mup mue M* V0 b0" << endl;
	//Adimensional temperature;
	if(temperature>0.)	temperature*=1./Mnucleon;
		
		//=== Loop over barionic density
		for(int irho=0; irho<iR; irho++){
			rhoB=(rhoBMax- (double)irho*dRho);

			qhd.setEOS_nucleons(rhoB, Yle, temperature);
			
			auto UN= qhd.getNucleonPotential();
			auto UH= qhd.getHyperonPotential();
			auto UD= qhd.getDeltaPotential();
			Unv.push_back(UN[0]); 
			Upv.push_back(UN[1]); 
			Ul0v.push_back(UH[0]); 
			Usmv.push_back(UH[1]);
			Us0v.push_back(UH[2]);
			Uspv.push_back(UH[3]);
			Uxmv.push_back(UH[4]); 
			Ux0v.push_back(UH[5]);
			Udmv.push_back(UD[0]);
			Ud0v.push_back(UD[1]); 
			Udpv.push_back(UD[2]); 
			Udppv.push_back(UD[3]);

			Energy= qhd.getEnergy()		;//+electron.energy; 		
			Pressure= qhd.getPressure();//+electron.pressure;
			Entropy= qhd.getEntropy()	;//+electron.entropy;
			FreeEn= Energy -temperature*Entropy;
			
			rhobv.push_back(rhoB);
			enerv.push_back(FreeEn/rhoB);
			enerdensv.push_back(FreeEn);
			pressv.push_back(Pressure);

			outFile << rhoB*pow(Mnucleon/hc, 3) << " "
					<< Pressure*Mnucleon*pow(Mnucleon/hc, 3) << " " 
					<< (FreeEn/rhoB - 1.)*Mnucleon  << " " 
					<< qhd.rhoS*pow(Mnucleon/hc, 3)  << " " 
					<< qhd.proton.chemPot  << " "
					<< qhd.neutron.chemPot << " " 
					<< qhd.proton.mass_eff << " " 
					<< qhd.neutron.mass_eff << " " 
					<< qhd.phi0 << " " 
					<< qhd.V0 << " " 
					<< qhd.b0 << " " 
					<< qhd.muB << " " 
					<< qhd.muQ << " " 
					<< qhd.proton.chemPot_eff  << " " << qhd.neutron.chemPot_eff + qhd.getRearrangementEnergy() << " " 
					<< qhd.getRearrangementEnergy() << " " << qhd.getDerivativeCoupling_sigma(rhoB) << " "
				 << qhd.getDerivativeCoupling_omega(rhoB) << " " << qhd.getDerivativeCoupling_rho(rhoB) 
					<< endl;


			double h=1e-2;
			double e0, em1, em2, ep1, ep2;

			qhd.setEOS_nucleons(rhoB, Yle, temperature);
			e0= qhd.getEnergy()/rhoB;

			qhd.setEOS_nucleons(rhoB, Yle-h, temperature);
			em1= qhd.getEnergy()/rhoB;

			qhd.setEOS_nucleons(rhoB, Yle-2.*h, temperature);
			em2= qhd.getEnergy()/rhoB;

			qhd.setEOS_nucleons(rhoB, Yle+h, temperature);
			ep1= qhd.getEnergy()/rhoB;

			qhd.setEOS_nucleons(rhoB, Yle+2.*h, temperature);
			ep2= qhd.getEnergy()/rhoB;

			esymv.push_back( (-ep2+ 16.*ep1 -30.*e0 + 16.*em1 - em2)/(8.*12.*h*h) );	


		}
		
	outFile.close();


  reverse(rhobv.begin(), 		rhobv.end());
	reverse(enerv.begin(), 		enerv.end());
	reverse(enerdensv.begin(), enerdensv.end());
	reverse(pressv.begin(), 		pressv.end());
	reverse(esymv.begin(), 		esymv.end());
	reverse(Unv.begin(),		Unv.end());
	reverse(Upv.begin(), 	Upv.end());
	reverse(Ul0v.begin(), 	Ul0v.end());
	reverse(Usmv.begin(), 	Usmv.end());
	reverse(Us0v.begin(), 	Us0v.end());
	reverse(Uspv.begin(), 	Uspv.end());
	reverse(Uxmv.begin(), 	Uxmv.end());
	reverse(Ux0v.begin(), 	Ux0v.end());
	reverse(Udmv.begin(), 	Udmv.end());
	reverse(Ud0v.begin(), 	Ud0v.end());
	reverse(Udpv.begin(), 	Udpv.end());
	reverse(Udppv.begin(),	Udppv.end()); 

	vector<double> enerd1v, enerd2v,enerd3v, esymd1v, esymd2v, esymd3v; 
	for(int irho=0; irho<iR; irho++){
		rhoB=	rhobv[0]+irho*dRho;
		esymd1v.push_back(deriv_func(rhoB, esymv, rhobv));
		// esymd2v.push_back(deriv2_func(rhoB, esymv, rhobv));
		enerd1v.push_back(deriv_func(rhoB, enerv, rhobv));
		// enerd2v.push_back(deriv_func(rhoB, enerv, rhobv));
	}
	for(int irho=0; irho<iR; irho++){
		rhoB=	rhobv[0]+irho*dRho;
		// lsymv.push_back(3.*rhoB*esymd1v[irho]);		
		enerd2v.push_back(deriv_func(rhoB, enerd1v, rhobv));
		esymd2v.push_back(deriv_func(rhoB, esymd1v, rhobv));
	}
	for(int irho=0; irho<iR; irho++){
		rhoB=	rhobv[0]+irho*dRho;
		// lsymv.push_back(3.*rhoB*esymd1v[irho]);		
		enerd3v.push_back(deriv_func(rhoB, enerd2v, rhobv));
		esymd3v.push_back(deriv_func(rhoB, esymd2v, rhobv));
	}
	vector<double> lsymv, ksymv, k0v, Kv, q0v, cs2v;
	for(int irho=0; irho<iR; irho++){
		rhoB=	rhobv[0]+irho*dRho;
		lsymv.push_back(3.*rhoB*esymd1v[irho]);
		ksymv.push_back(pow(3.*rhoB, 2.)*esymd2v[irho]);
		k0v.push_back(pow(3.*rhoB, 2.)*enerd2v[irho]);
		Kv.push_back(9.*deriv_func(rhoB, pressv, rhobv));
		q0v.push_back(pow(3.*rhoB, 3.)*enerd3v[irho]);
		cs2v.push_back(deriv_func(rhoB, pressv, rhobv)/deriv_func(rhoB, enerdensv, rhobv));
	}

	ofstream outBulk("data/bulk_"+parametrization+"_sym.txt");
	ofstream outUpar("data/U_delta_sym_"+parametrization+".txt");
	for(int irho=0; irho<iR; irho++){
		rhoB=	rhobv[0]+irho*dRho;

		outBulk	<<  rhoB*pow(Mnucleon/hc, 3) << " " 
				<< interpolation_func(rhoB, pressv, rhobv)*Mnucleon*pow(Mnucleon/hc, 3) << " "
				<< interpolation_func(rhoB, enerdensv, rhobv)*Mnucleon*pow(Mnucleon/hc, 3) << " "
				<< (interpolation_func(rhoB, enerv, rhobv)-1.)*Mnucleon << " "
			 	<< interpolation_func(rhoB, Kv, rhobv)*Mnucleon << " "
				<< interpolation_func(rhoB, k0v, rhobv)*Mnucleon << " "
				<< fabs(interpolation_func(rhoB, cs2v, rhobv)) << " "
				<< interpolation_func(rhoB, esymv, rhobv)*Mnucleon << " "
				<< interpolation_func(rhoB, lsymv, rhobv)*Mnucleon << " "
				<< endl;

		outUpar << rhoB*pow(Mnucleon/hc, 3) << " " 
					<< interpolation_func(rhoB, Unv, rhobv)*Mnucleon		<< " "
					<< 	interpolation_func(rhoB, Ud0v, rhobv)*Mnucleon 	
					<< endl;
					// << 	interpolation_func(rhoB, Ul0v, rhobv)*Mnucleon 	<< " "
					// << 	interpolation_func(rhoB, Us0v, rhobv)*Mnucleon 	<< " "
					// << 	interpolation_func(rhoB, Ux0v, rhobv)*Mnucleon 	<< " " 
		}

	outBulk.close();
	outUpar.close();

// 
// //Print bulk properties at saturation: 
// //assures code is working properly for infinite symmetric matter.

	qhd.setEOS_nucleons(qhd.rho0, 0.5, 0.0);
			  
	double rho0= qhd.rho0;
	Bind0En	= interpolation_func(rho0, enerv, rhobv) -1.;
	Esym		= interpolation_func(rho0, esymv, rhobv);
	Lsym		= interpolation_func(rho0, lsymv, rhobv); //3.*rho0*deriv_func(rho0, esymv, rhobv);
	K0			= interpolation_func(rho0, k0v, 	rhobv); //pow(3.*rho0, 2.)*deriv2_func(rho0, enerv, rhobv);
	double K 	= interpolation_func(rho0, Kv, 	rhobv);
	cs2			= interpolation_func(rho0, cs2v, 	rhobv); // deriv_func(rho0, pressv, rhobv)/deriv_func(rho0, enerdensv, rhobv);
	double Un= 	interpolation_func(rho0, Unv, rhobv);
	double Ul0= 	interpolation_func(rho0, Ul0v, rhobv);
	double Us0= 	interpolation_func(rho0, Us0v, rhobv);
	double Ux0= 	interpolation_func(rho0, Ux0v, rhobv);
	double Ud0= 	interpolation_func(rho0, Ud0v, rhobv);

	
	cout << " rho_0" << setw(25)<< "B/A (MeV)" << setw(25) <<"K_0 (MeV)" 
					  <<setw(25) << "J(MeV)" << setw(25) <<  "L(MeV)" << setw(25) <<"m*/m" << setw(25) << "cs2"
		  			<< endl;	

	cout <<rho0*pow(Mnucleon/hc, 3) << setw(25)
			 << Bind0En*Mnucleon << setw(25)
			 << K0*Mnucleon << setw(25)
			 << Esym*Mnucleon << setw(25)
			 << Lsym*Mnucleon << setw(25)
 			 << qhd.proton.mass_eff << setw(25) 
			 << fabs(cs2) << setw(25) 
			 << endl;

	cout 	<< "Un"								<< setw(25)
							<< "Ul0" 							<< setw(25)
							<< "Us0" 							<< setw(25)
							<< "Ux0" 							<< setw(25) 
							<< "Ud0" 							<< setw(25)<< endl 
							<< Un*Mnucleon		<< setw(25) 
							<< Ul0*Mnucleon 	<< setw(25)
							<< Us0*Mnucleon 	<< setw(25)
							<< Ux0*Mnucleon 	<< setw(25) 
							<< Ud0*Mnucleon 	<< setw(25)
							<< endl;

  return 0;
 }

