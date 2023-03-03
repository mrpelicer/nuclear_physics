// Pasta with non linear Walecka model using Mean Field Theory
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
#include "../../include/pasta.h"
#include "../../include/interpolator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

//===========================main code =============================

int main(){

	std::string parametrization;
	double rhoB, Yp, temperature;

	cout << "Specify the parametrization and temperature: " << endl;
	cin >> parametrization >> temperature;
	cout << "You chose: " 
				<< parametrization  << " parametrization " << endl 
				<< "T = " << temperature << endl;

//Declare nuclear matter: pasta and gas
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);

	double rhoBMax= 0.093*pow(hc/Mnucleon, 3);
	// if(temperature>4.){rhoBMax=}
	int iRhoMax=100;
	double dRho= rhoBMax/iRhoMax;
		
	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= Me/Mnucleon;
	electron.spin=1./2.;
	electron.Q=-1.;


	//Pasta values:
	Eigen::MatrixXd PressureM(3,2), BulkEnM(3,2), GibbsEnM(3,2), FreeEnM(3,2), EnergyM(3,2),
					EntropyM(3,2), coulEnM(3,2), surfEnM(3,2), RdM(3,2), RwM(3,2), VcM(3,2), VwM(3,2),
					AeM(3, 2), ZeM(3, 2);
	
	double Pressure, BulkEn, FreeEn, GibbsEn, Entropy, coulEn, surfEn, f, Rd, Rw, sigma;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	

	double Ae=0., Ze=0.;
	
	std::cout << "Yp, T= " << Yp << " " << temperature << std::endl;
	std::ofstream outGlobal("data/cpa_"+parametrization+"_betaEq"+
												   +"_T"+std::to_string(temperature)+".txt");

	std::ofstream outFree("data/freeEn_"+parametrization+"_betaEq"
													+"_T"+std::to_string(temperature)+".txt");

	std::ofstream outSol("data/solution_cpa_"+parametrization+"_betaEq"+
												  +"_T"+std::to_string(temperature)+".txt");
	
	// if(Yp<0.5){
	// 	firstGuess={0.70, 0.73, 0.65,0.92, 0.97, 0.95};
	// }
	temperature*=1./Mnucleon;
	
	electron.temperature=temperature;

	for(int irho=0; irho<iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
		pasta.solveCPA_betaEq(rhoB, temperature, electron);
		
		f=pasta.f;
		Yp=pasta.YpG;
		
		if(pasta.f>=0. && pasta.f<=1.){
	
			double alpha, dim;
		
			for(int iDim=0; iDim<=2; iDim++){
			dim = (double) (iDim+1);
			
			for(int iType=0; iType<=1; iType++){
			if(iType==0){alpha=f;} /*droplets*/
			else if(iType==1){alpha=(1.-f);}/*bubbles*/

			sigma= getSurfaceTension(cluster, Yp, temperature);
			
			RdM(iDim, iType)= getRadiusD(dim, alpha, Yp, cluster, gas);
			
			RwM(iDim, iType)= RdM(iDim, iType)/pow(alpha, 1./dim);
			surfEnM(iDim, iType)= sigma*dim/RdM(iDim, iType);
			coulEnM(iDim, iType)= surfEnM(iDim, iType)/2.;

			PressureM(iDim, iType)= gas.getPressure()+ electron.pressure;
			ZeM(iDim, iType)= (cluster.proton.density-gas.proton.density)
													*4*M_PI*pow( getRadiusD(3., alpha, Yp, cluster, gas), 3.)/3.;
													
			AeM(iDim, iType)= (cluster.rhoB-gas.rhoB)
													*4*M_PI*pow( getRadiusD(3., alpha, Yp, cluster, gas), 3.)/3.;
													
			EnergyM(iDim, iType)= f*cluster.getEnergy()+(1.-f)*gas.getEnergy()+3.*alpha*coulEnM(iDim, iType) + electron.energy;
			EntropyM(iDim, iType) = f*cluster.getEntropy() +(1.-f)*gas.getEntropy() + electron.entropy;
			FreeEnM(iDim, iType) = (EnergyM(iDim, iType) - temperature*EntropyM(iDim, iType))/rhoB;
			GibbsEnM(iDim, iType)= FreeEnM(iDim, iType) - f*(cluster.proton.chemPot + cluster.neutron.chemPot)
									- (1.-f)*(gas.proton.chemPot + gas.neutron.chemPot);
			BulkEnM(iDim, iType)= (f*cluster.getEnergy()+(1.-f)*gas.getEnergy())/rhoB - 1.;
		    }
			}

			//Get phase that minimizes energy and fix dimensions:
			Eigen::MatrixXd::Index minRow, minCol;
			FreeEn= FreeEnM.minCoeff(&minRow, &minCol);
			FreeEn*=Mnucleon;

			Pressure=PressureM(minRow, minCol)*Mnucleon*pow(Mnucleon/hc, 3.);
			BulkEn=BulkEnM(minRow, minCol)*Mnucleon;	
			GibbsEn=GibbsEnM(minRow, minCol)*Mnucleon;
			Entropy=EntropyM(minRow, minCol);
			coulEn=coulEnM(minRow, minCol)*Mnucleon/rhoB;
			surfEn=surfEnM(minRow, minCol)*Mnucleon/rhoB;
			Rd=RdM(minRow, minCol)*(hc/Mnucleon) ;
			Rw=RwM(minRow, minCol)*(hc/Mnucleon) ;
			Ae= AeM(minRow, minCol);
			Ze= ZeM(minRow, minCol);
			iDimension=minRow+1;

		if(minCol==0){
			if(iDimension==1){iPlot=3;}
			if(iDimension==2){iPlot=2;}
			if(iDimension==3){iPlot=1;}
		}else if(minCol==1){
			if(iDimension==1){iPlot=3;}
			if(iDimension==2){iPlot=4;}
			if(iDimension==3){iPlot=5;}
		}
	

		outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " << electron.kf*Mnucleon/hc << " " << Pressure << " " 
				  << FreeEn - Mnucleon << " " << GibbsEn - Mnucleon << " " 
				  << BulkEn << " "   << Entropy*temperature*Mnucleon/rhoB << " " 
				  << coulEn << " " << surfEn << " "
				  << f << " " << Rd << " " << Rw << " " << iPlot << " " 
					<< Ae << " " << Ze
				  << std::endl;

		outFree << rhoB*pow(Mnucleon/hc, 3.) << " " 
				<< FreeEn - Mnucleon 			<< " " << (FreeEnM(2, 0)-1.)*Mnucleon << " " 
				<< (FreeEnM(1, 0)-1.)*Mnucleon << " " << (FreeEnM(0, 0)-1.)*Mnucleon << " " 
				<< (FreeEnM(1, 1)-1.)*Mnucleon << " " << (FreeEnM(2, 1)-1.)*Mnucleon 
			  	<< std::endl;

		outSol  << rhoB*pow(Mnucleon/hc, 3.) << " " << cluster.Yp << " " 
				<< cluster.proton.chemPot_eff << " " << cluster.neutron.chemPot_eff << " "
				<< cluster.proton.mass_eff << " "
				<< gas.proton.chemPot_eff << " " << gas.neutron.chemPot_eff << " " 
				<< gas.proton.mass_eff  << " "
				<< cluster.proton.chemPot << " " << cluster.neutron.chemPot << " "
				<< gas.proton.chemPot << " " << gas.neutron.chemPot << " "
				<< cluster.proton.density*pow(Mnucleon/hc, 3.) << " " 
				<< cluster.neutron.density*pow(Mnucleon/hc, 3.) << " " 
				<< gas.proton.density*pow(Mnucleon/hc, 3.) << " " 
				<< gas.neutron.density*pow(Mnucleon/hc, 3.) << " " 
			    << std::endl;

		}
		else{
			std::cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << std::endl;
			//setInitialGibbs(nup1, nun1, Mef1, nup2, nun2, Mef2);

		}
	}

	outGlobal.close();
	outFree.close();
	outSol.close();
  return 0;
}