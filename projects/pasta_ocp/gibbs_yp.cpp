// Pasta with non linear Walecka model using Mean Field Theory

#include "../../include/rmf_walecka.h"
#include "../../include/pasta.h"
#include "../../include/interpolator.h"
#include "../../include/constant.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

//===========================main code =============================

int main(){

	std::string parametrization;
	double rhoB, Yp, temperature;

	cout << "Specify the parametrization, proton fraction and temperature: " << endl;
	cin >> parametrization >> Yp >> temperature;
	cout << "You chose: " 
				<< parametrization  << " parametrization " << endl 
				<< "Yp= "<< Yp << endl  
				<< "T = " << temperature << endl;

//Declare nuclear matter: pasta and gas
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);

	nlwm_class hmg_matter(parametrization);

	double rhoBMax;

	// if(Yp>0.3)			rhoBMax= 0.16;
	// else if(Yp>0.2)	rhoBMax= 0.14;
	// else if(Yp>0.15)	rhoBMax= 0.1;
	// else if(Yp>0.1)	rhoBMax= 0.08;
	// else 						rhoBMax= 0.1;
	rhoBMax=0.15;
	rhoBMax*=pow(hc/Mnucleon, 3);

	int iRhoMax=200;
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
	double FreeEn_hmg;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	

	double Ae=0., Ze=0.;
	
	std::cout << "Yp, T= " << Yp << " " << temperature << std::endl;
	std::ofstream outGlobal("data/cpa_"+parametrization+"_yp"+std::to_string(Yp)+
												   +"_T"+std::to_string(temperature)+".txt");

	std::ofstream outFree("data/freeEn_"+parametrization+"_yp"+std::to_string(Yp)
													+"_T"+std::to_string(temperature)+".txt");

	std::ofstream outSol("data/solution_cpa_"+parametrization+"_yp"+std::to_string(Yp)+
												  +"_T"+std::to_string(temperature)+".txt");
	
	// if(Yp<0.5){
	// 	firstGuess={0.70, 0.73, 0.65,0.92, 0.97, 0.95};
	// }
	temperature*=1./Mnucleon;
	
	electron.temperature=temperature;

	for(int irho=0; irho<iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
		
		electron.density=Yp*rhoB;
		electron.kf=pow(3.*pi2*electron.density, 1./3.);
		electron.solveChemPotEff();
		electron.chemPot = electron.chemPot_eff;
		electron.calculateProperties();
		// std::cout << electron.chemPot << std::endl;


		//homogeneous matter:
		hmg_matter.setEOS_nucleons(rhoB, Yp, temperature);
		FreeEn_hmg	= (hmg_matter.getEnergy() + electron.energy -temperature*(hmg_matter.getEntropy()+electron.entropy))/rhoB;
		
		//pasta:
		pasta.solveCPA(rhoB, Yp, temperature);
		f=pasta.f;
		
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

		if(FreeEn< FreeEn_hmg*Mnucleon){
			if(minCol==0){
				if(iDimension==1){iPlot=3;}
				if(iDimension==2){iPlot=2;}
				if(iDimension==3){iPlot=1;}
			}else if(minCol==1){
				if(iDimension==1){iPlot=3;}
				if(iDimension==2){iPlot=4;}
				if(iDimension==3){iPlot=5;}
			}
		}else{
			iPlot=6;
		}

		outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " << Pressure << " " 
				  << FreeEn - Mnucleon << " " << GibbsEn - Mnucleon << " " 
				  << BulkEn << " "   << Entropy*temperature*Mnucleon/rhoB << " " 
				  << coulEn << " " << surfEn << " "
				  << f << " " << Rd << " " << Rw << " " << iPlot << " " 
					<< Ae << " " << Ze << " " << (FreeEn_hmg -1.)*Mnucleon  << " " 
				<< cluster.neutron.chemPot*Mnucleon << " "<< cluster.proton.chemPot*Mnucleon << " "
				  << gas.neutron.chemPot*Mnucleon << " "<< gas.proton.chemPot*Mnucleon
				  << std::endl;

		outFree << rhoB*pow(Mnucleon/hc, 3.) << " " 
				<< FreeEn - FreeEn_hmg*Mnucleon 			<< " " << (FreeEnM(2, 0)-FreeEn_hmg)*Mnucleon << " " 
				<< (FreeEnM(1, 0)-FreeEn_hmg)*Mnucleon << " " << (FreeEnM(0, 0)-FreeEn_hmg)*Mnucleon << " " 
				<< (FreeEnM(1, 1)-FreeEn_hmg)*Mnucleon << " " << (FreeEnM(2, 1)-FreeEn_hmg)*Mnucleon 
			  	<< std::endl;

		outSol  << rhoB*pow(Mnucleon/hc, 3.) << " " << cluster.Yp << " " << gas.Yp << " " 
				<< cluster.proton.chemPot_eff*Mnucleon << " " << cluster.neutron.chemPot_eff*Mnucleon << " "
				<< cluster.proton.mass_eff*Mnucleon << " "
				<< gas.proton.chemPot_eff*Mnucleon << " " << gas.neutron.chemPot_eff*Mnucleon << " " 
				<< gas.proton.mass_eff*Mnucleon  << " "
				<< cluster.proton.chemPot*Mnucleon << " " << cluster.neutron.chemPot*Mnucleon << " "
				<< gas.proton.chemPot*Mnucleon << " " << gas.neutron.chemPot*Mnucleon << " "
				<< cluster.proton.density*pow(Mnucleon/hc, 3.) << " " 
				<< cluster.neutron.density*pow(Mnucleon/hc, 3.) << " " 
				<< gas.proton.density*pow(Mnucleon/hc, 3.) << " " 
				<< gas.neutron.density*pow(Mnucleon/hc, 3.) << " " 
			    << std::endl;

		}
		else{
			std::cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << std::endl;
			//setInitialGibbs(nup1, nun1, Mef1, nup2, nun2, Mef2);
			outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " << 1./0. << " " 
				  << 1./0.<< " " << 1./0. << " " 
				  << 1./0. << " "   << 1./0. << " " 
				  << 1./0. << " " << 1./0. << " "
				  << 1./0 << " " << 1./0. << " " << 1./0. << " " << 1./0. << " " 
					<< 1./0. << " " << 1./0. << " " << (FreeEn_hmg -1.)*Mnucleon 
				  << std::endl;
		}
	}

	outGlobal.close();
	outFree.close();
	outSol.close();
  return 0;
}