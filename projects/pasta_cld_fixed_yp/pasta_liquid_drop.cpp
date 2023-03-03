// Pasta with non linear Walecka model using Mean Field Theory
#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
#include "../../include/pasta.h"
#include "../../include/interpolator.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

// typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

double getZeff(double Z3, int idim_, int iAxis, double R3, double Ri);

using namespace std;

//===========================main code =============================

int main(){

	std::string parametrization;
	
	//Declare nuclear matter: pasta and gas:
	
	double rhoB, Yp, temperature;

	cout << "Specify the parametrization, proton fraction and temperature: " << endl;
	cin >> parametrization >> Yp >> temperature;
	cout << "You chose: " 
				<< parametrization  << " parametrization " << endl 
				<< "Yp= "<< Yp << endl  
				<< "T = " << temperature << endl;

	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);
	nlwm_class hmg(parametrization);
	
	double rhoBMax=0.13*pow(hc/Mnucleon, 3);
	if(Yp<0.4 && Yp >0.2)rhoBMax=0.11*pow(hc/Mnucleon, 3);
	else if(Yp<0.2)rhoBMax=0.11*pow(hc/Mnucleon, 3);

	int iRhoMax=400;
	double dRho= rhoBMax/iRhoMax;
		
	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= Me/Mnucleon;

	//Pasta values:
	Eigen::MatrixXd pressM(3,2), bulk_enM(3,2), gibbs_enM(3,2), free_enM(3,2), enerM(3,2),
					entropyM(3,2), cou_enM(3,2), surf_enM(3,2), RdM(3,2), RwM(3,2), vol_cM(3,2), vol_wM(3,2),
					AeM(3, 2), ZeM(3, 2), alphaM(3,2), munV(3,2), mupV(3,2);
	
	double press, bulk_en, free_en, gibbs_en, entropy, cou_en, surf_en,
					f, Rd, Rw, sigma, vol_c, vol_w, mun, mup;

	double free_en_hmg;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	

	double Ae=0., Ze=0.;
	
	std::cout << "Yp, T= " << Yp << " " << temperature << std::endl;
	std::ofstream outGlobal("data/cpa_"+parametrization+"_yp"+std::to_string(Yp)+
												   +"_T"+std::to_string(temperature)+".txt");

	std::ofstream outFree("data/free_en_"+parametrization+"_yp"+std::to_string(Yp)
													+"_T"+std::to_string(temperature)+".txt");

	std::ofstream outSol("data/solution_cpa_"+parametrization+"_yp"+std::to_string(Yp)+
												  +"_T"+std::to_string(temperature)+".txt");
	
	temperature*=1./Mnucleon;
	
	electron.temperature=temperature;

	for(int irho=0; irho<iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
		// rhoB=(rhoBMax+(double)irho*dRho);
		electron.density=Yp*rhoB;
		electron.kf=pow(3.*pi2*electron.density, 1./3.);
		electron.solveChemPotEff();
		electron.chemPot = electron.chemPot_eff;
		electron.calculateProperties();
	
		double alpha, dim;
	
		hmg.setEOS_nucleons(rhoB, Yp, temperature);
		free_en_hmg= (hmg.getEnergy() + electron.energy - temperature*(hmg.getEntropy() + electron.entropy))/rhoB;
		for(int iDim=0; iDim<=2; iDim++){
		dim = (double) (iDim+1);
		
		for(int iType=0; iType<=1; iType++){
		
			if (iType==1 && iDim ==0){
				free_enM(iDim, iType) = 1./0.;
				break;
			}

			
			pasta.solveCLD(rhoB, Yp, temperature, dim, iType);
		
			f=pasta.f;
			if(iType==0)			alpha=f; /*droplets*/
			else if(iType==1)	alpha=(1.-f);/*bubbles*/
			sigma= getSurfaceTension(cluster, Yp, temperature);
			
			alphaM(iDim, iType)= alpha;
			RdM(iDim, iType)= getRadiusD(dim, alpha, Yp, cluster, gas);
			RwM(iDim, iType)= RdM(iDim, iType)/pow(alpha, 1./dim);
			surf_enM(iDim, iType)= sigma*dim/RdM(iDim, iType);
			cou_enM(iDim, iType)= surf_enM(iDim, iType)/2.;
			double fsc_ = surf_enM(iDim, iType) +cou_enM(iDim, iType);

			pressM(iDim, iType)= gas.getPressure()+ electron.pressure;
			ZeM(iDim, iType)= (cluster.proton.density-gas.proton.density)
													*4*M_PI*pow( getRadiusD(3., alpha, Yp, cluster, gas), 3.)/3.;
													
			AeM(iDim, iType)= (cluster.rhoB-gas.rhoB)
													*4*M_PI*pow( getRadiusD(3., alpha, Yp, cluster, gas), 3.)/3.;
													
			enerM(iDim, iType)= f*cluster.getEnergy()+(1.-f)*gas.getEnergy()+3.*alpha*cou_enM(iDim, iType) + electron.energy;
			entropyM(iDim, iType) = f*cluster.getEntropy() +(1.-f)*gas.getEntropy() + electron.entropy;
			free_enM(iDim, iType) = (enerM(iDim, iType) - temperature*entropyM(iDim, iType))/rhoB;
			gibbs_enM(iDim, iType)= free_enM(iDim, iType) - f*(cluster.proton.chemPot + cluster.neutron.chemPot)
									- (1.-f)*(gas.proton.chemPot + gas.neutron.chemPot);
			bulk_enM(iDim, iType)= (f*cluster.getEnergy()+(1.-f)*gas.getEnergy())/rhoB - 1.;
			
			vol_cM(iDim,iType) = 4.*M_PI*pow(getRadiusD(3., alpha, Yp, cluster, gas), 3.)/3.;
			vol_wM(iDim,iType) = vol_cM(iDim,iType)/alpha;

			mupV(iDim,iType)  = gas.proton.chemPot
							+	(alpha/(1.-f))*
				(-2.*fsc_/(3.*(cluster.proton.density- gas.proton.density))
				+dim*getSurfaceTensionDerivative(cluster, Yp, temperature)
						*(1.-f)*(1.-Yp)/(rhoB*RdM(iDim,iType))
				);
								 
			munV(iDim,iType)  = gas.neutron.chemPot		-alpha/(1.-f)
						*dim*getSurfaceTensionDerivative(cluster, Yp, temperature)*(1.-f)
						*Yp/(rhoB*RdM(iDim,iType));	
		
			
		}
		}

	//Get phase that minimizes energy and fix dimensions:
		Eigen::MatrixXd::Index minRow, minCol;
		free_en		= free_enM.minCoeff(&minRow, &minCol);
		press		=pressM(minRow, minCol);
		bulk_en		=bulk_enM(minRow, minCol);	
		gibbs_en	=gibbs_enM(minRow, minCol);
		entropy		=entropyM(minRow, minCol);
		cou_en		=cou_enM(minRow, minCol)*Mnucleon/rhoB;
		surf_en		=surf_enM(minRow, minCol)*Mnucleon/rhoB;
		f			= alphaM(minRow,minCol);
		Rd			=RdM(minRow, minCol);
		Rw			=RwM(minRow, minCol);
		Ae			= AeM(minRow, minCol);
		Ze			= ZeM(minRow, minCol);
		vol_c		= vol_cM(minRow, minCol);
		vol_w		= vol_wM(minRow, minCol);
		mun			= munV(minRow,minCol);
		mup			= mupV(minRow,minCol);
		iDimension=minRow+1;

	if(free_en< free_en_hmg){
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

		double L=1./0.;
		if(iDimension==1) L = sqrt(vol_cM(2,minCol)/(2.*Rd));
		if(iDimension==2) L = vol_cM(2,minCol)/(M_PI*Rd*Rd);
		
		outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " << press*Mnucleon*pow(Mnucleon/hc, 3.) << " " 
				  << (free_en -1.)* Mnucleon << " " << (gibbs_en -1.)*Mnucleon << " " 
				  << bulk_en*Mnucleon << " "   << entropy/rhoB << " " 
				  << cou_en << " " << surf_en << " "
				  << f << " " << Rd*(hc/Mnucleon)  << " " << Rw*(hc/Mnucleon) << " " << iPlot << " " << iDimension << " "
				  << Ae << " " << Ze << " " << L*(hc/Mnucleon) << " "<< vol_c*pow(hc/Mnucleon, 3.) << " " 
				  << vol_w*pow(hc/Mnucleon, 3.) << " " << mun*Mnucleon  << " " << mup*Mnucleon << " " 
				  << (free_en_hmg-1.)*Mnucleon << " " 
				  << cluster.neutron.chemPot*Mnucleon << " "<< cluster.proton.chemPot*Mnucleon << " "
				  << gas.neutron.chemPot*Mnucleon << " "<< gas.proton.chemPot*Mnucleon
				  << std::endl;

		outFree << rhoB*pow(Mnucleon/hc, 3.) << " " 
				<< (free_en - 1.)*Mnucleon 			<< " " << (free_enM(2, 0)-1.)*Mnucleon << " " 
				<< (free_enM(1, 0)-1.)*Mnucleon << " " << (free_enM(0, 0)-1.)*Mnucleon << " " 
				<< (free_enM(1, 1)-1.)*Mnucleon << " " << (free_enM(2, 1)-1.)*Mnucleon 
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

	outGlobal.close();
	outFree.close();
	outSol.close();
  return 0;
}