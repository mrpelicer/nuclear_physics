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
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

double getFrequency(double Z_,  double ni_, double coul_int_, particle electron);
double Hfunc(double nua_, double nup_, double omega_);
double sigmaAverage_Parallel			(double nua_, double nup_, double omega, particle electron);
double sigmaAverage_Perpendicular (double nua_, double nup_, double omega, particle electron);
double sigmaAverage_Hall 					(double nua_, double nup_, double omega, particle electron);
//===========================main code =============================

int main(){

	string parametrization;
	double rhoB, temperature;

	cout << "Specify the parametrization, density (1/fm3), temperature (MeV): " << endl;
	cin >> parametrization >> rhoB >> temperature ;
	cout << "You chose: " << endl
				<< parametrization  << " parametrization " << endl 
				<< "nB = " << rhoB << "fm^-3" <<  endl
				<< "T = " << temperature << " MeV" <<  endl;


	rhoB*=pow(hc/Mnucleon, 3);
	temperature*=1./Mnucleon;
	
	//Declare nuclear matter: pasta and gas:
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);

	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= electron.mass;
	electron.Q=-1;
	electron.temperature=temperature;

	//Pasta values:
	MatrixXd FreeEnM(3,2), RdM(3,2); 
					// BulkEnM(3,2), GibbsEnM(3,2),  
	
	double FreeEn, f, Rd;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
		
	pasta.solveCPA_betaEq(rhoB, temperature, electron);

	f=pasta.f;
	double Yp=pasta.YpG;
	double sigma= getSurfaceTension(cluster, Yp, temperature);
	
	if(pasta.f>=0. && pasta.f<=1. && sigma>0){

		double alpha, dim;
		for(int iDim=0; iDim<=2; iDim++){
		dim = (double) (iDim+1);
		for(int iType=0; iType<=1; iType++){
		if(iType==0){alpha=f;} /*droplets*/
		else if(iType==1){alpha=(1.-f);} /*bubbles*/
		
			RdM(iDim, iType)= getRadiusD(dim, alpha, Yp, cluster, gas);
			double surf= sigma*dim/RdM(iDim, iType);
			double coul= surf/2.;

			double energy_= f*cluster.getEnergy()+(1.-f)*gas.getEnergy()
																	+3.*alpha*coul + electron.energy;
			double entropy_= f*cluster.getEntropy() +(1.-f)*gas.getEntropy() + electron.entropy;
			FreeEnM(iDim, iType) = (energy_ - temperature*entropy_)/rhoB;
		}
		}

		//Get phase that minimizes energy and fix dimensions:
		MatrixXd::Index minRow, minCol;
		FreeEn= FreeEnM.minCoeff(&minRow, &minCol);
		(void) FreeEn;
		Rd=RdM(minRow, minCol) ;
		iDimension=minRow+1;
		if(minCol==0){
			if(iDimension==1){iPlot=3;}
			if(iDimension==2){iPlot=2;}
			if(iDimension==3){iPlot=1;}
		}else if(minCol==1){
			if(iDimension==1){iPlot=3;}
			if(iDimension==2){iPlot=2;}
			if(iDimension==3){iPlot=1;}
		}

 		//double ri= pow(3./(4.*M_PI*cluster.rhoB), 1./3. );
		//double Gamma= pow(Ze*eGS, 2.)/(ri*temperature);

		cDroplet  droplet(electron);
		cRod      rod(electron);
		cSlab     slab(electron);
		
		rod.setLim_CosTheta(-1., 1.);
		slab.setLim_CosTheta(-1., 1.);
		slab.setLim_Phi(0., 2.*M_PI);

		droplet.setRadius(RdM(2, 0));
		double rad_rods=  RdM(1, 0); 
		double rad_slab=  RdM(0, 0);

		double q;
		vector<double> qv, F3v2;
		double qmin=0.;
		double qmax=2.;
		int iqmax= 200;
		cout << "Form factor calculation" << endl;


//always calculate for 3d
		double Ze=(cluster.proton.density - gas.proton.density)*droplet.getVolume();

		for(int iq=0; iq<=iqmax; iq++){
			q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
			droplet.setMomentum(q);
			double f3_= droplet.getStructureFunction(q);
			qv.push_back(q);
			F3v2.push_back(f3_*f3_);
		}
		
    reverse(qv.begin(), qv.end());
    reverse(F3v2.begin() , F3v2.end());

		droplet.setMomentumVec(qv);
		droplet.setFormFactor2Vec(F3v2);
    // double coulInt3d=  integrate_coulomb(coulomb_gsl, &droplet); // getCoulombIntegral(droplet);
		double coulInt3d=getCoulombIntegral(droplet);
 		double nu3d			=  getFrequency(Ze, cluster.rhoB, coulInt3d, electron);
		
		double coulIa=0., coulIp=0., nua=0., nup=0., nu_avg=0., nu_avg_inverse=0.;
		double L;

		int ilmax=4;
		double Lmin=Rd;
		double Lmax= 10000.*Rd;
		double logLmin= log10(Lmin);
		double logLmax= log10(Lmax);
		ofstream outTransport("data/transport_Ld_T"+to_string(temperature*Mnucleon)+"_"
					+to_string(iDimension)+".txt");

		for( int il=0; il<=ilmax; il++){
			vector<double> Fav2, Fpv2;
//			L= Lmin + (Lmax-Lmin)*il/ilmax;
			L= pow(10., logLmax- il*(logLmax-logLmin)/ilmax);
			ofstream outform("data/form_factor"+to_string(L/Rd)+"_"
					+to_string(iDimension)+".txt");

		if(iDimension==2){	
			double L_rod	=L;
			rod.setLengths(rad_rods,L_rod);
			Ze= (cluster.proton.density - gas.proton.density)*rod.getVolume();

			for(int iq=0; iq<=iqmax; iq++){
				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
				rod.setMomentum(q);
				double fa_= rod.getStructureFunction2_Axial(q);
				double fp_= rod.getStructureFunction2_Trans(q);
				Fav2.push_back(fa_);
				Fpv2.push_back(fp_);
				outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ << endl;
			}

    	reverse(Fav2.begin() , Fav2.end());
			reverse(Fpv2.begin() , Fpv2.end());
			rod.setMomentumVec(qv);

			rod.setFormFactor2Vec(Fav2);
			// coulIa=   		 integrate_coulomb(coulomb_gsl, &rod);
			coulIa=   		 getCoulombIntegral(rod);

			rod.setFormFactor2Vec(Fpv2);
			// coulIp=   		 integrate_coulomb(coulomb_gsl, &rod);
			coulIp=   		 getCoulombIntegral(rod);


			nua	= 3.*getFrequency(Ze, cluster.rhoB, coulIa, electron);
			nup	= 3.*getFrequency(Ze, cluster.rhoB, coulIp, electron);

			nu_avg= (nua + 2.*nup)/3.;
			nu_avg_inverse= (2./nup + 1./nua)/3.;
		}else if(iDimension==1){
			// double L_slab	= sqrt(droplet.getVolume()/(2.*rad_slab));
			double L_slab=L;
			slab.setLengths(L_slab, L_slab, 2.*rad_slab);
			Ze= (cluster.proton.density - gas.proton.density)*slab.getVolume();
			slab.setMomentum(q);
			for(int iq=0; iq<=iqmax; iq++){
				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
				rod.setMomentum(q);
				double fa_= slab.getStructureFunction2_Axial(q);
				double fp_= slab.getStructureFunction2_Trans(q);
				Fav2.push_back(fa_);
				Fpv2.push_back(fp_);
				outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ << endl;

			}

    	reverse(Fav2.begin() , Fav2.end());
			reverse(Fpv2.begin() , Fpv2.end());
			slab.setMomentumVec(qv);

			slab.setFormFactor2Vec(Fav2);
			// coulIa=   		 integrate_coulomb(coulomb_gsl, &slab);
			coulIa=   		 getCoulombIntegral(slab);

			slab.setFormFactor2Vec(Fpv2);
			// coulIp=   		 integrate_coulomb(coulomb_gsl, &slab);
			coulIp=   		 getCoulombIntegral(slab);


			nua	= 3.*getFrequency(Ze, cluster.rhoB, coulIa, electron);
			nup	= 3.*getFrequency(Ze, cluster.rhoB, coulIp, electron);

			nu_avg				= (nua + 2.*nup)/3.;
			nu_avg_inverse= (1./nua + 2./nup)/3.;

		}else{
			L=RdM(minCol, minCol);
			coulIa=   		 coulInt3d;
			coulIp=   		 coulInt3d;
			nua	= nu3d;
			nup	= nu3d;
			nu_avg= nu3d;
			nu_avg_inverse= nu3d;

		}
		

    outTransport << L/RdM(minRow, minCol) << " " 
 						<< nu3d*Mnucleon/MeVto_Sec << " " //s-1
 						<< nua*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup/nua << " "
						<< nu_avg*Mnucleon/MeVto_Sec << " " << nu_avg_inverse*MeVto_Sec/Mnucleon << " " 
						<< coulIp/coulIa << " "
						<< coulIa << " " << coulIp << " " 
						<< Rd*(hc/Mnucleon) << " " << rhoB*pow(Mnucleon/hc, 3.) << " " << Ze << " " 
						<< iDimension << " " << iPlot
						<< endl;
		outform.close();
	}
	  outTransport.close();



		}else{
			cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << endl;
		}

  return 0;
}

double getFrequency(double Z_, double ni_, double coul_int_, particle electron){

	return 4.*M_PI*pow(Z_, 2.)*pow(eGS, 4.)*ni_*hypot(electron.mass, electron.kf)
		 																							*coul_int_	/pow(electron.kf, 3.);

}

double sigmaAverage_Parallel(double nua, double nup, double omega, particle electron){
	double s_= (omega*omega+nup*nup)*(omega*omega+nup*nua)*Hfunc(nua, nup, omega) - nup;
	return pow(eGS, 2.)*electron.density*s_/(omega*omega*hypot(electron.kf, electron.mass)								);
 }

double sigmaAverage_Perpendicular(double nua, double nup, double omega, particle electron){
	double s_= ( nup*nua*(omega*omega - nup*nup)*Hfunc(nua, nup, omega) + nup)/2.;
	return pow(eGS, 2.)*electron.density*s_/(omega*omega*hypot(electron.kf, electron.mass));
}

double sigmaAverage_Hall (double nua, double nup, double omega, particle electron){
	double s_= omega*(1. - nua*nup*nup*Hfunc(nua, nup, omega));
	return pow(eGS, 2.)*electron.density*s_/(omega*omega*hypot(electron.kf, electron.mass));
}


double Hfunc(double nua, double nup, double omega){
	double h;
	double omega2=omega*omega;
	double r= sqrt(nup*( omega2 + nup*nua) );
	double s= omega*fabs(nup-nua);
	
	if(nua>nup) 			h=atan(s/r)/(s*r);
	else if(nua<nup)	h=atanh(s/r)/(s*r);
	else 							h=1./(pow(nua, 3.)+omega2*nua);

	return h;
}