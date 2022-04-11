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
	
	
	double rhoB, Yp, temperature, Bfield, thetab=0.;

	cout << "Specify the parametrization, temperature (MeV) and magnetic field (G): " << endl;
	cin >> parametrization >> temperature >> Bfield ;
	cout << "You chose: " << endl
				<< parametrization  << " parametrization " << endl 
				<< "T = " << temperature << " MeV" <<  endl
				<< "B = " << Bfield << " G" <<  endl;

	//Declare nuclear matter: pasta and gas:
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);

	double rhoBMax=0.09*pow(hc/Mnucleon, 3);
	double rhoBMin=(5e-5)*pow(hc/Mnucleon, 3);

	int iRhoMax=10;
	double dRho= (rhoBMax-rhoBMin)/iRhoMax;

	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= electron.mass;
	electron.Q=-1;

	Bfield*=Gauss_to_Mev2_GS/pow(Mnucleon, 2.);
	//Pasta values:
	MatrixXd PressureM(3,2),  FreeEnM(3,2), EnergyM(3,2), EntropyM(3,2), 
					coulEnM(3,2), surfEnM(3,2), 
					RdM(3,2), RwM(3,2), VcM(3,2), VwM(3,2), AeM(3, 2), ZeM(3, 2); 
					// BulkEnM(3,2), GibbsEnM(3,2),  
	
	double Pressure, FreeEn, f, Rd, Rw, sigma;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	

	double Ae=0., Ze=0.;
	
	ofstream outGlobal("data/cpa_"+parametrization+
																			"_betaEq_T"+to_string(temperature)+".txt");

	// ofstream outFree("data/freeEn_"+parametrization+"_yp"+to_string(Yp)
	// 												+"_T"+to_string(temperature)+".txt");

	// ofstream outSol("data/solution_cpa_"+parametrization+"_yp"+to_string(Yp)+
	// 											  +"_T"+to_string(temperature)+".txt");
	ofstream outTransport("data/transport_T"+to_string(temperature)+".txt");
	// ofstream outTransport("data/transportdominant_T"+to_string(temperature)+".txt");


	// if(Yp<0.5){
	// 	firstGuess={0.70, 0.73, 0.65,0.92, 0.97, 0.95};
	// }
	temperature*=1./Mnucleon;
	
	electron.temperature=temperature;

	for(int irho=0; irho<=iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
		
		// cout << electron.chemPot << endl;
		// pasta.solveCPA(rhoB, Yp, temperature);		
		pasta.solveCPA_betaEq(rhoB, temperature, electron);

		f=pasta.f;
		Yp=pasta.YpG;
		sigma= getSurfaceTension(cluster, Yp, temperature);

		if(pasta.f>=0. && pasta.f<=1. && sigma>0){
	
			double alpha, dim;
			for(int iDim=0; iDim<=2; iDim++){
			dim = (double) (iDim+1);
			
			for(int iType=0; iType<=1; iType++){
			if(iType==0){alpha=f;} /*droplets*/
			else if(iType==1){alpha=(1.-f);} /*bubbles*/
			
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
		    }
			}

			//Get phase that minimizes energy and fix dimensions:
			MatrixXd::Index minRow, minCol;
			FreeEn= FreeEnM.minCoeff(&minRow, &minCol);

			Pressure=PressureM(minRow, minCol);
			// BulkEn=BulkEnM(minRow, minCol)*Mnucleon;	
			// GibbsEn=GibbsEnM(minRow, minCol)*Mnucleon;
			// Entropy=EntropyM(minRow, minCol);
			// coulEn=coulEnM(minRow, minCol)*Mnucleon/rhoB;
			// surfEn=surfEnM(minRow, minCol)*Mnucleon/rhoB;
			Rd=RdM(minRow, minCol) ;
			Rw=RwM(minRow, minCol) ;
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

		double ri= pow(3./(4.*M_PI*cluster.rhoB), 1./3. );
		double Gamma= pow(Ze*eGS, 2.)/(ri*temperature);

		cDroplet  droplet(electron);
		cRod      rod(electron);
		cSlab     slab(electron);
		
		rod.setLim_CosTheta(-1., 1.);
		slab.setLim_CosTheta(-1., 1.);
		slab.setLim_Phi(0., 2.*M_PI);

		droplet.setRadius(RdM(2));
		double rad_rods= RdM(1); 
		double rad_slab= RdM(0);
		double Vcl= droplet.getVolume();			
		double L_rod	=	Vcl/(M_PI*pow(rad_rods, 2.));
		double L_slab	= sqrt(Vcl/(2.*rad_slab));
			
		outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " 
			<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.) << " " 
		  << (FreeEn - 1.)*Mnucleon << " " 
			<< Ze << " " << " " << Ae  << " " 
		  << Gamma << " " << Yp << " "<< f << " " 
			<< Rd*(hc/Mnucleon) << " " << Rw*(hc/Mnucleon) << " " 
			<< L_rod*(hc/Mnucleon) << " " << L_slab*(hc/Mnucleon) << " "
			<< iPlot << " " << sigma
		  << endl;
	
		rod.setLengths(rad_rods,L_rod);
		slab.setLengths(L_slab, L_slab, 2.*rad_slab);
		

		vector<double> qv, F3v2, F2pv2, F1pv2, F2av2, F1av2;
		double qmin=0.;
		double qmax=2.;
		int iqmax= 300;
		double q;
		ofstream outform("data/form_factor"+to_string(irho)+".txt");
		cout << "Form factor calculation" << endl;
		
		for(int iq=0; iq<=iqmax; iq++){
			q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
			droplet.setMomentum(q);
			rod.setMomentum(q);
			slab.setMomentum(q);
			double f3_= droplet.getStructureFunction(q);
			double f2fa_= rod.getStructureFunction(q, cos(0.));
			double f2fp_= rod.getStructureFunction(q, cos(M_PI/2.));
			double f1fa_= slab.getStructureFunction(q, cos(0.), 0.);
			double f1fp_= slab.getStructureFunction(q, cos(M_PI/2.), M_PI/2.);
			
			double f2a_2= rod.getStructureFunction2_Axial(q);
			double f2p_2= rod.getStructureFunction2_Trans(q);
			double f1a_2= slab.getStructureFunction2_Axial(q);
			double f1p_2= slab.getStructureFunction2_Trans(q);

			outform << q/electron.kf << " " 
						<< pow(f3_, 2.) << " " 	
						<< pow(f2fa_, 2.) << " " 	<< pow(f2fp_, 2.) << " " 
						<< pow(f1fa_, 2.) << " " 	<< pow(f1fp_, 2.) << " " 
						<< f2a_2 << " " 	<< f2p_2 << " " 
						<< f1a_2 << " " 	<< f1p_2 << endl;

			cout << q/electron.kf << " " 
						<< pow(f3_, 2.) << " " 	
						<< pow(f2fa_, 2.) << " " 	<< pow(f2fp_, 2.) << " " 
						<< pow(f1fa_, 2.) << " " 	<< pow(f1fp_, 2.) << " " 
						<< f2a_2 << " " 	<< f2p_2 << " " 
						<< f1a_2 << " " 	<< f1p_2 << endl;


				qv.push_back(q);
				F3v2.push_back(f3_*f3_);
				F2av2.push_back(f2a_2);
				F2pv2.push_back(f2p_2);
				F1av2.push_back(f1a_2);
				F1pv2.push_back(f1p_2);
		}
	  outform.close();

    reverse(qv.begin(), qv.end());
    reverse(F3v2.begin() , F3v2.end());
		reverse(F2av2.begin(), F2av2.end());
		reverse(F2pv2.begin(), F2pv2.end());
		reverse(F1av2.begin(), F1av2.end());
		reverse(F1pv2.begin(), F1pv2.end());

    droplet.setMomentumVec(qv);
		rod.setMomentumVec(qv);
		slab.setMomentumVec(qv);
		// droplet.setFormFactor2Vec(F3v2);
		// rod.setFormFactorA2Vec(F2av2);
		// rod.setFormFactorT2Vec(F2pv2);
		// slab.setFormFactorA2Vec(F1av2);
		// slab.setFormFactorT2Vec(F1pv2);

		cout << "Coulomb Integrals" << endl;

		droplet.setFormFactor2Vec(F3v2);
    double coulInt3d =   		 integrate_coulomb(coulomb_gsl, &droplet);
    // double coulInt3d=   getCoulombIntegral(droplet);
		
		rod.setFormFactor2Vec(F2av2);
		double coulInta2d=   		 integrate_coulomb(coulomb_gsl, &rod);
		rod.setFormFactor2Vec(F2pv2);
		double coulIntp2d=   		 integrate_coulomb(coulomb_gsl, &rod);

		slab.setFormFactor2Vec(F1av2);
		double coulInta1d=   		 integrate_coulomb(coulomb_gsl, &slab);
		slab.setFormFactor2Vec(F1pv2);
		double coulIntp1d=   		 integrate_coulomb(coulomb_gsl, &slab);
		
		// rod.setFormFactor2Vec(F2av2);
		// double coulInta2d=  getCoulombIntegral(rod);
		// rod.setFormFactor2Vec(F2pv2);
		// double coulIntp2d=  getCoulombIntegral(rod);
		// slab.setFormFactor2Vec(F1av2);
		// double coulInta1d=  getCoulombIntegral(slab);
		// slab.setFormFactor2Vec(F1pv2);
		// double coulIntp1d=  getCoulombIntegral(slab);


 		double nu3d		=  	 getFrequency(ZeM(2), cluster.rhoB, coulInt3d, electron);
		
		double nua2d	= 3.*getFrequency(ZeM(1), cluster.rhoB, coulInta2d, electron);
		double nup2d	= 3.*getFrequency(ZeM(1), cluster.rhoB, coulIntp2d, electron);
						
		double nua1d	= 3.*getFrequency(ZeM(0), cluster.rhoB, coulInta1d, electron);
 		double nup1d	= 3.*getFrequency(ZeM(0), cluster.rhoB, coulIntp1d, electron);

		double nu2d_avg= (nua2d + 2.*nup2d)/3.;
		double nu1d_avg= (nua1d + 2.*nup1d)/3.; 

		double nu2d_avg_inverse= (2./nup2d + 1./nua2d)/3.;
		double nu1d_avg_inverse= (2./nup1d + 1./nua1d)/3.;

		double sigma0_3d = pow(eGS, 2.)*electron.density/(nu3d*hypot(electron.kf, electron.mass));
		
		double sigma0_a2d= pow(eGS, 2.)*electron.density/(nua2d*hypot(electron.kf, electron.mass));
		double sigma0_p2d= pow(eGS, 2.)*electron.density/(nup2d*hypot(electron.kf, electron.mass));
		double sigma0_2d = pow(eGS, 2.)*electron.density*nu2d_avg_inverse/hypot(electron.kf, electron.mass);	

		double sigma0_a1d= pow(eGS, 2.)*electron.density/(nua1d*hypot(electron.kf, electron.mass));
		double sigma0_p1d= pow(eGS, 2.)*electron.density/(nup1d*hypot(electron.kf, electron.mass));
		double sigma0_1d = pow(eGS, 2.)*electron.density*nu1d_avg_inverse/hypot(electron.kf, electron.mass);	
		
	// Calculation with magnetic field:
		Matrix3d sigmaB_2d, sigmaB_1d;
		double bx= cos(thetab), bz=sin(thetab);
		double omega= eGS*Bfield/hypot(electron.kf, electron.mass);		

		double omega2=omega*omega, bx2= bx*bx, bz2=bz*bz;
		 
		double Delta2d= nua2d*pow(nup2d, 2.)+ omega2*bx2*nup2d+omega2*bz2*nua2d;
		double Delta1d= nua1d*pow(nup1d, 2.)+ omega2*bx2*nup1d+omega2*bz2*nua1d;

		sigmaB_2d(0,0) = nua2d*nup2d+omega2*bx2;
		sigmaB_2d(0,1) =-omega*bz*nua2d;
		sigmaB_2d(0,2) = omega2*bx*bz;
		sigmaB_2d(1,0) = omega*bz*nua2d;
		sigmaB_2d(1,1) = nua2d*nup2d;
		sigmaB_2d(1,2) = -omega*bx*nup2d;
		sigmaB_2d(2,0) = omega2*bx*bz;
		sigmaB_2d(2,1) = omega*bx*nup2d;
		sigmaB_2d(2,2) = nup2d*nup2d+omega2*bz2;
		sigmaB_2d*= pow(eGS, 2.)*electron.density/(Delta2d*hypot(electron.kf, electron.mass));

		sigmaB_1d(0,0) = nua1d*nup1d+omega2*bx2;
		sigmaB_1d(0,1) =-omega*bz*nua1d;
		sigmaB_1d(0,2) = omega2*bx*bz;
		sigmaB_1d(1,0) = omega*bz*nua1d;
		sigmaB_1d(1,1) = nua1d*nup1d;
		sigmaB_1d(1,2) = -omega*bx*nup1d;
		sigmaB_1d(2,0) = omega2*bx*bz;
		sigmaB_1d(2,1) = omega*bx*nup1d;
		sigmaB_1d(2,2) = nup1d*nup1d+omega2*bz2;
		sigmaB_1d*= pow(eGS, 2.)*electron.density/(Delta1d*hypot(electron.kf, electron.mass));

//Average over domains:		
// filled with ( parallel, perpendicular, Hall)
		double x3d= omega/nu3d;
		Vector3d sigma_avg_3d(1., 1./(1+x3d*x3d), x3d/(1.+x3d*x3d) );
		sigma_avg_3d*=sigma0_3d;
		
		Vector3d sigma_avg_2d(sigmaAverage_Parallel(nua2d, nup2d, omega, electron),
													sigmaAverage_Perpendicular(nua2d, nup2d, omega, electron),
													sigmaAverage_Hall(nua2d, nup2d, omega, electron)); 
		
		Vector3d sigma_avg_1d(sigmaAverage_Parallel(nua1d, nup1d, omega, electron),
											sigmaAverage_Perpendicular(nua1d, nup1d, omega, electron),
											sigmaAverage_Hall(nua1d, nup1d, omega, electron)); 

		double xa2d, xp2d, xa1d, xp1d;
		xa2d=omega/nua2d;
		xp2d=omega/nup2d;	
		xa1d=omega/nua1d;
		xp1d=omega/nup1d;
		// const IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

		double ratio_nu2d= max(nua2d, nup2d)/min(nua2d, nup2d);
		double ratio_nu1d= max(nua1d, nup1d)/min(nua1d, nup1d);

	
    outTransport << rhoB*pow(Mnucleon/hc, 3.) << " " 
						<< x3d << " " << xa2d << " " << xp2d << " " << xa1d << " " << xp1d << " "
						<< ratio_nu2d << " " << ratio_nu1d << "  " 
 						<< nu3d*Mnucleon/MeVto_Sec << " " //s-1
 						<< nua2d*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup2d*Mnucleon/MeVto_Sec << " " //s-1
 						<< nua1d*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup1d*Mnucleon/MeVto_Sec << " " //s-1
						<< nu2d_avg*Mnucleon/MeVto_Sec << " " 
						<< nu1d_avg*Mnucleon/MeVto_Sec << " " 
						<< nu2d_avg_inverse*Mnucleon/MeVto_Sec << " " 
						<< nu1d_avg_inverse*Mnucleon/MeVto_Sec << " " 
						<< sigma0_3d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_2d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_1d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_a2d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_p2d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_a1d*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_p1d*Mnucleon/MeVto_Sec << " " //s-1
						<< endl;


		}else{
			cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << endl;
		}
	}

	outGlobal.close();
	// outFree.close();
	// outSol.close();
	outTransport.close();
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
	else 								h=1./(pow(nua, 3.)+omega2*nua);

	return h;
}