// Pasta with non linear Walecka model using Mean Field Theory
#include "../../include/constant.hpp"
#include "../../include/particles.hpp"
#include "../../include/rmf_walecka.hpp"
#include "../../include/pasta.hpp"
#include "../../include/interpolator.hpp"
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
	
	
	double rhoB, Yp;

	cout << "Specify the parametrization and density (1/fm3): " << endl;
	cin >> parametrization >> rhoB;
	cout << "You chose: " << endl
				<< parametrization  << " parametrization " << endl 
				<< "nB = " << rhoB << "fm^-3" <<  endl;

	// double Lratio;
	// cout << "Specify L_d/R_W: " << endl;
	// cin >> Lratio ;	

	//Declare nuclear matter: pasta and gas:
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);


	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= electron.mass;
	electron.Q=-1;

	//Pasta values:
	MatrixXd PressureM(3,2),  FreeEnM(3,2), EnergyM(3,2), EntropyM(3,2), 
					coulEnM(3,2), surfEnM(3,2), 
					RdM(3,2), RwM(3,2), VcM(3,2), VwM(3,2), AeM(3, 2), ZeM(3, 2); 
					// BulkEnM(3,2), GibbsEnM(3,2),  
	
	double  f, Rd, Rw, sigma;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	

	double  Ze=0.;
	
	ofstream outTransport("data/transport_Tloop_nB"+to_string(rhoB)+".txt");

	rhoB*=pow(hc/Mnucleon, 3);
	double tempMin=0.1/Mnucleon;
  	double tempMax=5./Mnucleon;
  	int itMax=10;
  	double dt=  (tempMax-tempMin)/itMax;


	for(int it=0; it<=itMax; it++){
 		double temperature=tempMin+ it*dt;
		electron.temperature=temperature;

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
			double FreeEn= FreeEnM.minCoeff(&minRow, &minCol);
			(void) FreeEn;
			// Pressure=PressureM(minRow, minCol);
			Rd=RdM(minRow, minCol) ;
			Rw=RwM(minRow, minCol) ;
			// Ae= AeM(minRow, minCol);
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

		cDroplet  droplet(electron);
		cRod      rod(electron);
		cSlab     slab(electron);
	
		rod.setLim_CosTheta(-1., 1.);
		slab.setLim_CosTheta(-1., 1.);
		slab.setLim_Phi(0., 2.*M_PI);
		droplet.setRadius(RdM(2, minCol));
		double rad_rods= RdM(1, minCol); 
		double rad_slab= RdM(0, minCol);


		double L2w= droplet.getVolume()/(M_PI*rad_rods*rad_rods);	
		double L1w= sqrt(droplet.getVolume()/(2.*rad_slab));
		// double Vcl= droplet.getVolume();	

		double L_rod	=sqrt(2)*L1w;//*RwM(1,minCol);
		double L_slab	=sqrt(2)*L1w;//3.*RwM(0, minCol);

		rod.setLengths(rad_rods,L_rod);
		slab.setLengths(L_slab, L_slab, 2.*rad_slab);

		double qmin=0.;
		double qmax=2.;
		int iqmax= 300;
		double q;

		vector<double> qv, F3v2, Fav2, Fpv2;
		double coulIa=0., coulIp=0., nua=0., nup=0., nu_avg=0., nu_avg_inverse=0.;
		double nu3d, sigma0_3d;

		for(int iq=0; iq<=iqmax; iq++){
			q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
			qv.push_back(q);

		}


	  reverse(qv.begin(), qv.end());

		if(iDimension==2){	
			cout << " form factor calculation" << endl;
			Ze= (cluster.proton.density - gas.proton.density)*rod.getVolume();
			for(int iq=0; iq<=iqmax; iq++){
				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
				rod.setMomentum(q);
				double fa_= rod.getStructureFunction2_Axial(q);
				double fp_= rod.getStructureFunction2_Trans(q);
				Fav2.push_back(fa_);
				Fpv2.push_back(fp_);

				droplet.setMomentum(q);
				double f3_= droplet.getStructureFunction(q);
				F3v2.push_back(f3_*f3_);

			}

	   	reverse(Fav2.begin() , Fav2.end());
			reverse(Fpv2.begin() , Fpv2.end());
			rod.setMomentumVec(qv);
			rod.setFormFactor2Vec(Fav2);
			cout << " coulomb log calculation" << endl;
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
			cout << " form factor calculation" << endl;
			Ze= (cluster.proton.density - gas.proton.density)*slab.getVolume();
			for(int iq=0; iq<=iqmax; iq++){
				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
				slab.setMomentum(q);
				double fa_= slab.getStructureFunction2_Axial(q);
				double fp_= slab.getStructureFunction2_Trans(q);
				Fav2.push_back(fa_);
				Fpv2.push_back(fp_);

				droplet.setMomentum(q);
				double f3_= droplet.getStructureFunction(q);
				F3v2.push_back(f3_*f3_);

			}
	   	reverse(Fav2.begin() , Fav2.end());
			reverse(Fpv2.begin() , Fpv2.end());
			slab.setMomentumVec(qv);
			cout << " coulomb log calculation" << endl;
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
		cout << "3d form factor calculation" << endl;

			Ze= (cluster.proton.density - gas.proton.density)*droplet.getVolume();
			for(int iq=0; iq<=iqmax; iq++){
				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
				droplet.setMomentum(q);
				double f3_= droplet.getStructureFunction(q);
				F3v2.push_back(f3_*f3_);

			}

	  	reverse(F3v2.begin() , F3v2.end());

			droplet.setMomentumVec(qv);
			droplet.setFormFactor2Vec(F3v2);
	  	 // double coulInt3d=  integrate_coulomb(coulomb_gsl, &droplet);

			cout << "3d coulomb log calculation" << endl;

			double coulInt3d=getCoulombIntegral(droplet);
	 		nu3d			=  getFrequency(Ze, cluster.rhoB, coulInt3d, electron);
			sigma0_3d=  pow(eGS, 2.)*electron.density/(nu3d*hypot(electron.kf, electron.mass));

			coulIa=   		 coulInt3d;
			coulIp=   		 coulInt3d;
			nua	= nu3d;
			nup	= nu3d;
			nu_avg= nu3d;
			nu_avg_inverse= nu3d;

		}
		(void) nu_avg;
		(void) sigma0_3d;

		double sigma0_a= pow(eGS, 2.)*electron.density/(nua*hypot(electron.kf, electron.mass));
		double sigma0_p= pow(eGS, 2.)*electron.density/(nup*hypot(electron.kf, electron.mass));	
		double sigma0_avg = pow(eGS, 2.)*electron.density*nu_avg_inverse/hypot(electron.kf, electron.mass);		

		vector<double> c0v_={6.*alpha*coulEnM(0, minCol), 
										(3./2. + 2.*pow(10., 2.1*(alpha-0.3)))*alpha*coulEnM(1, minCol)};
		vector<double> av_= {2.*RwM(0, minCol), sqrt(2.*M_PI/sqrt(3.))*RwM(1, minCol)};
		vector<double> q0v_={2*M_PI/av_[0], 2*M_PI/av_[1]};
		vector<double> qcv_={4*M_PI/av_[0], 4*M_PI/av_[1]};
		vector<double> lambdav_={sqrt((1.+2.*alpha-2.*alpha*alpha)*RwM(0, minCol)*RwM(0, minCol)/45.),
								sqrt( 0.131*alpha*coulEnM(1, minCol)*RwM(1, minCol)*RwM(1, minCol)/c0v_[1])};

		double xi= 8.*pi2*pow(lambdav_[1], 2)*c0v_[1]/(pow(q0v_[1], 2.)*pow(qcv_[1], 2.)*temperature);

		double eta= pow(q0v_[0], 2.)*temperature/(8.*M_PI*lambdav_[0]*c0v_[0]);

    outTransport << temperature*Mnucleon << " " 
 						<< Rw*(hc/Mnucleon)  << " " //s-1
 						<< nua*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup*Mnucleon/MeVto_Sec << " " //s-1
 						<< nup/nua << " "
						<< nu_avg*Mnucleon/MeVto_Sec << " " << nu_avg_inverse*MeVto_Sec/Mnucleon << " " 
						<< coulIp/coulIa << " "
						<< coulIa << " " << coulIp << " " 
						<< Rd*(hc/Mnucleon) << " "
						<< sigma0_a*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_p*Mnucleon/MeVto_Sec << " " //s-1
						<< sigma0_avg*Mnucleon/MeVto_Sec << " " //s-1
						<< iDimension << " " << iPlot <<  " "
						<< xi << " " <<eta << " "  <<pow(2, -eta) << " " 
						<< c0v_[0]*Mnucleon*pow(Mnucleon/hc, 3.) << " " << c0v_[1]*Mnucleon*pow(Mnucleon/hc, 3.) << " "
						<< av_[0]*(hc/Mnucleon)  << " "<<  av_[1]*(hc/Mnucleon)  << " "
						<< lambdav_[0]*(hc/Mnucleon)  << " " << lambdav_[1]*(hc/Mnucleon) 
						<< endl;

	}else{
		cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << endl;
	}

}
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
	double hfunc_=0.;
	double omega2=omega*omega;
	double r= sqrt(nup*( omega2 + nup*nua) );
	double s= omega*sqrt(fabs(nup-nua));
	// double sr= s/r;

	if(nua>nup){
		// hfunc_= sr>1e-10 ? atan(s/r) : 0.;
		// hfunc_*=1./(s*r);
		hfunc_= atan(s/r)/(s*r);
	}else if(nua<nup){
		// hfunc_= sr>1e-10 ? atanh(s/r) : 0.;
		// hfunc_*=1./(s*r);
		hfunc_= atanh(s/r)/(s*r);
	}else hfunc_=1./(pow(nua, 3.)+omega2*nua);

	return hfunc_;
}