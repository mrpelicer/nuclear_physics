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
	
	
	double rhoB, Yp, temperature;

	cout << "Specify the parametrization, density (1/fm3) and temperature (MeV): " << endl;
	cin >> parametrization >> rhoB >> temperature;
	cout << "You chose: " << endl
				<< parametrization  << " parametrization " << endl 
				<< "nB = " << rhoB << "fm^-3" <<  endl
				<< "T = " << temperature << " MeV" <<  endl;

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
	
	double f, sigma, Ze=0., Rw;
	int iDimension, iPlot=0;
	//iPlot= 1, 2, 3, 4, 5 (spheres, rods, slabs, tubes, bubbles)	
	
// ofstream outGlobal("data/cpa_"+parametrization+
	// 																		"_betaEq_T"+to_string(temperature)+".txt");

	rhoB*=pow(hc/Mnucleon, 3);
	temperature*=1./Mnucleon;
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
		// Rd=RdM(minRow, minCol) ;
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
	// double ri= pow(3./(4.*M_PI*cluster.rhoB), 1./3. );
	// double Gamma= pow(Ze*eGS, 2.)/(ri*temperature);
	cDroplet  droplet(electron);
	cRod      rod(electron);
	cSlab     slab(electron);
	
	rod.setLim_CosTheta(-1., 1.);
	slab.setLim_CosTheta(-1., 1.);
	slab.setLim_Phi(0., 2.*M_PI);
	droplet.setRadius(RdM(2, minCol));
	double rad_rods= RdM(1, minCol); 
	double rad_slab= RdM(0, minCol);
	// double Vcl= droplet.getVolume();	
			vector<double> c0v_={6.*alpha*coulEnM(0, minCol), 
											(3./2. + 2.*pow(10., 2.1*(alpha-0.3)))*alpha*coulEnM(1, minCol)};
			vector<double> av_= {2.*RwM(0, minCol), sqrt(2.*M_PI/sqrt(3.))*RwM(1, minCol)};
			vector<double> q0v_={2*M_PI/av_[0], 2*M_PI/av_[1]};
			vector<double> qcv_={4*M_PI/av_[0], 4*M_PI/av_[1]};
			vector<double> lambdav_={sqrt((1.+2.*alpha-2.*alpha*alpha)*RwM(0, minCol)*RwM(0, minCol)/45.),
									sqrt( 0.131*alpha*coulEnM(1, minCol)*RwM(1, minCol)*RwM(1, minCol)/c0v_[1])};

			double xi= 8.*pi2*pow(lambdav_[1], 2)*c0v_[1]/(pow(q0v_[1], 2.)*pow(qcv_[1], 2.)*temperature);
			double eta= pow(q0v_[0], 2.)*temperature/(8.*M_PI*lambdav_[0]*c0v_[0]);

			double L1w= sqrt(droplet.getVolume()/(2.*rad_slab));
			double epsilon = pow(2, -eta);
			
			double L_slab	=L1w/pow(epsilon, 1/(2*eta));
			double L_rod	=L_slab;

			double Lsize=0;

			if(iDimension==1){Lsize=L_slab;}
			if(iDimension==2){Lsize=L_rod;}
			if(iDimension==3){Lsize=NAN;}

			double Lratio = Lsize/Rw;
	// double L_rod	=Lratio*RwM(1, minCol);
	// double L_slab	=Lratio*RwM(0, minCol);	
	// double L_rod	=	Vcl/(M_PI*pow(rad_rods, 2.));
	// double L_slab	= sqrt(Vcl/(2.*rad_slab));
	
	// outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " 
	// 	<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.) << " " 
	//   << (FreeEn - 1.)*Mnucleon << " " 
	// 	<< Ze << " " << " " << Ae  << " " 
	//   << Gamma << " " << Yp << " "<< f << " " 
	// 	<< Rd*(hc/Mnucleon) << " " << Rw*(hc/Mnucleon) << " " 
	// 	<< L_rod*(hc/Mnucleon) << " " << L_slab*(hc/Mnucleon) << " "
	// 	<< iPlot << " " << sigma
	//   << endl;

	rod.setLengths(rad_rods,L_rod);
	slab.setLengths(L_slab, L_slab, 2.*rad_slab);
	
	double qmin=0.;
	double qmax=2.;
	int iqmax= 300;
	double q;
	// ofstream outform("data/form_factor"+to_string(rhoB*pow(Mnucleon/hc, 3.) )+"_L"+to_string(Lratio)+".txt");
	ofstream outform("data/form_factor"+to_string(rhoB*pow(Mnucleon/hc, 3.) )+".txt");
	// ofstream outform_ani("data/form_factor"+to_string(rhoB*pow(Mnucleon/hc, 3.) )+"_ani_L"+to_string(Lratio)+".txt");
	ofstream outform_ani("data/form_factor"+to_string(rhoB*pow(Mnucleon/hc, 3.) )+"_ani.txt");
	
	for(int it=0; it<=50; it++){
		double theta_= M_PI*(1. - it/50.);
		outform_ani << theta_ << " " 
								<< rod.getStructureFunction(electron.kf/2., cos(theta_)) << " "
								<< rod.getStructureFunction(electron.kf, cos(theta_)) << " "
								<< rod.getStructureFunction(2.*electron.kf, cos(theta_)) << " "
								<< slab.getStructureFunction(electron.kf/2., cos(theta_), 0.) << " "
								<< slab.getStructureFunction(electron.kf, cos(theta_), 0.) << " "
								<< slab.getStructureFunction(2.*electron.kf, cos(theta_), 0.) << " "
								<< slab.getStructureFunction(electron.kf/2., cos(theta_), M_PI/4.) << " " 
								<< slab.getStructureFunction(electron.kf, cos(theta_), M_PI/4.) << " " 
								<< slab.getStructureFunction(2.*electron.kf, cos(theta_), M_PI/4.)
								<<endl;
	}
	outform_ani.close();
	vector<double> qv, F3v2, Fav2, Fpv2;
	double coulIa=0., coulIp=0., nua=0., nup=0., nu_avg=0., nu_avg_inverse=0.;
	double nu3d, sigma0_3d=0.;

	for(int iq=0; iq<=iqmax; iq++){
		q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
		qv.push_back(q);
		
	}
	
		
  reverse(qv.begin(), qv.end());
	
	if(iDimension==2){	
		cout << " form factor calculation" << endl;
		Ze= lround((cluster.proton.density - gas.proton.density)*rod.getVolume());
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

			outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ 
						<< " " << f3_<< endl;
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
		Ze= lround((cluster.proton.density - gas.proton.density)*slab.getVolume());
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

			outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ 
					<< " " << f3_ << endl;
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

		Ze= lround((cluster.proton.density - gas.proton.density)*droplet.getVolume());
		for(int iq=0; iq<=iqmax; iq++){
			q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
			droplet.setMomentum(q);
			double f3_= droplet.getStructureFunction(q);
			F3v2.push_back(f3_*f3_);
			outform << q/electron.kf << " " << f3_*f3_ << " " << f3_*f3_ << " " << f3_*f3_ << endl;

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
	outform.close();
	double sigma0_a= pow(eGS, 2.)*electron.density/(nua*hypot(electron.kf, electron.mass));
	double sigma0_p= pow(eGS, 2.)*electron.density/(nup*hypot(electron.kf, electron.mass));	
	double sigma0_avg = pow(eGS, 2.)*electron.density*nu_avg_inverse/hypot(electron.kf, electron.mass);		
	
		//to do

	// Calculation with magnetic field:

	double	Bfield;
	double BgMin=1.e12;
	double lBgMin=log10(BgMin);
	double BgMax=1.e20; // Gauss
	double lBgMax=log10(BgMax);
	int ibMax=100;
	double dlogB= (lBgMax-lBgMin)/ibMax;
	vector<double> thetavec={0., M_PI/4., M_PI/2.};

	//Bfield *= 1.95e-14/pow(Mnucleon, 2.);  // Conversion factor Gauss to MeV^2 to adim
	for(auto thetab : thetavec){
		ofstream outTransport("data/transport_t"+to_string(thetab*180./M_PI)+
						"_L"+to_string(Lratio)+"_T"+to_string(temperature*Mnucleon)+"_"
						+to_string(iDimension)+".txt");
	
		for(int ib=0; ib<ibMax; ib++){
			Bfield= pow(10., lBgMax- ib*dlogB);
			// Bfield= BgMax - ib*(BgMax - BgMin)/ibMax; 
			cout << Bfield << endl;
			Bfield*=Gauss_to_Mev2_GS/pow(Mnucleon, 2.);
		cout << Bfield << endl;
		cout <<"--------------------------------" << endl;


			Matrix3d sigma3d, sigma_pasta;
			double bx= sin(thetab), bz=cos(thetab);
			double omega= eGS*Bfield/hypot(electron.kf, electron.mass);		
			double omega2=omega*omega, bx2= bx*bx, bz2=bz*bz;

			double Delta= nua*pow(nup, 2.)+ omega2*(bx2*nup+bz2*nua);

			sigma_pasta(0,0) = nua*nup+omega2*bx2;
			sigma_pasta(0,1) =-omega*bz*nua;
			sigma_pasta(0,2) = omega2*bx*bz;
			sigma_pasta(1,0) = omega*bz*nua;
			sigma_pasta(1,1) = nua*nup;
			sigma_pasta(1,2) = -omega*bx*nup;
			sigma_pasta(2,0) = omega2*bx*bz;
			sigma_pasta(2,1) = omega*bx*nup;
			sigma_pasta(2,2) = nup*nup+omega2*bz2;

			sigma_pasta*= pow(eGS, 2.)*electron.density/(Delta*hypot(electron.kf, electron.mass));

//	Average over domains:		
//	 filled with ( parallel, perpendicular, Hall)

			Vector3d sigma_avg_pasta(sigmaAverage_Parallel(nua, nup, omega, electron),
															sigmaAverage_Perpendicular(nua, nup, omega, electron),
															sigmaAverage_Hall(nua, nup, omega, electron));

			double xa, xp;
			xa=omega/nua;
			xp=omega/nup;	

			double x3d= omega/nu3d;
			Vector3d sigma_avg_3d(1., 1./(1+x3d*x3d), x3d/(1.+x3d*x3d) );
			sigma_avg_3d*=sigma0_3d;

			if(iDimension==2 && xp<1e-7){
				sigma_avg_pasta(0)= sigma0_avg;
				sigma_avg_pasta(1)= sigma0_avg;
				sigma_avg_pasta(2)= 0;
			}else if(iDimension==1 && xa<1e-7){
				sigma_avg_pasta(0)= sigma0_avg;
				sigma_avg_pasta(1)= sigma0_avg;
				sigma_avg_pasta(2)= 0;
			}

			const IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
  	  // cout << nua << " " << nup << " "  << omega << " " << xa << " " << xp << " "
			// 				<< Bfield*pow(Mnucleon, 2.)/Gauss_to_Mev2_GS << " "
			// 				<< Hfunc(nua, nup, omega) << " " 
			// 				<< (omega*omega+nup*nup)*(omega*omega+nup*nua)*Hfunc(nua, nup, omega) << " "  
			// 				<< - nup << " " 
			// 				<< sigma_avg_pasta(0)*Mnucleon/MeVto_Sec  << " " 
			// 				<< ( nup*nua*(omega*omega - nup*nup)*Hfunc(nua, nup, omega) + nup)/2. << " " 
			// 				<< sigma_avg_pasta(1)*Mnucleon/MeVto_Sec  << " " 
			// 				<< omega*(1. - nua*nup*nup*Hfunc(nua, nup, omega)) << " "
			// 				<< sigma_avg_pasta(2)*Mnucleon/MeVto_Sec  << " " 
			// 				<< nup/nua 
			// 				<< endl;

  	  outTransport << xa << " " << xp << " " 
							<< Bfield*pow(Mnucleon, 2.)/Gauss_to_Mev2_GS << " " << omega*Mnucleon << " " 
							<< sigma0_3d*Mnucleon/MeVto_Sec << " " //s-1
							<< sigma0_a*Mnucleon/MeVto_Sec << " " //s-1
							<< sigma0_p*Mnucleon/MeVto_Sec << " " //s-1
							<< sigma0_avg*Mnucleon/MeVto_Sec << " " //s-1
							<< sigma_pasta(0,0)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(0,1)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(0,2)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(1,0)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(1,1)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(1,2)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(2,0)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(2,1)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_pasta(2,2)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_avg_pasta(0)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_avg_pasta(1)*Mnucleon/MeVto_Sec  << " " 
							<< sigma_avg_pasta(2)*Mnucleon/MeVto_Sec  << " " 
							<< nup/nua << " " 
							<< iDimension << " " << iPlot
							<< endl;

							// 							<< rhoB*pow(Mnucleon/hc, 3.) << " " 
							// << Lratio << " " << RdM(2, minCol) << " " << RdM(minRow, minCol) << " "
							// << coulInt3d << " " << coulIa << " " << coulIp << " "
 							// << nu3d*Mnucleon/MeVto_Sec << " " //s-1
 							// << nua*Mnucleon/MeVto_Sec << " " //s-1
 							// << nup*Mnucleon/MeVto_Sec << " " //s-1
 							// << nup/nua << " " 

		}
	outTransport.close();

	}
		}else{
			cout << "no pasta: " << rhoB*pow(Mnucleon/hc, 3.) << " " << f << endl;
		}

	// outGlobal.close();
	// outFree.close();
	// outSol.close();
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