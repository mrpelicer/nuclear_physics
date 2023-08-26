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
	
	double rhoB, Yp, temperature;
	
	//, Bfield;, thetab=0.;and magnetic field (G):

	cout << "Specify the parametrization and temperature (MeV)  " << endl;
	cin >> parametrization >> temperature ;
	cout << "You chose: " << endl
				<< parametrization  << " parametrization " << endl 
				<< "T = " << temperature << " MeV" <<  endl;
				// << "B = " << Bfield << " G" <<  endl;

	double Lratio;
	// cout << "Specify L_d/R_d: " << endl;
	// cin >> Lratio ;

	//Declare nuclear matter: pasta and gas:
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class pasta(cluster, gas);

	// nlwm_class hmg(parametrization);
	double rhoBMax=0.09*pow(hc/Mnucleon, 3);
	double rhoBMin=0.04*pow(hc/Mnucleon, 3);

	int iRhoMax=300;
	double dRho= (rhoBMax-rhoBMin)/iRhoMax;

	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff= electron.mass;
	electron.Q=-1;
 
	// Bfield*=G auss_to_Mev2_GS/pow(Mnucleon, 2.);
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

	ofstream outTransport("data/transport_L"+to_string(Lratio)+"_T"+to_string(temperature)+".txt");


	temperature*=1./Mnucleon;
	electron.temperature=temperature;
	double alpha, dim;

	for(int irho=0; irho<=iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
		
		// pasta.solveCPA(rhoB, Yp, temperature);		
		pasta.solveCPA_betaEq(rhoB, temperature, electron);
		//solve hmg eos to see if pasta is in this density 
		// particle electron_hmg=electron;
		// hmg.setEOS_betaEq(rhoB, temperature, electron);

		f=pasta.f;
		Yp=pasta.YpG;
		sigma= getSurfaceTension(cluster, Yp, temperature);

		if(pasta.f>=0. && pasta.f<=1. && sigma>0){
	
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

		// double freeEn_hmg= hmg.getEnergy() + electron_hmg.energy 
		// 					  - temperature*(hmg.getEntropy() + electron_hmg.entropy) ;
		// 	cout << "freeen: " << FreeEn << " " << freeEn_hmg << endl;
		// if(FreeEn <= freeEn_hmg){
			double ri= pow(3./(4.*M_PI*cluster.rhoB), 1./3. );
			double Gamma= pow(Ze*eGS, 2.)/(ri*temperature);

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
		
	//Newton's calculation
			// double b2d= 1.5*alpha*coulEnM(1, minCol);
			// double c2d= pow(10., 2.1*(alpha-0.3))*alpha*coulEnM(1, minCol);
			// double k3 = 0.0655*alpha*coulEnM(1, minCol)*pow(RwM(1, minCol), 2.); 
			// double lambda= sqrt(2.*k3/(b2d+2.*c2d));
			//double L_rod=2.*RwM(1, minCol)
							// *sqrt( (b2d + 2.*c2d)*(RwM(1, minCol) - RdM(1, minCol))*sqrt(M_PI*lambda*a_)
											// /(temperature)	);
			// Lratio=0;

			// double r_ = 1e3*Mnucleon/hc;
			// double b1d=6.*alpha*coulEnM(0, minCol);
			// double k1=2.*alpha*coulEnM(0, minCol)*(1.+2.*alpha-2.*alpha*alpha)*RwM(0, minCol)*RwM(0, minCol)/15.;
			// double L_slab= 2.*RwM(0, minCol)
			// 				*sqrt(4.*M_PI*RwM(0, minCol) - RdM(0, minCol)*sqrt(b1d*k1)/(temperature*log(r_/(2.*RwM(0, minCol))) ) );


	//Fixed by cell size
			// double L_rod	=2.*RwM(1,minCol);
			// double L_slab	=1.*RwM(0, minCol);


	//deGenne's estimate of thermal coherence:
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
			double L2w= droplet.getVolume()/(M_PI*rad_rods*rad_rods);	
			double epsilon = pow(2, -eta);
			
			double L_slab	=L1w/pow(epsilon, 1/(2*eta));
			double L_rod	=L_slab;

			double Lsize=0;

			if(iDimension==1){Lsize=L_slab;}
			if(iDimension==2){Lsize=L_rod;}
			if(iDimension==3){Lsize=NAN;}


			double Z10=0., Z100=0., Z1000=0., Zg=0.;
			if(iDimension==3){
				Z10	=   (cluster.proton.density - gas.proton.density)*4.*M_PI*pow(RdM(2, minCol), 3)/3.;
				Z100	= lround(Z10);
				Z1000 = lround(Z10);
				Ze		= lround(Z10);
				Zg		= lround(gas.proton.density*4.*M_PI*pow(RwM(2, minCol), 3)/3.);
			}else  if(iDimension==2){
				Z10		=  lround((cluster.proton.density - gas.proton.density)*M_PI*pow(rad_rods, 2.)*10.*rad_rods);
				Z100	=  lround((cluster.proton.density - gas.proton.density)*M_PI*pow(rad_rods, 2.)*100.*rad_rods);
				Z1000	=  lround((cluster.proton.density - gas.proton.density)*M_PI*pow(rad_rods, 2.)*1000.*rad_rods);
				Ze		=  lround((cluster.proton.density - gas.proton.density)*M_PI*pow(rad_rods, 2.)*L_rod);
				Zg		=  lround(gas.proton.density*M_PI*pow(RwM(1, minCol), 2.)*L_rod);
			}else if(iDimension==1){
				Z10		= lround( (cluster.proton.density - gas.proton.density)*2.*rad_slab*pow(10.*rad_slab, 2.));
				Z100	= lround( (cluster.proton.density - gas.proton.density)*2.*rad_slab*pow(100.*rad_slab, 2.));
				Z1000	= lround( (cluster.proton.density - gas.proton.density)*2.*rad_slab*pow(1000.*rad_slab, 2.));
				Ze		= lround( (cluster.proton.density - gas.proton.density)*2.*rad_slab*pow(L_slab, 2.));
				Zg		= lround( gas.proton.density*2.*RwM(0, minCol)*pow(L_slab, 2.));
			}else{
				cout << "error in the dimensionality" << endl;
				return 1;
			}

			outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " " 
				<< Pressure*Mnucleon*pow(Mnucleon/hc, 3.) << " " 
			  << (FreeEn - 1.)*Mnucleon << " " 
				<< Ze << " " << Ae  << " " 
			  << Gamma << " " << Yp << " "<< f << " " 
				<< Rd*(hc/Mnucleon) << " " << Rw*(hc/Mnucleon) << " " 
				<< L_rod*(hc/Mnucleon) << " " << L_slab*(hc/Mnucleon) << " "
				<< iPlot << " " << sigma << " " 
				<< Z10 << " " << Z100 << " " << Z1000 << " " << Lsize*(hc/Mnucleon) << " "
				<< xi*(hc/Mnucleon) << " " << eta  << " " << epsilon << " " << pow(epsilon, 1/(2*eta)) << " "
				<< c0v_[0]*Mnucleon*pow(Mnucleon/hc, 3.) << " " << c0v_[1]*Mnucleon*pow(Mnucleon/hc, 3.) << " "
				<< av_[0]*(hc/Mnucleon)  << " "<<  av_[1]*(hc/Mnucleon)  << " "
				<< lambdav_[0]*(hc/Mnucleon)  << " " << lambdav_[1]*(hc/Mnucleon) << " "
				<<  Lsize/Rd << " " << Lsize/Rw << " " << Lsize/L1w << " " 
				<< Lsize/L2w << " " << RdM(2,0) << " " <<  RdM(1,0) << " " << RdM(0,0) 
				<< endl;


	// 		rod.setLengths(rad_rods,L_rod);
	// 		slab.setLengths(L_slab, L_slab, 2.*rad_slab);


	// 		double qmin=0.;
	// 		double qmax=2.;
	// 		int iqmax= 200;
	// 		double q;
	// 		vector<double> qv, Fav2, Fpv2;

	// 		for(int iq=0; iq<=iqmax; iq++){
	// 			q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
	// 			qv.push_back(q);
	// 		}

	//     reverse(qv.begin(), qv.end());
	// 		double coulIa=0., coulIp=0., nua=0., nup=0., nu_avg=0., nu_avg_inverse=0.;

	// 		if(iDimension==2){	
	// 			cout << "2d form factor calculation" << endl;
	// 			Ze= lround((cluster.proton.density - gas.proton.density)*rod.getVolume());

	// 			for(int iq=0; iq<=iqmax; iq++){
	// 				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
	// 				rod.setMomentum(q);
	// 				double fa_= rod.getStructureFunction2_Axial(q);
	// 				double fp_= rod.getStructureFunction2_Trans(q);
	// 				Fav2.push_back(fa_);
	// 				Fpv2.push_back(fp_);
	// 				// outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ << endl;
	// 			}

	//     	reverse(Fav2.begin() , Fav2.end());
	// 			reverse(Fpv2.begin() , Fpv2.end());
	// 			rod.setMomentumVec(qv);
	// 			rod.setFormFactor2Vec(Fav2);

	// 			cout << "2d coulomb log calculation" << endl;

	// 			// coulIa=   		 integrate_coulomb(coulomb_gsl, &rod);
	// 			coulIa=   		 getCoulombIntegral(rod);

	// 			rod.setFormFactor2Vec(Fpv2);
	// 			// coulIp=   		 integrate_coulomb(coulomb_gsl, &rod);
	// 			coulIp=   		 getCoulombIntegral(rod);


	// 			nua	= 3.*getFrequency(Ze, cluster.rhoB, coulIa, electron);
	// 			nup	= 3.*getFrequency(Ze, cluster.rhoB, coulIp, electron);

	// 			nu_avg= (nua + 2.*nup)/3.;
	// 			nu_avg_inverse= (2./nup + 1./nua)/3.;

	// 		}else if(iDimension==1){
	// 			cout << "1d form factor calculation" << endl;

	// 			Ze= lround((cluster.proton.density - gas.proton.density)*slab.getVolume());

	// 			for(int iq=0; iq<=iqmax; iq++){
	// 				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
	// 				slab.setMomentum(q);
	// 				double fa_= slab.getStructureFunction2_Axial(q);
	// 				double fp_= slab.getStructureFunction2_Trans(q);
	// 				Fav2.push_back(fa_);
	// 				Fpv2.push_back(fp_);
	// 				// outform << q/electron.kf << " " << fa_ << " " << fp_ << " " << fa_ +2.*fp_ << endl;
	// 			}

	//     	reverse(Fav2.begin() , Fav2.end());
	// 			reverse(Fpv2.begin() , Fpv2.end());
	// 			slab.setMomentumVec(qv);
	// 			cout << "1d coulomb log calculation" << endl;

	// 			slab.setFormFactor2Vec(Fav2);
	// 			// coulIa=   		 integrate_coulomb(coulomb_gsl, &slab);
	// 			coulIa=   		 getCoulombIntegral(slab);

	// 			slab.setFormFactor2Vec(Fpv2);
	// 			// coulIp=   		 integrate_coulomb(coulomb_gsl, &slab);
	// 			coulIp=   		 getCoulombIntegral(slab);


	// 			nua	= 3.*getFrequency(Ze, cluster.rhoB, coulIa, electron);
	// 			nup	= 3.*getFrequency(Ze, cluster.rhoB, coulIp, electron);

	// 			nu_avg				= (nua + 2.*nup)/3.;
	// 			nu_avg_inverse= (1./nua + 2./nup)/3.;

	// 		}else{

	// 			vector<double> F3v2;
	// 			Ze= lround((cluster.proton.density - gas.proton.density)*droplet.getVolume());

	// 			for(int iq=0; iq<=iqmax; iq++){
	// 				q= (qmax - (qmax-qmin)*iq/iqmax)*electron.kf;
	// 				droplet.setMomentum(q);
	// 				double f3_= droplet.getStructureFunction(q);
	// 				F3v2.push_back(f3_*f3_);
	// 			}


	// 	    reverse(F3v2.begin() , F3v2.end());
	// 			droplet.setMomentumVec(qv);
	// 			droplet.setFormFactor2Vec(F3v2);
	//     	// double coulInt3d=  integrate_coulomb(coulomb_gsl, &droplet);
	// 			cout << "3d coulomb log calculation" << endl;
	// 			double coulInt3d=getCoulombIntegral(droplet);
	//  			double nu3d			=  getFrequency(Ze, cluster.rhoB, coulInt3d, electron);

	// 			// double sigma0_3d=  pow(eGS, 2.)*electron.density/(nu3d*hypot(electron.kf, electron.mass));


	// 			coulIa=   		 coulInt3d/3.;
	// 			coulIp=   		 coulInt3d/3.;
	// 			nua	= nu3d;
	// 			nup	= nu3d;
	// 			nu_avg= nu3d;
	// 			nu_avg_inverse= 1./nu3d;

	// 		}



	// 		double sigma0_a= pow(eGS, 2.)*electron.density/(nua*hypot(electron.kf, electron.mass));
	// 		double sigma0_p= pow(eGS, 2.)*electron.density/(nup*hypot(electron.kf, electron.mass));	
	// 		double sigma0_avg = pow(eGS, 2.)*electron.density*nu_avg_inverse/hypot(electron.kf, electron.mass);		

	// 		(void) nu_avg;

	// 		double kappa0_a= sigma0_a*pi2*temperature/(3.*pow(eGS, 2.));
	// 		double kappa0_p= sigma0_p*pi2*temperature/(3.*pow(eGS, 2.));
	// 		double kappa0_avg = sigma0_avg*pi2*temperature/(3.*pow(eGS, 2.));


	// // 	// Calculation with magnetic field:
	// // 		Matrix3d sigmaB_2d, sigmaB_1d;
	// // 		double bx= sin(thetab), bz=cos(thetab);
	// // 		double omega= eGS*Bfield/hypot(electron.kf, electron.mass);		

	// // 		double omega2=omega*omega, bx2= bx*bx, bz2=bz*bz;

	// // 		double Delta2d= nua2d*pow(nup2d, 2.)+ omega2*bx2*nup2d+omega2*bz2*nua2d;
	// // 		double Delta1d= nua1d*pow(nup1d, 2.)+ omega2*bx2*nup1d+omega2*bz2*nua1d;

	// // 		sigmaB_2d(0,0) = nua2d*nup2d+omega2*bx2;
	// // 		sigmaB_2d(0,1) =-omega*bz*nua2d;
	// // 		sigmaB_2d(0,2) = omega2*bx*bz;
	// // 		sigmaB_2d(1,0) = omega*bz*nua2d;
	// // 		sigmaB_2d(1,1) = nua2d*nup2d;
	// // 		sigmaB_2d(1,2) = -omega*bx*nup2d;
	// // 		sigmaB_2d(2,0) = omega2*bx*bz;
	// // 		sigmaB_2d(2,1) = omega*bx*nup2d;
	// // 		sigmaB_2d(2,2) = nup2d*nup2d+omega2*bz2;
	// // 		sigmaB_2d*= pow(eGS, 2.)*electron.density/(Delta2d*hypot(electron.kf, electron.mass));

	// // 		sigmaB_1d(0,0) = nua1d*nup1d+omega2*bx2;
	// // 		sigmaB_1d(0,1) =-omega*bz*nua1d;
	// // 		sigmaB_1d(0,2) = omega2*bx*bz;
	// // 		sigmaB_1d(1,0) = omega*bz*nua1d;
	// // 		sigmaB_1d(1,1) = nua1d*nup1d;
	// // 		sigmaB_1d(1,2) = -omega*bx*nup1d;
	// // 		sigmaB_1d(2,0) = omega2*bx*bz;
	// // 		sigmaB_1d(2,1) = omega*bx*nup1d;
	// // 		sigmaB_1d(2,2) = nup1d*nup1d+omega2*bz2;
	// // 		sigmaB_1d*= pow(eGS, 2.)*electron.density/(Delta1d*hypot(electron.kf, electron.mass));

	// // //Average over domains:		
	// // // filled with ( parallel, perpendicular, Hall)
	// // 		double x3d= omega/nu3d;
	// // 		Vector3d sigma_avg_3d(1., 1./(1+x3d*x3d), x3d/(1.+x3d*x3d) );
	// // 		sigma_avg_3d*=sigma0_3d;

	// // 		Vector3d sigma_avg_2d(sigmaAverage_Parallel(nua2d, nup2d, omega, electron),
	// // 													sigmaAverage_Perpendicular(nua2d, nup2d, omega, electron),
	// // 													sigmaAverage_Hall(nua2d, nup2d, omega, electron)); 

	// // 		Vector3d sigma_avg_1d(sigmaAverage_Parallel(nua1d, nup1d, omega, electron),
	// // 											sigmaAverage_Perpendicular(nua1d, nup1d, omega, electron),
	// // 											sigmaAverage_Hall(nua1d, nup1d, omega, electron)); 

	// // 		double xa2d, xp2d, xa1d, xp1d;
	// // 		xa2d=omega/nua2d;
	// // 		xp2d=omega/nup2d;	
	// // 		xa1d=omega/nua1d;
	// // 		xp1d=omega/nup1d;
	// 		// const IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

	//     outTransport << rhoB*pow(Mnucleon/hc, 3.) << " " << Lratio << " " 
	// 						<< RdM(2, minCol)*hc/Mnucleon << " " << RdM(minRow, minCol)*hc/Mnucleon << " "
	// 					  << RwM(2, minCol)*hc/Mnucleon << " " << RwM(minRow, minCol)*hc/Mnucleon << " "
	// 						<< Lsize*hc/Mnucleon << " "
	// 						<< coulIa << " " << coulIp << " "
	//  						<< nua*Mnucleon/MeVto_Sec << " " //s-1
	//  						<< nup*Mnucleon/MeVto_Sec << " " //s-1
	//  						<< nup/nua << " " 
	// 						<< sigma0_a*Mnucleon/MeVto_Sec << " " //s-1
	// 						<< sigma0_p*Mnucleon/MeVto_Sec << " " //s-1
	// 						<< sigma0_avg*Mnucleon/MeVto_Sec << " " //s-1
	// 						<< iDimension << " " << iPlot << " " 
	// 						<< kappa0_a*pow(Mnucleon, 2.)*JouletoErg* //erg s-1 cm-1 K-1
	// 																								(kBoltz/(MeVto_Sec*MeVto_Cm))	<< " " 
	// 						<< kappa0_p*pow(Mnucleon, 2.)*JouletoErg* //erg s-1 cm-1 K-1
	// 																								(kBoltz/(MeVto_Sec*MeVto_Cm))	<< " "
	// 						<< kappa0_avg*pow(Mnucleon, 2.)*JouletoErg* //erg s-1 cm-1 K-1
	// 																								(kBoltz/(MeVto_Sec*MeVto_Cm)) 
	// 						<< endl;
		

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