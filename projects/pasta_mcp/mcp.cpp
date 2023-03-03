// Pasta with non linear Walecka model using Mean Field Theory

#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_walecka.h"
#include "../../include/pasta.h"
#include "../../include/interpolator.h"
#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

double getZeff(double Z3, int idim_, int iAxis, double R3, double Ri);
//===========================main code =============================

int main(){

	std::string parametrization= "iufsu";

	//Declare nuclear matter: pasta and gas:
	nlwm_class cluster(parametrization);
	nlwm_class gas(parametrization);
	pasta_class sna(cluster, gas);
		
	double rhoB;
	double rhoBMax=0.06*pow(hc/Mnucleon, 3);
	int iRhoMax=1;
	double dRho= rhoBMax/iRhoMax;
	double Yp=0.1;
	double temperature=3.;
	
	//Declare variables for SNA:
	Eigen::VectorXd freeEn_ocp(3), freeEn_gas(3), Rd_ocp(3), Rw_ocp(3), 
									f_ocp(3), interfaceEn_ocp(3), Vn_ocp(3), munv(3), mupv(3), rearv(3),
									rhop_c(3),rhon_c(3), omega0v(3), rhop_g(3), rhon_g(3);
	
	double vn_ocp, vw_ocp, A_ocp, Z_ocp, omega0_ocp, omega_ocp;
	
	int id_ocp, iPlot=0;	
	
	double dim_ocp;
	
	particle electron;
	electron.mass= Me/Mnucleon;


	// saveSol=firstGuess;
	std::cout << "Yp, T= " << Yp << " " << temperature << std::endl;
	std::ofstream outGlobal("data/flu_"+parametrization+"_yp"+std::to_string(Yp)+
												   +"_T"+std::to_string(temperature)+".txt");

	outGlobal << "#rho_b freeEn Rd Rw A Z IP" << std::endl;


	temperature*=1./Mnucleon;
	electron.temperature=temperature;
	
	for(int irho=0; irho<iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
	
		electron.density=Yp*rhoB;
		electron.kf=pow(3.*pi2*electron.density, 1./3.);
		electron.solveChemPotEff();
		electron.chemPot = electron.chemPot_eff;
		electron.calculateProperties();
		
		double freeEnEl = electron.energy - temperature*electron.entropy;
		
		double dim;

		for(int iDim=2; iDim>=0; iDim--){
		dim = (double) (iDim)+1;
						
			sna.solveCLD(rhoB, Yp, temperature, dim, 0);
			f_ocp(iDim)= sna.f;
	
			double sigma= getSurfaceTension(cluster, Yp, temperature);
			Rd_ocp(iDim)= getRadiusD(dim, f_ocp(iDim), Yp, cluster, gas);
			Rw_ocp(iDim)= Rd_ocp(iDim)/pow(f_ocp(iDim), 1./dim);
			Vn_ocp(iDim) = 4.*M_PI*pow(getRadiusD(3., f_ocp(iDim), Yp, cluster, gas), 3.)/3.;
			
			freeEn_gas(iDim)= gas.getEnergy() - temperature*gas.getEntropy();
			interfaceEn_ocp(iDim)= 3.*sigma*dim/(2.*Rd_ocp(iDim));
			
			double energy = f_ocp(iDim)*cluster.getEnergy()  + (1.-f_ocp(iDim))*gas.getEnergy()
										+ f_ocp(iDim)*interfaceEn_ocp(iDim) + electron.energy;
										
			double entropy= f_ocp(iDim)*cluster.getEntropy() +(1.-f_ocp(iDim))*gas.getEntropy() 
											+ electron.entropy;
			
			freeEn_ocp(iDim) = energy-temperature*entropy;
			
			mupv(iDim) = gas.proton.chemPot+ (f_ocp(iDim)/(1.-f_ocp(iDim)))*	
				(-2.*interfaceEn_ocp(iDim)/(3.*(cluster.proton.density- gas.proton.density))
				+dim*getSurfaceTensionDerivative(cluster, Yp, temperature)
						*(1.-f_ocp(iDim))*(1.-Yp)/(rhoB*Rd_ocp(iDim))
				);

			munv(iDim) = gas.neutron.chemPot -f_ocp(iDim)/(1.-f_ocp(iDim))
			*dim*getSurfaceTensionDerivative(cluster, Yp, temperature)*(1.-f_ocp(iDim))
						*Yp/(rhoB*Rd_ocp(iDim));

			rearv(iDim)=f_ocp(iDim)*(interfaceEn_ocp(iDim)*getPhiFuncDerivative(dim,f_ocp(iDim))
			/(3.*getPhiFunc(dim, f_ocp(iDim)))
			+ dim*getSurfaceTensionDerivative(cluster, Yp, temperature)*(
				(cluster.proton.density - gas.proton.density)
				-Yp*(cluster.rhoB - gas.rhoB))/(rhoB*Rd_ocp(iDim)))
						/(cluster.proton.density - gas.proton.density);
			rhop_c(iDim) = cluster.proton.density;
			rhon_c(iDim) = cluster.neutron.density;
			rhop_g(iDim) = gas.proton.density;
			rhon_g(iDim) = gas.neutron.density;
			omega0v(iDim)=Vn_ocp(iDim)*(cluster.getEnergy() -gas.getEnergy() + interfaceEn_ocp(iDim)
										- temperature*(cluster.getEntropy() -gas.getEntropy() )
										- mupv(iDim)*(cluster.proton.density - gas.proton.density)
										- munv(iDim)*(cluster.neutron.density- gas.neutron.density));
		// 			std::cout << "Test00: "
		// << rhoB*pow(Mnucleon/hc, 3.) << " " << cluster.neutron.density*pow(Mnucleon/hc, 3.) << " " << cluster.proton.density*pow(Mnucleon/hc, 3.)
		// << std::endl;
		
//						std::cout <<rhoB*pow(Mnucleon/hc, 3.) << " " << cluster.Ye << " "
//						<< (freeEn_ocp(iDim)/rhoB - 1.)*Mnucleon << std::endl;
						
		}	
		

	//Get phase that minimizes energy:
		Eigen::MatrixXd::Index minRow;
		
		double freeEn_c= freeEn_ocp.minCoeff(&minRow);
		id_ocp=minRow;
		
		//double freeEn_ocp= freeEn_ocp(2);
		//id_ocp=0;
		dim_ocp=(double)(id_ocp+1);
					
		// saveSol[0]=nup1_ocp(id_ocp);
		// saveSol[1]=nun1_ocp(id_ocp);
		// saveSol[2]=Mef1_ocp(id_ocp);
		// saveSol[3]=nup2_ocp(id_ocp);
		// saveSol[4]=nun2_ocp(id_ocp);
		// saveSol[5]=Mef2_ocp(id_ocp);
		
		// cluster.setEOS_coexistence(nup1_ocp(id_ocp), nun1_ocp(id_ocp), Mef1_ocp(id_ocp));
		// gas.setEOS_coexistence(nup2_ocp(id_ocp), nun2_ocp(id_ocp), Mef2_ocp(id_ocp));
		
		vn_ocp= Vn_ocp(id_ocp);
		vw_ocp= vn_ocp/f_ocp(id_ocp);
		
		A_ocp= vn_ocp*(rhop_c(id_ocp) + rhon_c(id_ocp) - gas.rhoB);
		Z_ocp= vn_ocp*(rhop_c(id_ocp) - gas.proton.density);
				
//---------- FLUCTUATION ---------------------:
	
//External fields and rearrangement:
		double LambdaP = mupv(id_ocp);
		// gas.proton.chemPot+ (f_ocp(id_ocp)/(1.-f_ocp(id_ocp)))*
		// 		(-2.*interfaceEn_ocp(id_ocp)/(3.*(cluster.proton.density- gas.proton.density))
		// 		+dim_ocp*getSurfaceTensionDerivative(cluster, Yp, temperature)
		// 				*(1.-f_ocp(id_ocp))*(1.-Yp)/(rhoB*Rd_ocp(id_ocp))
		// 		);
								 
		double LambdaN = munv(id_ocp);
		// gas.neutron.chemPot
		// -f_ocp(id_ocp)/(1.-f_ocp(id_ocp))
		// 	*dim_ocp*getSurfaceTensionDerivative(cluster, Yp, temperature)*(1.-f_ocp(id_ocp))
		// 				*Yp/(rhoB*Rd_ocp(id_ocp));	
		
		double rear_avg=rearv(id_ocp);
		omega0_ocp=  omega0v(id_ocp);
										
		omega_ocp= omega0_ocp+rear_avg*Z_ocp;
		// cout <<dim_ocp << " " << A_ocp << " " << Z_ocp << " " << omega0_ocp << " " << omega_ocp << endl;
		double Vtot=0.;
		
		if(id_ocp==0){iPlot=3;}
		if(id_ocp==1){iPlot=2;}
		if(id_ocp==2){iPlot=1;}
		
		std::cout << "Test: "
		<< rhoB*pow(Mnucleon/hc, 3.) << " " 
		<< A_ocp << " "  << Z_ocp << " "
		<< dim_ocp << " " << omega_ocp << " " << omega0_ocp << " " 
		<< LambdaP << " " << LambdaN << " " 
		<< rhop_c(id_ocp)*pow(Mnucleon/hc, 3.) << " " << rhon_c(id_ocp)*pow(Mnucleon/hc, 3.)
		<< " " << rear_avg << " " << gas.proton.density*pow(Mnucleon/hc, 3.) << " "
		 << gas.proton.density*pow(Mnucleon/hc, 3.) <<" "
		 << rhop_g(id_ocp)*pow(Mnucleon/hc, 3.) << " " << rhon_g(id_ocp)*pow(Mnucleon/hc, 3.)
		<< std::endl;
		double rhob_g=  rhop_g(id_ocp)+ rhon_g(id_ocp);

		if(isnan(LambdaP) ){
			cout << "lambdap problem:";
			return 1;
		}
		if(isnan(LambdaN) ){
			cout << "lambdan problem:";
			return 1;
		}
		if(isnan(rear_avg) ){
			cout << "rear problem:";
			return 1;
		}
	if( (f_ocp(id_ocp) >0.) && (f_ocp(id_ocp)<1.) ){
		
		nlwm_class clusterN(parametrization);
		
		double Yp_cell, uN, sigmaN;

		double fracpP=1.2;	
		double fracpM=0.8;
		int const np=10;

		double fracnP=1.2;
		double fracnM=0.8;
		int const nn=10;


		Eigen::MatrixXd rhobNT(nn, np), rhonNT(nn, np), rhopNT(nn, np), omega_mpd(nn, np);
			
		Eigen::TensorFixedSize<double, Eigen::Sizes<nn, np, 3>> 
																	RdNT, RwNT, VNT, uNT,
																	F0NT, FNT, MuNT, rearNT,
																	GNT0, GNT, ANT, ZNT,					
																	PressureNT, G0densNT, GdensNT,
																	NN, N0N, pN, p0N, nN, n0N,
																	Zeff_xNT, Zeff_yNT, Zeff_zNT;			
		
		double norma=0., norma0=0.;
		VectorXd distDimN(3), distDimN0(3), normaV0(3), normaV(3);


		for(int id=0; id<3; id++){
			normaV(id)=0.;
			distDimN(id)=0.;
			normaV0(id)=0.;
			distDimN0(id)=0.;
		}

		for(int in=nn-1; in>=0; in--){
			for(int ip=np-1; ip>=0; ip--){
			
		
			rhonNT(in, ip) = (fracnP + (fracnM - fracnP)*in/nn)*rhon_c(id_ocp);
			rhopNT(in, ip) = (fracpP + (fracpM - fracpP)*ip/np)*rhop_c(id_ocp);


			rhobNT(in, ip) = rhonNT(in, ip) + rhopNT(in, ip);

			clusterN.setEOS_nucleons(rhobNT(in, ip), rhopNT(in, ip)/rhobNT(in, ip), temperature);
			
			uN= (Yp*rhoB -  rhop_g(id_ocp))/(clusterN.proton.density -   rhop_g(id_ocp));
			
			Yp_cell=  ( uN*(rhopNT(in, ip) -  rhop_g(id_ocp)) +   rhop_g(id_ocp) )
											/( uN*(rhobNT(in, ip) - rhob_g) +rhob_g);
			
			sigmaN	= getSurfaceTension(clusterN, Yp_cell, temperature);


			double dimN;	
			for(int iDim=2; iDim>=0; iDim--){
				dimN = (double) (iDim)+1; 
				
				uNT(in, ip, iDim)= uN;
				RdNT(in, ip, iDim)= getRadiusD(dimN, uN, Yp_cell, clusterN, gas);
				RwNT(in, ip, iDim)= RdNT(in, ip, iDim)/pow(uN, 1./dimN);
				
				VNT(in, ip, iDim)= 4.*M_PI*pow(getRadiusD(3., uN, Yp_cell, clusterN, gas), 3.)/3.;
								
				ANT(in, ip, iDim) = (clusterN.rhoB - 	rhob_g)*VNT(in, ip, iDim);
				ZNT(in, ip, iDim) = (clusterN.proton.density -  rhop_g(id_ocp))*VNT(in, ip, iDim);
				
				rearNT(in, ip, iDim) =ZNT(in, ip, iDim)*rear_avg;
				
				Zeff_xNT(in, ip, iDim)= getZeff(ZNT(in, ip, iDim), iDim, 1, 
																			getRadiusD(3., uN, Yp_cell, clusterN, gas), 
																			RdNT(in, ip, iDim) );
																			
				Zeff_yNT(in, ip, iDim)= getZeff(ZNT(in, ip, iDim), iDim, 2, 
																			getRadiusD(3., uN, Yp_cell, clusterN, gas), 
																			RdNT(in, ip, iDim) );		
				
				Zeff_zNT(in, ip, iDim)= getZeff(ZNT(in, ip, iDim), iDim, 3, 
																			getRadiusD(3., uN, Yp_cell, clusterN, gas), 
																			RdNT(in, ip, iDim) );
			
				double interface=	3.*sigmaN*dimN/(2.*RdNT(in, ip, iDim));
				
				F0NT(in, ip, iDim) =  VNT(in, ip, iDim)*
															(clusterN.getEnergy() - temperature*clusterN.getEntropy()
															+interface
															- freeEn_gas(id_ocp) );  
				
				FNT(in, ip, iDim)= F0NT(in, ip, iDim) + rearNT(in, ip, iDim);
				
				MuNT(in, ip, iDim) = VNT(in, ip, iDim)*
													( LambdaP*(clusterN.proton.density -   rhop_g(id_ocp))
													+ LambdaN*(clusterN.neutron.density-   rhon_g(id_ocp)));
								
				GNT0(in, ip, iDim)= F0NT(in, ip, iDim) - MuNT(in, ip, iDim);
				GNT(in, ip, iDim)= 	FNT(in, ip, iDim)  - MuNT(in, ip, iDim);


				double Cdist0= omega0_ocp/temperature;
				double Cdist= omega_ocp/temperature;
								cout << "i---" << endl;

				cout << rhobNT(in, ip)*pow(Mnucleon/hc, 3.) << " " <<  rhonNT(in, ip)*pow(Mnucleon/hc, 3.) << " " 
					<<  rhopNT(in, ip)*pow(Mnucleon/hc, 3.) << " " << Yp_cell << " " <<  MuNT(in, ip, iDim) << " " <<F0NT(in, ip, iDim)
					<< endl;
					cout << VNT(in, ip, iDim) << " " <<  clusterN.getEnergy()  <<  " " 
						<<  temperature*clusterN.getEntropy() << " "
						<< interface << " " <<  freeEn_gas(id_ocp)  << endl;
				cout << exp(- GNT0(in, ip, iDim)/temperature + Cdist0) << " " <<  GNT0(in, ip, iDim)  << " " << Cdist0 << endl;
				cout << exp(- GNT(in, ip, iDim)/temperature + Cdist) << " " <<  GNT(in, ip, iDim)  << " "	<< Cdist << endl;
				cout << "f---" << endl;
				if(uN>0 && uN<=1. ){
					N0N(in, ip, iDim)= exp(- GNT0(in, ip, iDim)/temperature + Cdist0)  ;
					NN(in, ip, iDim) = exp(- GNT(in, ip, iDim)/temperature  + Cdist)  ;
				}else{
					N0N(in, ip, iDim)= 0.;
					NN(in, ip, iDim) = 0.;
				}

				Vtot+=NN(in,ip, iDim)*VNT(in, ip, iDim)/uN;
				norma0+=N0N(in, ip, iDim);
				norma += NN(in, ip, iDim);

				normaV0(iDim)+=N0N(in, ip, iDim);
				normaV(iDim) += NN(in, ip, iDim);

				distDimN0(iDim)+=N0N(in, ip, iDim);
				distDimN(iDim)+=NN(in, ip, iDim);

			}

		}
		}	//end fluctuations
		
		
		double rhobbar=0., rhopbar= 0.;
		double Fbar=0., Abar=0., Zbar=0., Rdbar=0., Rwbar=0.;	
		double Zeff_xbar=0., Zeff_ybar=0., Zeff_zbar=0.;
		vw_ocp=Vtot/norma;	

		for(int in=0; in<nn; in++){
		for(int ip=0; ip<np; ip++){
		for(int iDim=2; iDim>=0; iDim--){
		
		p0N(in,ip,iDim)=N0N(in,ip,iDim)/norma0;		
		pN(in,ip,iDim)=NN(in,ip,iDim)/norma;
		nN(in,ip,iDim)=pN(in,ip,iDim)/vw_ocp;
		
		rhobbar += (rhobNT(in,ip) - rhob_g)*VNT(in, ip, iDim)*nN(in, ip, iDim);
		
		rhopbar += (rhopNT(in, ip) -   rhop_g(id_ocp))*VNT(in, ip, iDim)*nN(in, ip, iDim);
				
		Fbar+=nN(in, ip, iDim)*F0NT(in, ip, iDim);
		Rdbar+=pN(in, ip, iDim)*RdNT(in, ip, iDim);
		Rwbar+=pN(in, ip, iDim)*RwNT(in, ip, iDim);
		
		Abar+=pN(in, ip, iDim)*ANT(in, ip, iDim);
		Zbar+=pN(in, ip, iDim)*ZNT(in, ip, iDim);
		
		Zeff_xbar+=pN(in, ip, iDim)*Zeff_xNT(in, ip, iDim);
		Zeff_ybar+=pN(in, ip, iDim)*Zeff_yNT(in, ip, iDim);
		Zeff_zbar+=pN(in, ip, iDim)*Zeff_zNT(in, ip, iDim);
				
		}
		omega_mpd(in, ip)= GNT(in, ip, id_ocp);
		}
		}
		double Acell = Abar + vw_ocp*gas.rhoB;
		
		Eigen::MatrixXd::Index imp_n, imp_p;
		double omega_mpc= omega_mpd.minCoeff(&imp_n, &imp_p);
		(void)omega_mpc;
		std::cout << " Most probable cluster: " << imp_n << " " << imp_p << std::endl;
		double A_mp= ANT(imp_n, imp_p, id_ocp);
		double Z_mp= ZNT(imp_n, imp_p, id_ocp);
	
		
		rhobbar+=rhob_g;
		rhopbar+=  rhop_g(id_ocp);
		Fbar+=freeEn_gas(id_ocp)+freeEnEl;
	
		Eigen::VectorXd pastaProb(3);
		pastaProb(0)= distDimN(0)/norma;
		pastaProb(1)= distDimN(1)/norma;
		pastaProb(2)= distDimN(2)/norma;
		
		double p1= pastaProb(0);
		double p2= pastaProb(1);
		double p3= pastaProb(2);
	
		double Zvar=0.;
		double Zeff_xvar=0., Zeff_yvar=0., Zeff_zvar=0.;
		VectorXd ZvarD(3);
		ZvarD.setZero(3);
		
		std::cout << "Zero test: " << ZvarD(0) << " " << ZvarD(1) << " " << ZvarD(2) << std::endl;

		
		for(int iDim=2; iDim>=0; iDim--){
		for(int in=0; in<nn; in++){
		for(int ip=0; ip<np; ip++){
			Zvar+=pow( (ZNT(in, ip, iDim)-Zbar), 2.)*pN(in, ip, iDim);
			ZvarD(iDim) += pow( (ZNT(in, ip, iDim)-Zbar), 2.)*NN(in, ip, iDim);
			Zeff_xvar+=pow( (Zeff_xNT(in, ip, iDim) - Zeff_xbar), 2.)*pN(in, ip, iDim);
			Zeff_yvar+=pow( (Zeff_yNT(in, ip, iDim) - Zeff_ybar), 2.)*pN(in, ip, iDim);
			Zeff_zvar+=pow( (Zeff_zNT(in, ip, iDim) - Zeff_zbar), 2.)*pN(in, ip, iDim);
		}
		}
		}
		
		ZvarD(0)*= 1./distDimN(0);
		ZvarD(1)*= 1./distDimN(1);
		ZvarD(2)*= 1./distDimN(2);
		
		Eigen::MatrixXd::Index imaxPasta;
		double maxProb= pastaProb.maxCoeff(&imaxPasta);
		double Q_horowitz= (1.-maxProb)*Zbar*Zbar;
		
				
		std::cout << "(rho, maxP, A_mp, Z_mp)= "
							<< rhobbar*pow(Mnucleon/hc, 3.) << " " << maxProb << " "
							<< A_mp << " "  << Z_mp 
		<< std::endl;
					
		std::cout << "Z_bars's= " << Zbar << " " << Zeff_xbar  << " " 
							<< Zeff_ybar << " " << Zeff_zbar 
		<< std::endl;
		
		std::cout << "Z_vars's= " << Zvar << " " << Zeff_xvar  
															<< " " << Zeff_yvar << " " << Zeff_zvar 
		<< std::endl;

		
		
		for(int iDim=2; iDim>=0; iDim--){
	
			std::ofstream outdim("data/"+std::to_string(iDim+1)+"d_"+parametrization+"_yp"+std::to_string(Yp)+
																				 +"_T"+std::to_string(temperature*Mnucleon)+"_"
																				 +std::to_string(irho)+".txt");
			outdim << "#rhonN rhopN rhobN GN PN RdN VNT An Zn uN"<< std::endl;

				for(int in=0; in<nn; in++){
				for(int ip=0; ip<np; ip++){

				outdim << rhonNT(in, ip)/rhon_c(id_ocp) << " " 
							 << rhopNT(in, ip)/rhop_c(id_ocp) << " "
							 << rhobNT(in, ip)/(rhon_c(id_ocp)+rhop_c(id_ocp)) << " "
							 << GNT(in, ip, iDim)*Mnucleon << " "	<< pN(in, ip, iDim) << " " 
							 << GNT0(in, ip, iDim)*Mnucleon << " "	<< p0N(in,ip,iDim) << " "
							 << RdNT(in, ip, iDim)*(hc/Mnucleon) << " " << VNT(in, ip, iDim) <<  " "
							 << ANT(in, ip, iDim) << " " << ZNT(in, ip, iDim) << " " << uNT(in, ip, iDim) << " "
							 << nN(in, ip, iDim)
							 << std::endl;

				}
			}
					outdim.close();

		}

		std::cout << "Probs: "	<< p1 << " & " << p2 << " & " << p3 <<  std::endl;
		std::cout << "Probs0: "	<< distDimN0(0)/norma0 << " & " << distDimN0(1)/norma0 << " & " << distDimN0(2)/norma0 <<  std::endl;
	
		outGlobal << rhoB*pow(Mnucleon/hc, 3.) << " "
			  << (freeEn_c/rhoB - 1.)*Mnucleon << " "
			  << Rd_ocp(id_ocp)*(hc/Mnucleon) << " " 
				<< Rw_ocp(id_ocp)*(hc/Mnucleon) << " "
				<< A_ocp << " " << Z_ocp << " " << iPlot << " "
				<< p1  << "  " << p2  << "  " << p3 << " "
				<< rhobbar*pow(Mnucleon/hc, 3.) << " "  //col 11
				<< (Fbar/rhobbar - 1.)*Mnucleon << " "
				<< Rdbar*(hc/Mnucleon) << " " 
				<< Rwbar*(hc/Mnucleon) << " "
				<< Abar << " " << Zbar << " " //15, 16
				<< Acell << " " 
				<<   rhop_g(id_ocp)*pow(Mnucleon/hc, 3.) << " " <<   rhon_g(id_ocp)*pow(Mnucleon/hc, 3.)<< " "
				<< LambdaP << " " << LambdaN << " " 
				<< A_mp << " " << Z_mp << " " //22, 23
				<< gas.proton.chemPot << " " << gas.neutron.chemPot << " " 
				<< Zvar << " " << Q_horowitz << " " // 26, 27
				<< Zeff_xbar  << " " << Zeff_ybar << " " << Zeff_zbar << " "
				<< Zeff_xvar << " " << Zeff_yvar << " " << Zeff_zvar << " "
				<< ZvarD(0) << " " << ZvarD(1) << " " << ZvarD(2) << " " //34, 35, 36
			  << std::endl;


	}
	}

	outGlobal.close();

  return 0;
}

double getZeff(double Z3, int idim_, int iAxis, double R3, double Ri){

	double Si=0., Li=0.;
	double S3= 2.*M_PI*pow(R3, 2.);
	
	if(idim_==0){
		
		Li= sqrt(4.*M_PI*pow(R3, 3.)/(3.*Ri));
		if(iAxis==1){				//x-direction
			Si =Ri*Li;
		}else if(iAxis==2){ //y-direction
			Si =Ri*Li;
		}else if(iAxis==3){ //z-direction
			Si =Li*Li;
		}
	
	}else if(idim_==1){	
			
		Li= 4.*pow(R3, 3.)/(3.*pow(Ri, 2.));
		if(iAxis==1){				//x-direction
			Si =M_PI*pow(Ri, 2.);
		}else if(iAxis==2){ //y-direction
			Si =M_PI*Ri*Li;
		}else if(iAxis==3){ //z-direction
			Si =M_PI*Ri*Li;
		}
		
	}else if(idim_==2){
		Si=S3;
	}
	
	return Z3*Si/S3;
}