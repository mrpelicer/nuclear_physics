// Pasta with non linear Walecka model using Mean Field Theory

#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
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
	matter_class cluster(parametrization);
	matter_class gas(parametrization);
	phase_coexistence_class sna(cluster, gas);
		
	double rhoB;
	double rhoBMax=0.07*pow(hc/cluster.Mn, 3);
	int iRhoMax=35;
	double dRho= rhoBMax/iRhoMax;
	double Yp=0.1;
	double temperature=5.;
	
	//Declare variables for SNA:
	Eigen::VectorXd freeEn_sna(3), freeEn_gas(3), Rd_sna(3), Rw_sna(3), 
									f_sna(3), interfaceEn_sna(3), Vn_sna(3);
	Eigen::VectorXd nup1_sna(3), nun1_sna(3), nup2_sna(3), nun2_sna(3), Mef1_sna(3), Mef2_sna(3);
	
	double vn_ocp, vw_ocp, A_ocp, Z_ocp, omega0_ocp, omega_ocp;
	
	int id_ocp, iPlot=0;	
	
	double dim_ocp;
	
	//Guess solution for equilibrium variables: Eff_chemPots and Eff_masses
	std::vector<double> firstGuess={0.68, 0.68, 0.61, 0.98, 0.98, 0.99};
	std::vector<double> saveSol(6);

	if(Yp<=0.5 && Yp>0.4){
		firstGuess= { 0.682733, 0.68279, 0.626632, 0.989045, 0.980986, 0.999682};
	}else if(Yp<=0.4 && Yp>0.3){
		firstGuess= {0.680087, 0.694151, 0.632422, 0.97523, 0.992023, 0.999128};
		//firstGuess={0.669209, 0.70209, 0.629, 0.92, 0.97, 0.95};
		//firstGuess={0.70, 0.73, 0.65,0.92, 0.97, 0.95};
	}else if(Yp<=0.3){
		firstGuess={0.686883, 0.712638, 0.648993, 0.954912, 0.991327, 0.988764};
	}else if(Yp<=0.1){
		firstGuess={0.726012, 0.768844, 0.710641, 0.872514, 0.9208, 0.898116};
	}
	
	particle electron;
	electron.mass= Me/cluster.Mn;


	saveSol=firstGuess;
	std::cout << "Yp, T= " << Yp << " " << temperature << std::endl;
	std::ofstream outGlobal("../data/data_flu/flu_"+parametrization+"_yp"+std::to_string(Yp)+
												   +"_T"+std::to_string(temperature)+".txt");

	outGlobal << "#rho_b freeEn Rd Rw A Z IP" << std::endl;


	temperature*=1./cluster.Mn;
	electron.temperature=temperature;
	
	for(int idim=0; idim<3; idim++){
			nup1_sna(idim)= firstGuess[0];
			nun1_sna(idim)= firstGuess[1];
			Mef1_sna(idim)= firstGuess[2];
			nup2_sna(idim)= firstGuess[3];
			nun2_sna(idim)= firstGuess[4];
			Mef2_sna(idim)= firstGuess[5];
	}


	for(int irho=0; irho<iRhoMax; irho++){
		rhoB=(rhoBMax-(double)irho*dRho);
	
		electron.density=Yp*rhoB;
		electron.kf=pow(3.*pi2*electron.density, 1./3.);
		electron.solveChemPotEff();
		electron.chemPot = electron.chemPot_eff;
		electron.calculateProperties();
		
		double freeEnEl = electron.energy - temperature*electron.entropy;
		
		double dim;
		//sna.solveCPA(rhoB, Yp, temperature, saveSol);

		for(int iDim=2; iDim>=0; iDim--){
		dim = (double) (iDim)+1;
						
			sna.solveSNA(rhoB, Yp, temperature, dim, 0, saveSol);
			f_sna(iDim)= sna.f;
	
			nup1_sna(iDim)=cluster.proton.chemPot_eff;
			nun1_sna(iDim)=cluster.neutron.chemPot_eff;
			Mef1_sna(iDim)=cluster.neutron.mass;
			nup2_sna(iDim)=gas.proton.chemPot_eff;
			nun2_sna(iDim)=gas.neutron.chemPot_eff;
			Mef2_sna(iDim)=gas.neutron.mass;
			
			double sigma= getSurfaceTension(cluster, Yp, temperature);
			Rd_sna(iDim)= getRadiusD(dim, f_sna(iDim), Yp, cluster, gas);
			Rw_sna(iDim)= Rd_sna(iDim)/pow(f_sna(iDim), 1./dim);
			Vn_sna(iDim) = 4.*M_PI*pow(getRadiusD(3., f_sna(iDim), Yp, cluster, gas), 3.)/3.;
			
			freeEn_gas(iDim)= gas.getEnergy() - temperature*gas.getEntropy();
			interfaceEn_sna(iDim)= 3.*sigma*dim/(2.*Rd_sna(iDim));
			
			double energy = f_sna(iDim)*cluster.getEnergy()  + (1.-f_sna(iDim))*gas.getEnergy()
										+ f_sna(iDim)*interfaceEn_sna(iDim) + electron.energy;
										
			double entropy= f_sna(iDim)*cluster.getEntropy() +(1.-f_sna(iDim))*gas.getEntropy() 
											+ electron.entropy;
			
			freeEn_sna(iDim) = energy-temperature*entropy;
			
//						std::cout <<rhoB*pow(cluster.Mn/hc, 3.) << " " << cluster.Ye << " "
//						<< (freeEn_sna(iDim)/rhoB - 1.)*cluster.Mn << std::endl;
						
		}	
		

	//Get phase that minimizes energy:
		Eigen::MatrixXd::Index minRow;
		
		double freeEn_ocp= freeEn_sna.minCoeff(&minRow);
		id_ocp=minRow;
		
		//double freeEn_ocp= freeEn_sna(2);
		//id_ocp=0;
		dim_ocp=(double)(id_ocp+1);
					
		saveSol[0]=nup1_sna(id_ocp);
		saveSol[1]=nun1_sna(id_ocp);
		saveSol[2]=Mef1_sna(id_ocp);
		saveSol[3]=nup2_sna(id_ocp);
		saveSol[4]=nun2_sna(id_ocp);
		saveSol[5]=Mef2_sna(id_ocp);
		
		cluster.setEOS_coexistence(nup1_sna(id_ocp), nun1_sna(id_ocp), Mef1_sna(id_ocp));
		gas.setEOS_coexistence(nup2_sna(id_ocp), nun2_sna(id_ocp), Mef2_sna(id_ocp));
		
		vn_ocp= Vn_sna(id_ocp);
		vw_ocp= vn_ocp/f_sna(id_ocp);
		
		A_ocp= vn_ocp*(cluster.rhoB - gas.rhoB);
		Z_ocp= vn_ocp*(cluster.proton.density- gas.proton.density);
		
//---------- FLUCTUATION ---------------------:
	
//External fields and rearrangement:
		double LambdaP = gas.proton.chemPot
							+	(f_sna(id_ocp)/(1.-f_sna(id_ocp)))*
				(-2.*interfaceEn_sna(id_ocp)/(3.*(cluster.proton.density- gas.proton.density))
				+dim_ocp*getSurfaceTensionDerivative(cluster, Yp, temperature)
						*(1.-f_sna(id_ocp))*(1.-Yp)/(rhoB*Rd_sna(id_ocp))
				);
								 
		double LambdaN = gas.neutron.chemPot
		-f_sna(id_ocp)/(1.-f_sna(id_ocp))
			*dim_ocp*getSurfaceTensionDerivative(cluster, Yp, temperature)*(1.-f_sna(id_ocp))
						*Yp/(rhoB*Rd_sna(id_ocp));	
		
		double rear_avg=f_sna(id_ocp)*(
			interfaceEn_sna(id_ocp)*getPhiFuncDerivative(dim_ocp,f_sna(id_ocp))
			/(3.*getPhiFunc(dim_ocp, f_sna(id_ocp)))
			+ dim_ocp*getSurfaceTensionDerivative(cluster, Yp, temperature)*(
				(cluster.proton.density - gas.proton.density)
				-Yp*(cluster.rhoB - gas.rhoB))/(rhoB*Rd_sna(id_ocp)))
						/(cluster.proton.density - gas.proton.density);
					
		
		omega0_ocp= vn_ocp*(cluster.getEnergy() -gas.getEnergy() + interfaceEn_sna(id_ocp)
										- temperature*(cluster.getEntropy() -gas.getEntropy() )
										- LambdaP*(cluster.proton.density - gas.proton.density)
										- LambdaN*(cluster.neutron.density- gas.neutron.density));
										
		omega_ocp= omega0_ocp+rear_avg*Z_ocp;
		
		double Vtot=0.;
		
		if(id_ocp==0){iPlot=3;}
		if(id_ocp==1){iPlot=2;}
		if(id_ocp==2){iPlot=1;}
		
		std::cout 
		<< rhoB*pow(cluster.Mn/hc, 3.) << " " 
		<< A_ocp << " "  << Z_ocp << " "
		<< dim_ocp << " " << omega_ocp << " " 
		<< omega0_ocp
		<< std::endl;
		
	if( (f_sna(id_ocp) >0.) && (f_sna(id_ocp)<1.) ){
		
		matter_class clusterN(parametrization);
		
		double Yp_cell, uN, sigmaN;

		double fracpP=1.35;	
		double fracpM=0.65;
		int const np=1;

		double fracnP=1.35;
		double fracnM=0.65;
		int const nn=1;


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
			
		
			rhonNT(in, ip) = (fracnP + (fracnM - fracnP)*in/nn)*cluster.neutron.density;
			rhopNT(in, ip) = (fracpP + (fracpM - fracpP)*ip/np)*cluster.proton.density;

			rhobNT(in, ip) = rhonNT(in, ip) + rhopNT(in, ip);
			
			clusterN.setEOS_nucleons(rhobNT(in, ip), rhopNT(in, ip)/rhobNT(in, ip), temperature);
			
			uN= (Yp*rhoB - gas.proton.density)/(clusterN.proton.density - gas.proton.density);
			
			Yp_cell=  ( uN*(rhopNT(in, ip) - gas.proton.density) + gas.proton.density )
											/( uN*(rhobNT(in, ip) - gas.rhoB) + gas.rhoB);
			
			sigmaN	= getSurfaceTension(clusterN, Yp_cell, temperature);

			double dimN;	
			for(int iDim=2; iDim>=0; iDim--){
				dimN = (double) (iDim)+1; 
				
				uNT(in, ip, iDim)= uN;
				RdNT(in, ip, iDim)= getRadiusD(dimN, uN, Yp_cell, clusterN, gas);
				RwNT(in, ip, iDim)= RdNT(in, ip, iDim)/pow(uN, 1./dimN);
				
				VNT(in, ip, iDim)= 4.*M_PI*pow(getRadiusD(3., uN, Yp_cell, clusterN, gas), 3.)/3.;
								
				ANT(in, ip, iDim) = (clusterN.rhoB - 	gas.rhoB)*VNT(in, ip, iDim);
				ZNT(in, ip, iDim) = (clusterN.proton.density -gas.proton.density)*VNT(in, ip, iDim);
				
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
													( LambdaP*(clusterN.proton.density - gas.proton.density)
													+ LambdaN*(clusterN.neutron.density- gas.neutron.density));
								
				GNT0(in, ip, iDim)= F0NT(in, ip, iDim) - MuNT(in, ip, iDim);
				GNT(in, ip, iDim)= 	FNT(in, ip, iDim)  - MuNT(in, ip, iDim);

				double Cdist0= omega0_ocp/temperature;
				double Cdist= omega_ocp/temperature;
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
		
		rhobbar += (rhobNT(in,ip) - gas.rhoB)*VNT(in, ip, iDim)*nN(in, ip, iDim);
		
		rhopbar += (rhopNT(in, ip) - gas.proton.density)*VNT(in, ip, iDim)*nN(in, ip, iDim);
				
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
	
		
		rhobbar+=gas.rhoB;
		rhopbar+=gas.proton.density;
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
							<< rhobbar*pow(cluster.Mn/hc, 3.) << " " << maxProb << " "
							<< A_mp << " "  << Z_mp 
		<< std::endl;
					
		std::cout << "Z_bars's= " << Zbar << " " << Zeff_xbar  << " " 
							<< Zeff_ybar << " " << Zeff_zbar 
		<< std::endl;
		
		std::cout << "Z_vars's= " << Zvar << " " << Zeff_xvar  
															<< " " << Zeff_yvar << " " << Zeff_zvar 
		<< std::endl;

		
		
		for(int iDim=2; iDim>=0; iDim--){
	
			std::ofstream outdim("../data/data_flu_rho/tmp_"+std::to_string(iDim+1)+"d_"+parametrization+"_yp"+std::to_string(Yp)+
																				 +"_T"+std::to_string(temperature*cluster.Mn)+"_"
																				 +std::to_string(irho)+".txt");
			outdim << "#rhonN rhopN rhobN GN PN RdN VNT An Zn uN"<< std::endl;

				for(int in=0; in<nn; in++){
				for(int ip=0; ip<np; ip++){

				outdim << rhonNT(in, ip)/cluster.neutron.density << " " 
							 << rhopNT(in, ip)/cluster.proton.density << " "
							 << rhobNT(in, ip)/cluster.rhoB << " "
							 << GNT(in, ip, iDim)*cluster.Mn << " "	<< pN(in, ip, iDim) << " " 
							 << GNT0(in, ip, iDim)*cluster.Mn << " "	<< p0N(in,ip,iDim) << " "
							 << RdNT(in, ip, iDim)*(hc/cluster.Mn) << " " << VNT(in, ip, iDim) <<  " "
							 << ANT(in, ip, iDim) << " " << ZNT(in, ip, iDim) << " " << uNT(in, ip, iDim) << " "
							 << nN(in, ip, iDim)
							 << std::endl;

				}
			}
					outdim.close();

		}

		std::cout << "Probs: "	<< p1 << " & " << p2 << " & " << p3 <<  std::endl;
		std::cout << "Probs0: "	<< distDimN0(0)/norma0 << " & " << distDimN0(1)/norma0 << " & " << distDimN0(2)/norma0 <<  std::endl;
	
		outGlobal << rhoB*pow(cluster.Mn/hc, 3.) << " "
			  << (freeEn_ocp/rhoB - 1.)*cluster.Mn << " "
			  << Rd_sna(id_ocp)*(hc/cluster.Mn) << " " 
				<< Rw_sna(id_ocp)*(hc/cluster.Mn) << " "
				<< A_ocp << " " << Z_ocp << " " << iPlot << " "
				<< p1  << "  " << p2  << "  " << p3 << " "
				<< rhobbar*pow(cluster.Mn/hc, 3.) << " "  //col 11
				<< (Fbar/rhobbar - 1.)*cluster.Mn << " "
				<< Rdbar*(hc/cluster.Mn) << " " 
				<< Rwbar*(hc/cluster.Mn) << " "
				<< Abar << " " << Zbar << " " //15, 16
				<< Acell << " " 
				<< gas.proton.density*pow(cluster.Mn/hc, 3.) << " " << gas.neutron.density*pow(cluster.Mn/hc, 3.)<< " "
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