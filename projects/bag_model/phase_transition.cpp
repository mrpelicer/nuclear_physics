#include "../../include/constant.hpp"
#include "../../include/particles.hpp"
#include "../../include/rmf_walecka.hpp"
#include "../../include/bag_model.hpp"
#include "../../include/interpolator.hpp"

#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_interp.h>	
#include <gsl/gsl_spline.h>

using namespace std;

vector<double> findTransition(vector<double> mubHv_, vector<double> mubQv_, 
							vector<double> pressHv_, vector<double> pressQv_,
							vector<double> rhobHv_, vector<double> rhobQv_);

struct TransitionFunctor{
public:
	TransitionFunctor(vector<double> mubHv_, vector<double> mubQv_, 
							vector<double> pressHv_, vector<double> pressQv_, 
							vector<double> rhobHv_, vector<double> rhobQv_);

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    vector<double> mubHv, mubQv,pressHv,pressQv, rhobHv, rhobQv;
};

int main(){
//Choose pararametrization
	std::string parametrization= 	"nl3wr*";
	std::string hyperon_params =	"nl3wr*"; 

	nlwm_class hadrons(parametrization);

	bool doHyperons	=	true;
	cout << "dohyp " << doHyperons << endl;
	// cout << "DoHyperons? (0 = cin >> "
	// cin >> doHyperons; 
	if(doHyperons) hadrons.includeHyperons(doHyperons, hyperon_params);

	bag_model_class quarks;


	double Bag=pow(165./Mnucleon, 4);	//MeV
	double Xv=1; 											//adim
	double Gv=0.*pow(Mnucleon/hc, 2); // fm^2
	double xsi=0.;
	// double tcrit=0/Mnucleon;					//MeV

	cout << "Choose Bag^1/4 (Mev),  Gv (fm2),  Xv, xsi" << endl;
	cin >> Bag >> Gv >> Xv >> xsi;

	// cout << "Choose Tcrit (MeV) -- >0 if you want B(T) and 0 for cte B" << endl;
	// cin >> tcrit;
	Bag=pow(Bag/Mnucleon, 4);
	Gv*=pow(Mnucleon/hc, 2);
	// tcrit*=1./Mnucleon;

	quarks.setParameters(Bag, Gv, xsi, Xv, 0);
	int iflavor= (doHyperons== true) ? 3 : 2;
	quarks.setFlavorNumber(iflavor);


	particle electron;
	electron.mass= Me/Mnucleon;
	electron.mass_eff=Me/Mnucleon;
	electron.Q=-1.;
	electron.gamma=2.;

	particle muon;
	muon.mass= Mm/Mnucleon;
	muon.mass_eff=Mm/Mnucleon;
	muon.Q=-1.;
	muon.gamma=2.;

	
	bool doBfield		= false;
	double Bg=0; //3.e18; // Gauss
	cout << "Choose the magnetic field (G)" << endl;
	cin >> Bg;
	if(Bg>1e15) doBfield= true;
	double Bc=pow(electron.mass_eff, 2.)/eHL;  
	double Bfield=Bg*Bc/4.41e13;
	
	cout << "You chose B= " << Bg << "G \n";
	if(doBfield){
		hadrons.setBfield(	doBfield, Bfield);
		quarks.setBfield(		doBfield, Bfield);
		electron.setBfield(	doBfield, Bfield);
	 	muon.setBfield(			doBfield, Bfield);
	}

//Set system variables
	// double temperature;
	double tempMin=0./Mnucleon;
  double tempMax=200./Mnucleon;
  int itMax=1;
  double dt=  (tempMax-tempMin)/itMax;

		//=== Loop over barionic density
	// for(int imu=0; imu<imuMax; imu++){
		// muB= muBMax- imu*dMu;
	double rhoB;
	double rhoBMin=0.002/pow(Mnucleon/hc, 3);
  double rhoBMax=1./pow(Mnucleon/hc, 3);///7.5*hadrons.rho0;
	//0.62/pow(Mnucleon/hc, 3); fsu2h c amm ou b
  int iR=300;
  double dRho=  (rhoBMax-rhoBMin)/iR;
	double Yu, Yd, Ys, Ye, Ym;

	std::string filename1;
	filename1="data/diagram_H_" +parametrization+
												"_B"+to_string(pow(Bag, 1/4.)*Mnucleon)+
												".txt";

	std::ofstream outDiag(filename1);

		string file_densH, file_densQ;
		if(!doBfield){
			file_densH= "data/densH_" +parametrization+"_noB.txt";
			file_densQ= "data/densQ_" +parametrization+"_noB.txt";
		}else{
			file_densH= "data/densH_" +parametrization+"_wtB.txt";
			file_densQ= "data/densQ_" +parametrization+"_wtB.txt";
		}
		ofstream outDensH(file_densH);
		ofstream outDensQ(file_densQ);

	for(int it=0; it<itMax; it++){
 		double temperature=tempMin+ it*dt;
  //  for (double temperature : tempv){
		//  temperature*=1./Mnucleon;

		std::string filenameH, filenameQ;
		vector<double>mubHv, mubQv, pressHv, pressQv, rhobHv, rhobQv, enerHv, enerQv;

		if(!doBfield){
			filenameH="data/press_H_" +parametrization+
									"_T"+to_string(temperature*Mnucleon)+"_noB.txt";
			filenameQ="data/press_Q_B"+to_string(pow(Bag, 1/4.)*Mnucleon)+
									"_Gv" +to_string(Gv*pow(hc/Mnucleon, 2))+
									"_xsi"+to_string(xsi)+
									"_T"	+to_string(temperature*Mnucleon)+"_noB.txt";
		}else{
			filenameH="data/press_H_" +parametrization+
									"_T"+to_string(temperature*Mnucleon)+"_wtB.txt";
			filenameQ="data/press_Q_B"+to_string(pow(Bag, 1/4.)*Mnucleon)+
									"_Gv" +to_string(Gv*pow(hc/Mnucleon, 2))+
									"_xsi"+to_string(xsi)+
									"_T"	+to_string(temperature*Mnucleon)+"_wtB.txt";			
		}

		std::ofstream outHadron(filenameH);
		std::ofstream outQuark (filenameQ);

		electron.temperature=temperature;
		muon.temperature=temperature;
		cout << temperature*Mnucleon << endl;
		// int ip=0;
		// do{
		// // for(int ip=0; ip<=iP; ip++){
		// 	Press=(PressMax- (double)ip*dPress);
		for(int irho=0; irho<=iR; irho++){
			rhoB=rhoBMax- (double)irho*dRho;
			
			//Solve self-consistently:
			hadrons.setEOS_betaEq(rhoB, temperature, electron, muon);
			
			electron.pressure	=electron.chemPot*electron.density 	- electron.energy + temperature*electron.entropy;
			muon.pressure			=muon.chemPot*muon.density 					- muon.energy			+ temperature*muon.entropy;
			//Set variables:
			// double Energy		= hadrons.getEnergy() 		+ electron.energy	 + muon.energy;
			// double Entropy 	= hadrons.getEntropy() 		+ electron.entropy + muon.entropy;
			double PressureH	=hadrons.getPressure() + electron.pressure + muon.pressure;
			double EnergyH  	=hadrons.getEnergy() +  electron.energy + muon.energy;

			Yu= (2.*hadrons.proton.density + hadrons.neutron.density
								+hadrons.lambda0.density  + 2.*hadrons.sigmap.density+ hadrons.sigma0.density
								+hadrons.xi0.density)/(3.*hadrons.getBaryonDens()); 

			Yd= (hadrons.proton.density + 2.*hadrons.neutron.density
					+hadrons.lambda0.density  +  hadrons.sigma0.density + 2.*hadrons.sigmam.density
					+hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							

			Ys= (hadrons.lambda0.density  
					+  hadrons.sigmap.density + hadrons.sigma0.density
					+		hadrons.sigmam.density +2.*hadrons.xi0.density
					+		2.*hadrons.xim.density)/(3.*hadrons.getBaryonDens()); 							

			Ye=electron.density/rhoB;
			Ym=muon.density/rhoB;
			// cout << Yu << " " << Yd << " " << Ys << endl;

			quarks.setEoSFlavor_muBFixed(hadrons.muB, temperature, electron, muon,
														Yu, Yd, Ys);			

			double edensh= electron.density;
			double mdensh= muon.density;
			// quarks.setEoSFlavor_muBFixed2(hadrons.muB, temperature, electron, muon,
													// Yu, Yd, Ys, Ye, Ym);			

			// quarks.setEoSFlavorFixed(hadrons.getBaryonDens(), temperature, Yu, Yd, Ys);


			double PressureQ= quarks.getPressure()+electron.pressure+muon.pressure;
			double EnergyQ  =	quarks.getEnergy() +  electron.energy + muon.energy;
			
			outHadron << hadrons.muB*Mnucleon << " " 
								<< PressureH*Mnucleon*pow(Mnucleon/hc, 3) << " "
								<< hadrons.getBaryonDens()*pow(Mnucleon/hc, 3) << " " 
								<< EnergyH*Mnucleon*pow(Mnucleon/hc, 3)
								<< std::endl;
		
			outQuark << quarks.muB*Mnucleon << " " 
							 << PressureQ*Mnucleon*pow(Mnucleon/hc, 3) << " "
							 << quarks.rhoB*pow(Mnucleon/hc, 3) << " " 
							 << EnergyQ*Mnucleon*pow(Mnucleon/hc, 3) 
							 <<endl;

			outDensH << rhoB							*pow(Mnucleon/hc, 3)	<< " " // /hadrons.rho0 			 << " " // *pow(Mnucleon/hc, 3) << " " 
				<< hadrons.proton.density		*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.neutron.density	*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< edensh										*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< mdensh										*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.lambda0.density	*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.sigmap.density		*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.sigma0.density		*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.sigmam.density		*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB	 << " "
				<< hadrons.xi0.density			*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.xim.density			*pow(Mnucleon/hc, 3) 	<< " "  // /rhoB  << " "
				<< hadrons.muB*Mnucleon
				<< std::endl;

			outDensQ << quarks.getBaryonDens()*pow(Mnucleon/hc, 3) << " " // /hadrons.rho0 			 << " " // *pow(Mnucleon/hc, 3) << " " 
							 << quarks.qu.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
							 << quarks.qd.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
							 << electron.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
							 << muon.density					*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
							 << quarks.qs.density			*pow(Mnucleon/hc, 3) << " "  // /rhoB  << " "
			 				 << quarks.muB*Mnucleon
							 << std::endl;


			if(PressureH>=0. && 
						(irho==0 || (hadrons.rhoB < rhobHv.back() && quarks.rhoB < rhobQv.back()))){
				mubHv.push_back(hadrons.muB);
				mubQv.push_back(quarks.muB);
				pressHv.push_back(PressureH);
				pressQv.push_back(PressureQ);
				rhobHv.push_back(hadrons.rhoB);
				rhobQv.push_back(quarks.rhoB);
				enerHv.push_back(EnergyH);
				enerQv.push_back(EnergyQ);
			}
			// ip++;
		// }while((hadrons.getBaryonDens()*pow(Mnucleon/hc, 3))>0.05);
		}
		// }
		hadrons.firstRun=true;
		quarks.firstRun=true;
		reverse(mubHv.begin(), mubHv.end());
		reverse(mubQv.begin(), mubQv.end());
		reverse(pressHv.begin(), pressHv.end());
		reverse(pressQv.begin(), pressQv.end());
		reverse(rhobHv.begin(), rhobHv.end());
		reverse(rhobQv.begin(), rhobQv.end());
		reverse(enerHv.begin(), enerHv.end());
		reverse(enerQv.begin(), enerQv.end());



		vector<double> trans_point= findTransition(mubHv, mubQv, pressHv, pressQv, rhobHv, rhobQv);
		double rhoht= trans_point[0];
		double rhoqt= trans_point[1];
							
		cout << "transition: T= " << temperature*Mnucleon << " "  
			  << rhoht*pow(Mnucleon/hc, 3.)  << " " << rhoqt*pow(Mnucleon/hc, 3.)  << " " 
				<< interpolation_func(rhoht, pressHv, rhobHv)*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
				<< interpolation_func(rhoqt, pressQv, rhobQv)*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
				<< interpolation_func(rhoht, mubHv, rhobHv)*Mnucleon  << " " 
				<< interpolation_func(rhoqt, mubQv, rhobQv)*Mnucleon  << " "
				 << endl;

			double Le=interpolation_func(rhoht, pressHv, rhobHv)
						*(interpolation_func(rhoqt, enerQv, rhobQv)    
							- interpolation_func(rhoht, enerHv, rhobHv))
							/(interpolation_func(rhoqt, enerQv, rhobQv)*interpolation_func(rhoht, enerHv, rhobHv));

			cout << "$P_0 = " << interpolation_func(rhoqt, pressQv, rhobQv)*Mnucleon*pow(Mnucleon/hc, 3.) << " " 
			<< "$ \\\\ $ \\mu_0 =" << interpolation_func(rhoqt, mubQv, rhobQv)*Mnucleon 
			<< "$ \\\\ $ L \\big{|}_\\varepsilon =" << Le << "$" << endl;

			double mu0= interpolation_func(rhoht, mubHv, rhobHv);
			double p0	= interpolation_func(rhoht, pressHv, rhobHv);
			cout << pow(Bag, 1/4.)*Mnucleon << " " << xsi << " "
			<< p0*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
			<< mu0*Mnucleon  << " " 
			<< rhoht*pow(Mnucleon/hc, 3.)  << " " << rhoqt*pow(Mnucleon/hc, 3.)  << " " 
			<< interpolation_func(rhoht, enerHv, rhobHv)*Mnucleon*pow(Mnucleon/hc, 3.) << " "   
			<< interpolation_func(rhoqt, enerQv, rhobQv)*Mnucleon*pow(Mnucleon/hc, 3.) 
			<< endl;

			cout << Bg << " " << Le << endl;
		if(interpolation_func(rhoht, pressHv, rhobHv)-interpolation_func(rhoqt, pressQv, rhobQv)>1e-8 
		|| interpolation_func(rhoht, mubHv, rhobHv)-interpolation_func(rhoqt, mubQv, rhobQv) >1e-8){
			cout << "TRANSITION WRONGLY CALCULATED!!!" << endl;
		}else{

		double mub_= interpolation_func(rhoht, mubHv, rhobHv);
		double pressure_=interpolation_func(rhoht, pressHv, rhobHv);

		outDiag << mub_*Mnucleon  << " " 
						<< temperature*Mnucleon << " " 
						<<	pressure_*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
						<< rhoht*pow(Mnucleon/hc, 3.)  << " " << rhoqt*pow(Mnucleon/hc, 3.) << " "
						<< Yu << " " << Yd << " " << Ys
						<< endl;
		}
			outHadron.close();
			outQuark.close();

	}
		
	outDiag.close();
  return 0;
}

vector<double> findTransition(vector<double> mubHv_, vector<double> mubQv_, 
							vector<double> pressHv_, vector<double> pressQv_, 
							vector<double> rhobHv_, vector<double> rhobQv_){

	// double rhobh_=1.*pow(hc/Mnucleon, 3);
	// double rhobq_=1.2*pow(hc/Mnucleon, 3);
	// double rhobh_= rhobHv_[floor(rhobHv_.size()/2)];
	// double rhobq_= rhobQv_[floor(rhobQv_.size()/2)];

	/*bag~150 				-270
		bag~155--175 		-210	
		bag~180--195		-120
		bag~200					-120
	*/
	double rhobh_= rhobHv_[rhobHv_.size()-160];
	double rhobq_= rhobQv_[rhobQv_.size()-160];
	double x[]={rhobh_, rhobq_};
	 
	Problem p;
	CostFunction* cost= 
							new NumericDiffCostFunction<TransitionFunctor,ceres::CENTRAL, 2, 2>
							(new TransitionFunctor(mubHv_, mubQv_, pressHv_, pressQv_, rhobHv_ ,rhobQv_));
	p.AddResidualBlock(cost, NULL, x);
	//  cout<< rhobHv_.front()*pow(Mnucleon/hc, 3) << " " << rhobHv_.back()*pow(Mnucleon/hc, 3) << " " 
	//  			<< rhobQv_.front()*pow(Mnucleon/hc, 3) << " " << rhobQv_.back()*pow(Mnucleon/hc, 3) << endl;
  p.SetParameterLowerBound(x, 0, rhobHv_.front());
  p.SetParameterUpperBound(x, 0, rhobHv_.back());
	p.SetParameterLowerBound(x, 1, rhobQv_.front());
	p.SetParameterUpperBound(x, 1, rhobQv_.back());
	Solver::Options options;
	options.parameter_tolerance = 1e-10;
	options.function_tolerance = 1e-12;
	options.gradient_tolerance=1e-13;
	// options.line_search_direction_type= ceres::STEEPEST_DESCENT;
	options.use_nonmonotonic_steps= true;
	options.dense_linear_algebra_library_type=ceres::LAPACK;
	options.linear_solver_type= ceres::DENSE_QR;
  options.update_state_every_iteration = true;
	//options.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
	//options.linear_solver_type= ceres::DENSE_QR;
	//options.trust_region_strategy_type = ceres::DOGLEG;
	//options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
  options.minimizer_progress_to_stdout = true;
		 
	Solver::Summary summary;
	options.max_num_iterations=2e3;	
	//Run
	Solve(options, &p, &summary);
	//Print if convergence was achieved.
	 std::cout << summary.FullReport() <<  "\n" << summary.IsSolutionUsable() << "\n"; 
	 std::cout  << rhobh_ << " " << rhobq_  << " ---> " << x[0] << " " << x[1] 
	            << std::endl;
	
	vector<double> sol(2);
	sol[0] = summary.IsSolutionUsable() == true ? x[0] : NAN;
	sol[1] = summary.IsSolutionUsable() == true ? x[1] : NAN;
	return sol;
}

TransitionFunctor::TransitionFunctor(vector<double> mubHv_, vector<double> mubQv_, 
							vector<double> pressHv_, vector<double> pressQv_, 
							vector<double> rhobHv_, vector<double> rhobQv_){
	mubHv = mubHv_;
	mubQv = mubQv_;
	pressHv = pressHv_;
	pressQv = pressQv_;
	rhobHv= rhobHv_;
	rhobQv= rhobQv_;
}
template <typename T>
bool TransitionFunctor::operator()(const T* x, T* residuals) const{

	double mubh_, mubq_, pressh_, pressq_;
	
	 mubh_= interpolation_func(x[0], pressHv, rhobHv);
	 mubq_= interpolation_func(x[1], pressQv, rhobQv);
	 pressh_= interpolation_func(x[0], mubHv, rhobHv);
	 pressq_= interpolation_func(x[1], mubQv, rhobQv);

  residuals[0] = mubh_ - mubq_;
  residuals[1] = pressh_ - pressq_;

	return true;
	
}