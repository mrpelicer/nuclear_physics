#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/rmf_non_linear_walecka.h"
#include "../../include/quark_model.h"
#include "../../include/interpolator.h"

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
	std::string parametrization= 	"l3wr";
	std::string hyperon_params =	"l3wr3"; 

	nlwm_class hadrons(parametrization);
	quarks_class quarks;

	bool doHyperons	=	true;

	if(doHyperons) hadrons.includeHyperons(doHyperons, hyperon_params);
	int iflavor= (doHyperons== true) ? 3 : 2;
	quarks.setFlavorNumber(iflavor);

	double tcrit=170./Mnucleon;
	double C, sqrtD;
	cout << "Choose C and sqrt(D)" << endl;
	cin >> C >> sqrtD; 
	// double C= 0.2;
	double D= pow(sqrtD, 2.)/pow(Mnucleon, 2);
	quarks.setParameters(C, D, tcrit);


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
//Set system variables
	// double temperature;
	double tempMin=0./Mnucleon;
  double tempMax=150./Mnucleon;
  int itMax=50;
  double dt=  (tempMax-tempMin)/itMax;

	double rhoB;
	double rhoBMin=0.2/pow(Mnucleon/hc, 3);
  double rhoBMax=1.1/pow(Mnucleon/hc, 3);///7.5*hadrons.rho0;
	//0.62/pow(Mnucleon/hc, 3); fsu2h c amm ou b
  int iR=120;
  double dRho=  (rhoBMax-rhoBMin)/iR;
	double Yu, Yd, Ys, Ye, Ym;

	std::string filename1;
	filename1="data/diagram_H_" +parametrization+
												"_Q_"+to_string(C)+"_"+to_string(sqrt(D)*Mnucleon)+".txt";

		std::ofstream outDiag(filename1);

	for(int it=0; it<itMax; it++){
 		double temperature=tempMin+ it*dt;

		std::string filenameH, filenameQ;
		vector<double>mubHv, mubQv, pressHv, pressQv, rhobHv, rhobQv;

		filenameH="data/press_H_" +parametrization+
														"_T"+to_string(temperature*Mnucleon)+".txt";
		filenameQ="data/press_Q_"+to_string(C)+"_"+to_string(sqrt(D)*Mnucleon)+
			 											"_T"+to_string(temperature*Mnucleon)+".txt";
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
			double EnergyH		= hadrons.getEnergy() 		+ electron.energy	 + muon.energy;
			double EntropyH 	= hadrons.getEntropy() 		+ electron.entropy + muon.entropy;
			double PressureH	=hadrons.getPressure() + electron.pressure + muon.pressure;

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

			// quarks.setEoSFlavor_muBFixed(hadrons.muB, temperature, electron, muon,
														// Yu, Yd, Ys);			
			quarks.setEoSFlavor_muBFixed2(hadrons.muB, temperature, electron, muon,
																		Yu, Yd, Ys, Ye, Ym);			


			// quarks.setEoSFlavor_PressFixed(PressureH, temperature, electron, muon,
																			// Yu, Yd, Ys);			
			

			// quarks.setEoSFlavorFixed(hadrons.getBaryonDens(), temperature, Yu, Yd, Ys);
			
			double PressureQ= quarks.getPressure()+electron.pressure+muon.pressure;
			
			quarks.muB= (quarks.getEnergy() +electron.energy + muon.energy
      -temperature*(quarks.getEntropy()+electron.entropy + muon.entropy)
      + PressureQ)/quarks.getBaryonDens();

				outHadron << hadrons.muB*Mnucleon << " " 
									<< PressureH*Mnucleon*pow(Mnucleon/hc, 3) << " "
									<< hadrons.getBaryonDens()*pow(Mnucleon/hc, 3) << " " 
									<< Yu << " " << Yd << " " << Ys
									<< std::endl;
			
				outQuark << quarks.muB*Mnucleon << " " 
								 << PressureQ*Mnucleon*pow(Mnucleon/hc, 3) << " "
								 << quarks.rhoB*pow(Mnucleon/hc, 3)
								 <<endl;

			if(PressureH>=0. && 
						(irho==0 || (hadrons.rhoB < rhobHv.back() && quarks.rhoB < rhobQv.back()))){
				mubHv.push_back(hadrons.muB);
				mubQv.push_back(quarks.muB);
				pressHv.push_back(PressureH);
				pressQv.push_back(PressureQ);
				rhobHv.push_back(hadrons.rhoB);
				rhobQv.push_back(quarks.rhoB);
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

		if(mubHv[0] == mubHv[mubHv.size()-1]){
			break;
		}else{
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

		double mub_= interpolation_func(rhoht, mubHv, rhobHv);
		double pressure_=interpolation_func(rhoht, pressHv, rhobHv);
		if(interpolation_func(rhoht, pressHv, rhobHv)-interpolation_func(rhoqt, pressQv, rhobQv)>1e-6 
		|| interpolation_func(rhoht, mubHv, rhobHv)-interpolation_func(rhoqt, mubQv, rhobQv) >1e-6){
			mub_=NAN;
			pressure_=NAN;
		}
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

	double rhobh_=1.*pow(hc/Mnucleon, 3);
	double rhobq_=1.2*pow(hc/Mnucleon, 3);
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
	options.parameter_tolerance = 1e-8;
	options.function_tolerance = 1e-8;
	options.gradient_tolerance=1e-10;
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
	options.max_num_iterations=1e4;	
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