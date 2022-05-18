#include "../../include/constant.h"
#include "../../include/particles.h"
#include "../../include/bag_model.h"
#include "../../include/interpolator.h"

#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_interp.h>	
#include <gsl/gsl_spline.h>

using namespace std;

vector<double> findTransition(vector<double> mubQv_, 
															vector<double> pressQv_,
															vector<double> rhobQv_);

struct TransitionFunctor{
public:
	TransitionFunctor(vector<double> mubQv_, 
										vector<double> pressQv_, 
										vector<double> rhobQv_);

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    vector<double> mubHv, mubQv,pressHv,pressQv, rhobHv, rhobQv;
};

int main(){

	bool doHyperons	=	false;

	bag_model_class quarks;

	double Bag=pow(205./Mnucleon, 4);	//MeV
	double Gv=0.8*pow(Mnucleon/hc, 2); // fm^2
	double xsi=6.*20.;
	double Xv=1; 											//adim
	double tcrit=0/Mnucleon;					//MeV

	quarks.setParameters(Bag, Gv, xsi, Xv, tcrit);
	int iflavor= (doHyperons== true) ? 3 : 2;
	cout << iflavor << endl;
	quarks.setFlavorNumber(iflavor);

//Set system variables
	double tempMin=0./Mnucleon;
  double tempMax=200./Mnucleon;
  int itMax=100;
  double dt=  (tempMax-tempMin)/itMax;

	double rhoB;
	double rhoBMin=0.002/pow(Mnucleon/hc, 3);
  double rhoBMax=1.5/pow(Mnucleon/hc, 3);///7.5*hadrons.rho0;
	//0.62/pow(Mnucleon/hc, 3); fsu2h c amm ou b
  int iR=300;
  double dRho=  (rhoBMax-rhoBMin)/iR;
	
	std::string filename="data/phase_diag_B"+to_string(pow(Bag, 1/4.)*Mnucleon)+
																				"_Gv"+to_string(Gv*pow(hc/Mnucleon, 2))+
																				"_xsi"+to_string(xsi)+".txt";

	std::ofstream outDiag(filename);

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


	for(int it=0; it<itMax; it++){
 		double temperature=tempMin+ it*dt;
		electron.temperature=temperature;
		muon.temperature=temperature;
		
		std::string filenameH, filenameQ;
		vector<double> mubQv, pressQv, rhobQv;
		
		filenameQ="data/eos_sym_B"+to_string(pow(Bag, 1/4.)*Mnucleon)+
			 											"_T"+to_string(temperature*Mnucleon)+".txt";
		std::ofstream outQuark (filenameQ);

		// electron.temperature=temperature;
		// muon.temperature=temperature;
		cout << temperature*Mnucleon << endl;
		// int ip=0;
		// do{
		// // for(int ip=0; ip<=iP; ip++){
		// 	Press=(PressMax- (double)ip*dPress);
		for(int irho=0; irho<=iR; irho++){
			rhoB=rhoBMax- (double)irho*dRho;

			// quarks.setEOS_symmetric(rhoB, temperature);
			quarks.setEOS_betaEq(rhoB, temperature, electron, muon);

			double PressureQ= quarks.getPressure() + electron.pressure + muon.pressure;
		
				outQuark << quarks.muB*Mnucleon << " " 
								 << PressureQ*Mnucleon*pow(Mnucleon/hc, 3) << " "
								 << quarks.rhoB*pow(Mnucleon/hc, 3)
								 <<endl;

				mubQv.push_back(quarks.muB);
				pressQv.push_back(PressureQ);
				rhobQv.push_back(quarks.rhoB);
		}
		quarks.firstRun=true;
		reverse(mubQv.begin(), mubQv.end());
		reverse(pressQv.begin(), pressQv.end());
		reverse(rhobQv.begin(), rhobQv.end());

		vector<double> trans_point= findTransition(mubQv, pressQv, rhobQv);
		double rhoqt= trans_point[0];

			
		double mub_= interpolation_func(rhoqt, mubQv, rhobQv);
		double pressure_=interpolation_func(rhoqt, pressQv, rhobQv);


		cout << "transition: T= " << temperature*Mnucleon << " "  
			  << rhoqt*pow(Mnucleon/hc, 3.)  << " " 
				<< pressure_*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
				<< mub_*Mnucleon  << " "
				 << endl;

		outDiag << mub_*Mnucleon  << " " 
						<< temperature*Mnucleon << " " 
						<<	pressure_*Mnucleon*pow(Mnucleon/hc, 3.)  << " " 
						<< rhoqt*pow(Mnucleon/hc, 3.) << " "
						<< endl;
		}
		outDiag.close();
		
  return 0;
}

vector<double> findTransition(vector<double> mubQv_, 
															vector<double> pressQv_, 
															vector<double> rhobQv_){

	double rhobq_=.6*pow(hc/Mnucleon, 3);
	double x[]={rhobq_};
	 
	Problem p;
	CostFunction* cost= 
							new NumericDiffCostFunction<TransitionFunctor,ceres::CENTRAL, 1, 1>
							(new TransitionFunctor(mubQv_, pressQv_, rhobQv_));
	p.AddResidualBlock(cost, NULL, x);
	//  cout<< rhobHv_.front()*pow(Mnucleon/hc, 3) << " " << rhobHv_.back()*pow(Mnucleon/hc, 3) << " " 
	//  			<< rhobQv_.front()*pow(Mnucleon/hc, 3) << " " << rhobQv_.back()*pow(Mnucleon/hc, 3) << endl;
	p.SetParameterLowerBound(x, 0, rhobQv_.front());
	p.SetParameterUpperBound(x, 0, rhobQv_.back());
	Solver::Options options;
	options.parameter_tolerance = 1e-12;
	options.function_tolerance = 1e-12;
	options.gradient_tolerance=1e-14;
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
	 std::cout  << rhobq_  << " ---> " << x[0]
	            << std::endl;
	
	vector<double> sol(1);
	sol[0] = summary.IsSolutionUsable() == true ? x[0] : NAN;
	return sol;
}

TransitionFunctor::TransitionFunctor(vector<double> mubQv_, 
																		 vector<double> pressQv_, 
																		 vector<double> rhobQv_){
	mubQv = mubQv_;
	pressQv = pressQv_;
	rhobQv= rhobQv_;
}

template <typename T>
bool TransitionFunctor::operator()(const T* x, T* residuals) const{

	double pressq_;
	pressq_= interpolation_func(x[0], pressQv, rhobQv);

  residuals[0] = pressq_;

	return true;
	
}