#include <iostream>
#include <iterator>
#include <limits>
#include <cmath>
#include <cstdint>
#include <utility>
#include <iomanip>
#include <algorithm>
#include <vector> 
#include <fstream>

#include "../../include/constant.h"
#include "../../include/interpolator.h"
#include "../../include/tov_solver.h"

#include "ceres/ceres.h"
#include "glog/logging.h"
using ceres::NumericDiffCostFunction;
//using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/array.hpp>
// #include <boost/numeric/odeint.hpp>

using namespace std;

vector<double> findTransition(vector<double> mubhv_, vector<double> mubqv_, 
							vector<double> preshv_, vector<double> presqv_,
							vector<double> rhobhv_, vector<double> rhobqv_);

struct TransitionFunctor{
public:
	TransitionFunctor(vector<double> mubhv_, vector<double> mubqv_, 
							vector<double> preshv_, vector<double> presqv_, 
							vector<double> rhobhv_, vector<double> rhobqv_);

template <typename T>
  bool operator()(const T* arg, T* residuals) const;

private:
    vector<double> mubHv, mubQv,pressHv,pressQv, rhobHv, rhobQv;
};

int main(){

	//define hadron eos to be read
	string Lstr;
	cout << "Define L of the file (44, 60, 76, 92, 100, 108, 116)" << endl;
	cin >> Lstr;
	string heos_file = "input/EOS"+Lstr+".bin"; 
    ifstream heos_data(heos_file);
	
	//define quark eos to be read
	string gvstr;
	cout << "Define gv of the quark file (38, 40, 42)" << endl;
	cin >> gvstr;
	string qeos_file = "input/G"+gvstr+".bin"; 
	ifstream qeos_data(qeos_file);
    
    if(heos_data.fail()){ // checks to see if file opended
       cout << "Error in opening hadron EoS file " << heos_file << endl;
        return 1;
    }

	if(qeos_data.fail()){ // checks to see if file opended
       cout << "Error in opening quark EoS file " << qeos_file << endl;
        return 1;
    }

	double rhobh_, enerh_, mubh_, presh_, enerhmev_, preshmev_;
	vector<double> rhobhv, mubhv, enerhv, preshv;
    while (heos_data  >> rhobh_ >> enerh_ >> presh_ >> mubh_ >> enerhmev_ >> preshmev_){
        rhobhv.push_back(rhobh_);
		mubhv.push_back(mubh_);
        enerhv.push_back(enerhmev_);
        preshv.push_back(preshmev_);
		// enerhv.push_back(enerh_);
		// preshv.push_back(presh_);
    }

    heos_data.close();

	double rhobq_, mubq_, enerq_, presq_, enerqmev_, presqmev_;
	vector<double> rhobqv, mubqv, enerqv, presqv;
    while (qeos_data  >> rhobq_ >> enerq_ >> presq_ >> mubq_ >> enerqmev_ >> presqmev_){
        rhobqv.push_back(rhobq_);
		mubqv.push_back(mubq_);
        enerqv.push_back(enerqmev_);
        presqv.push_back(presqmev_);
    }

    qeos_data.close();

//find the phase transition:
	vector<double> trans_point= findTransition(mubhv, mubqv, preshv, presqv, rhobhv, rhobqv);

//save densities where it occurs:
	double rhoht= trans_point[0];
	double rhoqt= trans_point[1];

//find points from the density solution:
	double mu0	= interpolation_func(rhoht, mubhv, rhobhv);
	double p0	= interpolation_func(rhoht, preshv, rhobhv);					
	double eh0	= interpolation_func(rhoht, enerhv, rhobhv);
	double eq0	= interpolation_func(rhoqt, enerqv, rhobqv);					

//write the transition:
		cout << "check transition: " << "Densities :" <<  rhoht  << " " << rhoqt  << endl
				<< "press: " << interpolation_func(rhoht, preshv, rhobhv)  << " " 
				<< interpolation_func(rhoqt, presqv, rhobqv)  << endl
				<< "chempot: " << interpolation_func(rhoht, mubhv, rhobhv)  << " " 
				<< interpolation_func(rhoqt, mubqv, rhobqv)  << " "
				 << endl;

			double Le=interpolation_func(rhoht, preshv, rhobhv)
						*(interpolation_func(rhoqt, enerqv, rhobqv)    
							- interpolation_func(rhoht, enerhv, rhobhv))
							/(interpolation_func(rhoqt, enerqv, rhobqv)*interpolation_func(rhoht, enerhv, rhobhv));

			cout << "$P_0 = " << p0 << " " 
			<< "$ \\\\ $ \\mu_0 =" << mu0
			<< "$ \\\\ $ L \\big{|}_\\varepsilon =" << Le << "$" << endl;

			
		cout << Lstr << " "
				<< mu0  << " " <<	p0  << " " 
				<< eh0 << " " << eq0 << " "
				<< rhoht << " " << rhoqt << " " 
				<< Le
				<< endl;


	int iener= 1000;
    double en_max= enerhv.front();
    double en_min= enerhv.back();
    double de= (log10(en_max)- log10(en_min))/iener;

    double pr_max= preshv.front();
    double pr_min= preshv.back();
    // double dp= (pr_max -pr_min)/iener;
	double dp= (log10(pr_max) -log10(pr_min))/iener;

    ofstream eos_test("data/eos_L"+Lstr+"_gv"+gvstr+".txt");

    for(int ie=1; ie<iener; ie++){

		double p_ = pow(10., log10(pr_max) - ie*dp);
		vector<double> pv_= p_>p0 ? presqv : preshv;
		vector<double> ev_= p_>p0 ? enerqv : enerhv;
		vector<double> nv_= p_>p0 ? rhobqv : rhobhv;
		vector<double> mv_= p_>p0 ? mubqv : mubhv;


        double e_= interpolation_func(p_, ev_, pv_);
		double dens_= interpolation_func(p_, nv_, pv_);
		double mu_= interpolation_func(p_, mv_, pv_);


        eos_test << e_ << " " <<  p_<< endl;
    }
	eos_test.close();
// transform to mev/fm3

	// std::for_each(preshv.begin(), preshv.end(), [](double &el){el *= hc; });
    // std::for_each(enerhv.begin(), enerhv.end(), [](double &el){el *= hc; });
	// std::for_each(presqv.begin(), presqv.end(), [](double &el){el *= hc; });
	// std::for_each(enerqv.begin(), enerqv.end(), [](double &el){el *= hc; });
	// p0*=hc;
	// eh0*=hc;
	// eq0*=hc;

 
	// int ipres= 150;
    // double pq_max= presqv.back();
	// double pq_min= presqv.front() > 0. ? presqv.front()  : 1e-11;
	// double ph_max= preshv.back();
    // double ph_min= preshv.front();
    // cout <<  pq_max << " "
	// 	 <<  pq_min << " "
	// 	 <<  ph_max << " "
    // 	 <<  ph_min << " " << endl;
    
	// double dp_h=  (log10(ph_max)- log10(ph_min))/ipres;
	// double dp_q=  (log10(pq_max)- log10(pq_min))/ipres;
	// double dp_qh= (log10(pq_max)- log10(ph_min))/ipres;

	// ofstream tov_h_file("data/mr_h_"+Lstr+".txt");
	// ofstream tov_q_file("data/mr_q_"+gvstr+".txt");
	// ofstream tov_qh_file("data/mr_"+Lstr+"_"+gvstr+".txt");

	// tov_class tov_hybrid(enerhv, preshv, enerqv, presqv, p0, eh0, eq0);
	// tov_class tov_hadron(enerhv, preshv);
	// tov_class tov_quarks(enerqv, presqv);
	
	// tov_hadron.do_crust();
	// tov_hybrid.do_crust();
	// double mass_qh, radius_qh,  mass_h, radius_h,  mass_q, radius_q;
	// int ipmax= ipres-40;
    // for(int ip=1; ip<ipmax; ip++){
      
    //     // double p_qh_center=pq_max-ip*dp_qh;
	// 	double p_qh_center = pow(10., log10(pq_max) - ip*dp_qh);
	// 	// double p_h_center= ph_max-ip*dp_h;
 	// 	double p_h_center = pow(10., log10(ph_max) - ip*dp_h);
	// 	// double p_q_center= pq_max-ip*dp_q;
	// 	double p_q_center = pow(10., log10(pq_max) - ip*dp_q);
		

	// 	tov_hadron.solve_tov_euler_pcentral(p_h_center, 1e2);
	// 	tov_quarks.solve_tov_euler_pcentral(p_q_center, 1e2);
	// 	tov_hybrid.solve_tov_euler_pcentral(p_qh_center, 1e2);

	// 	mass_h=tov_hadron.getMass(); 	
	// 	radius_h=tov_hadron.getRadius();

	// 	mass_q=tov_quarks.getMass(); 	
	// 	radius_q=tov_quarks.getRadius();

	// 	mass_qh=tov_hybrid.getMass(); 	
	// 	radius_qh=tov_hybrid.getRadius();

    //     cout << ip << "/" << ipmax << " " << endl;

 	// 	tov_h_file  << radius_h << " " << mass_h << " " << p_h_center << endl;
	// 	tov_q_file  << radius_q << " " << mass_q << " " << p_q_center << endl;
    //     tov_qh_file << radius_qh << " " << mass_qh << " " << p_qh_center << endl;
	// 	// first=false;
	// }	
	// tov_h_file.close();
	// tov_q_file.close();
	// tov_qh_file.close();
  	return 0;
}

vector<double> findTransition(vector<double> mubhv_, vector<double> mubqv_, 
							vector<double> preshv_, vector<double> presqv_, 
							vector<double> rhobhv_, vector<double> rhobqv_){

	
	double rhobh_= rhobhv_[rhobhv_.size()-2];
	double rhobq_= rhobqv_[rhobqv_.size()-2];
	double x[]={rhobh_, rhobq_};
	 
	Problem p;
	CostFunction* cost= 
							new NumericDiffCostFunction<TransitionFunctor,ceres::CENTRAL, 2, 2>
							(new TransitionFunctor(mubhv_, mubqv_, preshv_, presqv_, rhobhv_ ,rhobqv_));
	p.AddResidualBlock(cost, NULL, x);
	//  cout<< rhobhv_.front()*pow(Mnucleon/hc, 3) << " " << rhobhv_.back()*pow(Mnucleon/hc, 3) << " " 
	//  			<< rhobqv_.front()*pow(Mnucleon/hc, 3) << " " << rhobqv_.back()*pow(Mnucleon/hc, 3) << endl;
  	p.SetParameterLowerBound(x, 0, rhobhv_.front());
  	p.SetParameterUpperBound(x, 0, rhobhv_.back());
	p.SetParameterLowerBound(x, 1, rhobqv_.front());
	p.SetParameterUpperBound(x, 1, rhobqv_.back());
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

//Class of the phase transition:
TransitionFunctor::TransitionFunctor(vector<double> mubhv_, vector<double> mubqv_, 
							vector<double> preshv_, vector<double> presqv_, 
							vector<double> rhobhv_, vector<double> rhobqv_){
	mubHv = mubhv_;
	mubQv = mubqv_;
	pressHv = preshv_;
	pressQv = presqv_;
	rhobHv= rhobhv_;
	rhobQv= rhobqv_;
}

template <typename T>
bool TransitionFunctor::operator()(const T* x, T* residuals) const{

	double mubh_, mubq_, pressh_, pressq_;
	
	 pressh_= interpolation_func(x[0], pressHv, rhobHv);
	 pressq_= interpolation_func(x[1], pressQv, rhobQv);
	 mubh_= interpolation_func(x[0], mubHv, rhobHv);
	 mubq_= interpolation_func(x[1], mubQv, rhobQv);

	residuals[0] = mubh_ - mubq_;
	residuals[1] = pressh_ - pressq_;

	return true;	
}



// 	int iener= 1000;
//     double en_max= enerhv.front();
//     double en_min= enerhv.back();
//     double de= (log10(en_max)- log10(en_min))/iener;

//     double pr_max= preshv.front();
//     double pr_min= preshv.back();
//     // double dp= (pr_max -pr_min)/iener;
// 	double dp= (log10(pr_max) -log10(pr_min))/iener;

//     ofstream eos_test("eos.txt");

//     for(int ie=0; ie<iener; ie++){
//      // cout << e_center << endl;
        
//         // double e_center= en_max-ie*de;
// 					//Bg= pow(10., lBgMax- ib*dlogB);
// 		double e_center = pow(10., log10(en_max) - ie*de);
//         double p_center=  interpolation_func(e_center, preshv, enerhv);// interpolatePres(e_center);
        
// 		double p_ = pow(10., log10(pr_max) - ie*dp);
//         // double p_ = pr_max- ie*dp;
//         double e_= interpolation_func(p_, enerhv, preshv);
//         // double e_= interpolateEner(p_);

//         eos_test << e_center <<  " " <<  p_center << " " << e_ << " " <<  p_<< endl;
//     }
// 	eos_test.close();