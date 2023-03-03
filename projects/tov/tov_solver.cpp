
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
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
using namespace std;


void skip(istream & in, size_t n , char delim)
{
   size_t i = 0;
   while ( i++ < n)
      in.ignore(1e12, delim);
// ignores up to 1e12 chars but stops ignoring after delim
}


int main(){
    string eos_file = "eos_iufsu_noB.txt";
    // string eos_file = "eos_iufsu_noB.txt"; 
    ifstream eos_data("input/"+eos_file);
    ofstream eos_test("data/eos.txt");
    
    if(eos_data.fail()){ // checks to see if file opended
       cout << "Error in opening EoS file " << eos_file << endl;
        return 1;
    }

    double rhob_, ener_, pres_;
    vector<double> rhobv, enerv, presv;

    // skip(eos_data, 1, '\n');
    while (eos_data  >> rhob_ >> ener_ >> pres_){
        rhobv.push_back(rhob_);
        enerv.push_back(ener_);
        presv.push_back(pres_);
    }

    eos_data.close();

    if(rhobv.back() < rhobv.front()){
        std::reverse(rhobv.begin(), rhobv.end());
    	std::reverse(enerv.begin(), enerv.end());
    	std::reverse(presv.begin(), presv.end());
    }


    //change units from mev/fm3 to 1/fm4
    //  std::for_each(presv.begin(), presv.end(), [](double &el){el *= hc; });
    //  std::for_each(enerv.begin(), enerv.end(), [](double &el){el *= hc; });
    ofstream outMR("data/mr_"+eos_file);


    //check the EOS:


    tov_class tov(enerv, presv);
    double dr=1e2;

    tov.do_crust();
    // int iener= 1000;
    // int iener=tov.e_eosv.size();
    // double en_max= tov.e_eosv.back();
    // double en_min= tov.e_eosv.front();
    // double de= (en_max- en_min)/iener;

    // double pr_max= tov.p_eosv.back();
    // double pr_min= tov.p_eosv.front();
    // double dp= (pr_max -pr_min)/iener;

    // for(int ie=0; ie<iener; ie++){
        
    //     double e_center= en_max-ie*de;
            //  cout << e_center << endl;

        // double p_center=  interpolation_func(e_center, tov.p_eosv, tov.e_eosv);// interpolatePres(e_center);
        
        // double p_ = pr_max- ie*dp;
        // double e_= interpolation_func(p_, tov.e_eosv, tov.p_eosv);
        // double e_= interpolateEner(p_);
        // eos_test << tov.e_eosv[ie] << " " << tov.p_eosv[ie] << endl;
        // eos_test << e_center <<  " " <<  p_center << " " << e_ << " " <<  p_<< endl;
    // }
    // double en_max= enerv.front();
    // double en_min= enerv.back();
    // double de= (en_max- en_min)/iener;

    // double pr_max= presv.front();
    // double pr_min= presv.back();
    // double dp= (pr_max -pr_min)/iener;
    // for(int ie=0; ie<iener; ie++){
    //  // cout << e_center << endl;
        
    //     double e_center= en_max-ie*de;
    //     double p_center=  interpolation_func(e_center, presv, enerv);// interpolatePres(e_center);
        
    //     double p_ = pr_max- ie*dp;
    //     double e_= interpolation_func(p_, enerv, presv);
    //     // double e_= interpolateEner(p_);

    //     eos_test << e_center <<  " " <<  p_center << " " << e_ << " " <<  p_<< endl;
    // }

    // eos_test.close();

    int iener= 100;
    double en_max= enerv.back();
    double en_min= enerv.front();
    double de= (en_max- en_min)/iener;
    

    cout << "Calculating Mass-Radius diagram" << endl;

    vector<double> massv, radiusv;
    //cout << en_max << " " << en_min << endl;
    for(int ie=0; ie<iener; ie++){

        double e_center= en_max-ie*de;

        tov.solve_tov_euler(e_center, dr);
        // if(ie>0 && ie < iener)tov.solve_tidal_euler(e_center, dr);
        double mass=tov.getMass();
        double radius=tov.getRadius();
        double compactness=  tov.getCompactness();
        // double  y= tov.gety();
        // double  beta= tov.getBeta();
        // double  h= tov.getH();
        // double Lambda= tov.getLambda();
        cout << radius << " " << mass << " " << e_center  << " " <<endl; 
        //<< Lambda << " "            << compactness <<  " " << y << " " << beta << " " << h 

        outMR << radius << " " << mass << " " << e_center << " " << endl; 
        massv.push_back(mass);
        radiusv.push_back(radius);
        //<< Lambda << " "            << compactness << " " <<  y << " " << beta << " " << h << endl;
    }

    double Mmax=  *max_element(massv.begin(), massv.end()) ;
    // auto ir= find_if(massv.begin(), massv.end(), [Mmax](double mmax_) { return abs(Mmax -mmax_) < 1e-3; });
    auto ir= find(massv.begin(), massv.end(), Mmax );
    cout << "Mmax: " << *max_element(massv.begin(), massv.end()) << " " << ir-massv.begin() << endl;
    cout << "R( Mmax) : " << radiusv[ir-massv.begin()] << endl;
    outMR.close();

    // ofstream outMR2("data/mr_diag2.txt");

    // tov_class_rk4 tovrk4(enerv, presv);

    
    //    for(int ie=0; ie<iener; ie++){

    //     double e_center= en_max-ie*de;
    //     double p_center=interpolation_func(e_center, presv, enerv);

    //     double mass=0., radius=0., nu=0.;
    //     double dr=1e2; //fm
    //     state_type x = { mass, p_center }; // initial conditions

        
    //     runge_kutta4< state_type > stepper;
    //     radius=dr;
    //      do{
            
    //         stepper.do_step( tovrk4 , x , radius , dr );
    //         radius+=dr;
    //     }while(x[1]> presv[0]   );

    //     mass=x[0];

    //     double compactness= mass/radius;

    //     radius*= fm_to_cm*1e-5; //fm to km

    //     mass*= hc*1.73e-30; // 1/fm to MeV and MeV to kg

    //     outMR2 << radius << " " << mass/Msun << " " << e_center << " " << p_center  << " " << compactness << endl;
    // }
    // outMR2.close();
    return 0;
}
