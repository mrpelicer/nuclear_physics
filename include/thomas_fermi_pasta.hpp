#pragma once
#include <iostream>
#include <string>
#include "rmf_walecka.hpp"
#include "pasta.hpp"
#include <gsl/gsl_sf_gamma.h>
#include <Eigen/Dense>



struct harmonic_integral_parameters{
    bool use_r_jacobian=true;
    std::vector<double> function;
    std::vector<double> r;
    int n;
};

class tf_walecka_class : public nlwm_class{
public:

    tf_walecka_class(std::string parametrization_);

    void solve_Thomas_Fermi(double rhoB, double Yp, double temperature, double radius_WS_init, int idim);

    void setInitialConditions(double Z_, double N_);

    double Phi_Oscilator(int n,  double x_, double b_);

    // std::vector<double> calculate_sigma_field(std::vector<double> coef_meson);
    // std::vector<double> calculate_delta_field(std::vector<double> coef_meson);
    // std::vector<double> calculate_omega_field(std::vector<double> coef_meson);
    // std::vector<double> calculate_rho_field(std::vector<double> coef_meson);

    std::vector<double> calculate_kf(double mu_, particle &particle_);
    std::vector<double> calculate_meson_field_coef(std::vector<double> coef_meson, double m2, double g, double g3, double g4);
    
    void calculate_meson_fields();
    Eigen::MatrixXd get_H_matrix(double m2);


    double Laguerre(int n, double alpha, double x);
    std::vector<double> sigma_meson_tf, delta_meson_tf, omega_meson_tf, rho_meson_tf, phi_meson_tf;

    double integrate_tf(double (func)(double, void *), void *parametersPointer);
    std::vector<double> get_function(std::vector<double> coefficients);
    int NL_points=100;
    
    double b_oscilator=3.*Mnucleon/hc;
    bool firstRun=true;
    // std::vector<double> r;
    std::vector<double> r={2.6093471482799946e-2,
  0.13493663331100003,
  0.32059043170099999,
  0.56660460587080008,
  0.85112566101840004,
   1.1488743389816001,
   1.4333953941291999,
   1.6794095682990000,
   1.8650633666890000,
   1.9739065285172002,
   2.0260934714827998,
   2.1349366333109998,
   2.3205904317010000,
   2.5666046058708001,
   2.8511256610184001,
   3.1488743389815999,
   3.4333953941291999,
   3.6794095682990000,
   3.8650633666890002,
   3.9739065285172002,
   4.0260934714827998,
   4.1349366333109998,
   4.3205904317009995,
   4.5666046058707996,
   4.8511256610183997,
   5.1488743389816003,
   5.4333953941292004,
   5.6794095682990005,
   5.8650633666890002,
   5.9739065285172002,
   6.0260934714827998,
   6.1349366333109998,
   6.3205904317009995,
   6.5666046058707996,
   6.8511256610183997,
   7.1488743389816003,
   7.4333953941292004,
   7.6794095682990005,
   7.8650633666890002,
   7.9739065285172002,
   8.0260934714828007,
   8.1349366333109998,
   8.3205904317009995,
   8.5666046058708005,
   8.8511256610184006,
   9.1488743389815994,
   9.4333953941291995,
   9.6794095682990005,
   9.8650633666890002,
   9.9739065285171993,
   10.026093471482801,
   10.134936633311000,
   10.320590431701000,
   10.566604605870801,
   10.851125661018401,
   11.148874338981599,
   11.433395394129199,
   11.679409568299000,
   11.865063366689000,
   11.973906528517199,
   12.026093471482801,
   12.134936633311000,
   12.320590431701000,
   12.566604605870801,
   12.851125661018401,
   13.148874338981599,
   13.433395394129199,
   13.679409568299000,
   13.865063366689000,
   13.973906528517199,
   14.026093471482801,
   14.134936633311000,
   14.320590431701000,
   14.566604605870801,
   14.851125661018401,
   15.148874338981599,
   15.433395394129199,
   15.679409568299000,
   15.865063366689000,
   15.973906528517199,
   16.026093471482799,
   16.134936633311000,
   16.320590431701000,
   16.566604605870801,
   16.851125661018401,
   17.148874338981599,
   17.433395394129199,
   17.679409568299000,
   17.865063366689000,
   17.973906528517201,
   18.026093471482799,
   18.134936633311000,
   18.320590431701000,
   18.566604605870801,
   18.851125661018401,
   19.148874338981599,
   19.433395394129199,
   19.679409568299000,
   19.865063366689000,
   19.973906528517201
   };

    int Nr_points=r.size(); //150.

    double N, Z; //change to int
    int idim;
    double volume;
    double radius;
    double jacobian;


    //These have size Nr_points
    std::vector<double> sigmav, deltav, omegav, rhov;
    std::vector<double> rhoBv, rhoSv, rho3v, rhoS3v, rhopv, rhonv, rhoev;
    std::vector<double> Energyv, Pressurev, Entropyv;
    std::vector<double> kfnv, kfpv;

    //These have size NL_points

    
    harmonic_integral_parameters integral_parameters;

};

double harmonic_integration_function(double r_, void * params);

double _integration_function(double r_, void * params); 

// double zriddr(double (*func)(double, void *), double x1, double x2, double xacc, void *p, int MAXIT= 60, double UNUSED= -1.11E-30);

// double diff_N(double x, void *params);

struct NFunctor{
public:
	NFunctor(tf_walecka_class & tf_):tf(tf_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		tf_walecka_class & tf;
};

struct PFunctor{
public:
	PFunctor(tf_walecka_class & tf_):tf(tf_)
	{}

  template <typename T>
  bool operator()(const T* x, T* residuals) const;
	
  private:
		tf_walecka_class & tf;
};