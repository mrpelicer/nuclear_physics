#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H 

#include <gsl/gsl_interp.h>	
#include <gsl/gsl_spline.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double interpolation_func(double x_, std::vector<double> yv_,	std::vector<double> xv_);

double deriv_func(double x_, std::vector<double> yv_,	std::vector<double> xv_);

double deriv2_func(double x_, std::vector<double> yv_,	std::vector<double> xv_);

double interpolate(double x_, const std::vector<double> &yv_, const std::vector<double> &xv_);

#endif