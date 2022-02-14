#include "interpolator.h"

double interpolation_func(double x_, std::vector<double> yv_,	std::vector<double> xv_)
{
    double res=0.;
		assert(yv_.size()==xv_.size());

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,	xv_.size());

    gsl_interp_init(interpolation, xv_.data(), yv_.data(), xv_.size()); //vector.data() transforms std::vector<double> in an array!

    res = gsl_interp_eval(interpolation, xv_.data(), yv_.data(), x_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;
}

double deriv_func(double x_, std::vector<double> yv_,	std::vector<double> xv_)
{
    double res=0.;
		assert(yv_.size()==xv_.size());

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,	xv_.size());

    gsl_interp_init(interpolation, xv_.data(), yv_.data(), xv_.size()); //vector.data() transforms std::vector<double> in an array!

    res = gsl_interp_eval_deriv(interpolation, xv_.data(), yv_.data(), x_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;
}

double deriv2_func(double x_, std::vector<double> yv_,	std::vector<double> xv_)
{
    double res=0.;
		assert(yv_.size()==xv_.size());

    gsl_interp_accel *acc =  gsl_interp_accel_alloc();
    gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,	xv_.size());

    gsl_interp_init(interpolation, xv_.data(), yv_.data(), xv_.size()); //vector.data() transforms std::vector<double> in an array!

    res = gsl_interp_eval_deriv2(interpolation, xv_.data(), yv_.data(), x_, acc);

    gsl_interp_free(interpolation);
    gsl_interp_accel_free (acc);

    return res;
}