#include "interpolator.h"


double interpolation_func(double x_, std::vector<double> yv_,	std::vector<double> xv_)
{
    double res=0.;
		assert(yv_.size()==xv_.size());

    if(x_<xv_.front() || x_>xv_.back()){
      res=NAN;
      cout << "INTERPOLATOR ERROR:  " << x_ << " " << xv_.back() << " " << xv_.front() << endl;
    }else{
      gsl_interp_accel *acc =  gsl_interp_accel_alloc();
      gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_cspline,	xv_.size());

      gsl_interp_init(interpolation, xv_.data(), yv_.data(), xv_.size()); //vector.data() transforms std::vector<double> in an array!

      res = gsl_interp_eval(interpolation, xv_.data(), yv_.data(), x_, acc);

      gsl_interp_free(interpolation);
      gsl_interp_accel_free (acc);
    }
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

double interpolate(double x_, const std::vector<double> &yv_, const std::vector<double> &xv_)
{
  //Local variables
  double y_;
  int k,l;

//Find the index of the element of xv_ that is nearest-above to x
auto i = lower_bound(xv_.begin(), xv_.end(), x_); 
//If the vector values are in decreasing order use:
//auto i = lower_bound(xv_.rbegin(), xv_.rend(), x);
k = i - xv_.begin(); //Nearest index
if (i == xv_.end())
  --k;  // extrapolating above
else if (*i == x_)
  return yv_[k];

l = k? k - 1: 1; //nearest-below index, except when extrapolating downward

//Interpolation:
  if(xv_[k]<xv_[l]) 
      y_ = yv_[k]+(x_-xv_[k])*(yv_[l]-yv_[k])/(xv_[l]-xv_[k]);
  else 
      y_ = yv_[l]+(x_-xv_[l])*(yv_[k]-yv_[l])/(xv_[k]-xv_[l]);

  return y_;
}

