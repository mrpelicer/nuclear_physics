#include "thomas_fermi_pasta.hpp"

tf_walecka_class::tf_walecka_class(std::string parametrization_) : 
    nlwm_class(parametrization_)
{}

void tf_walecka_class::solve_Thomas_Fermi(double rhoB_, double Yp_, double temperature_, double total_radius_, int idim_){

    rhoB=rhoB_;
    Yp=Yp_;
    setTemperature(temperature_);
    radius= total_radius_;
    idim=idim_;

    if(idim==3 && NL_points>170)
        NL_points=170; // maximum possible value before gamma function overflows

    if(idim==1 && NL_points>86)
    NL_points=86; // maximum possible value before gamma function overflows    
    // else if(idim==2 && NL_points>100)

    // if(idim==1){
    //     if(NL_points>80 )NL_points=50;
    // }

    for(auto &element : r)
        element *= Mnucleon/hc;

    // double r_min=1e-5*Mnucleon/hc;
    // // double r_min=0.0261*Mnucleon/hc;
    // for(int ir=0; ir<Nr_points; ir++){
    //     r.push_back(r_min+ (radius-r_min)*ir/Nr_points);
    // } 
    
    std::cout << r.size() << " " << Nr_points << endl;    
    integral_parameters.r= r;

    double xL= 10.*Mnucleon/hc; //results should be independent of this

    
    if(idim==3){
        volume=4.*M_PI*pow(radius, 3)/3.;
        cout << "volume= " << volume << " " << radius <<  " " << pow(radius, 3) << M_PI << endl;
        jacobian=4.*M_PI;
    }else if(idim==2){
        volume= M_PI*pow(radius, 2.)*xL;   
        jacobian=2.*M_PI*xL;
    }else if(idim==1){
        volume= 2.*radius*xL*xL ;
        jacobian=xL*xL;
    }else{
        cout << "idim not defined in solve_Thomas_Fermi()!" << endl;
        exit(1);
    }

    N=(1.-Yp)*rhoB*volume;
    Z=Yp*rhoB*volume;

    // cout << "n, y, V= " << rhoB << " " << Yp << " " << volume << " " << radius<< endl;
    cout << "N,Z= " << N << " " << Z << endl;
    setInitialConditions(Z, N);

    FILE *frhotest;
    frhotest=fopen("rhotest.dat", "w");


    for(auto i=0; i<Nr_points; i++){
        fprintf(frhotest, "%.12E %.12E %.12E %.12E %.12E %.12E\n", r[i]*hc/Mnucleon, rhoBv[i]*pow(Mnucleon/hc, 3.),
                    rho3v[i]*pow(Mnucleon/hc, 3.),rhoSv[i]*pow(Mnucleon/hc, 3.), rhoS3v[i]*pow(Mnucleon/hc, 3.),
                    rhoev[i]*pow(Mnucleon/hc, 3.));
    }
    fclose(frhotest);

    //Calculate the initial meson fields
    calculate_meson_fields();

    FILE *fmeson;
    fmeson=fopen("meson_test.dat", "w");


    for(int i=0; i<Nr_points; i++){
        // fprintf(fmeson, "%.12E %.12E %.12E %.12E %.12E\n", r[i]*hc/Mnucleon, sigmav[i]*hc,
            //  omegav[i]/rhoBv[i], rhov[i]*hc, deltav[i]*hc);
            fprintf(fmeson, "%.12E %.12E %.12E\n", r[i]*hc/Mnucleon,
     omegav[i]*Mnucleon/hc, rhov[i]*hc);
    }

    fclose(fmeson);


    double xi2=1.;
    double del=1e-14;
    std::vector<double> rhoBv_previous=rhoBv;
    std::vector<double> rhoBv_updated=rhoBv;


    do{
        xi2=0;

	    Problem NProblem;
        Problem PProblem;
        if(firstRun){
            proton.chemPot= 4.719*Mnucleon/hc;
            neutron.chemPot=4.719*Mnucleon/hc;
            firstRun=false;
        }

        double xn[]={neutron.chemPot};
        double xp[]={proton.chemPot};
        CostFunction* NCost =	new NumericDiffCostFunction<NFunctor, ceres::CENTRAL, 1, 1>
	    													(new NFunctor(*this) );
        CostFunction* PCost =	new NumericDiffCostFunction<PFunctor, ceres::CENTRAL, 1, 1>
	    													(new PFunctor(*this) );                                            


        NProblem.AddResidualBlock(NCost, NULL, xn);
        NProblem.SetParameterLowerBound(xn, 0,0.);

        PProblem.AddResidualBlock(PCost, NULL, xp);
        PProblem.SetParameterLowerBound(xp, 0,0.);

        Solver::Options options;
        options.parameter_tolerance = 1e-10;
        options.function_tolerance = 1e-10;
        options.gradient_tolerance=1e-10;
        //options.trust_region_strategy_type = ceres::DOGLEG;
        //options.dogleg_type = ceres::TRADITIONAL_DOGLEG;
		options.dense_linear_algebra_library_type=ceres::LAPACK;
		options.linear_solver_type= ceres::DENSE_QR;
        options.line_search_direction_type = ceres::LBFGS;


        options.minimizer_progress_to_stdout = true;

        Solver::Summary summary;
        Solve(options, &NProblem, &summary);
        // std::cout <<summary.BriefReport() << std::endl;
        Solve(options, &PProblem, &summary);
        // std::cout <<summary.BriefReport() << std::endl;

        std::cout << "mun= " <<  xn[0]*Mnucleon/hc << std::endl;
        std::cout << "mup= " <<  xp[0]*Mnucleon/hc << std::endl;


        kfnv=calculate_kf(xn[0], neutron);
        kfpv=calculate_kf(xp[0], proton);
// 
        // for(auto kf_: kfnv)
            // std::cout << "kf_n= " << kf_*Mnucleon/hc << std::endl;
// 
        // for(auto kf_: kfpv)
            // std::cout << "kf_p= " << kf_*Mnucleon/hc << std::endl;

        // exit(1);


        //weight given by the new and old densities when updating : 0.5= avg over them
        double MIX_=0.5;

        neutron.chemPot= xn[0];
        proton.chemPot= xp[0];

        for(int i=0; i<Nr_points; i++){


            neutron.kf= kfnv[i];
            proton.kf= kfpv[i];

            neutron.mass_eff= neutron.mass - gs*sigmav[i] - gd*2.*neutron.I3*deltav[i];
            proton.mass_eff = proton.mass - gs*sigmav[i]  - gd*2.*proton.I3*deltav[i];

            neutron.chemPot_eff= xn[0] - gv*omegav[i] - neutron.I3*gr*rhov[i];
            proton.chemPot_eff= xp[0] - gv*omegav[i] - proton.I3*gr*rhov[i];

            neutron.calculateDensity();
            proton.calculateDensity();

            neutron.calculateCondensate();
            proton.calculateCondensate();

            // cout << "rhoS= " << i << " " << (neutron.condensate + proton.condensate)*pow(Mnucleon/hc, 3.) << endl;


            rhoSv[i]= (1.-MIX_)*rhoSv[i] + MIX_*(neutron.condensate + proton.condensate);
            rhoS3v[i]= (1.-MIX_)*rhoS3v[i] + MIX_*2.*(neutron.I3*neutron.condensate + proton.I3*proton.condensate);

            rhoBv[i]= (1.- MIX_)*rhoBv[i] + MIX_*(neutron.density + proton.density);
            rho3v[i]= (1.- MIX_)*rho3v[i] + MIX_*(neutron.I3*neutron.density + proton.I3*proton.density);

            rhopv[i]= (1.- MIX_)*rhopv[i] + MIX_*proton.density;
            rhonv[i]= rhoBv[i] - rhopv[i];//(1.- MIX_)*rhonv[i] + MIX_*neutron.density;
            // rhoev[i]= (1.- MIX_)*rhoev[i] + MIX_*electron.density;
        }

        rhoBv_updated=rhoBv;

        // integral_parameters.function= rhonv;
        // double N_new= jacobian*integrate_tf(_integration_function, this);

        // integral_parameters.function= rhopv;
        // double Z_new= jacobian*integrate_tf(_integration_function, this);

    //     cout << "rhon, rhop:" << endl; 
    //     for(int i=0; i<Nr_points; i++){
    //         cout << rhonv[i]*pow(Mnucleon/hc, 3) << " " << rhopv[i]*pow(Mnucleon/hc, 3) << endl;
    //     }
    //     cout << N_new << " " << Z_new << endl;
        
    // std::cout << "Press Enter to continue..." << std::endl;
    // std::cin.get();        
    calculate_meson_fields();

        for(int i=0; i<Nr_points; i++){
           xi2+=pow( (rhoBv_updated[i]-rhoBv_previous[i]) , 2.);
           rhoBv_previous[i] = rhoBv_updated[i];
        }

        cout << "xi2= " << xi2 << endl;

    }while(xi2>del);


    FILE* fconverged;
    fconverged= fopen("converged.dat", "w");
    std::vector<double> sigma_integrand;
    std::vector<double> sigma_Centelles_integrand;

        for(int i=0; i<Nr_points; i++){

            neutron.kf= kfnv[i];
            proton.kf= kfpv[i];

            sigma_meson= sigmav[i];
            delta_meson=deltav[i];
            omega_meson=omegav[i];
            rho_meson=rhov[i];

            neutron.mass_eff= neutron.mass - gs*sigmav[i] - gd*2.*neutron.I3*deltav[i];
            proton.mass_eff = proton.mass - gs*sigmav[i]  - gd*2.*proton.I3*deltav[i];

            neutron.chemPot_eff= neutron.chemPot - gv*omegav[i] - neutron.I3*gr*rhov[i];
            proton.chemPot_eff = proton.chemPot - gv*omegav[i] - proton.I3*gr*rhov[i];

            neutron.calculateDensity();
            proton.calculateDensity();

            neutron.calculateCondensate();
            proton.calculateCondensate();

            neutron.calculateProperties();
            proton.calculateProperties();

            Energyv.push_back(getEnergy());
            Pressurev.push_back(getPressure_thermodynamic_relation());

            // double e2_= proton.energy+neutron.energy+ 0.5*gv*omegav[i]*rhoBv[i]
            //                                         + 0.5*gr*rhov[i]*rho3v[i]
            //                                         + 0.5*gs*sigmav[i]*rhoSv[i]
            //                                         + 0.5*gd*deltav[i]*rhoS3v[i]
            //                                         +gv4*pow(gv*omegav[i], 4.)/24.
            //                                         -gs3*pow(sigmav[i], 3.)/12.
            //                                         -gs4*pow(sigmav[i], 4.)/24.; // same as getEnergy()

            double dsigmadr = deriv_func(r[i], sigmav, r);
            double ddeltadr = deriv_func(r[i], deltav, r);
            double domegadr = deriv_func(r[i], omegav, r);
            double drhodr   = deriv_func(r[i], rhov, r);
           
            sigma_integrand.push_back(  dsigmadr*dsigmadr +ddeltadr*ddeltadr 
                                    - domegadr*domegadr - drhodr*drhodr );



            fprintf(fconverged,"%.12E %.12E %.12E %.12E %.12E %.12E %.12E %.12E\n", 
                    r[i]*hc/Mnucleon, rhoBv[i]*pow(Mnucleon/hc, 3), rhoSv[i]*pow(Mnucleon/hc, 3), 
                    rhopv[i]*pow(Mnucleon/hc, 3), rhonv[i]*pow(Mnucleon/hc, 3),
                    Energyv[i]*Mnucleon*pow(Mnucleon/hc, 3),
                    Pressurev[i]*Mnucleon*pow(Mnucleon/hc, 3)   ,
                    sigma_integrand[i]
                    );

        }

        for(int i=0; i<Nr_points; i++){            
            sigma_Centelles_integrand.push_back(Pressurev[i] - Pressurev[Nr_points-1] );
            //Getting reasonable results only if I assume the opposite sign.
        }

    for(int i=0; i<Nr_points; i++){
        cout << r[i]*hc/Mnucleon << " " << rhoBv[i]*pow(Mnucleon/hc, 3.) << " " << rhonv[i]*pow(Mnucleon/hc, 3.) << " " << rhopv[i]*pow(Mnucleon/hc, 3.) << endl;
    }

    fclose(fconverged);


    integral_parameters.function= rhonv;
    double N_= jacobian*integrate_tf(_integration_function, this);

    integral_parameters.function= rhopv;
    double Z_= jacobian*integrate_tf(_integration_function, this);

    double A_ = N_+Z_;
    integral_parameters.function= Energyv;
    double total_ener= jacobian*integrate_tf(_integration_function, this);
    std::cout <<" E= " << ((total_ener-A_*Mn/Mnucleon)/A_)*Mnucleon << endl;

    //PhysRevC.78.015802
    integral_parameters.use_r_jacobian=false;
    integral_parameters.function= sigma_integrand;
    double sigma_=integrate_tf(_integration_function, this);

    //Nucl Phys. A635(1998)193
    integral_parameters.function= sigma_Centelles_integrand;
    double sigma_centelles= integrate_tf(_integration_function, this);
    integral_parameters.use_r_jacobian=true;

    std::cout << "sigma= " << sigma_*Mnucleon*pow(Mnucleon/hc, 2) << std::endl;
    std::cout << "sigma_c= " << sigma_centelles*Mnucleon*pow(Mnucleon/hc, 2)/2. << std::endl;
    std::cout << "Converged: " << Z << " " << N << " " << Z_ <<  " " << N_ << endl;


    std::cout << " mu's: " << neutron.chemPot*Mnucleon/hc << " " << proton.chemPot*Mnucleon/hc << endl;




    cout << "fim " << endl;


}


template <typename T>
bool NFunctor::operator()(const T* x, T* residuals) const{

    double mun=x[0];
    std::vector<double> kfv=tf.calculate_kf(mun, tf.neutron);
    std::vector<double> rhonv_;

    for(auto kf_: kfv)
        rhonv_.push_back(tf.neutron.gamma*kf_*kf_*kf_/6./M_PI/M_PI);

    double N_cell=0.;
    tf.integral_parameters.function= rhonv_;
    N_cell= tf.jacobian*tf.integrate_tf(_integration_function, &tf);

    residuals[0]=tf.N - N_cell;

    return true;
}

template <typename T>
bool PFunctor::operator()(const T* x, T* residuals) const{
   
    double mup=x[0];
    std::vector<double> kfv=tf.calculate_kf(mup, tf.proton);
    std::vector<double> rhopv_;

    for(auto kf_: kfv)
        rhopv_.push_back(tf.proton.gamma*kf_*kf_*kf_/6./M_PI/M_PI);

    double Z_cell=0.;
    tf.integral_parameters.function= rhopv_;        
    Z_cell= tf.jacobian*tf.integrate_tf(_integration_function, &tf);

    residuals[0]=tf.Z - Z_cell;

    return true;
}

void tf_walecka_class::setInitialConditions(double Z, double N){
    
    for(int i=0; i<Nr_points; i++){
        double r0= 1.2*Mnucleon/hc ;
        double R_= r0*pow(N+Z, 1./3.); //r0* A^1/3 
        double a_=0.6*Mnucleon/hc;
        
        double rhoc=0.1/pow(Mnucleon/hc, 3.);
        double rhob_=rhoc/(1.+exp( (r[i]-R_)/a_)) ;
        rhoBv.push_back(rhob_);
        rhopv.push_back(Yp*rhob_);
        rhonv.push_back((1.-Yp)*rhob_);

        rho3v.push_back(proton.I3*Yp*rhob_ + neutron.I3*(1.-Yp)*rhob_);// /2
        rhoSv.push_back(0.95*rhob_);
        rhoS3v.push_back(0.95*(1.-2.*Yp)*rhob_);//2 ???

        rhoev.push_back(3.*Z/(4.*M_PI*pow(radius, 3.)) );
    }

}

std::vector<double> tf_walecka_class::calculate_kf(double mu_, particle &particle_){

    std::vector<double> kfv;
    for(int i=0; i<Nr_points; i++){
        double mass_eff= particle_.mass - gs*sigmav[i] - 2.*particle_.I3*gd*deltav[i];
        double mu_eff= mu_ - gv*omegav[i] - particle_.I3*gr*rhov[i];

        double kf2= mu_eff*mu_eff - mass_eff*mass_eff;

        double kf_= kf2>0.? sqrt(kf2) : 0.;
        kfv.push_back(kf_);
    }

    return kfv;

}


double tf_walecka_class::Phi_Oscilator(int n, double x_, double b_){
    
    if(n<1) return 0.;
    double Phi_n=0.;
    if(idim==3){

        double alpha= 1./2.;
        double norm_n= sqrt(2.*gsl_sf_gamma(n)/gsl_sf_gamma((double)n + .5));
        Phi_n=norm_n*Laguerre(n-1, alpha, x_*x_/b_/b_)*exp(-x_*x_/2./b_/b_)/pow(b_, 3./2.);

    }else if(idim==2){
        
        double alpha=0.;
        double norm_n= sqrt(2);
        Phi_n=norm_n*Laguerre(n-1, alpha, x_*x_/b_/b_)*exp(-x_*x_/2./b_/b_)/b_;
    
    }else if(idim==1){

        double alpha=-1./2.;
        double norm_n= pow(-2, n-1)*gsl_sf_gamma((double)n)/(sqrt(gsl_sf_gamma(2*n-1) )*pow(M_PI, 1./4.));
        Phi_n=norm_n*Laguerre(n-1, alpha, x_*x_/b_/b_)*exp(-x_*x_/2./b_/b_)/sqrt(b_);

    }else{
        cout << "idim not defined in Phi_Oscilator()!" << endl;
        exit(1);
    }
    return Phi_n;
}


std::vector<double> tf_walecka_class::get_function(std::vector<double> coefficients){
    std::vector<double> result;

    for(int i=0; i<Nr_points; i++){
        double function_=0.;
        for(auto n=0; n<NL_points; n++){
            function_+= coefficients[n]*Phi_Oscilator(n+1, r[i], b_oscilator);
        }
        result.push_back(function_);
    }

    return result;
}


double tf_walecka_class::Laguerre(int n, double alpha_,  double x_){
        
    //for n=0, for any alpha: 
    double L0=1.;
    double L1= 1.+alpha_-x_;

    if(n==0) return L0;
    else if(n==1) return L1;
    else{
        double L_=0;    
        double Lmm  =L0;
        double Lm   =L1;

        for(auto i=2; i<=n; i++){
            L_= ( (2.*(i-1.)+1. + alpha_ - x_)*Lm - (i-1.+alpha_)*Lmm)/((double)i);
            Lmm=Lm;
            Lm=L_;
        }
        return L_;
    }
}


void tf_walecka_class::calculate_meson_fields(){
    
    std::vector<double> coef_sigmav_source, coef_deltav_source, coef_omegav_source, coef_rhov_source;
    std::vector<double> coef_sigmav, coef_deltav, coef_omegav, coef_rhov;
    // sigma:
    // cout << " //sigma" << endl;

    //define source:
    integral_parameters.function= rhoSv;
    for (double& element :integral_parameters.function)
        element *= gs;
    
    //calculate source coefficients
    for(auto n=0; n<NL_points; n++){
        integral_parameters.n = n;
        coef_sigmav_source.push_back(integrate_tf(harmonic_integration_function, this ));
    }

    //calculate field coefficients
    coef_sigmav= calculate_meson_field_coef(coef_sigmav_source, Ms*Ms, gs, gs3 , gs4);

    //calculate field as a function of r
    sigmav= get_function(coef_sigmav);
 

 // hence and repeat:

    //delta:
    // cout << " //delta" << endl;
    integral_parameters.function= rhoS3v;

    for (double& element :integral_parameters.function)
        element *= gd;
    
    for(auto n=0; n<NL_points; n++){
        integral_parameters.n = n;
        coef_deltav_source.push_back(integrate_tf(harmonic_integration_function, this ));
    }
    coef_deltav= calculate_meson_field_coef(coef_deltav_source, Md*Md, gd, 0., 0.);
    //calculate delta(r)
    deltav=get_function(coef_deltav);
    
    
    // vector mesons:
    if(fabs(Lvr)<1e-12){ //calculate them separately:

        // cout << " //omega" << endl;
        integral_parameters.function= rhoBv;
        for (double& element :integral_parameters.function)
            element *= gv; // multiply the density by the coupling in source function
        
        for(int n=0; n<NL_points; n++){
            integral_parameters.n = n;
            coef_omegav_source.push_back(integrate_tf(harmonic_integration_function, this ));
        }
       
        coef_omegav=calculate_meson_field_coef(coef_omegav_source, Mv*Mv, gv, 0.0 , 0.0);

        //calculate omega(r)
        omegav= get_function(coef_omegav);
        
        
        // cout << " //rho" << endl;
        integral_parameters.function= rho3v;
        for (double& element :integral_parameters.function)
            element *= gr;

        
        for(auto n=0; n<NL_points; n++){
            integral_parameters.n = n;
            coef_rhov_source.push_back(integrate_tf(harmonic_integration_function, this ));
        }
        coef_rhov= calculate_meson_field_coef(coef_rhov_source, Mr*Mr, gr, 0.0 , 0.0);

        //calculate rho(r)
        rhov=get_function(coef_rhov);        
    }else{
        std::vector<std::vector<double>> vector_meson_coefs;
        
    }


}

Eigen::MatrixXd tf_walecka_class::get_H_matrix(double m2){
    
    Eigen::MatrixXd H_matrix = Eigen::MatrixXd::Zero(NL_points,NL_points);
    //H*meson = source in linear case

    if(idim==3){
        //Define first and last columns of matrix H:
        H_matrix(0,0)=3./2./b_oscilator/b_oscilator + m2;
        H_matrix(0,1)=sqrt(3./2.)/b_oscilator/b_oscilator;

        H_matrix(NL_points-1, NL_points-2)= sqrt( (NL_points-1.)*(NL_points-1./2.))/b_oscilator/b_oscilator;
        H_matrix(NL_points-1, NL_points-1)= (2.*(NL_points-1.)+3./2.)/b_oscilator/b_oscilator + m2;

        for(int i=1; i<NL_points-1; i++){
            double n= (double) i+1.;
            H_matrix(i, i-1)= sqrt((n-1.)*(n-1./2.))/b_oscilator/b_oscilator;
            H_matrix(i, i)= (2.*(n-1.)+3./2.)/b_oscilator/b_oscilator + m2;
            H_matrix(i, i+1)= sqrt( n*(n+1./2.))/b_oscilator/b_oscilator;

            // cout << H_matrix(i, i-1)*pow(<< " " << H_matrix(i, i) << " " << H_matrix(i, i+1) << endl;
        }
    }else if(idim==2){
        H_matrix(0,0)=1./b_oscilator/b_oscilator + m2;
        H_matrix(0,1)=1./b_oscilator/b_oscilator;

        H_matrix(NL_points-1, NL_points-2)= ( NL_points-1. )/b_oscilator/b_oscilator;
        H_matrix(NL_points-1, NL_points-1)= (2.*NL_points-1.)/b_oscilator/b_oscilator + m2;

        for(int i=1; i<NL_points-1; i++){
            double n= (double) i+1.;
            H_matrix(i, i-1)= (n-1.)/b_oscilator/b_oscilator;
            H_matrix(i, i)= (2.*n-1.)/b_oscilator/b_oscilator + m2;
            H_matrix(i, i+1)= ((double)n)/b_oscilator/b_oscilator;

            // cout << H_matrix(i, i-1)*pow(<< " " << H_matrix(i, i) << " " << H_matrix(i, i+1) << endl;
        }

    }else if(idim==1){
        H_matrix(0,0)=1./2./b_oscilator/b_oscilator + m2;
        H_matrix(0,1)=-sqrt(2.)/2./b_oscilator/b_oscilator;

        H_matrix(NL_points-1, NL_points-2)= - sqrt((2.*NL_points -2.)*(2.*NL_points -3.))/2./b_oscilator/b_oscilator;
        H_matrix(NL_points-1, NL_points-1)= (4.*NL_points-3.)/2./b_oscilator/b_oscilator + m2;

        for(int i=1; i<NL_points-1; i++){
            double n= (double) i+1.;
            H_matrix(i, i-1)= -sqrt(2.*(n-1.)*(2.*n-3.))/2./b_oscilator/b_oscilator;
            H_matrix(i, i)= (4.*n-3.)/2./b_oscilator/b_oscilator + m2;
            H_matrix(i, i+1)= -sqrt(2.*n*(2.*n-1.))/2./b_oscilator/b_oscilator;

            // cout << H_matrix(i, i-1)*pow(<< " " << H_matrix(i, i) << " " << H_matrix(i, i+1) << endl;
        }
    }

    return H_matrix;
}
std::vector<double> tf_walecka_class::calculate_meson_field_coef(std::vector<double> source_coeffs,double m2, 
                                                                double g, double g3, double g4 ){

    //copy std::vector to eigen::vector
    Eigen::VectorXd coefficients = Eigen::Map<Eigen::VectorXd>(source_coeffs.data(), source_coeffs.size());

    //coeff.size() is =  NL_points

    // //Define matrix H, where H*meson_coeffs= source_coeffs
    Eigen::MatrixXd H_matrix = get_H_matrix(m2);
    Eigen::MatrixXd H_inverse= H_matrix.inverse();

    //get meson field coefficients
    Eigen::VectorXd meson_eigen_coefs= H_inverse*coefficients;

    std::vector<double> mesonv_coefs(meson_eigen_coefs.data(), meson_eigen_coefs.data() + meson_eigen_coefs.size());
    // cout << "mesonv_coefs[n] "<< endl;
    // for(int n=0; n<NL_points; n++){
    //     cout << n+1 << " " << mesonv_coefs[n]*pow(hc/Mnucleon, 1./2.) << endl;
    // }

    // cout << "source_coefs[n] "<< endl;
    // for(int n=0; n<NL_points; n++){
    // std::cout << n+1 << " " <<  source_coeffs[n]*pow(Mnucleon/hc, 3./2.)<< endl;
    // }
    // std::vector<double> mesonv= get_function(mesonv_coefs);
    
    // for(int i=0; i<Nr_points; i++){
    //     cout << r[i]*hc/Mnucleon << " " << mesonv[i]*Mnucleon/hc << endl;
    // }
    // cout << endl<< endl<< endl<< endl;

    int icounter=0;

    if(fabs(g3)<1e-12 &&  fabs(g4)<1e-12)   return mesonv_coefs;
    else{
        
        std::vector<double> source_coeffs_nonlinear= source_coeffs;
        
        cout << "//calculate non-linear terms" << endl;
        //calculate non-linear terms
        double xi2=1.;
        double del=1.e-10;

        //make ceres solver? - not worth it!
        // std::vector<double> mesonv_nonlinear;
        do{

            xi2=0.;
            
            //calculate meson^2 and meson^3 at r[i]:
            std::vector<double> mesonv_= get_function(mesonv_coefs);
            std::vector<double> meson2, meson3;
            for(double meson_r: mesonv_){
                meson2.push_back(meson_r*meson_r);
                meson3.push_back(meson_r*meson_r*meson_r);
                // cout << meson_r*Mnucleon/hc << " " 
                // << meson_r*meson_r*Mnucleon/hc*Mnucleon/hc <<" " 
                // << meson_r*meson_r*meson_r*Mnucleon/hc*Mnucleon/hc*Mnucleon/hc << endl;
            }

//ok!
            double meson2_integrated=0.;
            double meson3_integrated=0.;

            for(int n=0; n<NL_points; n++){
                integral_parameters.n=n;
                integral_parameters.function= meson2;
                //integrate in r:
                meson2_integrated=integrate_tf(harmonic_integration_function, this );

                integral_parameters.function= meson3;
                meson3_integrated=integrate_tf(harmonic_integration_function, this );
            
                source_coeffs_nonlinear[n]= source_coeffs[n]-  g3*meson2_integrated/2. - g4*meson3_integrated/6.;                
                // std::cout << "int:" << n+1 << " " <<  source_coeffs_nonlinear[n]*pow(Mnucleon/hc, 3./2.) 
                // << " " << meson2_integrated*Mnucleon/hc << " " << meson3_integrated*Mnucleon*Mnucleon/hc/hc << endl;
            
            }

                //when source_coeffs_nonlinear is changed, so is Eigen vector coefficients
                //check:
                            //H does not change
                Eigen::VectorXd coefficients_nonlinear = Eigen::Map<Eigen::VectorXd>(source_coeffs_nonlinear.data(), source_coeffs_nonlinear.size());                            
                Eigen::VectorXd meson_coefs_eigen_updated= H_inverse*coefficients_nonlinear;
                std::vector<double> mesonv_coefs_updated(meson_coefs_eigen_updated.data(), 
                                                            meson_coefs_eigen_updated.data() + meson_coefs_eigen_updated.size());

                for(int n=0; n<NL_points; n++){
                    xi2+= pow(mesonv_coefs_updated[n]-mesonv_coefs[n], 2.);
                    mesonv_coefs[n]=mesonv_coefs_updated[n];
                }

    
        icounter++;

        // if(icounter>1e6) break;
        }while(xi2>del);
        // cout << "final?" << xi2 << " " << icounter << endl;
        
    }//end if g3 or g4 not zero
        std::vector<double> meson_new= get_function(mesonv_coefs);
    // cout << "mesonv_field[r] "<< endl;

    // for(int i=0; i<Nr_points; i++){
    //     cout << r[i]*hc/Mnucleon << " " << meson_new[i]*Mnucleon/hc << endl;
    // }
    // if(icounter>1e6){ cout << "no convergence in calculate_mesons" << endl; exit(1);};
    return mesonv_coefs;
}

double tf_walecka_class::integrate_tf(double (func)(double, void *), void *parametersPointer){

    tf_walecka_class &tf= *reinterpret_cast<tf_walecka_class*>(parametersPointer);
    double result = 0e0;
    double er_res = 0e0;
    double err_abs = 1e-8; //1e-13
    double err_rel = 1e-10; //1e-10;
    size_t max_steps = 1e8;
	
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(max_steps);

	gsl_function My_function;
    My_function.function = func;
    My_function.params = parametersPointer;

    // printf("lower limit: %e\n", tf.r.front());
    // printf("upper limit: %e\n", tf.r.back());
	gsl_integration_qag(&My_function, tf.r.front(), tf.r.back(), err_abs, err_rel, max_steps, 3,
                                                           w, &result, &er_res);


	gsl_integration_workspace_free(w);
    return result;
    
}



double harmonic_integration_function(double r_, void * params){
    tf_walecka_class &tf= *reinterpret_cast<tf_walecka_class*>(params);

    double weight_=1.;
    if(tf.idim==3)weight_=r_*r_;
    else if(tf.idim==2) weight_=r_;
    else if(tf.idim==1) weight_=2.;
    else{
        cout << "idim not defined in harmonic_integration_function()!" << endl;
        exit(1);
    }

    double function_integrated= weight_*interpolation_func(r_, tf.integral_parameters.function, tf.integral_parameters.r)
            *tf.Phi_Oscilator(tf.integral_parameters.n+1, r_, tf.b_oscilator);
        
    return function_integrated;

}



double _integration_function(double r_, void * params){
    tf_walecka_class &tf= *reinterpret_cast<tf_walecka_class*>(params);

    double weight_=1.;
    if(tf.integral_parameters.use_r_jacobian){
        if(tf.idim==3)      weight_=r_*r_;
        else if(tf.idim==2) weight_=r_;
        else if(tf.idim==1) weight_=2.;
        else{
            cout << "idim not defined in _integration_function()!" << endl;
            exit(1);
        }
    }
    double function_integrated= weight_*interpolation_func(r_, tf.integral_parameters.function, tf.integral_parameters.r);
        
    return function_integrated;

}


// double diff_N(double x, void *params) {
//     // Define your function here
//     tf_walecka_class &tf= *reinterpret_cast<tf_walecka_class*>(params);
//     double mun=x;

//     std::vector<double> rhonv_, kfv;
//     for(int i=0; i<tf.Nr_points; i++){
//         double masseff= tf.neutron.mass - tf.gs*tf.sigmav[i] - 2.*tf.neutron.I3*tf.gd*tf.deltav[i];
//         double muneff_r= mun - tf.gv*tf.omegav[i] - tf.neutron.I3*tf.gr*tf.rhov[i];

//         double kf2= muneff_r*muneff_r - masseff*masseff;

//         double kf_= kf2>0.? sqrt(kf2) : 0.;
//         std::cout << i << " " << kf_*Mnucleon/hc << endl;
//         kfv.push_back(kf_);
//         rhonv_.push_back(tf.neutron.gamma*kf_*kf_*kf_/6./M_PI/M_PI);
//     }

//     // double sumkf=0.;
//     // for(auto kf_ : kfv)
//     //     sumkf+=kf_;

//     double N_cell=0.;
//     // if(sumkf>0.){
//         tf.integral_parameters.function= rhonv_;
//         N_cell= tf.jacobian*tf.integrate_tf(_integration_function, &tf);
//     // }
//     std::cout << "N: " << N_cell << " " << tf.N <<  endl;

//     return tf.N - N_cell;
// }



// double zriddr(double (*func)(double, void *), double x1, double x2, double xacc, void *p, int MAXIT, double UNUSED) {
//     double fl = func(x1,p);
//     double fh = func(x2,p);

//     std::cout << "fs: " <<  fl*Mnucleon/hc << " " << fh*Mnucleon/hc << endl;
//     if ((fl > 0 && fh < 0) || (fl < 0 && fh > 0)) {
//         double xl = x1;
//         double xh = x2;
//         double zriddr_ = UNUSED;

//         for (int j = 0; j < MAXIT; ++j) {
//             double xm = 0.5 * (xl + xh);
//             double fm = func(xm, p);
//             double s = std::sqrt(fm * fm - fl * fh);

//             if (s == 0.0) return zriddr_;

//             double xnew = xm + (xm - xl) * ((fl - fh > 0) ? fm/s : -fm/s);

//             if (std::abs(xnew - zriddr_) <= xacc) return xnew;

//             zriddr_ = xnew;
//             double fnew = func(zriddr_,p);

//             if (fnew == 0.0) return zriddr_;

//             if ((fm * fnew > 0 && fm > 0) || (fm * fnew < 0 && fm < 0)) {
//                 xl = xm;
//                 fl = fm;
//                 xh = zriddr_;
//                 fh = fnew;
//             } else if ((fl * fnew > 0 && fl > 0) || (fl * fnew < 0 && fl < 0)) {
//                 xh = zriddr_;
//                 fh = fnew;
//             } else if ((fh * fnew > 0 && fh > 0) || (fh * fnew < 0 && fh < 0)) {
//                 xl = zriddr_;
//                 fl = fnew;
//             } else {
//                 std::cout << "Error: Shouldn't reach here in zriddr" << std::endl;
//                 return UNUSED;
//             }

//             if (std::abs(xh - xl) <= xacc) return xh;
//         }

//         std::cout << "Error: zriddr exceeded maximum iterations" << std::endl;
//         return UNUSED;
//     } else if (fl == 0.0) {
//         return x1;
//     } else if (fh == 0.0) {
//         return x2;
//     } else {
//         std::cout << "Error: Root must be bracketed in zriddr" << std::endl;
//         return UNUSED;
//     }
// }
