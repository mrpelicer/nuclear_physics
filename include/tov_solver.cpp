#include "tov_solver.hpp"



void tov_class::solve_tov_euler(double e_center, double dr){
        
        // massv.clear();
        // radiusv.clear();
        // pv.clear();
        // ev.clear();

        e_center*=MeV_fm3_to_pa_cgs/(pow(c_vel, 2.));
        
        double dm, dp;
        vector<double> presv_=p_eosv;
        vector<double> enerv_=e_eosv;

        // double p_center=interpolation_func(e_center, presv_, enerv_);
        double p_center=interpolate(e_center, presv_, enerv_);

        double mass_=0.;
        double radius_ =dr;
        double ener_=e_center;
        double pres_= p_center;
        

        set_stop_condition();

       do{     
        // cout << pres_ << " " << ener_ << " " << dp << endl;
            dm= dr*4.*M_PI*radius_*radius_*ener_;
            dp= -G_cte*dr*(ener_+pres_/pow(c_vel, 2.))*( mass_+4.*M_PI*pow(radius_, 3.)*pres_/pow(c_vel, 2.) )
                                /(radius_*(radius_ -2.*G_cte*mass_/pow(c_vel, 2.)) );


            
            massv.push_back(mass_);
            radiusv.push_back(radius_);
            pv.push_back(pres_);
            ev.push_back(ener_);

            mass_+=dm;
            radius_+=dr;
            pres_+=dp;
            if(pres_> presv_[0] && pres_>0.) ener_= interpolate(pres_, enerv_, presv_);
            // ener_= interpolation_func(pres_, enerv_, presv_);

    
            
            
           }while(pres_> presv_[0] && pres_>0. && ener_>0.);
//&& dm > dmrel
    mass= mass_;
    radius=radius_;
    
}



void tov_class::solve_tov_euler_pcentral(double p_center, double dr){

        massv.clear();
        radiusv.clear();
        pv.clear();
        ev.clear();
        p_center*=MeV_fm3_to_pa_cgs;
        double dm, dp;

        if(doqh_trans) set_eos_phase(p_center);
        
        vector<double> presv_=p_eosv;
        vector<double> enerv_=e_eosv;

        double e_center=interpolation_func(p_center, enerv_, presv_);

        double mass_=0.;
        double radius_=dr;
        double ener_=e_center;
        double pres_= p_center;
        
        set_stop_condition();

        do{

           
            dm= dr*4.*M_PI*radius_*radius_*ener_;
            dp= -G_cte*dr*(ener_+pres_/pow(c_vel, 2.))*( mass_+4.*M_PI*pow(radius_, 3.)*pres_/pow(c_vel, 2.) )
                                /(radius_*(radius_ -2.*G_cte*mass_/pow(c_vel, 2.)) );

            massv.push_back(mass_);
            radiusv.push_back(radius_);
            pv.push_back(pres_);
            ev.push_back(ener_);

            mass_+=dm;
            radius_+=dr;

            pres_+=dp;
            if(doqh_trans) set_eos_phase(pres_);
            presv_=p_eosv;
            enerv_=e_eosv;
            if(pres_> presv_[0] && pres_>0.) ener_= interpolation_func(pres_, enerv_, presv_);

        }while(pres_> p_stop && pres_>0.);
//&& dm>dmrel
    mass= mass_;
    radius=radius_;
        
}


void  tov_class::set_eos_phase(double pres_){
    if(pres_>p_trans){
        p_eosv=pq_eosv;
        e_eosv=eq_eosv;
    }else{
        p_eosv=ph_eosv;
        e_eosv=eh_eosv;
    }
    
}


double tov_class::getMass(){
    return mass/Msun;
}

double tov_class::getRadius(){
    return radius/1e5;// get in km
}


void tov_class::set_stop_condition(){
    if(docrust) p_stop= p_crust.front();
    else if(!doqh_trans) p_stop= p_eosv.front();
    else{p_stop= ph_eosv.front();}
}

void tov_class::do_crust(){

    /// BPS MODEL
    // cout << "This will not work with  pure quark stars!! Comment calling do_crust()" << endl;
    docrust=true;
    
    // in units of 1/fm4
    p_crust={1.2120000000000001e-11,
                8.2360000000000000e-11,
                2.7640000000000000e-10,
                5.1520000000000004e-10,
                1.5930000000000000e-09,
                4.0229999999999999e-09,
                1.3799999999999998e-08,
                3.3150000000000002e-08,
                1.0770000000000000e-07,
                2.5590000000000001e-07,
                3.4789999999999999e-07,
                4.7290000000000001e-07,
                6.4300000000000003e-07,
                9.1470000000000000e-07,
                1.0410000000000001e-06,
                1.8400000000000000e-06,
                2.4690000000000000e-06,
                2.4959999999999999e-06,
                2.6419999999999999e-06,
                2.8779999999999998e-06,
                3.1099999999999999e-06,
                3.4249999999999998e-06,
                3.8519999999999997e-06,
                4.4250000000000000e-06,
                5.1809999999999992e-06,
                6.1679999999999992e-06,
                8.1980000000000001e-06,
                1.1090000000000001e-05,
                1.5090000000000000e-05,
                2.0500000000000000e-05,
                2.7670000000000001e-05,
                3.7010000000000000e-05,
                5.3609999999999997e-05};

    e_crust={9.3870000000000006e-08,
                3.7380000000000003e-07,
                9.3920000000000005e-07,
                1.4890000000000001e-06,
                3.7409999999999998e-06,
                7.4650000000000005e-06,
                1.8770000000000002e-05,
                3.7469999999999999e-05,
                9.4179999999999996e-05,
                1.8809999999999999e-04,
                2.3690000000000001e-04,
                2.9819999999999998e-04,
                3.7579999999999997e-04,
                5.2419999999999995e-04,
                5.9579999999999995e-04,
                9.4519999999999999e-04,
                1.2220000000000000e-03,
                1.2680000000000000e-03,
                1.4859999999999999e-03,
                1.8790000000000000e-03,
                2.2640000000000000e-03,
                2.7650000000000001e-03,
                3.3999999999999998e-03,
                4.1820000000000000e-03,
                5.1310000000000001e-03,
                6.2599999999999991e-03,
                8.3289999999999996e-03,
                1.0900000000000000e-02,
                1.4019999999999999e-02,
                1.7760000000000001e-02,
                2.2180000000000002e-02,
                2.73200000000000010e-02,
                3.54200000000000000e-02};

 
    
    std::for_each(p_crust.begin(), p_crust.end(), [](double &el){el *= hc*MeV_fm3_to_pa_cgs; });
	std::for_each(e_crust.begin(), e_crust.end(), [](double &el){el *= hc*MeV_fm3_to_pa_cgs/(pow(c_vel, 2.)); });


    p_init_crust= p_crust.back();
    e_init_crust= e_crust.back();
    // cout << e_init_crust << " " << p_init_crust << endl;
    // cout << e_eosv.front() << " "  << p_eosv.front() << endl;
   
    if(doqh_trans){
        assert(p_init_crust< ph_eosv.front());
        assert(e_init_crust< eh_eosv.front());    
        ph_eosv.insert(ph_eosv.begin(), p_crust.begin(), p_crust.end()); 
        eh_eosv.insert(eh_eosv.begin(), e_crust.begin(), e_crust.end()); 

    }else{
        
        if( p_init_crust>= p_eosv.front() || e_init_crust>= e_eosv.front()){
            do{
                cout << "here?" <<endl;
                p_eosv.erase(p_eosv.begin());
                e_eosv.erase(e_eosv.begin());
                // p_init_crust= p_crust.back();
                // e_init_crust= e_crust.back();
            }while(p_init_crust>= p_eosv.front() || e_init_crust>= e_eosv.front());
            //  cout <<"finished?" << endl;
            //  cout << e_init_crust << " " << p_init_crust << endl;
                // cout << e_eosv.front() << " "  << p_eosv.front() << endl;
        }

        p_eosv.insert(p_eosv.begin(), p_crust.begin(), p_crust.end()); 
        e_eosv.insert(e_eosv.begin(), e_crust.begin(), e_crust.end()); 
    }

       
    // vector<double>edens_Crust={1.9900000000000000e-08,
    //                         7.9240000000000000e-08,
    //                         1.9900000000000000e-07,
    //                         3.1549999999999999e-07,
    //                         7.9240000000000003e-07,
    //                         1.5810000000000000e-06,
    //                         3.9720000000000003e-06,
    //                         7.9240000000000007e-06,
    //                         1.9899999999999999e-05,
    //                         3.9719999999999999e-05,
    //                         5.0000000000000002e-05,
    //                         6.2940000000000004e-05,
    //                         7.9239999999999993e-05,
    //                         1.1050000000000000e-04,
    //                         1.2559999999999999e-04,
    //                         1.9900000000000001e-04,
    //                         2.5720000000000002e-04,
    //                         2.6699999999999998e-04,
    //                         3.1260000000000001e-04,
    //                         3.9510000000000001e-04,
    //                         4.7590000000000002e-04,
    //                         5.8120000000000003e-04,
    //                         7.1429999999999996e-04,
    //                         8.7860000000000000e-04,
    //                         1.0770000000000000e-03,
    //                         1.3140000000000001e-03,
    //                         1.7480000000000000e-03,
    //                         2.2870000000000000e-03,
    //                         2.9420000000000002e-03,
    //                         3.7260000000000001e-03,
    //                         4.6499999999999996e-03,
    //                         5.7279999999999996e-03,
    //                         7.4240000000000000e-03};
}