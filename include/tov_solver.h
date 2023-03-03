#ifndef tov_h
#define tov_h

#include "constant.h"
#include "interpolator.h"
#include <vector>

using namespace std;

class tov_class{
    public:

        tov_class(vector<double> e_eosv_, vector<double> p_eosv_){
            e_eosv = e_eosv_;
            p_eosv = p_eosv_;
            std::for_each(p_eosv.begin(), p_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs; });
	        std::for_each(e_eosv.begin(), e_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs/(pow(c_vel, 2.)); });

        };

        tov_class(vector<double> eh_eosv_, vector<double> ph_eosv_, vector<double> eq_eosv_, vector<double> pq_eosv_, 
                    double p_trans_, double eh_trans_, double eq_trans_){
            eh_eosv = eh_eosv_;
            ph_eosv = ph_eosv_;
            eq_eosv = eq_eosv_;
            pq_eosv = pq_eosv_;

            std::for_each(ph_eosv.begin(), ph_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs; });
            std::for_each(eh_eosv.begin(), eh_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs/(pow(c_vel, 2.)); });
            std::for_each(pq_eosv.begin(), pq_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs; });
            std::for_each(eq_eosv.begin(), eq_eosv.end(), [](double &el){el *= MeV_fm3_to_pa_cgs/(pow(c_vel, 2.)); });

            p_trans= p_trans_*MeV_fm3_to_pa_cgs;
            eh_trans= eh_trans_*MeV_fm3_to_pa_cgs/(pow(c_vel, 2.));
            eq_trans= eq_trans_*MeV_fm3_to_pa_cgs/(pow(c_vel, 2.));

            doqh_trans=true;
        };

        void solve_tov_euler(double e_center, double dr);
        void solve_tov_euler_pcentral(double pe_center, double dr);

        void solve_tidal_euler(double e_center, double dr);

        double getMass();
        double getRadius();
        double gety();
        double getCompactness();
        double getBeta();
        double getH();
        double getk2();
        double getLambda();

        void do_crust();
        vector<double> e_eosv, p_eosv;
    private:
        //EoS definitions
        
        vector<double> eh_eosv, ph_eosv, eq_eosv, pq_eosv;
        //Phase transition:
        double p_trans;
        double eh_trans, eq_trans;
        void  set_eos_phase(double pres_);
        bool doqh_trans=false;

        //Crust variables:
        double p_init_crust, e_init_crust;
        vector<double> e_crust, p_crust;
        int icrust=0;
        bool docrust=false;
        double p_stop;
        void set_stop_condition();

        //Star:
        double mass=0., radius=0., H=0., beta=0.;
        double y=0., k2=0., compactness=0.;
        double dmrel= 1e-12;
        double e_center, p_center;
        vector<double> massv, radiusv, pv, ev;
        double dedp(double r_);
};

#endif