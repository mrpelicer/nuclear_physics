#include "rmf_non_linear_walecka.h"

nlwm_class::nlwm_class(std::string parametrization_){
  setParametrization(parametrization_);
	proton.spin		= 1./2.;	proton.I3			= 1./2.;	proton.Q			= 1.;
	neutron.spin	= 1./2.;	neutron.I3		=-1./2.;	neutron.Q			= 0.;
}

//=============== Set RMF parameters: nucleons/meson couplings and masses  ===============
void nlwm_class::setParametrization(std::string parametrization_){
		parametrization=parametrization_;

	//===== GM1 =========
	if(parametrization_=="gm1")		//Phys. Rev. Lett. 67, 2414 (1991)
	{
		Mn=938.99; 										//	Mn=938.930;
		Ms=400./Mnucleon; 									//	Ms=400./Mnucleon;
		Mv=783./Mnucleon; 									//	Mv=783./Mnucleon;
		Mb=770./Mnucleon; 									//	Mb=770./Mnucleon;
		gs=Ms*sqrt(11.785)*Mnucleon/hc; 		//	gs=Ms*sqrt(11.785)*Mnucleon/hc;
		gv=Mv*sqrt(7.148)*Mnucleon/hc; 	 		//	gv=Mv*sqrt(7.148)*Mnucleon/hc;
		gb=Mb*sqrt(3.870)*Mnucleon/hc; 			//	gb=Mb*sqrt(3.870)*Mnucleon/hc;
		gs3=2.*pow(gs, 3.)*0.002947;
		gs4=-6.*pow(gs, 4.)*0.001071;
		rho0=0.153*pow(hc/Mnucleon, 3); 
	}
	//===== GM1wr =========
	else if(parametrization_=="gm1wr")		//Phys. Rev. Lett. 67, 2414 (1991)
	{
		Mn=938.930;
		Ms=400./Mnucleon;
		Mv=783./Mnucleon;
		Mb=770./Mnucleon;
		gs=Ms*sqrt(11.79)*Mnucleon/hc;
		gv=Mv*sqrt(7.149)*Mnucleon/hc;
		gb=Mb*sqrt(4.411)*Mnucleon/hc;
		gs3=2.*pow(gs, 3.)*0.002947;
		gs4=-6.*pow(gs, 4.)*0.001070;
		Lv=0.02015;
		rho0=0.153*pow(hc/Mnucleon, 3); 
	}
	//===== NL1 =====
	else if(parametrization_=="nl1") 	//see nl3 ref: 	Phys. Rev. C 55, 540 (1997)
	{
		Mn=938.000;
		Mstar=0.57;
		Ms=492.250/Mnucleon;
		Mv=783.000/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.138;
		gv=13.285;
		gb=2.*4.976;
		gs3=2.*(12.172)*hc/Mnucleon;
		gs4=-6.*(36.265);
		rho0=0.153*pow(hc/Mnucleon, 3);
	}
	// ===== NL3 =====
	else if(parametrization_=="nl3") 		//Phys. Rev. C 55, 540 (1997)
	{
		Mn=939.000;
		Mstar=0.60;
		Ms=508.194/Mnucleon;
		Mv=782.501/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.217;
		gv=12.868;
		gb=8.948;      
		gs3=4.384; 
		gs4=-173.31;
		rho0=0.148*pow(hc/Mnucleon, 3);
	}
	else if(parametrization_=="nl3-src") 		//Private communication O.Lourenço
	{
		Mn=939.000;
		Mstar=0.60;
		Ms=508.194/Mnucleon;
		Mv=782.501/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.495965168181097  ;
		gv=12.013019178168182 ;
		gb= 14.848527065665982;      
		gs3=2.*2.9964708544616308; 
		gs4= -6.*45.697449558726071;
		Lv=(1.9291481766492405E-004)/2.;
		rho0=0.148*pow(hc/Mnucleon, 3);		
	}
	// ==== NL3wr ====
	else if(parametrization_=="nl3wr")		//Phys. Rev. C 55, 540 (1997)
	{
		Mn=939.000;
		Mstar=0.60;
		Ms=508.194/Mnucleon;
		Mv=782.501/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.217;
		gv=12.868;
		gb=11.276624779336951;
		gs3=4.384; //1/MeV
		gs4=-173.31;
		Lv=0.03;
		rho0=0.1481*pow(hc/Mnucleon, 3); 
	}
	// ===FSUGold ====
	else if(parametrization_=="fsu") 		//Phys. Rev. Lett. 95, 122501 (2005)
	{
		Mn=939.000;
		Ms=491.5/Mnucleon;
		Mv=782.5/Mnucleon;
		Mb=763./Mnucleon;
		gs=10.592;
		gv=14.302;
		gb=11.767;
		gs3=1.7976;
		gs4=299.11;
		xsi=0.06;
		Lv=0.03;
		rho0=0.1484*pow(hc/Mnucleon, 3); 
	}

	// ===IU-FSU ====
	else if(parametrization_=="iufsu")			//Phys. Rev. C 82, 055803 (2010)
	{
		Mn=939.000;
		Ms=491.5/Mnucleon;
		Mv=782.5/Mnucleon;
		Mb=763./Mnucleon;
		gs=9.971;
		gv=13.032;
		gb=13.590;
		gs3= 3.5695;
		gs4=2.926;
		xsi=0.03;
		Lv=0.046;
		rho0=0.155*pow(hc/Mnucleon, 3); 
	}else if(parametrization_=="iufsu-src")			//Private communication
	{
		Mn=939.000;
		Ms=491.5/Mnucleon;
		Mv=782.5/Mnucleon;
		Mb=763./Mnucleon;
		gs=10.132389686488604;
		gv=11.867320024781414;
		gb=15.551063183598112;
		gs3= 2.*2.9556603254255340;
		gs4=-6.*29.880147636307864;
		xsi=0.03;
		Lv=(1.0941815626967822e-2)/2.;
		rho0=0.155*pow(hc/Mnucleon, 3); 
	}
	//===== FSU2R =========
	else if(parametrization_=="fsu2")		//Phys. Rev. C 90, 044305 (2014)
	{
		Mn=939.000;
		Mstar=0.593;
		Ms=497.479/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=sqrt(108.0943);
		gv=sqrt(183.7893);
		gb=sqrt(80.4656);      
		gs3=3.0029*pow(gs, 3.)/Mnucleon; //(1/Mn)
		gs4=-0.000533*pow(gs, 4.);
		xsi=0.0256;
		Lv=0.000823;
		rho0=0.1505*pow(hc/Mnucleon, 3);
	}

	//===== FSU2R =========
	else if(parametrization_=="fsu2r")		//Astron.Soc.Austral. 34 e065 
	{
		Mn=939.000;
		Mstar=0.593;
		Ms=497.479/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=sqrt(107.5751);
		gv=sqrt(182.3949);
		gb=sqrt(206.4260);      
		gs3=3.0911*pow(gs, 3.)/Mnucleon; //(1/Mnucleon)
		gs4=-0.001680*pow(gs, 4.);
		xsi=0.024;
		Lv=0.045;
		rho0=0.150*pow(hc/Mnucleon, 3);
	}
	else if(parametrization_=="fsu2r-src")		//PhysRevD.105.023008
	{
		Mn=939.000;
		Mstar=0.593;
		Ms=497.479/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.5174;
		gv=12.3648;
		gb=15.5988;      
		gs3=2.*2.9133; //(1/Mnucleon)
		gs4=-6.*32.4432;
		xsi=0.024;
		Lv=0.0093/2.;
		rho0=0.1505*pow(hc/Mnucleon, 3);
	}
	//===== FSU2h =========
	else if(parametrization_=="fsu2h")		//Astron.Soc.Austral. 34 e065 
	{
		Mn=939.;
		// Mstar=0.593;
		Ms=497.479/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		Mt=1020.00/Mnucleon;
		gs=sqrt(102.7200);
		gv=sqrt(169.5315);
		gb=sqrt(197.2692);      
		gs3=4.0014*pow(gs, 3.)/Mnucleon; //(1/Mnucleon)
		gs4=-0.013298*pow(gs, 4.);
		xsi=0.008;
		Lv=0.045;
		rho0=0.1505*pow(hc/Mnucleon, 3);
	}

	else if(parametrization_=="l3wr"){
  	Mn=939.;        
  	Ms=512./Mnucleon;
  	Mv=783./Mnucleon;
  	Mb=770./Mnucleon;
		Mt=1020.00/Mnucleon;
  	gs=sqrt(12.108)*Ms*Mnucleon/hc;
  	gv=sqrt(7.1320)*Mv*Mnucleon/hc;
  	gb=sqrt(4.8010)*Mb*Mnucleon/hc;
  	gs3=2.*0.004138*pow(gs, 3);
  	gs4=-6.*0.0039*pow(gs, 4.);
  	xsi=0.;
  	Lv=0.0185;
		rho0=0.1555*pow(hc/Mnucleon, 3);
	}

	else if(parametrization_=="el3wr"){
  	Mn=939.;        
  	Ms=512./Mnucleon;
  	Mv=783./Mnucleon;
  	Mb=770./Mnucleon;
		Mt=1020.00/Mnucleon;
  	gs=sqrt(12.108)*Ms*Mnucleon/hc;
  	gv=sqrt(7.1320)*Mv*Mnucleon/hc;
  	gb=sqrt(5.85)*Mb*Mnucleon/hc;
  	gs3=2.*0.004138*pow(gs, 3);
  	gs4=-6.*0.0039*pow(gs, 4.);
  	xsi=0.;
  	Lv=0.0283;
		rho0=0.156*pow(hc/Mnucleon, 3);
	}

	else if(parametrization=="iobp-I"){
		Mn=939.;        
		Ms=0.533;
		Mv=0.833;
		Mb=0.812;
		gs=10.417594;
		gv=13.354502;
		gb=11.121523;
		gs3=4.3545284;
		gs4=-88.51362;
		xsi=0.017014393;
		Lv=0.0146;
		rho0=0.149*pow(hc/Mnucleon, 3);
		
	}

	else if(parametrization=="nl3*"){ //https://arxiv.org/pdf/0909.1432.pdf : wrong gs3 in paper
		Mn=939.000;
		Ms=502.574/Mnucleon;
		Mv=782.600/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.0944;
		gv=12.8065 ;
		gb=2.*4.441;
		gs3=2.*10.8093*hc/Mnucleon; //1/MeV
		gs4=-6.*30.1486;
		Lv=0.0;
		rho0=0.150*pow(hc/Mnucleon, 3); 
	}

	else if(parametrization=="nl3wr*"){ //https://arxiv.org/pdf/2111.02247.pdf : wrong gs3 in paper
		Mn=939.000;
		Ms=502.574/Mnucleon;
		Mv=782.600/Mnucleon;
		Mb=763.000/Mnucleon;
		Mt=1020.00/Mnucleon;
		gs=10.0944;
		gv=12.8065 ;
		gb=14.441;
		gs3=2.*10.8093*hc/Mnucleon; //1/MeV
		gs4=-6.*30.1486;
		Lv=0.045;
		rho0=0.150*pow(hc/Mnucleon, 3); 
	}
	else if(parametrization=="ddme1"){ //Phys. Rev. C 66,024306 (2002)
		Mn=939.000;
		Ms=549.5255/Mnucleon;
		Mv=783.0000/Mnucleon;
		Mb=763.0000/Mnucleon;
		Mt=1020.00/Mnucleon;
		gs=10.4434;
		gv=12.8939;
		gb=2.*3.8053;

		as=1.3854;
		bs=0.9781;
		cs=1.5342;
		ds=0.4661;

		av=1.3879;
		bv=0.8525;
		cv=1.3566;
		dv=0.4957;

		ab=0.5008;

 		rho0=0.152*pow(hc/Mnucleon, 3); 
		dogRhob=true;
	}
	
	else if(parametrization=="ddme2"){ //Phys. Rev. C 71, 024312 (2005)
		Mn=939.000;
		Ms=550.1238/Mnucleon;
		Mv=783.0000/Mnucleon;
		Mb=763.0000/Mnucleon;
		Mt=1020.00/Mnucleon;
		gs=10.5396;
		gv=13.0189;
		gb=7.3672;
	
		as=1.3881;
		bs=1.0943;
		cs=1.7057;
		ds=0.4421;
		
		av=1.3892;
		bv=0.9240;
		cv=1.4620;
		dv=0.4775;

		ab=0.5647;

		rho0=0.152*pow(hc/Mnucleon, 3); 
		dogRhob=true;
	}else if(parametrization=="ffg"){ //PhysRevC.93.014619
		Mn=939.000;
		Ms=500.000/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		gs=10.9310;
		gv=14.5947;
		gb=5.9163;
		gs3=2.*0.0007473*pow(gs, 3.); //1/MeV
		gs4=6.*0.003882*pow(gs, 4.);
		Lv=0.2736/2.;
		xsi= 6.*0.01;
		rho0=0.150*pow(hc/Mnucleon, 3); 
	}else if(parametrization=="hmt"){ //PhysRevC.93.014619
		Mn=939.000;
		Ms=500.000/Mnucleon;
		Mv=782.500/Mnucleon;
		Mb=763.000/Mnucleon;
		gs= 10.8626;
		gv= 12.9185 ;
		gb= 7.8712;
		gs3=2.*0.0007473*pow(gs, 3.); //1/MeV
		gs4=6.*0.0005139*pow(gs, 4.);
		Lv=0.03740/2.;
		xsi= 6.*0.01;
		rho0=0.150*pow(hc/Mnucleon, 3); 
	}
  else{parametrization_= "";}
	std::cout << parametrization << " parametrization." << std::endl;
}


//=============== Set RMF parameters: hyperons couplings and masses  ===============
void nlwm_class::includeHyperons(bool do_, std::string parameters_){
	doHyperons=do_;
	std::cout << "do Hyperons! " << doHyperons << std::endl;

	// lambda0.mass=1116./Mnucleon;
	// sigmam.mass =1193./Mnucleon;
	// sigma0.mass =1193./Mnucleon;
	// sigmap.mass =1193./Mnucleon;
	// xim.mass		=1318./Mnucleon;
	// xi0.mass		=1318./Mnucleon;

	lambda0.mass=1115.7/Mnucleon;
	sigmam.mass =1197.5/Mnucleon;
	sigma0.mass =1192.6/Mnucleon;
	sigmap.mass =1189.4/Mnucleon;
	xim.mass		=1321.7/Mnucleon;
	xi0.mass		=1314.9/Mnucleon;

	lambda0.spin	= 1./2.;			 	lambda0.I3	= 0.;			lambda0.Q		= 0;  lambda0.stg =-1;
	sigmap.spin		=	1./2.;			 	sigmap.I3		=	1.;			sigmap.Q		= 1.;	sigmap.stg	=-1;
	sigma0.spin		=	1./2.;			 	sigma0.I3		=	0.;			sigma0.Q		= 0; 	sigma0.stg	=-1;
	sigmam.spin		= 1./2.;				sigmam.I3		=-1.;			sigmam.Q		=-1.;	sigmam.stg	=-1;
	xi0.spin			=	1./2.;				xi0.I3			=	1./2.;	xi0.Q				= 0; 	xi0.stg			=-2;
	xim.spin			= 1./2.;				xim.I3			=-1./2.;	xim.Q				=-1.;	xim.stg			=-2;

	parhyp= parameters_;

	if(parhyp=="gm"){ //Glendenning and Moszkowski PRL 67, 18
		if(parametrization=="gm1" || parametrization=="gm1wr" || parametrization=="gm3"  ){
			double xsH=0.7;
			double xvH=0.783;
			double xbH=0.783;
			xsl=xsH; xss=xsH; xsx=xsH;
			xvl=xvH; xvs=xvH; xvx=xvH;
			xbl=xbH; xbs=xbH; xbx=xbH;
		}else if(parametrization=="nl3" || parametrization=="nl3wr"){//2106.09515
			xsl=0.613; xss=0.460; xsx=0.317;
			xvl=0.667; xvs=0.667; xvx=0.333;
			xbl=0.; xbs=1.; xbx=1.;
		}else if(parametrization=="fsu" || parametrization=="iufsu"){//2106.09515
			xsl=0.611; xss=0.454; xsx=0.316;
			xvl=0.667; xvs=0.667; xvx=0.333;
			xbl=0.; xbs=1.; xbx=1.;
		}else{
			std::cout << "Unspecified hyperon parametrization." << std::endl;
			exit(1);
		}
	}
	else if(parhyp=="su3"){ //código do Kauan --
		xsl= 0.610; xss= 0.3957; xsx= 0.113;
		//xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xvl=.667; xvs=.667; xvx=.332;
		xbl=0.; xbs=1.; xbx=1.; //*** CHECK THIS xbs=2, WHICH MAYBE WAS ADDED IF ISOSPIN WAS NOT ACCOUNTED FOR!!!!
	}	
	else if(parhyp=="ddme2-a"){ //Fortin et al. Phys. Rev. C 95 065803
		xsl= 0.621; xss= 0.467; xsx= 0.321;
		xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xbl=0.; xbs=1.; xbx=1.;
		xtl=-sqrt(2.)/3.; xts=-sqrt(2.)/3.; xtx=-2.*sqrt(2.)/3.;
	}
	else if(parhyp=="ddme2-b"){ //Fortin et al. Phys. Rev. C 95 065803
		xsl= 0.896; xss= 0.467; xsx= 0.321;
		xvl=1.; xvs=2./3.; xvx=1./3.;
		xbl=0.; xbs=1.; xbx=1.;
		xtl=-sqrt(2.)/3.; xts=-sqrt(2.)/3.; xtx=-2.*sqrt(2.)/3.;
	}

	else if(parhyp=="fsu2h"){//Tolos et al. 1708.08681 -- canonical set
		xsl= 0.611; xss= 0.467; xsx= 0.316;
		xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xbl=0.; xbs=1.; xbx=1.;
		xtl=-sqrt(2.)/3.; xts=-sqrt(2.)/3.; xtx=-2.*sqrt(2.)/3.;
	}
	else if(parhyp=="fsu2h_1"){//Tolos et al. 1708.08681 -- most repulsive set
		xsl= 0.611; xss= 0.467; xsx= 0.271;
		xvl=.667; xvs=.667; xvx=.332;
		xbl=0.; xbs=1.; xbx=1.;
		xtl=-sqrt(2.)/3.; xts=-sqrt(2.)/3.; xtx=-2.*sqrt(2.)/3.;
	}
	else if(parhyp=="fsu2h_2"){//Tolos et al. 1708.08681 -- most attractive set
		xsl= 0.611; xss= 0.541; xsx= 0.316;
		xvl=.667; xvs=.667; xvx=.332;
		xbl=0.; xbs=1.; xbx=1.;
		xtl=-sqrt(2.)/3.; xts=-sqrt(2.)/3.; xtx=-2.*sqrt(2.)/3.;
	}
	else if(parhyp=="l3wr1"){
		double av_=1.;
		double as_=1.582;
	
		xsl= (10.+6.*as_)/(13.+12.*as_); 		 xss= (22.-6.*as_)/(13.+12.*as_); 			xsx=(13.-6.*as_)/(13.+12.*as_);
		//xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xvl=(4.+2.*av_)/(5+4.*av_);					 xvs=(8.-2.*av_)/(5+4.*av_); 						xvx=(5.-2.*av_)/(5+4.*av_);
		xbl=0.;											 			 xbs=2.*av_; 													xbx=-(1.-2.*av_);
		xtl=sqrt(2.)*(2.*av_-5.)/(5+4.*av_); xts=-sqrt(2.)*(2.*av_+1.)/(5+4.*av_); 	xtx=-sqrt(2.)*(2.*av_+4.)/(5+4.*av_);
	}
	else if(parhyp=="l3wr2"){
		double av_=.75;
		double as_=1.240;
		xsl= (10.+6.*as_)/(13.+12.*as_); 		 xss= (22.-6.*as_)/(13.+12.*as_); 			xsx=(13.-6.*as_)/(13.+12.*as_);
		//xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xvl=(4.+2.*av_)/(5+4.*av_);					 xvs=(8.-2.*av_)/(5+4.*av_); 						xvx=(5.-2.*av_)/(5+4.*av_);
		xbl=0.;											 			 xbs=2.*av_; 													xbx=-(1.-2.*av_);
		xtl=sqrt(2.)*(2.*av_-5.)/(5+4.*av_); xts=-sqrt(2.)*(2.*av_+1.)/(5+4.*av_); 	xtx=-sqrt(2.)*(2.*av_+4.)/(5+4.*av_);
	}
	else if(parhyp=="l3wr3"){
		double av_=.5;
		double as_=.911;
		xsl= (10.+6.*as_)/(13.+12.*as_); 		 xss= (22.-6.*as_)/(13.+12.*as_); 			xsx=(13.-6.*as_)/(13.+12.*as_);
		//xvl=2./3.; xvs=2./3.; xvx=1./3.;
		xvl=(4.+2.*av_)/(5+4.*av_);					 xvs=(8.-2.*av_)/(5+4.*av_); 						xvx=(5.-2.*av_)/(5+4.*av_);
		xbl=0.;											 			 xbs=2.*av_; 													xbx=-(1.-2.*av_);
		xtl=sqrt(2.)*(2.*av_-5.)/(5+4.*av_); xts=-sqrt(2.)*(2.*av_+1.)/(5+4.*av_); 	xtx=-sqrt(2.)*(2.*av_+4.)/(5+4.*av_);
	}
	else if(parametrization=="iufsu" && parhyp=="fw"){
		xsl=0.7098; xss=0.552; xsx=0.522;
		xvl=0.79; xvs=0.79; xvx=0.59;
		xbl=0.; xbs=1.; xbx=1.;
	}else if(parhyp=="iufsu_str"){// must fix!!! wrong potentials
		xsl=0.590;		xss=0.429; 		xsx=0.306;
		xvl=0.667;		xvs=0.667; 		xvx=0.333;
		xbl=0.; 			xbs=1.; 			xbx=1.;
		xtl=-0.471*sqrt(169.8349); 	xts=-0.471*sqrt(169.8349);		xtx=-0.943*sqrt(169.8349);


	}else if(parhyp=="nl3wr*"){ /* see tab II of https://arxiv.org/pdf/2111.02247.pdf
															too determine av_ and as_ */
		double av_=0.5;
		xsl= 0.651; 		 xss= 0.730; 			xsx=0.473; // 0.428 in the paper! wrong potential here
		xvl=(4.+2.*av_)/(5+4.*av_);					 xvs=(8.-2.*av_)/(5+4.*av_); 						xvx=(5.-2.*av_)/(5+4.*av_);
		xbl=0.;											 			 xbs=2.*av_; 													xbx=-(1.-2.*av_);
		xtl=sqrt(2.)*(2.*av_-5.)/(5+4.*av_); xts=-sqrt(2.)*(2.*av_+1.)/(5+4.*av_); 	xtx=-sqrt(2.)*(2.*av_+4.)/(5+4.*av_);
	}

	
	else{
		std::cout << "Unspecified hyperon parametrization." << std::endl;
		exit(1);
	}

}


//=============== Set RMF parameters: delta isobar couplings and masses  ===============
void nlwm_class::includeDeltas(bool do_, std::string parameters_){
	doDeltas=do_;
	std::cout << "do Deltas! " << doDeltas << std::endl;

	double mdl=1232./Mnucleon;
	deltapp.mass=mdl;
	deltap.mass	=mdl;
	delta0.mass	=mdl;
	deltam.mass	=mdl;

	deltapp.spin	= 3./2.;	deltapp.I3	= 3./2.;		deltapp.Q		= 2.;  deltapp.gamma=4.;
	deltap.spin		= 3./2.;	deltap.I3		=	1./2.;		deltap.Q		= 1.;  deltap.gamma =4.;
	delta0.spin		= 3./2.;	delta0.I3		=-1./2.;		delta0.Q		= 0.;  delta0.gamma =4.;
	deltam.spin		= 3./2.;	deltam.I3		=-3./2.;		deltam.Q		=-1.;  deltam.gamma =4.;
	


	pardelta= parameters_;
	if(pardelta=="su6"){ //SU(6) symmetry
		xsd=1.2; xvd=1.2; xbd=1.;
	}
	if(pardelta=="prd89_1"){ //PRD89, 043014
		xsd=1.25; xvd=1.; xbd=1.;
	}
	if(pardelta=="prd89_2"){ //PRD89, 043014
		xsd=1.15; xvd=0.9; xbd=1.;
	}
	if(pardelta=="mplA_1"){ //Modern Physics Letters A, Vol. 15, No. 24 (2000) 1529–1537
		xsd=1.1; xvd=1.; xbd=1.;
	}
	if(pardelta	=="mplA_2"){ //Modern Physics Letters A, Vol. 15, No. 24 (2000) 1529–1537
		xsd=1.2; xvd=1.; xbd=1.;
	}
}


//Print parameters of chosen set: 
void nlwm_class::printParameters(void){

  if(parametrization== ""){std::cout << "You have not chosen a valid parameter set!" 
																		 << std::endl;}
  else{
    std::cout << "Masses (Mn, Ms, Mv, Mb, Mstar): "  << endl
							<< Mn << " " << Ms << " " << Mv << " " << Mb << " " << Mstar <<
    std::endl << "Couplings (gs, gv, gb, gs3, gs4, xsi, Lv): "  << endl
							<< gs << " " << gv << " " << gb << " " 
              << gs3 << " " << gs4 << " " << xsi << " " << Lv 
		<< std::endl;

	if(doHyperons){
	    std::cout << "Hyperons: " << endl
							<< xsl << " " << xss << " " << xsx << std::endl
							<< xvl << " " << xvs << " " << xvx << std::endl
							<< xbl << " " << xbs << " " << xbx << std::endl
							<< xtl << " " << xts << " " << xtx << std::endl;
	}
	if(doDeltas){
	    std::cout << "Deltas: " << endl
							<< xsd << " " << xvd << " " << xbd << std::endl;
	}

  }
}


//=============== Calculate thermodynamic properties for all baryons: ===============
void nlwm_class::setThermodynamics(){
	proton.calculateProperties();
	neutron.calculateProperties();
	if(doHyperons){
		lambda0.calculateProperties();
		sigmap.calculateProperties();
		sigma0.calculateProperties();
		sigmam.calculateProperties();
		xi0.calculateProperties();
		xim.calculateProperties();
	}
	if(doDeltas){
		deltapp.calculateProperties();
		deltap.calculateProperties();
		delta0.calculateProperties();
		deltam.calculateProperties();
	}
}


//=============== Set temperature for all baryons: ===============
void nlwm_class::setTemperature(double temp_){
  temperature=temp_;
  proton.temperature 	=temperature; 
	neutron.temperature	=temperature;
	if(doHyperons){
		lambda0.temperature	=temperature;
		sigmap.temperature	=temperature;
		sigma0.temperature	=temperature;
		sigmam.temperature	=temperature;
		xi0.temperature			=temperature;
		xim.temperature			=temperature;
	}
	if(doDeltas){
		deltapp.temperature	=temperature;
		deltap.temperature	=temperature;
		delta0.temperature	=temperature;
		deltam.temperature	=temperature;
	}

}


//=============== Set the magnetic field for all baryons: ===============
void nlwm_class::setBfield(bool dob_, double B_){
	Bfield=B_;
		proton.setBfield(dob_, B_);
		neutron.setBfield(dob_, B_);
	if(doHyperons){
		lambda0.setBfield(dob_, B_);
		sigmap.setBfield(dob_, B_);
		sigma0.setBfield(dob_, B_);
		sigmam.setBfield(dob_, B_);
		xi0.setBfield(dob_, B_);
		xim.setBfield(dob_, B_);
	}
	if(doDeltas){
		deltapp.setBfield(dob_, B_);
		deltap.setBfield(dob_, B_);
		delta0.setBfield(dob_, B_);
		deltam.setBfield(dob_, B_);
	}
}


//=============== Set the magnetic moment for all baryons: ===============          
void nlwm_class::setAMM(bool doa_){							//MUB			Kb
		proton.setAMM(doa_, 1.79);									// 2.79     1.79 
		neutron.setAMM(doa_, -1.91);								//-1.91		 -1.91
		// proton.setAMM(doa_, 0.);									// 2.79     1.79 
		// neutron.setAMM(doa_,0.);								//-1.91		 -1.91

	if(doHyperons){ //PhysRevC.79.025803
		lambda0.setAMM(doa_, -0.61);								//-0.61    -0.61
		sigmap.setAMM(doa_, 1.67);									// 2.46  		1.67
		sigma0.setAMM(doa_, 1.61);									// 1.61			1.61
		sigmam.setAMM(doa_, -0.37);									//-1.16		 -0.37
		xi0.setAMM(doa_, -1.25);										//-1.25		 -1.25
		xim.setAMM(doa_, 0.06);											//-0.65			0.06
		// lambda0.setAMM(doa_,  	0.);
		// sigmap.setAMM(doa_, 		0.);
		// sigma0.setAMM(doa_, 		0.);
		// sigmam.setAMM(doa_,  	0.);
		// xi0.setAMM(doa_,  			0.);
		// xim.setAMM(doa_, 			0.);
	}
	if(doDeltas){//https://arxiv.org/pdf/hep-lat/0302008v2.pdf
		deltapp.setAMM(doa_, 3.47);										//4.99			3.47
		deltap.setAMM(doa_,  1.73);										//2.49  		1.73
		delta0.setAMM(doa_,  0.06);										//0.06			0.06
		deltam.setAMM(doa_, -1.69);										//-2.45		 -1.69
		// deltapp.setAMM(doa_,0.);
		// deltap.setAMM( doa_,0.);
		// delta0.setAMM( doa_,0.);
		// deltam.setAMM( doa_,0.);
	}	
}


//=============== Set nucleon EoS with fixed proton fraction and temperature  ===============
void nlwm_class::setEOS_nucleons(double rhoB_, double Yp_, double temp_){
  rhoB=rhoB_;
  Yp=Yp_;  
	setTemperature(temp_);

  proton.density=Yp*rhoB;
  proton.kf=pow(3.*pi2*proton.density, 1./3.);
  proton.kf2=pow(proton.kf, 2.);

  neutron.density=(1.-Yp)*rhoB;
  neutron.kf=pow(3.*pi2*neutron.density, 1./3.);
  neutron.kf2=pow(neutron.kf, 2.);

  rho3= proton.I3*proton.density+neutron.I3*neutron.density;
	
	setVectorMeanFields();
	setScalarMeanFields();			
	
  proton.solveChemPotEff();
  neutron.solveChemPotEff();

	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) 	: gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) 		: gb;

  proton.chemPot  =  proton.chemPot_eff  + gv_*V0 + gb_*b0*proton.I3;
  neutron.chemPot =  neutron.chemPot_eff + gv_*V0 + gb_*b0*neutron.I3;
	
	if(dogRhob){
		proton.chemPot += getRearrangementEnergy();
		neutron.chemPot += getRearrangementEnergy();
	}
	muB = neutron.chemPot;
	muQ = proton.chemPot - neutron.chemPot;

	proton.calculateCondensate();
	neutron.calculateCondensate();
	rhoS= proton.condensate + neutron.condensate;
	proton.calculateProperties();
  neutron.calculateProperties();

}


//=============== Set nucleon EoS with short range correlations 
//								with fixed proton fraction and temperature  ===============
void nlwm_class::setEOS_src_nucleons(double rhoB_, double Yp_, double temp_){
  rhoB=rhoB_;
  Yp=Yp_;  
	setTemperature(temp_);

	proton.dosrc= true;
	neutron.dosrc= true;

//set constants:
	double c0_=0.161;
	double c1_= -0.25;
	double phi0_= 2.38;
	double phi1_= -0.56;
	
	//test if recovers no src
	// double c0_=0.;
	// double c1_= 0.;
	// double phi0_= 1.; 
	// double phi1_= 0.;
	
	
	proton.phi_ = phi0_*(1.-phi1_*(1.-2.*Yp_));
	neutron.phi_= phi0_*(1.+phi1_*(1.-2.*Yp_));

	proton.c_= c0_*(1.-c1_*(1.-2.*Yp_));
	neutron.c_= c0_*(1.+c1_*(1.-2.*Yp_));

	proton.delta_ = 1.- 3.*proton.c_*(1.-1./proton.phi_);
	neutron.delta_= 1.- 3.*neutron.c_*(1.-1./neutron.phi_);

  proton.density=Yp*rhoB;
  proton.kf=pow(3.*pi2*proton.density, 1./3.);
  proton.kf2=pow(proton.kf, 2.);

  neutron.density=(1.-Yp)*rhoB;
  neutron.kf=pow(3.*pi2*neutron.density, 1./3.);
  neutron.kf2=pow(neutron.kf, 2.);

  rho3= proton.I3*proton.density+neutron.I3*neutron.density;
	
	setVectorMeanFields();
	setScalarMeanFields();			
	
  proton.solveChemPotEff();
  neutron.solveChemPotEff();

	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) 	: gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) 		: gb;

	double enerp_phi= sqrt(pow(proton.phi_*proton.kf, 2.) + pow(proton.mass_eff, 2.));
	double enern_phi= sqrt(pow(neutron.phi_*neutron.kf, 2.) + pow(neutron.mass_eff, 2.));


	double proton_chempot_src= 3.*proton.c_*(proton.chemPot_eff - enerp_phi/proton.phi_ )
						+ 4.*proton.c_*proton.kf
						*log( (proton.phi_*proton.kf + enerp_phi)/(proton.kf +proton.chemPot_eff));

	double neutron_chempot_src= 3.*neutron.c_*(neutron.chemPot_eff - enern_phi/neutron.phi_ )
					+ 4.*neutron.c_*neutron.kf
					*log( (neutron.phi_*neutron.kf + enern_phi)/(neutron.kf +neutron.chemPot_eff) );

  proton.chemPot  =  proton_chempot_src + proton.delta_*proton.chemPot_eff  + gv_*V0 + gb_*b0*proton.I3;
  neutron.chemPot =  neutron_chempot_src + neutron.delta_*neutron.chemPot_eff + gv_*V0 + gb_*b0*neutron.I3;

	
	if(dogRhob){
		proton.chemPot += getRearrangementEnergy();
		neutron.chemPot += getRearrangementEnergy();
	}

	muB = neutron.chemPot;
	muQ = proton.chemPot - neutron.chemPot;

	proton.calculateCondensate();
	neutron.calculateCondensate();

	rhoS= proton.condensate + neutron.condensate;
	
	proton.calculateProperties();
  neutron.calculateProperties();

}

//=============== Solve for the vector meson fields with Ceres library: ===============
void nlwm_class::setVectorMeanFields(){
	
	double v0_= getOmegaEffDens()*gv/pow(Mv , 2) ;       
	double b0_= getIsoEffDens()*gb/(pow(Mb, 2) );   

	if(dogRhob){
		v0_*=getCoupling_omega(rhoB);
		b0_*=getCoupling_rho(rhoB);
	}else{
		Problem pV;
    CostFunction* costV =	new NumericDiffCostFunction<VFunctor, ceres::CENTRAL, 1, 1>
														(new VFunctor(*this) );


    pV.AddResidualBlock(costV, NULL, &v0_);

   Solver::Options optionsV;
    // optionsV.parameter_tolerance = 1e-10;
    // optionsV.function_tolerance = 1e-8;
    // optionsV.gradient_tolerance=1e-12;
    //optionsV.trust_region_strategy_type = ceres::DOGLEG;
    //optionsV.dogleg_type = ceres::TRADITIONAL_DOGLEG;
		optionsV.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsV.linear_solver_type= ceres::DENSE_QR;

    optionsV.minimizer_progress_to_stdout = false;

    Solver::Summary summaryV;
    Solve(optionsV, &pV, &summaryV);
 //	std::cout << "v: " << summaryV.BriefReport() << "\n";

    b0_= (gb*getIsoEffDens())/( pow(Mb, 2.)+2.*Lv*pow(gv*gb*v0_, 2.));
	
	}
	V0=v0_;
	b0=b0_;
}
//=============== Functor Vector meson solver for Ceres: ===============
template <typename T>
bool VFunctor::operator()(const T* x, T* residuals) const{

	baryons.V0=x[0];
	baryons.b0=(baryons.gb*baryons.rho3)/( pow(baryons.Mb, 2.)
																+2.*baryons.Lv*pow(baryons.gv*baryons.gb*x[0], 2.));

  residuals[0] = baryons.omegaMeson_eom_residue(baryons.getOmegaEffDens());
	return true;
}


//===============  Solve for the scalar meson fields with Ceres library:  ===============
void nlwm_class::setScalarMeanFields(){
	double mef_=neutron.mass_eff;
	
	Problem pS;
	CostFunction* costS =	new NumericDiffCostFunction<SFunctor,ceres::CENTRAL, 1, 1>
														(new SFunctor(*this));

	pS.AddResidualBlock(costS, NULL, &mef_);
	pS.SetParameterLowerBound(&mef_, 0, 0.);
	pS.SetParameterUpperBound(&mef_, 0, 1.);

	// Set solver
	Solver::Options optionsS;
	optionsS.parameter_tolerance = 1e-10;
	optionsS.function_tolerance = 1e-10;
	optionsS.gradient_tolerance=1e-12;
	optionsS.dense_linear_algebra_library_type=ceres::LAPACK;
	optionsS.linear_solver_type= ceres::DENSE_QR;
	//optionsS.trust_region_strategy_type = ceres::DOGLEG;
	//optionsS.dogleg_type = ceres::TRADITIONAL_DOGLEG;
	optionsS.minimizer_progress_to_stdout = false;
	Solver::Summary summaryS;
	
	//Run
	Solve(optionsS, &pS, &summaryS);

	//Print if convergence was achieved.
	//std::cout << "M*: "  << summaryS.BriefReport() << "\n";
	//std::cout << "Meff: " << neutron.mass << " -> " << mef_ << std::endl;

	Mef= mef_;
	proton.mass_eff=Mef;
	neutron.mass_eff=Mef;

	double gs_= dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	phi0=(1.-Mef)/gs_;

}
//=============== Functor scalar meson solver for Ceres: ===============
template <typename T>
bool SFunctor::operator()(const T* x, T* residuals) const{

	baryons.proton.mass_eff=x[0];
	baryons.neutron.mass_eff=x[0];

	baryons.proton.solveChemPotEff();
	baryons.neutron.solveChemPotEff();

	baryons.proton.calculateCondensate();
	baryons.neutron.calculateCondensate();
	baryons.rhoS= baryons.proton.condensate + baryons.neutron.condensate;

	double gs_= baryons.dogRhob ? baryons.gs*baryons.getCoupling_sigma(baryons.rhoB) : baryons.gs;
	baryons.phi0= (1.-x[0])/gs_;
	double res=baryons.sigmaMeson_eom_residue(baryons.rhoS);
	residuals[0] = res;
	return true;
}


//=============== Set EoS for pure neutron matter  ===============
void nlwm_class::setEOS_neutrons(double rhoB_, double temp_, particle &electron_, particle &muon_){

	rhoB=rhoB_;
	
  setTemperature(temp_);
  
//	if(rhoB/rho0 < 1.) {doHyperons=false; doDeltas=false;}
	double mub_;
	double mue_;
	double phi0_ ;
	double v0_ ;       
	double b0_;
	
	if(firstRun){
		if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
		else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
	}else{
		mub_	= muB;
		mue_= electron_.chemPot;
		phi0_=phi0 ;
		v0_	=V0;       
		b0_	= b0;
	}

	// if(parametrization=="fsu2h"){

		double x[]={mub_, mue_, phi0_, v0_, b0_};
		

		Problem pBetaEq;
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<NeutronFunctor,ceres::CENTRAL, 5, 5>
								(new  NeutronFunctor(*this, electron_, muon_));

		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

		// Set solver
		Solver::Options optionsBetaEq;
	//if(parametrization!="iufsu"){
		optionsBetaEq.parameter_tolerance = 1e-8;
		optionsBetaEq.function_tolerance = 1e-10;
		optionsBetaEq.gradient_tolerance=1e-12;
	//}
		// optionsBetaEq.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
		// optionsBetaEq.linear_solver_type=ceres::SPARSE_NORMAL_CHOLESKY;
		
		
		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		//optionsBetaEq.update_state_every_iteration = true;
		//optionsBetaEq.use_explicit_schur_complement= true;
		
		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;
		optionsBetaEq.max_num_iterations=1e5;	

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ << 
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]
		<< std::endl << std::endl;
				
		mub_	 =x[0];
		mue_	 =x[1];
		phi0_	 =x[2];
		v0_		 =x[3];
		b0_		 =x[4];
		electron_.setLepton(mue_);
		//electron_.calculateProperties();
		muon_.setLepton(mue_);
		//muon_.calculateProperties();

		setDensities(mub_, -mue_,  phi0_,  v0_, b0_);
		setThermodynamics();
	//  std::cout << 
	//   rhoB 		-	getBaryonDens() << " "  <<
	//   getChargeDens() +	electron_.Qdens + muon_.Qdens << " "  <<
	//   sigmaMeson_eom_residue(	getSigmaEffDens()) << " "  <<
	//   omegaMeson_eom_residue(	getOmegaEffDens()) << " "  <<
	//   rhoMeson_eom_residue(		getIsoEffDens()) << " " << std::endl;
	 	firstRun=false;

}
//=============== Functor neutron matter solver for Ceres: ===============
template<typename T>
bool NeutronFunctor::operator()(const T* x, T* residuals) const{

	electron.setLepton(x[1]);
	muon.setLepton(x[1]);

	baryons.muB		=x[0];
	baryons.muQ		=x[1];
	baryons.phi0	=x[2];
	baryons.V0		=x[3];
	baryons.b0		=x[4];
	baryons.neutron.setBaryonEff(	x[0], -x[1], baryons.gs*x[2], baryons.gv*x[3], baryons.gb*x[4]);
	baryons.proton.setBaryonEff(		x[0], -x[1], baryons.gs*x[2], baryons.gv*x[3], baryons.gb*x[4]);

	residuals[0] = baryons.rhoB - baryons.getBaryonDens();
	residuals[1] = baryons.getChargeDens();//+	electron.Qdens + muon.Qdens;
	residuals[2] = baryons.sigmaMeson_eom_residue(	baryons.getSigmaEffDens());
	residuals[3] = baryons.omegaMeson_eom_residue(	baryons.getOmegaEffDens());
	residuals[4] = baryons.rhoMeson_eom_residue(		baryons.getIsoEffDens());
	//residuals[5] = baryons.thetaMeson_eom_residue(	baryons.getThetaEffDens()) ;
	return true;
}


//=============== Set EoS in beta equilibrium w/ electrons+muons  ===============
void nlwm_class::setEOS_betaEq(double rhoB_, double temp_, particle &electron_, particle &muon_){

	rhoB=rhoB_;
	setTemperature(temp_);
  
	double mub_;
	double mue_;
	double phi0_ ;
	double v0_ ;       
	double b0_;
	double theta0_=0.;
	double rear_=0.;

	if(dogRhob){ 
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
			theta0_= -0.028034;
			rear_= 	-0.0129762;
		}else{
			mub_	= muB;
			mue_= electron_.chemPot;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_=theta0;
			rear_=getRearrangementEnergy();
		}

		double x[]={mub_, mue_, phi0_, v0_, b0_, theta0_, rear_};
		Problem pBetaEq;
		
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<BetaEqFunctorDD,ceres::CENTRAL, 7, 7>
								(new  BetaEqFunctorDD(*this, electron_, muon_));
		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);
		pBetaEq.SetParameterLowerBound(x, 6, -0.35);
		pBetaEq.SetParameterUpperBound(x, 6, 0.35);


		Solver::Options optionsBetaEq;
	
		optionsBetaEq.parameter_tolerance = 1e-10;
		optionsBetaEq.function_tolerance = 1e-10;
		optionsBetaEq.gradient_tolerance=1e-12;

		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		optionsBetaEq.use_nonmonotonic_steps= true;
		optionsBetaEq.update_state_every_iteration = true;
		
		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;
		optionsBetaEq.max_num_iterations=1e5;	

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ <<  " " << theta0_ << " " << rear_ <<
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]  << " " << x[5] << " " << x[6]
		<< std::endl << std::endl;
				
		mub_	 =x[0];
		mue_	 =x[1];
		phi0_	 =x[2];
		v0_		 =x[3];
		b0_		 =x[4];
		theta0_=x[5];
		rear_=x[6];

		setDensitiesDD(mub_, -mue_,  phi0_,  v0_, b0_, theta0_, rear_);
		setThermodynamics();
		electron_.setLepton(mue_);
		electron_.calculateProperties();
		muon_.setLepton(mue_);
		muon_.calculateProperties();


	}else if( xtl!=0. || xts!=0. || xtx!=0.){ 
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_ = -0.0118828;
			if(parametrization=="fsu2h")theta0_ = -0.0252404;
			if(parametrization=="ddme2")theta0_= -0.028034;
		}else{
			mub_	= muB;
			mue_= electron_.chemPot;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_=theta0;
		}

		//must solve 4 meson equations + charge and mass equilibrium
		double x[]={mub_, mue_, phi0_, v0_, b0_, theta0_};

		Problem pBetaEq;
		
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<BetaEqFunctor2,ceres::CENTRAL, 6, 6>
								(new  BetaEqFunctor2(*this, electron_, muon_));
		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

		// if(temperature<Tmin_integration){
		//  	pBetaEq.SetParameterLowerBound(x, 0, 0.);
		// 	pBetaEq.SetParameterLowerBound(x, 1, electron_.mass_eff);
		//  	pBetaEq.SetParameterLowerBound(x, 2, 0.);
		//pBetaEq.SetParameterLowerBound(x, 3, 0.);
		// }
		// Set solver
		Solver::Options optionsBetaEq;
	//if(parametrization!="iufsu"){
		// if(Bfield==0){
			optionsBetaEq.parameter_tolerance = 1e-10;
			optionsBetaEq.function_tolerance = 1e-10;
			optionsBetaEq.gradient_tolerance=1e-12;
		
		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		optionsBetaEq.use_nonmonotonic_steps= true;
		optionsBetaEq.update_state_every_iteration = true;
		
		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;
		optionsBetaEq.max_num_iterations=1e5;	

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ <<  " " << theta0_ << 
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]  << " " << x[5]
		<< std::endl << std::endl;
				
		mub_	 =x[0];
		mue_	 =x[1];
		phi0_	 =x[2];
		v0_		 =x[3];
		b0_		 =x[4];
		theta0_=x[5];

		setDensities(mub_, -mue_,  phi0_,  v0_, b0_, theta0_);
		setThermodynamics();
		electron_.setLepton(mue_);
		electron_.calculateProperties();
		muon_.setLepton(mue_);
		muon_.calculateProperties();

	}else{
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
		}else{
			mub_	= muB;
			mue_= electron_.chemPot;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
		}

		double x[]={mub_, mue_, phi0_, v0_, b0_};

		Problem pBetaEq;
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<BetaEqFunctor,ceres::CENTRAL, 5, 5>
								(new  BetaEqFunctor(*this, electron_, muon_));

		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

		// if(temperature<Tmin_integration){
		// 	pBetaEq.SetParameterLowerBound(x, 0, 0.);
		//	pBetaEq.SetParameterLowerBound(x, 1, electron_.mass_eff);
		// 	pBetaEq.SetParameterLowerBound(x, 2, 0.);
		// 	pBetaEq.SetParameterLowerBound(x, 3, 0.);
		// }
		// Set solver
		Solver::Options optionsBetaEq;
	//if(parametrization!="iufsu"){
		optionsBetaEq.parameter_tolerance = 1e-10;
		optionsBetaEq.function_tolerance = 1e-10;
		optionsBetaEq.gradient_tolerance=1e-12;
		optionsBetaEq.max_num_iterations=1e6;	

		if(rhoB*pow(Mnucleon/hc, 3.) < (5.e-4) ){
			// optionsBetaEq.parameter_tolerance = 1e-22;
			// optionsBetaEq.function_tolerance = 1e-22;
			// optionsBetaEq.gradient_tolerance=1e-25;
			optionsBetaEq.parameter_tolerance = 1e-20;
			optionsBetaEq.function_tolerance = 1e-20;
			optionsBetaEq.gradient_tolerance=1e-23;
			optionsBetaEq.max_num_iterations=1e7;	

		}	
			if(rhoB*pow(Mnucleon/hc, 3.) < (5.e-7) ){
			// optionsBetaEq.parameter_tolerance = 1e-22;
			// optionsBetaEq.function_tolerance = 1e-22;
			// optionsBetaEq.gradient_tolerance=1e-25;
			optionsBetaEq.parameter_tolerance = 1e-35;
			optionsBetaEq.function_tolerance = 1e-35;
			optionsBetaEq.gradient_tolerance=1e-35;
			optionsBetaEq.max_num_iterations=1e8;	

		}	
	
	//}
		// optionsBetaEq.line_search_direction_type= ceres::STEEPEST_DESCENT;
		// optionsBetaEq.line_search_type=ceres::ARMIJO;

		// optionsBetaEq.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
		// optionsBetaEq.linear_solver_type=ceres::SPARSE_NORMAL_CHOLESKY;
		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		// optionsBetaEq.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		
		optionsBetaEq.use_nonmonotonic_steps= true;
		optionsBetaEq.update_state_every_iteration = true;
		// optionsBetaEq.use_explicit_schur_complement= true;
		

		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " << phi0_ << " " << v0_  << " " << b0_ << 
			"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4] 
		<< std::endl <<  std::endl;
				
		mub_	=x[0];
		mue_	=x[1];
		phi0_	=x[2];
		v0_		=x[3];
		b0_		=x[4];

		setDensities(mub_, -mue_,  phi0_,  v0_, b0_);
		
		setThermodynamics();
		electron_.setLepton(mue_);
		electron_.calculateProperties();
		muon_.setLepton(mue_);
		muon_.calculateProperties();

	}
	
 	firstRun=false;

}
//=============== Functor beta-equilibrium for Ceres w/ 3 mesons: ===============
template <typename T>
bool BetaEqFunctor::operator()(const T* x, T* residuals) const{

	electron.setLepton(x[1]);
	muon.setLepton(x[1]);
	baryons.setDensities(x[0], -x[1], x[2], x[3], x[4]);
											//mub, muq, phi0,  v0,    b0
	residuals[0] = baryons.rhoB 		-	baryons.getBaryonDens();
	residuals[1] = baryons.getChargeDens() +	electron.Qdens + muon.Qdens;
	residuals[2] = baryons.sigmaMeson_eom_residue(	baryons.getSigmaEffDens());
	residuals[3] = baryons.omegaMeson_eom_residue(	baryons.getOmegaEffDens());
	residuals[4] = baryons.rhoMeson_eom_residue(		baryons.getIsoEffDens());
	return true;
}
//=============== Functor beta-equilibrium for Ceres w/ 4 mesons: ===============
template <typename T>
bool BetaEqFunctor2::operator()(const T* x, T* residuals) const{

	electron.setLepton(x[1]);
	muon.setLepton(x[1]);
	baryons.setDensities(x[0], -x[1], x[2], x[3], x[4], x[5]);
											//mub, muq, phi0,  v0,    b0, theta0

	residuals[0] = baryons.rhoB - baryons.getBaryonDens();
	residuals[1] = baryons.getChargeDens() +	electron.Qdens + muon.Qdens;
	residuals[2] = baryons.sigmaMeson_eom_residue(	baryons.getSigmaEffDens());
	residuals[3] = baryons.omegaMeson_eom_residue(	baryons.getOmegaEffDens());
	residuals[4] = baryons.rhoMeson_eom_residue(		baryons.getIsoEffDens());
	residuals[5] = baryons.thetaMeson_eom_residue(	baryons.getThetaEffDens()) ;
	return true;
}

template <typename T>
bool BetaEqFunctorDD::operator()(const T* x, T* residuals) const{

	electron.setLepton(x[1]);
	muon.setLepton(x[1]);
	baryons.setDensitiesDD(x[0], -x[1], x[2], x[3], x[4], x[5], x[6]);
											//mub, muq, phi0,  v0,    b0, theta0

	residuals[0] = baryons.rhoB - baryons.getBaryonDens();
	residuals[1] = baryons.getChargeDens() +	electron.Qdens + muon.Qdens;
	residuals[2] = baryons.sigmaMeson_eom_residue(	baryons.getSigmaEffDens());
	residuals[3] = baryons.omegaMeson_eom_residue(	baryons.getOmegaEffDens());
	residuals[4] = baryons.rhoMeson_eom_residue(		baryons.getIsoEffDens());
	residuals[5] = baryons.thetaMeson_eom_residue(	baryons.getThetaEffDens()) ;
	residuals[6] = baryons.getRearrangementEnergy() - x[6];
	return true;
}


//=============== Set EoS in beta equilibrium w/ electrons (needed for low densities only)===============
void nlwm_class::setEOS_betaEq(double rhoB_, double temp_, particle &electron_){

	rhoB=rhoB_;
	setTemperature(temp_);
	
	double mue_= electron_.chemPot; 
	Problem pBetaEq;
	CostFunction* costBetaEq= 
							new NumericDiffCostFunction<BetaEqFunctor3,ceres::CENTRAL, 1, 1>
							(new  BetaEqFunctor3(*this, electron_));

	pBetaEq.AddResidualBlock(costBetaEq, NULL, &mue_);

	// Set solver
	Solver::Options optionsBetaEq;
	optionsBetaEq.parameter_tolerance = 1e-8;
	optionsBetaEq.function_tolerance = 1e-10;
	optionsBetaEq.gradient_tolerance=1e-12;
	optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
	optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
	//optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
	//optionsBetaEq.dogleg_type = ceres::TRADITIONAL_DOGLEG;
	optionsBetaEq.minimizer_progress_to_stdout = false;
	Solver::Summary summaryBetaEq;
	
	//Run
	Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

	//Print if convergence was achieved.
	std::cout << summaryBetaEq.BriefReport() << "\n";

	electron_.chemPot= mue_;
	electron_.chemPot_eff=mue_;
	electron_.calculateDensity();
	electron_.calculateProperties();

	Yp= electron_.density/rhoB;
  proton.density=electron_.density;
  proton.kf=pow(3.*pi2*proton.density, 1./3.);
  proton.kf2=pow(proton.kf, 2.);

  neutron.density=(1.-Yp)*rhoB;
  neutron.kf=pow(3.*pi2*neutron.density, 1./3.);
  neutron.kf2=pow(neutron.kf, 2.);

  rho3=proton.I3*proton.density+neutron.I3*neutron.density;
	
	setVectorMeanFields();
	setScalarMeanFields();			
	
  proton.solveChemPotEff();
  proton.calculateProperties();

  neutron.solveChemPotEff();
  neutron.calculateProperties();
	
  proton.chemPot  =  proton.chemPot_eff  + gv*V0 +proton.I3*gb*b0;
  neutron.chemPot =  neutron.chemPot_eff + gv*V0 +neutron.I3*gb*b0;
	
	proton.calculateCondensate();
	neutron.calculateCondensate();
	rhoS= proton.condensate+ neutron.condensate;
}
//=============== Functor beta-equilibrium w/ only electrons for Ceres (3 mesons only): ===============
template <typename T>
bool BetaEqFunctor3::operator()(const T* x, T* residuals) const{

	// double Yp_=x[0];
	electron.chemPot=x[0];
	electron.chemPot_eff=x[0];
	electron.kf2=pow(electron.chemPot, 2.)- pow(electron.mass_eff, 2.);	
  electron.kf2<=0. ? electron.kf=0. : electron.kf=sqrt(electron.kf2);

	electron.calculateDensity();
	// electron.solveChemPotEff();
	// electron.chemPot=electron.chemPot_eff;

	double Yp_=electron.density/baryons.rhoB;
  baryons.proton.density=Yp_*baryons.rhoB;
  baryons.proton.kf=pow(3.*pi2*baryons.proton.density, 1./3.);
  baryons.proton.kf2=pow(baryons.proton.kf, 2.);

  baryons.neutron.density=(1.-Yp_)*baryons.rhoB;
  baryons.neutron.kf=pow(3.*pi2*baryons.neutron.density, 1./3.);
  baryons.neutron.kf2=pow(baryons.neutron.kf, 2.);

  baryons.rho3=baryons.proton.I3*baryons.proton.density
						 +baryons.neutron.I3*baryons.neutron.density;
	
	baryons.setVectorMeanFields();
	baryons.setScalarMeanFields();			
	
  baryons.proton.solveChemPotEff();
  baryons.neutron.solveChemPotEff();
	
  baryons.proton.chemPot  =  baryons.proton.chemPot_eff  
														+baryons.gv*baryons.V0 + baryons.gb*baryons.b0/2.;
  baryons.neutron.chemPot =  baryons.neutron.chemPot_eff 
														+baryons.gv*baryons.V0 - baryons.gb*baryons.b0/2.;
	residuals[0] =baryons.neutron.chemPot - baryons.proton.chemPot- electron.chemPot;
	return true;
}

//=============== Functor beta-equilibrium w/ only electrons for Ceres (3 mesons only): ===============
void nlwm_class::setEOS_fixedYl(double rhoB_, double temp_, double Yle_, double Ylm_,	
										particle &electron_, particle &muon_, particle &ne_, particle &nm_){
	rhoB=rhoB_;
	setTemperature(temp_);
  
	double mub_;
	double mue_;
	double munue_;
	double munum_;
	double phi0_ ;
	double v0_ ;       
	double b0_;
	double theta0_=0.;
	double rear_=0.;

	if(dogRhob){ 
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);
			munue_=mue_/100.;
			munum_= munue_;
			theta0_= -0.028034;
			rear_= 	-0.0129762;
		}else{
			mub_	= muB;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_=theta0;
			rear_=getRearrangementEnergy();
		}


		double x[]={mub_, mue_, munue_, munum_, phi0_,  v0_, b0_, theta0_, rear_};
		Problem pYl;
		
		CostFunction* costYl= 
								new NumericDiffCostFunction<YlFunctorDD,ceres::CENTRAL, 9, 9>
								(new  YlFunctorDD(*this, electron_, muon_, ne_, nm_, Yle_, Ylm_));
		pYl.AddResidualBlock(costYl, NULL, x);
		pYl.SetParameterLowerBound(x, 8, -0.35);
		pYl.SetParameterUpperBound(x, 8, 0.35);


		Solver::Options optionsYl;
	
		optionsYl.parameter_tolerance = 1e-10;
		optionsYl.function_tolerance = 1e-10;
		optionsYl.gradient_tolerance=1e-12;

		optionsYl.linear_solver_type= ceres::DENSE_QR;
		optionsYl.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsYl.trust_region_strategy_type = ceres::DOGLEG;
		optionsYl.dogleg_type = ceres::SUBSPACE_DOGLEG;
		optionsYl.use_nonmonotonic_steps= true;
		optionsYl.update_state_every_iteration = true;
		
		optionsYl.minimizer_progress_to_stdout = false;
		Solver::Summary summaryYl;
		optionsYl.max_num_iterations=1e5;	

		//Run
		Solve(optionsYl, &pYl, &summaryYl);

		//Print if convergence was achieved.
		std::cout << summaryYl.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " <<  munue_ << " " <<  munum_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ <<  " " << theta0_ << " " << rear_ <<
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]  << " " 
					<< x[5] << " "<< x[6] << " " << x[7] << " " << x[8]
		<< std::endl << std::endl;
				
		mub_	 =x[0];
		mue_	 =x[1];
		munue_ =x[2];
		munum_ =x[3];
		phi0_	 =x[4];
		v0_		 =x[5];
		b0_		 =x[6];
		theta0_=x[7];
		rear_	 =x[8];

		double muq_= munue_- mue_;
		setDensitiesDD(mub_, muq_,  phi0_,  v0_, b0_, theta0_, rear_);
		setThermodynamics();

		double mum_= x[3]- muq_;

		electron_.setLepton(mue_);
		ne_.setLepton(munue_);
		muon_.setLepton(mum_);
		nm_.setLepton(munum_);

		electron_.calculateProperties();
		muon_.calculateProperties();
		ne_.calculateProperties();
		nm_.calculateProperties();

	}else if( xtl!=0. || xts!=0. || xtx!=0.){ 
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_ = -0.0118828;
			if(parametrization=="fsu2h")theta0_ = -0.0252404;
			if(parametrization=="ddme2")theta0_= -0.028034;
		}else{
			mub_	= muB;
			mue_= electron_.chemPot;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
			if(xtl!=0. || xts!=0. || xtx!=0.) theta0_=theta0;
		}

		//must solve 4 meson equations + charge and mass equilibrium
		double x[]={mub_, mue_, phi0_, v0_, b0_, theta0_};

		Problem pBetaEq;
		
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<BetaEqFunctor2,ceres::CENTRAL, 6, 6>
								(new  BetaEqFunctor2(*this, electron_, muon_));
		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

		// if(temperature<Tmin_integration){
		//  	pBetaEq.SetParameterLowerBound(x, 0, 0.);
		// 	pBetaEq.SetParameterLowerBound(x, 1, electron_.mass_eff);
		//  	pBetaEq.SetParameterLowerBound(x, 2, 0.);
		//pBetaEq.SetParameterLowerBound(x, 3, 0.);
		// }
		// Set solver
		Solver::Options optionsBetaEq;
	//if(parametrization!="iufsu"){
		// if(Bfield==0){
		optionsBetaEq.parameter_tolerance = 1e-10;
		optionsBetaEq.function_tolerance = 1e-10;
		optionsBetaEq.gradient_tolerance=1e-12;
		
		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		optionsBetaEq.use_nonmonotonic_steps= true;
		optionsBetaEq.update_state_every_iteration = true;
		
		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;
		optionsBetaEq.max_num_iterations=1e5;	

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " 
							<< phi0_ << " " << v0_  << " " << b0_ <<  " " << theta0_ << 
		"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4]  << " " << x[5]
		<< std::endl << std::endl;
				
		mub_	 =x[0];
		mue_	 =x[1];
		phi0_	 =x[2];
		v0_		 =x[3];
		b0_		 =x[4];
		theta0_=x[5];

		setDensities(mub_, -mue_,  phi0_,  v0_, b0_, theta0_);
		setThermodynamics();
		electron_.setLepton(mue_);
		electron_.calculateProperties();
		muon_.setLepton(mue_);
		muon_.calculateProperties();

	}else{
		if(firstRun){
			if(Bfield==0) setInitial_hd(mub_, mue_, phi0_, v0_, b0_);
			else					setInitial_hdb(mub_, mue_, phi0_, v0_, b0_);		
		}else{
			mub_	= muB;
			mue_= electron_.chemPot;
			phi0_=phi0 ;
			v0_	=V0;       
			b0_	= b0;
		}

		double x[]={mub_, mue_, phi0_, v0_, b0_};

		Problem pBetaEq;
		CostFunction* costBetaEq= 
								new NumericDiffCostFunction<BetaEqFunctor,ceres::CENTRAL, 5, 5>
								(new  BetaEqFunctor(*this, electron_, muon_));

		pBetaEq.AddResidualBlock(costBetaEq, NULL, x);

		// if(temperature<Tmin_integration){
		// 	pBetaEq.SetParameterLowerBound(x, 0, 0.);
		//	pBetaEq.SetParameterLowerBound(x, 1, electron_.mass_eff);
		// 	pBetaEq.SetParameterLowerBound(x, 2, 0.);
		// 	pBetaEq.SetParameterLowerBound(x, 3, 0.);
		// }
		// Set solver
		Solver::Options optionsBetaEq;
	//if(parametrization!="iufsu"){
		optionsBetaEq.parameter_tolerance = 1e-10;
		optionsBetaEq.function_tolerance = 1e-10;
		optionsBetaEq.gradient_tolerance=1e-12;
		optionsBetaEq.max_num_iterations=1e6;	

		if(rhoB*pow(Mnucleon/hc, 3.) < (5.e-4) ){
			// optionsBetaEq.parameter_tolerance = 1e-22;
			// optionsBetaEq.function_tolerance = 1e-22;
			// optionsBetaEq.gradient_tolerance=1e-25;
			optionsBetaEq.parameter_tolerance = 1e-20;
			optionsBetaEq.function_tolerance = 1e-20;
			optionsBetaEq.gradient_tolerance=1e-23;
			optionsBetaEq.max_num_iterations=1e7;	

		}	
			if(rhoB*pow(Mnucleon/hc, 3.) < (5.e-7) ){
			// optionsBetaEq.parameter_tolerance = 1e-22;
			// optionsBetaEq.function_tolerance = 1e-22;
			// optionsBetaEq.gradient_tolerance=1e-25;
			optionsBetaEq.parameter_tolerance = 1e-35;
			optionsBetaEq.function_tolerance = 1e-35;
			optionsBetaEq.gradient_tolerance=1e-35;
			optionsBetaEq.max_num_iterations=1e8;	

		}	
	
	//}
		// optionsBetaEq.line_search_direction_type= ceres::STEEPEST_DESCENT;
		// optionsBetaEq.line_search_type=ceres::ARMIJO;

		// optionsBetaEq.sparse_linear_algebra_library_type=ceres::SUITE_SPARSE;
		// optionsBetaEq.linear_solver_type=ceres::SPARSE_NORMAL_CHOLESKY;
		optionsBetaEq.linear_solver_type= ceres::DENSE_QR;
		optionsBetaEq.dense_linear_algebra_library_type=ceres::LAPACK;
		// optionsBetaEq.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		
		optionsBetaEq.trust_region_strategy_type = ceres::DOGLEG;
		optionsBetaEq.dogleg_type = ceres::SUBSPACE_DOGLEG;
		
		optionsBetaEq.use_nonmonotonic_steps= true;
		optionsBetaEq.update_state_every_iteration = true;
		// optionsBetaEq.use_explicit_schur_complement= true;
		

		optionsBetaEq.minimizer_progress_to_stdout = false;
		Solver::Summary summaryBetaEq;

		//Run
		Solve(optionsBetaEq, &pBetaEq, &summaryBetaEq);

		//Print if convergence was achieved.
		std::cout << summaryBetaEq.BriefReport() << "\n";
		std::cout << "rhob= " << rhoB*pow(Mnucleon/hc, 3) << std::endl;
		std::cout << mub_ << " " << mue_ << " " << phi0_ << " " << v0_  << " " << b0_ << 
			"---> "<< x[0] << " " << x[1] << " " << x[2]  << " " << x[3] << " " << x[4] 
		<< std::endl <<  std::endl;
				
		mub_	=x[0];
		mue_	=x[1];
		phi0_	=x[2];
		v0_		=x[3];
		b0_		=x[4];

		setDensities(mub_, -mue_,  phi0_,  v0_, b0_);
		
		setThermodynamics();
		electron_.setLepton(mue_);
		electron_.calculateProperties();
		muon_.setLepton(mue_);
		muon_.calculateProperties();

	}
	
 	firstRun=false;

}


template <typename T>
bool YlFunctorDD::operator()(const T* x, T* residuals) const{
//mub_, mue_, munue_, munum_, phi0_,  v0_, b0_, theta0_, rear_
	double muq_= x[2] - x[1];
	double mum_= x[3]- muq_;

	electron.setLepton(x[1]);
	ne.setLepton(x[2]);
	muon.setLepton(mum_);
	nm.setLepton(x[3]);

	baryons.setDensitiesDD(x[0], muq_, x[4], x[5], x[6], x[7], x[8]);
											//mub, muq, phi0,  v0,    b0, theta0

	residuals[0] = baryons.rhoB - baryons.getBaryonDens();
	residuals[1] = baryons.getChargeDens() +	electron.Qdens + muon.Qdens;
	residuals[2] = baryons.sigmaMeson_eom_residue(	baryons.getSigmaEffDens());
	residuals[3] = baryons.omegaMeson_eom_residue(	baryons.getOmegaEffDens());
	residuals[4] = baryons.rhoMeson_eom_residue(		baryons.getIsoEffDens());
	residuals[5] = baryons.thetaMeson_eom_residue(	baryons.getThetaEffDens()) ;
	residuals[6] = baryons.getRearrangementEnergy() - x[8];
	residuals[7] = Yle*baryons.getBaryonDens() - electron.density 	- ne.density;
	residuals[8] = Ylm*baryons.getBaryonDens() - muon.density 			- nm.density;
 	return true;
}

//=============== Set nucleon EoS inputing effective chemical potentials/mass ===============
void nlwm_class::setEOS_coexistence(double nup_, double nun_, double mef_){
	
  proton.chemPot_eff= nup_;
  proton.mass_eff= mef_;
  proton.kf2= proton.chemPot_eff*proton.chemPot_eff -proton.mass_eff*proton.mass_eff;
  proton.kf2<=0. ? proton.kf=0. : proton.kf=sqrt(proton.kf2);
	

  neutron.chemPot_eff= nun_;
  neutron.mass_eff= mef_;
  neutron.kf2= neutron.chemPot_eff*neutron.chemPot_eff -neutron.mass_eff*neutron.mass_eff;
  neutron.kf2<0. ? neutron.kf=0	: neutron.kf=sqrt(neutron.kf2);
	
	if(temperature<Tmin_integration){
		proton.kf==0. ? proton.density=0. : proton.density		=	proton.gamma*pow(proton.kf, 3.)/(6.*pi2); //integrate(densityFunc, &proton);
		neutron.kf==0. ? neutron.density=0. : neutron.density	= neutron.gamma*pow(neutron.kf, 3.)/(6.*pi2); //integrate(densityFunc, &neutron);

	}else{
		proton.density= integrate(densityFunc, &proton);
		neutron.density= integrate(densityFunc, &neutron);
	}
  proton.calculateProperties();
  neutron.calculateProperties();

  rhoB= proton.density + neutron.density;
  rho3= proton.I3*proton.density + neutron.I3*neutron.density;
	proton.calculateCondensate();
	neutron.calculateCondensate();
  rhoS= proton.condensate + neutron.condensate;
	Yp = proton.density/rhoB;	
	
	Mef=mef_;
	phi0=(1.-mef_	)/gs;
	setVectorMeanFields();
	
  proton.chemPot  = proton.chemPot_eff  + gv*V0 + proton.I3*gb*b0;
  neutron.chemPot = neutron.chemPot_eff + gv*V0 + neutron.I3*gb*b0;

}


//=============== Set nucleon EoS inputing effective chemical potentials/mass 
//																 					with short range correlations	 ===============
void nlwm_class::setEOS_coexistence_src(double nup_, double nun_, double mef_, double yp_){
	
  proton.chemPot_eff= nup_;
  proton.mass_eff= mef_;
  proton.kf2= proton.chemPot_eff*proton.chemPot_eff -proton.mass_eff*proton.mass_eff;
  proton.kf2<=0. ? proton.kf=0. : proton.kf=sqrt(proton.kf2);

  neutron.chemPot_eff= nun_;
  neutron.mass_eff= mef_;
  neutron.kf2= neutron.chemPot_eff*neutron.chemPot_eff -neutron.mass_eff*neutron.mass_eff;
  neutron.kf2<0. ? neutron.kf=0	: neutron.kf=sqrt(neutron.kf2);
	
	Yp= yp_;

	double c0_=0.161;
	double c1_= -0.25;
	double phi0_= 2.38;
	double phi1_= -0.56;
	
	// double phi0_= 1.; //test if recovers no src
	// double phi1_= 0.;
	// double c0_=0.;
	// double c1_=0.;
	proton.phi_ = phi0_*(1.-phi1_*(1.-2.*Yp));
	neutron.phi_= phi0_*(1.+phi1_*(1.-2.*Yp));

	proton.c_= c0_*(1.-c1_*(1.-2.*Yp));
	neutron.c_= c0_*(1.+c1_*(1.-2.*Yp));

	proton.delta_ = 1.- 3.*proton.c_*(1.-1./proton.phi_);
	neutron.delta_= 1.- 3.*neutron.c_*(1.-1./neutron.phi_);

	double enerp_phi= sqrt(pow(proton.phi_*proton.kf, 2.) + pow(proton.mass_eff, 2.));
	double enern_phi= sqrt(pow(neutron.phi_*neutron.kf, 2.) + pow(neutron.mass_eff, 2.));

	if(temperature<Tmin_integration){
		// proton.kf==0. ? proton.density=0. : proton.density		=	proton.gamma*pow(proton.kf, 3.)/(6.*pi2); //integrate(densityFunc, &proton);
		// neutron.kf==0. ? neutron.density=0. : neutron.density	= neutron.gamma*pow(neutron.kf, 3.)/(6.*pi2); //integrate(densityFunc, &neutron);

			proton.kf==0. ? proton.density=0. : 
					 	proton.density		=	proton.gamma*pow(proton.kf, 3.)/(6.*pi2);
			
			neutron.kf==0. ? neutron.density=0. : 
						neutron.density =	neutron.gamma*pow(neutron.kf, 3.)/(6.*pi2);

	}else{
		proton.density= integrate(densityFunc, &proton);
		neutron.density= integrate(densityFunc, &neutron);
	}
	
  proton.calculateProperties();
  neutron.calculateProperties();

	double proton_chempot_src= 3.*proton.c_*(proton.chemPot_eff - enerp_phi/proton.phi_ )
						+ 4.*proton.c_*proton.kf*log( (proton.phi_*proton.kf + enerp_phi)/(proton.kf +proton.chemPot_eff));

	double neutron_chempot_src= 3.*neutron.c_*(neutron.chemPot_eff - enern_phi/neutron.phi_ )
					+ 4.*neutron.c_*neutron.kf*log( (neutron.phi_*neutron.kf + enern_phi)/(neutron.kf +neutron.chemPot_eff) );

// cout << "neutron okay: " << neutron.density << " " << endl;
 
  rhoB= proton.density + neutron.density;
  rho3= proton.I3*proton.density + neutron.I3*neutron.density;
	proton.calculateCondensate();
	neutron.calculateCondensate();
  rhoS= proton.condensate + neutron.condensate;
	
	Mef=mef_;
	phi0=(1.-mef_	)/gs;
	setVectorMeanFields();
	
  proton.chemPot  =  proton_chempot_src + proton.delta_*proton.chemPot_eff  + gv*V0 + gb*b0*proton.I3;
  neutron.chemPot =  neutron_chempot_src + neutron.delta_*neutron.chemPot_eff + gv*V0 + gb*b0*neutron.I3;
}

//=============== Set baryon densities inputing chemical potentials and (3) meson fields ===============
void nlwm_class::setDensities(double mub_, double muq_, double phi0_, double v0_, double b0_){
	 
	muB=mub_;
	muQ=muq_;
	phi0=phi0_;
	V0=v0_;
	b0=b0_;
	double gs_= 	dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	double gv_= 	dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= 	dogRhob ? gb*getCoupling_rho(rhoB) : gb;
	double rear_= dogRhob ? getRearrangementEnergy() : 0.;	

	neutron.setBaryonEff(		muB, muQ, gs_*phi0, gv_*V0 + rear_, gb_*b0);
	proton.setBaryonEff(		muB, muQ, gs_*phi0, gv_*V0 + rear_, gb_*b0);

	if(doHyperons){ 
		lambda0.setBaryonEff(	muB, muQ, gs_*xsl*phi0, gv_*xvl*V0 + rear_, gb_*xbl*b0);
		sigmap.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*xvs*V0 + rear_, gb_*xbs*b0);
		sigma0.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*xvs*V0 + rear_, gb_*xbs*b0);
		sigmam.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*xvs*V0 + rear_, gb_*xbs*b0);
		xi0.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*xvx*V0 + rear_, gb_*xbx*b0);
		xim.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*xvx*V0 + rear_, gb_*xbx*b0);
	}
	
	if(doDeltas){
		deltam.setBaryonEff(	muB, muQ, xsd*gs_*phi0, xvd*gv_*V0 + rear_, xbd*gb_*b0);
		delta0.setBaryonEff(	muB, muQ, xsd*gs_*phi0, xvd*gv_*V0 + rear_, xbd*gb_*b0);
		deltap.setBaryonEff(	muB, muQ, xsd*gs_*phi0, xvd*gv_*V0 + rear_, xbd*gb_*b0);
		deltapp.setBaryonEff(	muB, muQ, xsd*gs_*phi0, xvd*gv_*V0 + rear_, xbd*gb_*b0);
	}

}


//=============== Set baryon densities inputing chemical potentials and (4) meson fields ===============
void nlwm_class::setDensities(double mub_, double muq_, double phi0_, double v0_, double b0_, double theta0_){
	 
	muB=mub_;
	muQ=muq_;
	phi0=phi0_;
	V0=v0_;
	b0=b0_;
	theta0=theta0_;

	double gs_= 	dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	double gv_= 	dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= 	dogRhob ? gb*getCoupling_rho(rhoB) 	 : gb;
	double rear_= dogRhob ? getRearrangementEnergy() 		 : 0.;	

	neutron.setBaryonEff(	muB, muQ, gs_*phi0, gv_*V0+rear_, gb_*b0);
	proton.setBaryonEff(	muB, muQ, gs_*phi0, gv_*V0+rear_, gb_*b0);

	if(doHyperons){ 
		lambda0.setBaryonEff(	muB, muQ, gs_*xsl*phi0, gv_*(xvl*V0+xtl*theta0)+rear_, gb_*xbl*b0);
		sigmap.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		sigma0.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		sigmam.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		xi0.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*(xvx*V0+xtx*theta0)+rear_, gb_*xbx*b0);
		xim.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*(xvx*V0+xtx*theta0)+rear_, gb_*xbx*b0);
	}
	
	if(doDeltas){
		deltam.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		delta0.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		deltap.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		deltapp.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
	}

	yN=(proton.density + neutron.density)/getBaryonDens();
	yH=(lambda0.density+ sigmap.density+ sigma0.density+ sigmam.density
													+ xi0.density+ xim.density)/getBaryonDens();
	yD=(deltapp.density + deltap.density+ delta0.density + deltam.density)/getBaryonDens();
}

void nlwm_class::setDensitiesDD(double mub_, double muq_, double phi0_, double v0_, double b0_, double theta0_, double rear_){
	 
	muB=mub_;
	muQ=muq_;
	phi0=phi0_;
	V0=v0_;
	b0=b0_;
	theta0=theta0_;

	double gs_= 	dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	double gv_= 	dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= 	dogRhob ? gb*getCoupling_rho(rhoB) 	 : gb;

	neutron.setBaryonEff(	muB, muQ, gs_*phi0, gv_*V0+rear_, gb_*b0);
	proton.setBaryonEff(	muB, muQ, gs_*phi0, gv_*V0+rear_, gb_*b0);

	if(doHyperons){ 
		lambda0.setBaryonEff(	muB, muQ, gs_*xsl*phi0, gv_*(xvl*V0+xtl*theta0)+rear_, gb_*xbl*b0);
		sigmap.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		sigma0.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		sigmam.setBaryonEff(	muB, muQ, gs_*xss*phi0, gv_*(xvs*V0+xts*theta0)+rear_, gb_*xbs*b0);
		xi0.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*(xvx*V0+xtx*theta0)+rear_, gb_*xbx*b0);
		xim.setBaryonEff(			muB, muQ, gs_*xsx*phi0, gv_*(xvx*V0+xtx*theta0)+rear_, gb_*xbx*b0);
	}
	
	if(doDeltas){
		deltam.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		delta0.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		deltap.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
		deltapp.setBaryonEff(	muB, muQ, gs_*xsd*phi0, gv_*xvd*V0+rear_, gb_*xbd*b0);
	}

	yN=(proton.density + neutron.density)/getBaryonDens();
	yH=(lambda0.density+ sigmap.density+ sigma0.density+ sigmam.density
													+ xi0.density+ xim.density)/getBaryonDens();
	yD=(deltapp.density + deltap.density+ delta0.density + deltam.density)/getBaryonDens();
}


//=============== Set initial guess for beta-equil. solver w/o magnetic field: ===============
void nlwm_class::setInitial_hd(double &mub_, double &mue_, double  &phi0_, double &v0_, double &b0_){

		// mub_	= 1.5;	//1.5
		// mue_=0.1;			//0.1								
	 	// phi0_=0.1; 		// 0.1 low densities		.2							
	 	// v0_	=0.07;  	// 0.1 low densities 			.15			
	 	// b0_	= -0.01;	// -1e-5 low dens	
	if(parametrization=="iufsu" && parhyp=="gm"){
		// mub_	= 1.41069;	//1.37237 0.0848789 0.0938846 0.0715153 -0.0025282
		// mue_=0.0619547;		
	 	// phi0_=0.0999732; 	
	 	// v0_	=0.0767674;  	
	 	// b0_	= -0.00229432;
		 mub_=1.40352;
		 mue_= 0.0668899;
		 phi0_= 0.0981881;
		 v0_= 0.0754749;
		 b0_= -0.00229697;
	}else if(parametrization=="iufsu" && parhyp=="su3"){
		mub_	= 1.5;	//1.5 1.42001 0.0749067 0.100246 0.0771754 -0.00261042
		mue_=0.13;			//0.1								
	 	phi0_=0.09; 		// 0.1 low densities		.2							
	 	v0_	=0.08;  	// 0.1 low densities 			.15			
	 	b0_	= -0.003;	// -1e-5 low dens				
	}else if(parametrization=="fsu2h"){
		// mub_=1.45583;//1.45583 0.0638229 0.0984247 0.0822707 -0.00164294 -0.00889494
		// mue_= 0.0638229;
		// phi0_= 0.0984247;
		// v0_= 0.0822707;
		// b0_= -0.00164294;	
		mub_	=1.72376;	//1.5 1.42001 0.0749067 0.100246 0.0771754 -0.00261042
		mue_=		 0.13766; 			//0.1								
	 	phi0_=		 0.0953097;// 0.1 low densities		.2							
	 	v0_	=		0.0952903;   	// 0.1 low densities 			.15			
	 	b0_	=-0.00105593 ;	// -1e-5 low dens				
	
	}else if (parametrization=="nl3"){
		mub_=1.40371;
		mue_= 0.125273;
		phi0_= 0.0899285;
		v0_= 0.0776314;
		b0_= -0.00972853;
	}else if(parametrization=="nl3wr"){
		mub_=1.50497;
		mue_= 0.0242015;
		phi0_= 0.0956183;
		v0_= 0.0901114;
		b0_= -0.000437992;
	}else if(parametrization=="l3wr"){
		mub_=1.67397;
		mue_=0.0795405;
		phi0_=0.110646;
		v0_=0.120147;
		b0_=-0.00466873;
	}else if(parametrization=="nl3wr*"){
		mub_=2.56661;
		mue_= 0.12955;
		phi0_= 0.0990648;
		v0_= 0.161139;
		b0_= -0.000221779	;
	}else if(parametrization=="ddme2"){
		mub_	= 1.82908;
		mue_=0.0883155;		
	 	phi0_=0.0966946; 	
	 	v0_	=0.12012;  	
	 	b0_	= -0.000751007;	
	}else{
		mub_	= 1.5;	//1.5
		mue_=0.1;			//0.1								
	 	phi0_=0.09; 		// 0.1 low densities		.2							
	 	v0_	=0.08;  	// 0.1 low densities 			.15			
	 	b0_	= -0.003;	// -1e-5 low dens				
	}
}
//=============== Set initial guess for beta-equil. solver W/ magnetic field: ===============
void nlwm_class::setInitial_hdb(double &mub_, double &mue_, double  &phi0_, double &v0_, double &b0_){

		// mub_	= 1.5;	//1.5
		// mue_=0.1;			//0.1								
	 	// phi0_=0.1; 		// 0.1 low densities		.2							
	 	// v0_	=0.07;  	// 0.1 low densities 			.15			
	 	// b0_	= -0.01;	// -1e-5 low dens	
	if(parametrization=="iufsu"){// && parhyp=="gm"
		if(pardelta=="su6"){ 
		// nyd, b>1
		// mub_=1.39907; 
		// mue_= 0.0679872;
		// phi0_= 0.0976974;
		// v0_= 0.0749879;
		// b0_= -0.0022934;
		//nyd, b<1
		// mub_=1.37436;
		// mue_= 0.0683921;
		// phi0_= 0.0981716;
		// v0_= 0.0736323;
		// b0_= -0.0029034;
		//nd,b>1
		mub_=1.26449;
		mue_=0.134712;
		phi0_=0.0800846;
		v0_=0.0576118;
		b0_= -0.00341266;
		//nd, b<1
		// mub_=1.40589;
		// mue_= 0.0913061;
		// phi0_= 0.100226;
		// v0_= 0.0750137;
		// b0_= -0.00402991;
		//nd, b<1, rhoinit=0.55
		// mub_=1.23423 ;  
		// mue_=0.142305;
		// phi0_=0.0779439;
		// v0_=0.0549942;
		// b0_=-0.00335551;
		 }else{//(pardelta=="mplA_1"){1.1.18884 0.103163 0.081786 0.0554884 -0.00334198
		 	mub_=1.18884;
		 	mue_= 0.103163;
		 	phi0_= 0.081786;
		 	v0_= 0.0554884;
		 	b0_= -0.00334198;
		 }

	}else if(parametrization=="iufsu" && parhyp=="su3"){
		// mub_	= 1.5;	//1.5 1.42001 0.0749067 0.100246 0.0771754 -0.00261042
		// mue_=0.13;			//0.1								
	 	// phi0_=0.09; 		// 0.1 low densities		.2							
	 	// v0_	=0.08;  	// 0.1 low densities 			.15			
	 	// b0_	= -0.003;	// -1e-5 low dens				
		mub_=1.21859;
		mue_= 0.135118;
		phi0_= 0.0790006;
		v0_= 0.0552466;
		b0_= -0.00307039;
	}else if(parametrization=="fsu2h"){
		// mub_=1.6868;
		// mue_= 0.141791;
		// phi0_= 0.0943894;
		// v0_= 0.0926038;
		// b0_= -0.00107964;		
		// mub_=1.6868;
		// mue_= 0.141791;
		// phi0_= 0.0943894;
		// v0_= 0.0926038;
		// b0_= -0.00107964;		
		mub_=1.41062;
		mue_=0.0811537;
		phi0_= 0.0945031;
		v0_= 0.077437;
		b0_= -0.00167563;
		if(doDeltas){
			// mub_=1.41062; //no amm
			// mue_=0.0811537;
			// phi0_= 0.0945031;
			// v0_= 0.077437;
			// b0_= -0.00167563;
			mub_=1.31723; // with amm, dens=3.5 rh0
			mue_=0.105105;
			phi0_=0.0878622;
			v0_=0.068489 ;
			b0_=-0.00188802;
			// mub_=1.41299;
			// mue_= 0.0906757;
			// phi0_= 0.0897645;
			// v0_= 0.0749307;
			// b0_= -0.00204853;
		}
	}else if (parametrization=="nl3"){
		mub_=1.35417; //1.48564 0.0434809 0.0972714 0.088873 -0.00439729
		mue_= 0.153066;
		phi0_= 0.0876225;
		v0_= 0.0723747;
		b0_= -0.0116196;
	}else if(parametrization=="nl3wr"){
	 
		mub_=	1.36911;
		mue_=0.0772368;
		phi0_=0.0902372; 
		v0_=0.0768166;
		b0_=-0.00231245;
		// mub_=1.50497;
		// mue_= 0.0242015;
		// phi0_= 0.0956183;
		// v0_= 0.0901114;
		// b0_= -0.000437992;
	}else if(parametrization=="l3wr"){
		// mub_=	1.29247; 
		// mue_=0.165258;
		// phi0_=0.0841176; 
		// v0_=0.0724188;
		// b0_=-0.00737542;   
		mub_=1.57618;
		mue_=0.107493;
		phi0_=0.108622;
		v0_=0.110745;
		b0_=-0.00511695;
		// mub_=1.89686;
		// mue_= 0.162052;
		// phi0_=0.102812;
		// v0_= 0.131473;
		// b0_=-0.00427271;
		}else if(parametrization=="nl3wr*"){
			mub_=2.56661;
			mue_= 0.12955;
			phi0_= 0.0990648;
			v0_= 0.161139;
			b0_= -0.000221779	;
		}else{		// mub_	= 1.31921;	//1.1.31921 0.228353 0.141197 0.0810486 -0.0103867
		// mue_=0.228353;			//0.1								
	 	// phi0_=0.141197; 		// 0.1 low densities		.2							
	 	// v0_	=0.0810486;  	// 0.1 low densities 			.15			
	 	// b0_	= -0.0108451;	// -1e-5 low dens				
		mub_	= 1.56849;	//1.15
		mue_=0.178599;			//0.1								
	 	phi0_=0.11565; 		// 0.1 low densities		.2							
	 	v0_	=0.100307;  	// 0.1 low densities 			.15			
	 	b0_	= -0.0120915;	// -1e-5 low dens				
	}
}


//=============== Calculate nucleon potential: ===============
std::vector<double> nlwm_class::getNucleonPotential(){
	double gs_= dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) : gb;
	
	double Un= gv_*V0 + gb_*neutron.I3*b0 - gs_*phi0;
	double Up= gv_*V0 + gb_*proton.I3*b0 	- gs_*phi0;

	return {Un, Up};
}

//=============== Calculate hyperon potential: ===============
std::vector<double> nlwm_class::getHyperonPotential(){
	double gs_= dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) : gb;

	double Ul0 = gv_*(xvl*V0+xtl*theta0) + gb_*xbl*lambda0.I3*b0		- gs_*xsl*phi0;
	double Usm = gv_*(xvs*V0+xts*theta0) + gb_*xbs*sigmam.I3*b0		 	- gs_*xss*phi0;
	double Us0 = gv_*(xvs*V0+xts*theta0) + gb_*xbs*sigma0.I3*b0		 	- gs_*xss*phi0;
	double Usp = gv_*(xvs*V0+xts*theta0) + gb_*xbs*sigmap.I3*b0		 	- gs_*xss*phi0;
	double Uxm = gv_*(xvx*V0+xtx*theta0) + gb_*xbx*xim.I3*b0		 		- gs_*xsx*phi0;
	double Ux0 = gv_*(xvx*V0+xtx*theta0) + gb_*xbx*xi0.I3*b0		 		- gs_*xsx*phi0;

	return{Ul0, Usm, Us0, Usp, Uxm, Ux0};
}

//=============== Calculate delta potential: ===============
std::vector<double> nlwm_class::getDeltaPotential(){

	double gs_= dogRhob ? gs*getCoupling_sigma(rhoB)	: gs;
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB)	: gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB)		: gb;

	double Udm= gv_*xvd*V0 + gb_*xbd*deltam.I3*b0 	- gs_*xsd*phi0;
	double Ud0= gv_*xvd*V0 + gb_*xbd*delta0.I3*b0 	- gs_*xsd*phi0;
	double Udp= gv_*xvd*V0 + gb_*xbd*deltap.I3*b0 	- gs_*xsd*phi0;
	double Udpp=gv_*xvd*V0 + gb_*xbd*deltapp.I3*b0 	- gs_*xsd*phi0;

	return {Udm, Ud0, Udp, Udpp};
}

//=============== Density dependent coupling for sigma  meson: ===============
double nlwm_class::getCoupling_sigma(double rhob_){
		double x_= rhob_/rho0;
	
	double ed2_=pow(x_+ds, 2.);
	return as*(1.+bs*ed2_)/(1.+cs*ed2_);
}

//=============== Density dependent coupling for sigma  meson: ===============
double nlwm_class::getCoupling_omega(double rhob_){
	double x_= rhob_/rho0;

	double ed2_=pow(x_+dv, 2.);
	return av*(1.+bv*ed2_)/(1.+cv*ed2_);
}

//=============== Density dependent coupling for rho meson: ===============
double nlwm_class::getCoupling_rho(double rhob_){
	double x_= rhob_/rho0;

	return exp(-ab*(x_-1.));
}

//=============== Derivative of density dependent coupling for sigma  meson: ===============
double nlwm_class::getDerivativeCoupling_sigma(double rhob_){
	double x_= rhob_/rho0;

	double ed2_=pow(x_+ds, 2.);
	return (2.*as/rho0)*(bs-cs)*(x_+ds)/pow((1.+cs*ed2_), 2.);
}

//=============== Derivative of density dependent coupling for sigma  meson: ===============
double nlwm_class::getDerivativeCoupling_omega(double rhob_){
	double x_= rhob_/rho0;

	double ed2_=pow(x_+dv, 2.);
	return (2.*av/rho0)*(bv-cv)*(x_+dv)/pow((1.+cv*ed2_), 2.);
}

//=============== Derivative of density dependent coupling for rho meson: ===============
double nlwm_class::getDerivativeCoupling_rho(double rhob_){
	double x_= rhob_/rho0;

	return -ab*exp(-ab*(x_-1.))/rho0;
}

//=============== rearrangement term due to density dependent coupling for rho meson: ===============
double nlwm_class::getRearrangementEnergy(void){

	double sigmar_=0.;
	double sigmar_s_=0., sigmar_v_=0., sigmar_b_=0., sigmar_t_=0.;
	
	if(dogRhob){
		sigmar_s_=getDerivativeCoupling_sigma(rhoB)*gs*phi0*(proton.condensate+neutron.condensate);
		sigmar_v_=getDerivativeCoupling_omega(rhoB)*gv*V0*(proton.density+neutron.density);
		sigmar_b_=getDerivativeCoupling_rho(rhoB)*gb*b0*(proton.I3*proton.density+
																											neutron.I3*neutron.density);

		if(doHyperons){
			sigmar_s_+=getDerivativeCoupling_sigma(rhoB)*gs*phi0*(xsl*lambda0.condensate
						+xss*(sigmap.condensate+sigma0.condensate+sigmam.condensate) 
						+xsx*(xi0.condensate+xim.condensate) );

			sigmar_v_+=getDerivativeCoupling_omega(rhoB)*gv*V0*(xvl*lambda0.density
						+xvs*(sigmap.density+sigma0.density+sigmam.density) 
						+xvx*(xi0.density+xim.density) );

			sigmar_b_+=getDerivativeCoupling_rho(rhoB)*gb*b0*(xbl*lambda0.I3*lambda0.density
					+xbs*(sigmap.I3*sigmap.density+sigma0.I3*sigma0.density+sigmam.I3*sigmam.density) 
					+xbx*(xi0.I3*xi0.density+xim.I3*xim.density) );

			sigmar_t_=getDerivativeCoupling_omega(rhoB)*gv*theta0*(xtl*lambda0.density
									+xts*(sigmap.density+sigma0.density+sigmam.density) 
									+xtx*(xi0.density+xim.density) );																									
		}																							
		if(doDeltas){
			sigmar_s_+=getDerivativeCoupling_sigma(rhoB)*gs*phi0*xsd*(deltapp.condensate
						+deltap.condensate+delta0.condensate+deltam.condensate);

			sigmar_v_+=getDerivativeCoupling_omega(rhoB)*gv*V0*xvd*(deltapp.density
						+deltap.density+delta0.density+deltam.density);

			sigmar_b_+=getDerivativeCoupling_rho(rhoB)*gb*b0*xbd*(deltapp.I3*deltapp.density
					+deltap.I3*deltap.density+delta0.I3*delta0.density+deltam.I3*deltam.density);
		}	

		sigmar_=sigmar_v_+ sigmar_b_ + sigmar_t_ - sigmar_s_;
	}

	return sigmar_;
}


//=============== Calculate total baryon density: ===============
double nlwm_class::getBaryonDens(){
	double densb_=	proton.density + neutron.density;
	
	if(doHyperons){densb_+= lambda0.density + sigmap.density + sigma0.density + sigmam.density 
												+ xi0.density + xim.density;
	}
	if(doDeltas){ densb_+= deltapp.density + deltap.density + delta0.density + deltam.density;}

	return densb_;
}

//=============== Calculate total isospin density: ===============
double nlwm_class::getIsoDens(){
	double dens3_=	proton.I3*proton.density + neutron.I3*neutron.density;
	
	if(doHyperons){
		dens3_+= lambda0.I3*lambda0.density
					+ sigmap.I3*sigmap.density + sigma0.I3*sigma0.density + sigmam.I3*sigmam.density 
					+ xi0.I3*xi0.density + xim.I3*xim.density;
	}
	if(doDeltas){ 
		dens3_+= deltapp.I3*deltapp.density + deltap.I3*deltap.density 
					+ delta0.I3*delta0.density + deltam.I3*deltam.density;
	}

	return dens3_;

}

//=============== Calculate total charge density: ===============
double nlwm_class::getChargeDens(){
	double charge_=	proton.Qdens + neutron.Qdens;
	
	if(doHyperons){charge_+= lambda0.Qdens + sigmap.Qdens + sigma0.Qdens + sigmam.Qdens 
												+ xi0.Qdens + xim.Qdens;
	}
	if(doDeltas){ charge_+= deltapp.Qdens + deltap.Qdens + delta0.Qdens + deltam.Qdens;}

	return charge_;
}


//=============== Calculate condensate density weighted by x_{sb} ratios: ===============
double nlwm_class::getSigmaEffDens(){
	double cond_=	proton.condensate + neutron.condensate;
	
	if(doHyperons){cond_+=	xsl*lambda0.condensate + 
												+	xss*(sigmap.condensate + sigma0.condensate + sigmam.condensate)
												+ xsx*(xi0.condensate + xim.condensate);
	}
	if(doDeltas){ cond_+= xsd*(deltapp.condensate + deltap.condensate 
													 + delta0.condensate + deltam.condensate);												 
	}
	return cond_;
}

//=============== Calculate baryon density weighted by x_{vb} ratios: ===============
double nlwm_class::getOmegaEffDens(){
	double densv_=	proton.density + neutron.density;
	
	if(doHyperons){
		densv_+= 	xvl*lambda0.density 
						+ xvs*(sigmap.density + sigma0.density + sigmam.density )
						+ xvx*(xi0.density + xim.density);
	}
	if(doDeltas){ 
		densv_+= xvd*(deltapp.density + deltap.density + delta0.density + deltam.density);
	}

	return densv_;
}

//=============== Calculate isospin density weighted by x_{rb} ratios: ===============
double nlwm_class::getIsoEffDens(){
	double dens3_=	proton.I3*proton.density + neutron.I3*neutron.density;
	
	if(doHyperons){
		dens3_+= xbl*lambda0.I3*lambda0.density
					+ xbs*(sigmap.I3*sigmap.density + sigma0.I3*sigma0.density + sigmam.I3*sigmam.density)
					+ xbx*(xi0.I3*xi0.density + xim.I3*xim.density);
	}
	if(doDeltas){ 
		dens3_+= xbd*(deltapp.I3*deltapp.density + deltap.I3*deltap.density 
					+ delta0.I3*delta0.density + deltam.I3*deltam.density);
	}

	return dens3_;

}

//=============== Calculate strange meson density weighted by x_{sb} ratios: ===============
double nlwm_class::getThetaEffDens(){
	double denst_=0.;
	
	if(doHyperons){
		denst_= 	xtl*lambda0.density 
						+ xts*(sigmap.density + sigma0.density + sigmam.density )
						+ xtx*(xi0.density + xim.density);
	}
	return denst_;
}


//=============== Calculate residue of scalar meson eom: ===============
double nlwm_class::sigmaMeson_eom_residue(double rhoS_){
	double gs_= dogRhob ? gs*getCoupling_sigma(rhoB) : gs;
	
  return gs_*rhoS_- ( pow(Ms, 2.)*phi0 +gs3*pow(phi0, 2.)/2. + gs4*pow(phi0, 3.)/6. );
}

//=============== Calculate residue of isoscalar-vector meson eom: ===============
double nlwm_class::omegaMeson_eom_residue(double rhoB_){
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) : gb;
  
	return gv_*rhoB_ -(V0*pow(Mv, 2.) 
								+xsi*pow(gv_, 4.)*pow(V0, 3.)/6.
                + 2.*Lv*V0*pow(gv_*gb_*b0, 2.) );
}

//=============== Calculate residue of isovector-vector meson eom: ===============
double nlwm_class::rhoMeson_eom_residue(double rho3_){
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	double gb_= dogRhob ? gb*getCoupling_rho(rhoB) 	 : gb;

  return   (gb_*rho3_) -b0*(pow(Mb, 2.)+2.*Lv*pow(gv_*gb_*V0, 2.));
}

//=============== Calculate residue of strange-scalar meson eom: ===============
double nlwm_class::thetaMeson_eom_residue(double rhoT_){
	double gv_= dogRhob ? gv*getCoupling_omega(rhoB) : gv;
	
	return gv_*rhoT_ - pow(Mt, 2.)*theta0;
}


//=============== Calculate total baryonic energy: ===============
double nlwm_class::getEnergy(void){
  double ener= proton.energy + neutron.energy ;
              
	double enerMeson= pow(Ms*phi0, 2)/2. + gs3*pow(phi0, 3.)/6.+ gs4*pow(phi0, 4)/24.
              +pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/8.+ pow(Mt*theta0, 2)/2.
              +pow(Mb*b0, 2)/2. 	+ 3.*Lv*pow(gv*gb*V0*b0, 2);//+ pow(Bfield, 2.)/2.;

	if(doHyperons){
		ener+= lambda0.energy+ sigmap.energy + sigma0.energy + sigmam.energy + xi0.energy + xim.energy;
	}
	if(doDeltas){
		ener+= deltapp.energy + deltap.energy + delta0.energy + deltam.energy;
	}

	return ener+ enerMeson;
}


//=============== Calculate total baryonic pressure: ===============
double nlwm_class::getPressure(void){
 
  // double press= proton.pressure+ neutron.pressure ;
	// if(dogRhob) press+=getRearrangementEnergy()*rhoB;

	// press+=
							// -pow(Ms*phi0, 2)/2. - gs3*pow(phi0, 3)/6.- gs4*pow(phi0, 4)/24.
          		// +pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/24. + pow(Mt*theta0, 2)/2.
          		// +pow(Mb*b0, 2)/2. 	+ Lv*pow(gv*gb*V0*b0, 2);

	double press= proton.chemPot*proton.density+ neutron.chemPot*neutron.density 
 					- getEnergy();
	if(temperature>Tmin_integration){press+=temperature*getEntropy();}
	
	if(doHyperons){
		press+= lambda0.chemPot*lambda0.density + sigmap.chemPot*sigmap.density 
					+ sigma0.chemPot*sigma0.density + sigmam.chemPot*sigmam.density 
					+ xi0.chemPot*xi0.density+ xim.chemPot*xim.density;
		// press+= lambda0.pressure + sigmap.pressure 
					// + sigma0.pressure + sigmam.pressure 
					// + xi0.pressure+ xim.pressure;
		}
	if(doDeltas){
		press+= deltapp.chemPot*deltapp.density + deltap.chemPot*deltap.density 
					+ delta0.chemPot*delta0.density + deltam.chemPot*deltam.density;
		// press+= deltapp.pressure + deltap.pressure 
					// + delta0.pressure + deltam.pressure;

	}
	return press;
}


//=============== Calculate baryonic pressure parallel to magnetic field: ===============
double nlwm_class::getPressureParallel(void){
	 double Lmeson = - pow(Ms*phi0, 2)/2. - gs3*pow(phi0, 3.)/6.- gs4*pow(phi0, 4)/24.
              +pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/24.+ pow(Mt*theta0, 2)/2.
              +pow(Mb*b0, 2)/2. 	+ Lv*pow(gv*gb*V0*b0, 2);

	double pressp_= proton.pressureParallel + neutron.pressureParallel;

	if(doHyperons){
		pressp_+= lambda0.pressureParallel+ sigmap.pressureParallel 
					+ sigma0.pressureParallel + sigmam.pressureParallel 
					+ xi0.pressureParallel + xim.pressureParallel;
	}
	if(doDeltas){
		pressp_+= deltapp.pressureParallel + deltap.pressureParallel 
					+ delta0.pressureParallel + deltam.pressureParallel;
	}

	return pressp_  + Lmeson;
}

//=============== Calculate baryonic pressure trasnsverse to magnetic field: ===============
double nlwm_class::getPressureTransverse(void){
	 double Lmeson = - pow(Ms*phi0, 2)/2. - gs3*pow(phi0, 3.)/6.- gs4*pow(phi0, 4)/24.
              +pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/24.+ pow(Mt*theta0, 2)/2.
              +pow(Mb*b0, 2)/2. 	+ Lv*pow(gv*gb*V0*b0, 2);

	double presst_= proton.pressureTransverse + neutron.pressureTransverse;
	
	 	if(doHyperons){
		presst_+= lambda0.pressureTransverse+ sigmap.pressureTransverse 
					+ sigma0.pressureTransverse + sigmam.pressureTransverse 
					+ xi0.pressureTransverse + xim.pressureTransverse;
	}
	if(doDeltas){
		presst_+= deltapp.pressureTransverse + deltap.pressureTransverse 
					+ delta0.pressureTransverse + deltam.pressureTransverse;
	}
	//double press_meson= -(+pow(Ms*phi0, 2)/2. + gs3*pow(phi0, 3.)/6.+ gs4*pow(phi0, 4)/24.
              //+pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/24.+ pow(Mt*theta0, 2)/2.
              //+pow(Mb*b0, 2)/2. 	+ 3.*Lv*pow(gv*gb*V0*b0, 2));
	return presst_  + Lmeson;

}

double nlwm_class::getMagnetization(void){
	// double emeson = pow(Ms*phi0, 2)/2. + gs3*pow(phi0, 3.)/6.+ gs4*pow(phi0, 4)/24.
              // +pow(Mv*V0, 2)/2.   + xsi*pow(gv*V0, 4)/8.+ pow(Mt*theta0, 2)/2.
              // +pow(Mb*b0, 2)/2. 	+ 3.*Lv*pow(gv*gb*V0*b0, 2);

	double mag_= proton.magnetization + neutron.magnetization;
	
	if(doHyperons){
		mag_+= lambda0.magnetization+ sigmap.magnetization 
				+ sigma0.magnetization + sigmam.magnetization 
				+ xi0.magnetization + xim.magnetization;
	}
	if(doDeltas){
		mag_+= deltapp.magnetization + deltap.magnetization 
				+ delta0.magnetization + deltam.magnetization;
	}
	return mag_ ;

}


//=============== Calculate total baryonic entropy: ===============
double nlwm_class::getEntropy(void){
  double entrp= proton.entropy + neutron.entropy;
	if(doHyperons){
		entrp+= lambda0.entropy + sigmap.entropy + sigma0.entropy + sigmam.entropy 
						+ xi0.entropy + xim.entropy;
	}
	if(doDeltas){
		entrp+= deltapp.entropy + deltap.entropy + delta0.entropy + deltam.entropy;
	}
	return entrp;
}

