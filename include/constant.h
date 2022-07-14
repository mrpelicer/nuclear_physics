#ifndef CONSTANT_H
#define CONSTANT_H

#include <iostream>
#include <cmath>
using namespace std;

//Set global constants:
const double pi2=M_PI*M_PI;                  //pi squared

//physical parameters
const double hc=197.326;                    // 1 fm^-1 = 197.326 MeV for hc=1
const double alphaEM= 1./137.;              //= 1.44 MeV fm/hc
const double Mnucleon= 939.;                // nucleon mass in MeV
const double Me=0.511;                      //electron mass in MeV
const double Mm=105.66;                     //muon mass in MeV
const double eHL = sqrt(4.*M_PI*alphaEM);   //Heaviside-Lorentz units
const double eGS = sqrt(alphaEM);           //Gaussian units
const double kBoltz   =  1.381e-23;         // J/K

//unit conversor  
const double MeVto_Sec = 6.582e-22;      
const double MeVto_Cm  = 1.973e-15;
const double JouletoErg=1e7;
const double MeVtoJoule= 1.60218e-13;
const double MeVdivfm3_to_gdivcm3		= 1.7827e12;
const double MeVdivfm3_to_dyndivcm2	= 1.6022e33;
const double Gauss_to_Mev2_HL = 1.95e-14;
const double Gauss_to_Mev2_GS = 6.91e-14;
const double _fm3_to_g_cm3 = 1.687e15;
// 
const double Tmin_integration=0.7/Mnucleon ;  // minimum temperature for fermi dirac integration (in MeV)
const double hdif= 5e-3;

#endif
