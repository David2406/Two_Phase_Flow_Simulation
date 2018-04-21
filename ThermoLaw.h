#ifndef THERMO
#define THERMO
#include <Eigen/Dense>
#include <string>
#include <math.h>
#include <iostream>
#include "Param.h"

using namespace std;
using Eigen::VectorXd;

class ThermoLaw
{

    public:

        //Type of the thermodynamics
        string ThermoLawType1_; 
        string ThermoLawType2_; 
        // Specific heat ratio
        double Gamma1_;
        double Gamma2_;
        // Stiffened Gas Reference Pressure
        double PiSG1_;  
        double PiSG2_;  
        //Isentropic EOS Constant Entropy
        double kBAR1_;  
        double kBAR2_;  

        //constructor:
        ThermoLaw(\
                string ThermoLawType1, string ThermoLawType2,\
                double Gamma1, double Gamma2,\
                double PiSG1 , double PiSG2,\
                double kBAR1 , double kBAR2\
                );

        //constructor by copy:
        ThermoLaw(ThermoLaw& Thermo);

        //constructor by default:
        ThermoLaw(); 

        //methods:
        string Get_ThermoLawType();
        double Get_Gamma();
        double Get_PiSG();
        double Get_kBAR();

};

//External methods:

//EOS of the specific internal energy function of rho and p, external to the class ThermoLaw, only taking the string parameter ThermoLawType
double Internal_Energy_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double rho, double p);

//EOS    mesh tab of the specific internal energy function of VectorXd rho and VectorXd p, external to the class ThermoLaw, only taking the string parameter ThermoLawType
VectorXd Internal_Energy_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p);

//Specific Entropy
double Entropy_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
        double rho, double p);

//Specific Entropy
VectorXd Entropy_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p);

//Sound Speed
double Sound_Speed_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double rho, double p); 

//EOS mesh tab of the sound speed function of VectorXd rho and VectorXd p
VectorXd Sound_Speed_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p); 

//EOS of pressure function of rho and e (internal energy per unit mass)
double Pressure_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    long double rho, long double e); 

//EOS mesh tab of the pressure function of VectorXd rho and VectorXd e
VectorXd Pressure_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& e); 

//EOS of density function of p and e (internal energy per unit mass)
double Density_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double e); 

//EOS mesh tab of the density function of VectorXd p and VectorXd e
VectorXd Density_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& p, VectorXd& e); 

/*************** CONVECTION SOURCE STABILITY: U-P RELAXATION ************/

//Gruneisen-like coefficient (in symmetrization variables p, s) 
double Gruneisen_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s); 

//Stiffness-like coefficient (in symmetrization variables p, s)
//(see note for more details)
double Stiffness_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s); 

//Sound Speed (in symmetrization variables p, s)
double Sound_Speed_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s); 

//Rho*Sound_Speed^2 (in symmetrization variables p, s)
double Rho_Sound_Speed_Squared_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s); 

//Density (in symmetrization variables p, s)
double Rho_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s); 

/*************** ISENTROPIC B-N RELAXATION ************/

double Rho_Sound_Speed_Squared_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double rho, double p); 

double Dp_Rho_Sound_Speed_Squared_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double rho, double p); 
#endif
