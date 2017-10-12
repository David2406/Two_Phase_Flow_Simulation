#ifndef SOL
#define SOL
#include <Eigen/Dense>
#include "Mesh.h"
#include "ThermoLaw.h"
#include <string>
#include <stdlib.h> //to have exit EXIT_FAILURE

using namespace std;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::MatrixXd;


//For the 7 eqs model
typedef Eigen::Matrix<double, 7, 1> Vector7d;

class Sol
{
    private:

    public:

	//String containing the world characterizing
	//the left boundary condition type, ex: 'transparent'	
	string LeftBCType_;

	//String containing the world characterizing
	//the right boundary condition type, ex: 'transparent'	
	string RightBCType_;

	//The type of theoretical solution, ex: PURE_CONTACT    
	string SolType_; 

//        //Number of conservative variables of the model
//	//int Nb_Variables_; 
//
//	int FieldsNumber_; //Number of stored flow fields (used for output treatment)
//
//    //Statistical fraction
//    int alpha_id1_;
//
//    //Phase 1 
//	int rho_id1_;
//	int u_id1_;
//	int p_id1_;
//	int c_id1_;
//	int M_id1_;
//	int internal_energy_id1_;
//	int s_id1_;
//
//    //Phase 2 
//	int rho_id2_;
//	int u_id2_;
//	int p_id2_;
//	int c_id2_;
//	int M_id2_;
//	int internal_energy_id2_;
//	int s_id2_;
//
//	//exact id
//
//    //Phase 1
//	int rhoe_id1_; 
//	int ue_id1_;
//	int pe_id1_;
//	int internal_energye_id1_;
//	int se_id1_;
//
//    //Phase 2
//	int rhoe_id2_; 
//	int ue_id2_;
//	int pe_id2_;
//	int internal_energye_id2_;
//	int se_id2_;

	//Saving the initial variable values for the theoretical solution update
	double x_0_;

    //initial conditions L , R

    // alpha_L, rho1_L, u1_L, p1_L, rho2_L, u2_L, p2_L
    Vector7d InitL_;

    // alpha_R, rho1_R, u1_R, p1_R, rho2_R, u2_R, p2_R
    Vector7d InitR_;

	//Error in L1 norm for alpha1, rho1, u1, p1, rho2, u2, p2
	Vector7d errL1_;
	//Norm L1 of the exact solutions (needed to understand CV curves)
	Vector7d normL1_exact_;

    //The class Sol contains an object ThermoLaw which will allow to
    //Update its SoundSpeed_ vector at each iteration of the simulation
    //FOR BOTH PHASES
    ThermoLaw SolTherm_;    

    //Variables array				
	MatrixXd ConsVar_;//Dynamic matrix of 'NcellExt' rows and 7 columns:
	                        // column 0:  alpha1             values
	                        // column 1:  m1 = alpha1*rho1   values
	                        // column 2:  m1*u1              values
	                        // column 3:  m1*e1              values
	                        // column 4:  m2 = alpha2*rho2   values
	                        // column 5:  m2*u2              values
	                        // column 6:  m2*e2              values
	MatrixXd NConsVar_;//Dynamic matrix of 'NcellExt' rows and 7 columns:
	                        // column 0: alpha1 values
	                        // column 1: rho1   values
	                        // column 2: u1     values
	                        // column 3: p1     values
	                        // column 4: rho2   values
	                        // column 5: u2     values
	                        // column 6: p2     values

	MatrixXd SoundSpeed_;//Dynamic Vector of 'NcellExt' rows and 2 columns:
	                        // column 0:  c1             values
	                        // column 1:  c2             values

	MatrixXd Mach_;//Dynamic Vector of 'NcellExt' rows and 2 columns:
	                        // column 0:  Mach1          values
	                        // column 1:  Mach2          values
				
    MatrixXd ConsVarFlux_;//Dynamic matrix of 'Nfaces' rows and 7 columns:
	                        // column 0:  0                       values
	                        // column 1:  m1*u1                   values
	                        // column 2:  m1*(u1**2) + alpha1*p1  values
	                        // column 3:  (m1*e1 + alpha1*p1)*u1  values
	                        // column 4:  m2*u2                   values
	                        // column 5:  m*(u2**2) + alpha2*p2   values
	                        // column 6:  (m2*e2 + alpha2*p2)*u2  values

    MatrixXd NConsVarFlux_;//Dynamic matrix of 'Nfaces' rows and 7 columns:
	                        // column 0:  u_I*div(alpha1)      values
	                        // column 1:  0                    values
	                        // column 2: -p_I*div(alpha1)      values
	                        // column 3: -p_I*u_I*div(alpha1)  values
	                        // column 4:  0                    values
	                        // column 5: +p_I*div(alpha1)      values
	                        // column 6: +p_I*u_I*div(alpha1)  values

    MatrixXd SourceTerms_;//Dynamic matrix of 'NcellExt' rows and 7 columns:
	                        // column 0:  tau_p*(p1-p2)                  values
	                        // column 1:  tau_mu*(mu1-mu2)               values
	                        // column 2:  tau_u*(u1-u2) + (u1*u2/2)*col0 values
	                        // column 3:  tau_T*(T1-T2) + u_I*...
	                        // column 1: -tau_mu*(mu1-mu2)               values
	                        // column 2: -tau_u*(u1-u2) + (u1*u2/2)*col0 values
	                        // column 3: -tau_T*(T1-T2) + u_I*...

	MatrixXd SolExact_; //Dynamic matrix of 'NcellExt' rows and 7 columns:
	                        // column 0: alpha1_exact values
	                        // column 1: rho1_exact   values
	                        // column 2: u1_exact     values
	                        // column 3: p1_exact     values
	                        // column 4: rho2_exact   values
	                        // column 5: u2_exact     values
	                        // column 6: p2_exact     values
	//constructor:

	Sol(Mesh& M,\
		ThermoLaw& Therm,\
        Vector7d& InitL, Vector7d& InitR,\
		string SolType,\
		string LeftBCType, string RightBCType,\
        double x0\
		);

	//constructor by copy:
	Sol( Sol&); 

	//constructor by default:
	Sol();

	//methods:

    //Initialization 

    void ConsVarInit(Vector7d& InitR, Vector7d& InitL,\
            int NcellExt, int Ncells, double Length);

    void NConsVarInit(Vector7d& InitR, Vector7d& InitL,\
            int NcellExt, int Ncells, double Length);

	//Update methods:
//    void SoundSpeed_Update();//Updates SoundSpeed_ vector of Sol using new values
//	void Mach_Update();//Updates MachMax_ , MachMin_ using the new Sol values

    //Once ConsVarFlux_ has been updated using the solver, this function updates the conservative variables matrix ConsVar_
//	void ConsVar_Update(  Mesh& mesh, double TimeStep, \
        string SchemeType\
        );

    //Once ConsVar_ has been updated, this function updates the non conservative variables matrix NConsVar_
//	void NConsVar_Update();

    //Theoretical solution update
//	void SolExact_Update(Mesh& mesh, double time); 

	//Error computation

    //compute the L1 errors of the variables rho, u, p
//	void Compute_L1_err(Mesh& mesh);  

};

//External functions

// returns the two values vector in order to initialize Riemann Problem
VectorXd Riemann_Field(Vector7d& val_L, Vector7d& val_R, double Ncells, double NcellExt, double Length, double x_0); 

double Euler_k_Resolvant(double& p, ThermoLaw& SolTherm, double& rho_k, double& u_k, double& p_k, double& c_k); //Return the local 'k' (k==L or k==R) resolvant function of p for euler equations

double Euler_Resolvant(double& p, string& SolType, ThermoLaw& SolTherm, double& rho_L, double& rho_R, double& u_L, double& u_R, double& p_L, double& p_R, double& c_L, double& c_R); //Return the resolvant function of p for euler equations

double P_Euler_Dichotomy(double& p_inf, double& p_sup, string& SolType, ThermoLaw& SolTherm, double& rho_L, double& rho_R, double& u_L, double& u_R, double& p_L, double& p_R, double& c_L, double& c_R, double& eps); // returns the p_star for euler equations

//Determines rho_star in the left or right intermediate areas
double Rho_k_star_Resolvant(double rho0_k, double p0_k, double p_star, ThermoLaw& SolTherm );

//Determines c_star in the left or right intermediate areas
double C_k_star_Resolvant(double c0_k, double p0_k, double p_star, ThermoLaw& SolTherm );

//Determines sigma_star the shock velocity in the left or right intermediate areas
double Sigma_k_star_Resolvant(double u0_k, double c0_k, \
	double p0_k, double p_star, \
	int sign, ThermoLaw& SolTherm );

//Compute admissible right density for an isolated 3-shock for Euler
double Rho_pure_shock(double rho0_L, double u0_L, double p0_L,\
	double sigma, ThermoLaw& SolTherm);

//Compute admissible right velocity for an isolated 3-shock for Euler
double U_pure_shock(double rho0_L, double u0_L, double p0_L,\
	double rho_R_shock, ThermoLaw& SolTherm);

//Compute admissible right pressure for an isolated 3-shock for Euler
double P_pure_shock(double rho0_L, double p0_L,\
	double rho_R_shock, ThermoLaw& SolTherm);

//Compute admissible right velocity for an isolated 1-rarefaction for Euler
double U_pure_rarefaction(double rho0_L, double u0_L, double p0_L,\
	double rho_R_rarefaction, ThermoLaw& SolTherm);

//Compute admissible right pressure for an isolated 1-rarefaction for Euler
double P_pure_rarefaction(double rho0_L, double p0_L,\
	double rho_R_rarefaction, ThermoLaw& SolTherm);

//Compute the update of Gamma for an isolated sinus impulse
double Gamma_Update_Impulse(double rho0_L, double u0_L, double p0_L, \
	double M0, double Cv, double T0, string ThermoType);

//Compute the update of PiSG for an isolated sinus impulse
double PiSG_Update_Impulse(double Gamma, double rho0_L, double u0_L, double p0_L,\
	double M0, string ThermoType);

Vector3d State_k_Fan(double& rho_k, double& u_k, double& p_k, double& c_k, ThermoLaw& SolTherm, int sign, double& x, double& time); //returns the solution of the fan for a k rarefaction wave

#endif
