#ifndef SOL_ISENTROPIC
#define SOL_ISENTROPIC
#include <Eigen/Dense>
#include "Mesh.h"
#include "ThermoLaw.h"
#include <string>
#include <stdlib.h> //to have exit EXIT_FAILURE

using namespace std;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::DiagonalMatrix;


//For the isentropic 7 eqs model
typedef Eigen::Matrix<double, 5, 1> Vector5d;
typedef Eigen::Matrix< double , 1, 5 > RowVector5d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;

class Sol_Isen
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

	//The type of theoretical variables, ex: ConsVar, NConsVar, EqRelaxVar    
	string VariableType_; 

	//Saving the initial variable values for the theoretical solution update
	double x_0_;

    //initial conditions L , R

    // alpha_L, rho1_L, u1_L, rho2_L, u2_L
    Vector5d InitL_;

    // alpha_R, rho1_R, u1_R, rho2_R, u2_R
    Vector5d InitR_;

    //Frozen Jacobian matrix, Frozen Rusanov radius
    //built with (InitR_ + InitL_)/2
    Matrix5d JacConvFrozen_;
    Vector5d DiracMassFrozen_;
    double uI_Frozen_;
    double   EigenvaluesFrozen_;
    Matrix5d TimeMat_;
    Vector5d W_eq_;
    Vector5d W_ref_;
    Matrix5d EigenVectorBasis_;
    Matrix5d EigenVectorBasisInv_;

    //Relaxation Times
    //tauRelax_(0) = pressure relaxation time
    //tauRelax_(1) = velocity relaxation time
    //tauRelax_(3) = mass     relaxation time
    Vector3d tauRelax_;

    //etaRelax_(0) = tau = min(tauP, tauU)
    //etaRelax_(1) = tau/tauP
    //etaRelax_(3) = tau/tauU
    Vector3d etaRelax_; 

    //reference pressure and density
    double pRef_;
    double mRef_;

	//Error in L1 norm for alpha1, p1, u1, p2, u2
	Vector5d errL1_;

	//Norm L1 of the exact solutions (needed to understand CV curves)
	Vector5d normL1_exact_;

    //The class Sol contains an object ThermoLaw which will allow to
    //Update its SoundSpeed_ vector at each iteration of the simulation
    //FOR BOTH PHASES

    ThermoLaw SolTherm_;    

    //Variables array				
	MatrixXd ConsVar_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0:  alpha1             values
	                        // column 1:  m1 = alpha1*rho1   values
	                        // column 2:  m1*u1              values
	                        // column 4:  m2 = alpha2*rho2   values
	                        // column 5:  m2*u2              values

	MatrixXd ConsVarOld_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0:  alpha1             values
	                        // column 1:  m1 = alpha1*rho1   values
	                        // column 2:  m1*u1              values
	                        // column 4:  m2 = alpha2*rho2   values
	                        // column 5:  m2*u2              values

	MatrixXd NConsVar_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0: alpha1 values
	                        // column 1: p1     values 
	                        // column 2: u1     values 
	                        // column 4: p2     values
	                        // column 5: u2     values

    //STAGGERED MESH
	MatrixXd ConsVarDual_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0:  alpha1             values
	                        // column 1:  m1 = alpha1*rho1   values
	                        // column 2:  m1*u1              values
	                        // column 4:  m2 = alpha2*rho2   values
	                        // column 5:  m2*u2              values

	MatrixXd NConsVarDual_;//Dynamic matrix of 'NcellExt' rows and 5 columns: 
	                        // column 0: alpha1 values
	                        // column 1: p1     values 
	                        // column 2: u1     values 
	                        // column 4: p2     values
	                        // column 5: u2     values

	MatrixXd NConsVarEqRelax_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0: alpha1 values
	                        // column 1: U      values sum (m_k/m) u_k
	                        // column 2: P      values sum alpha_k p_k
	                        // column 4: du     values
	                        // column 5: dp     values

	MatrixXd NConsVarEqRelaxDual_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0: alpha1 values
	                        // column 1: U      values sum (m_k/m) u_k
	                        // column 2: P      values sum alpha_k p_k
	                        // column 4: du     values
	                        // column 5: dp     values

	MatrixXd SoundSpeed_;//Dynamic Vector of 'NcellExt' rows and 2 columns:
	                        // column 0:  c1             values
	                        // column 1:  c2             values

	MatrixXd Mach_;//Dynamic Vector of 'NcellExt' rows and 2 columns:
	                        // column 0:  Mach1          values
	                        // column 1:  Mach2          values
				
    MatrixXd ConsVarFlux_;//Dynamic matrix of 'Nfaces' rows and 5 columns:
	                        // column 0:  0                       values
	                        // column 1:  m1*u1                   values
	                        // column 2:  m1*(u1**2) + alpha1*p1  values
	                        // column 4:  m2*u2                   values
	                        // column 5:  m*(u2**2) + alpha2*p2   values

    MatrixXd NConsVarFlux_;//Dynamic matrix of 'Nfaces' rows and 5 columns:
	                        // column 0:  u_I*div(alpha1)      values
	                        // column 1:  0                    values
	                        // column 2: -p_I*div(alpha1)      values
	                        // column 4:  0                    values
	                        // column 5: +p_I*div(alpha1)      values

    MatrixXd NConsVarFluxR_;//Dynamic matrix of 'Nfaces' rows and 5 columns:
	                        // column 0:  u_I*div(alpha1)      values
	                        // column 1:  0                    values
	                        // column 2: -p_I*div(alpha1)      values
	                        // column 4:  0                    values
	                        // column 5: +p_I*div(alpha1)      values

    MatrixXd NConsVarFluxL_;//Dynamic matrix of 'Nfaces' rows and 5 columns:
	                        // column 0:  u_I*div(alpha1)      values
	                        // column 1:  0                    values
	                        // column 2: -p_I*div(alpha1)      values
	                        // column 4:  0                    values
	                        // column 5: +p_I*div(alpha1)      values

    MatrixXd NConsVarFace_;//Dynamic matrix of 'Nfaces' rows and 1 column:
                                // column 0: alpha1 face values
                                // used for the non-conservative terms discretization

    MatrixXd SourceTerms_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0:  Kp*(p1-p2)                  values
	                        // column 1:  Kmu*(mu1-mu2)               values
	                        // column 2:  Ku*(u1-u2) + (u1+u2)/2*col1 values
	                        // column 1: -Kmu*(mu1-mu2)               values
	                        // column 2: -Ku*(u1-u2) + (u1+u2)/2*col1 values

    MatrixXd SourceTermsEqRelax_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0:  -Kp*dp                          values
	                        // column 1:  0                               values
	                        // column 2:  Kp*(dp-(C2-C1))*dp              values
	                        // column 1: -Ku*(m/(m1*m2))*du               values
	                        // column 2: -Kp*(C1/alpha_1 + C2/alpha_2)*dp values

	MatrixXd SolExact_; //Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0: alpha1_exact values
	                        // column 1: p1_exact     values
	                        // column 2: u1_exact     values
	                        // column 4: p2_exact     values
	                        // column 5: u2_exact     values

	MatrixXd SolExactEqRelax_;//Dynamic matrix of 'NcellExt' rows and 5 columns:
	                        // column 0: alpha1_exact values
	                        // column 1: U_exact      values sum (m_k/m) u_k
	                        // column 2: P_exact      values sum alpha_k p_k
	                        // column 4: du_exact     values
	                        // column 5: dp_exact     values
	//constructor:

	Sol_Isen(Mesh& M,\
		ThermoLaw& Therm,\
        double pRef, double mRef,\
        Vector3d tauRelax,\
        string VariableType,\
        Vector5d& InitL, Vector5d& InitR,\
		string SolType,\
		string LeftBCType, string RightBCType,\
        double x0\
		);

	//constructor by copy:
	Sol_Isen( Sol_Isen&); 

	//constructor by default:
	Sol_Isen();

	//methods:

    //Initialization 

    void ConsVarInit(Vector5d& InitR, Vector5d& InitL,\
            Vector3d tauRelax, double pRef, double mRef,\
            int NcellExt, int Ncells, double Length);

    void NConsVarInit(Vector5d& InitR, Vector5d& InitL,\
            Vector3d tauRelax, double pRef, double mRef,\
            int NcellExt, int Ncells, double Length);

	//Update methods:

    //Updates SoundSpeed_  of Sol using new values
    void SoundSpeed_Update();

    //Updates MachMax_ , MachMin_ using the new Sol values
	void Mach_Update();

    //Update NConsVar_ using the matrix ConsVar_
    void ConsVarToNConsVar(string MeshTag);

    //Update NConsVarEqRelax_ using the matrix NConsVar_
    void NConsVarToEqRelax(string MeshTag);

    //Update NConsVar_ using the matrix NConsVarEqRelax_
    void EqRelaxToNConsVar(string MeshTag);

    //Exact solutions update methods:
    void SolExact_Update(Mesh& mesh, double time);

    //Compute the L1 error and the L1 norms of the exact solutions
    void Compute_L1_err(Mesh& mesh);

    //Checking methods
    //Print the integral over the 1D domain of the conversative variable nVar
    void Check_Conservativity(int nVar, Mesh& mesh);
};

//External functions

// Functions in front the source terms
double Kp_Coeff(double p0, double alpha1);
double Ku_Coeff(double m0, double alpha1);
double Da1_Kp_Coeff(double p0, double alpha1);
double Da1_Ku_Coeff(double m0, double alpha1);

/************************************************/
/*************  CHANGE OF VARIABLE  *************/
/************************************************/

// Local change of variable W_state = [a1, p1, u1, p2, u2] --> U_state =  [a1, m1, m1*u1, m2, m2*u2]
Vector5d NConsVarToConsVarLoc(Vector5d& W_state,\
        ThermoLaw& Therm\
        );

// Local change of variable U_state = [a1, m1, m1u1, m2, m2u2] --> W_state =  [a1, p1, u1, p2, u2]
Vector5d ConsVarToNConsVarLoc(Vector5d& U_state,\
        ThermoLaw& Therm\
        );

// Local change of variable W_state = [a1, p1, u1, p2, u2] --> V_state =  [a1, U, P, du, dp]
Vector5d NConsVarToEqRelaxLoc(Vector5d& W_state,\
        ThermoLaw& Therm\
        );

// Local change of variable V_state = [a1, U, P, du, dp] --> W_state =  [a1, p1, u1, p2, u2]
Vector5d EqRelaxToNConsVarLoc(Vector5d& V_state,\
        ThermoLaw& Therm\
        );

//Returns the averaged state between W_state_L and W_state_R
Vector5d NConsVarAveraged(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        double weight_L);

/************************************************/
/**********  BOUNDARY LAYER RESOLUTION  *********/
/************************************************/

/*%%%%%%%%%%  All Variables Functions %%%%%%%%%%*/

//Returns the source term according to the choice of variables ConsVar, NConsVar, EqRelaxVar
Vector5d SourceTerm(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        );

//Returns the gradient of the source term according to the choice of variables ConsVar, NConsVar, EqRelaxVar
Matrix5d LinearizedSourceTerm(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        );

//Returns the jacobian matrix according to the choice of variables ConsVar, NConsVar, EqRelaxVar
Matrix5d LinearizedJacobian(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm\
        );

/*%%%%%%%%%%  EqRelax Variables %%%%%%%%%%*/

// Return the gradient of the SourceTermsEqRelax column, taken at Equilibrium du = dp = 0
Matrix5d LinearizedSourceTermsEq(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        );

// Return the gradient of the SourceTermsEqRelax column, contribution out of equilibrium du = dp = 0
Matrix5d LinearizedSourceTermsOutOfEq(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        );

// Return the Jacobian matrix related to the isentropic baer-nunziato flux expressed in EqRelax variables
Matrix5d LinearizedJacobianEqRelax(Vector5d& W_state,\
        ThermoLaw& Therm\
        );

/*%%%%%%%%%%  Conservative Variables %%%%%%%%%%*/

// Return the gradient of the SourceTerms column
Matrix5d LinearizedSourceTermsCons(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        );

// Return the Jacobian matrix related to the isentropic baer-nunziato flux expressed in ConsVar_ variables
Matrix5d LinearizedJacobianCons(Vector5d& W_state,\
        ThermoLaw& Therm\
        );

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/************************************************/
/**************  CONSERVATIVE FLUX  *************/
/************************************************/

// Return the conservative of the isentropic baer-nunziato system
Vector5d ConsVarFlux(Vector5d& W_state, ThermoLaw& Therm);

/************************************************/
/********  ANALYTICAL SOLUTION FUNCIONS  ********/
/************************************************/

//Returns the function related to the entropy jump condition accross the contact discontinuity u_I
double EntropyJumpFunction(double mL, double alpha1_in, double rho1, ThermoLaw& Therm);

//Returns the density for which the EntropyJumpFunction changes variations
double Rho1UpperBoundEntropyJumpFunction(double mL, double alpha1_in, ThermoLaw& Therm);

//Returns the density rho1_R involved in the entropy jump condition through the contact u_I
double EntropyJumpConditionDichotomy(\
        double& rho_inf, double& rho_sup, ThermoLaw& Therm,\
        Vector5d& W_state_L, double alpha1_R,\
        double& eps);

//Returns the velocity u1_R involved in the mass1 jump condition accros the contact discontinuity u_I
double VelocityJumpFunction(\
        double alpha1_L, double rho1_L, double u1_L, double u2_L,\
        double alpha1_R, double rho1_R\
        );

//Returns the velocity p2_R involved in the momentum mixture jump condition accros the contact discontinuity u_I
double PressureJumpFunction(\
        double mL, double alpha1_L, double p1_L, double u1_L, double p2_L,\
        double alpha1_R, double p1_R, double u1_R\
        );

//Provided a left state W_L and a right alpha1_R, solves the 4 equations system based on the jump conditions of the u_I-wave
Vector5d WstateContactResolution(Vector5d W_state_L, double alpha1_R,\
        ThermoLaw& Therm, double eps);

//Returns the specific volume for which the PressureJumpFunction changes variations
double TauUpperBoundPressureJumpFunction(int Phase_Id,\
        double Jshock,\
        ThermoLaw& Therm);

//Returns the function related to the R-H momentum jump relations
double  PressureRHJumpFunction(\
        int Phase_Id,\
        double Jshock,\
        double tau_R,\
        ThermoLaw& Therm);

//Returns the density rho1_L involved in the Rankine-hugoniot jump conditions for the Euler
//isentropic system
double RHJumpConditionDichotomy(\
        int& Phase_Id,\
        double& tau_inf, double& tau_sup, ThermoLaw& Therm,\
        Vector5d& W_state_R,\
        double& Jshock,\
        double& eps);

//Returns the velocity u1_L involved in the momentum RH relation
double VelocityRHJumpFunction(\
        double Jshock,\
        double tau_L,\
        double tau_R, double u_R\
        );

//Provided a right state W_R , solves the 2 equations system based on R-H
Vector5d WstateRHJumpCondition(\
        Vector5d& W_state_R,\
        double sigma1Shock, double sigma2Shock,\
        ThermoLaw& Therm, double eps);


//Returns the eigenvectors of the isentropic Baer Nunziato system
Vector5d IsentropicBN_Eigenvalues(\
                Vector5d& W_state_avr,\
                ThermoLaw& Therm\
                );

//Returns the eigenvectors of the isentropic Baer Nunziato system according to the 
//variable type
Matrix5d IsentropicBN_Eigenvectors(\
        string VariableType,\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        );

//Returns the eigenvectors of the isentropic Baer Nunziato system in conservative variables
Matrix5d IsentropicBN_EigenvectorsCons(\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        );

//Returns the inverse of the eigen vector matrix for the isentropic Baer Nunziato system in conservative variables
Matrix5d IsentropicBN_EigenvectorsInv(\
        string VariableType,\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        );

//Returns the inverse of the eigen vector matrix for the isentropic Baer Nunziato system in conservative variables
Matrix5d IsentropicBN_EigenvectorsInvCons(\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        );

//Returns the coordinates of X_R - X_L projected into the eigenvectors base
//according to the variable type
Vector5d IsentropicBN_EigenvectorsBaseProjection(\
        string VariableType,\
        Vector5d& W_state_avr,\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        );

//Returns the coordinates of U_R - U_L projected into the eigenvectors base
Vector5d IsentropicBN_EigenvectorsBaseProjectionCons(\
        Vector5d& W_state_avr,\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        );

//The Heaviside funtion
double Heaviside(double x, double time, double lambda);

//Problem linear-conv linear-relax dirac
double AkRelaxDirac(double x, double t,\
        double uI_Frozen, double lambda_k,\
        double ak, double ak_eq,\
        double dirac_weight, double tauvec);

//Problem linear-conv linear-relax dirac
double AkRelaxDiracTest(double x, double t,\
        double uI_Frozen, double lambda_k,\
        double dirac_weight, double tauvec);

//Returns the relaxation time matrix weighted by reference state (Cubic relaxation)
Matrix5d TimeRelaxMat(Vector5d& U_state_ref, Vector5d& TimeWeight);

//Returns the diagonal matrix of the vector R0^-1 * (U - Uinfty) taken to the square (Cubic relaxation)
Matrix5d QuadraticDistance(Vector5d& A_coords);

//Returns the inverse of the diagonal matrix of the vector R0^-1 * (U - Uinfty) taken to the square (Cubic relaxation)
Matrix5d QuadraticDistanceInv(Matrix5d& TimeMat, Vector5d& A_coords, double dtRelax);

//Returns the vector of the vector R0^-1 * (U - Uinfty) taken to the square (Cubic relaxation)
Vector5d CubicVectorDistance(Vector5d& A_coords);


#endif
