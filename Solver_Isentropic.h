#ifndef SOLVER_ISENTROPIC
#define SOLVER_ISENTROPIC
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <vector>
#include "Mesh.h"
#include "ThermoLaw.h"
#include "Sol_Isentropic.h"
#include <string>
#include <iostream>
#include <stdio.h>
#include <time.h> //Measure the CPU time
#include <fstream>  //To write inside a file
#include <algorithm> //For the min/max functions

using namespace std;
using Eigen::Array;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::MatrixXcd;
using Eigen::ComplexEigenSolver;
using Eigen::Vector3d; //For the function Compute_k_Star_State
using Eigen::MatrixXd;
using Eigen::SparseLU; //For acoustic implicitation
using Eigen::DiagonalMatrix;

//Defining a type for sparse matrix used for acoustic implicitation
typedef Eigen::SparseMatrix<double> SpMat; 

//Defining a type for matrix double 4 by 4
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 2, 1> Vector2d;

//Defining a type for triplets used for acoustic implicitation
typedef Eigen::Triplet<double> T;

//: public Sol_Isen //Heritage: the Solver class heritates of the methods of the Sol 
	                 //class for the sake of simplicity
class Solver_Isen
{
    private:


    public:
 
        string FileOutputFormat_;
        bool FractionalStep_;
        string SourceTermType_;
        string TimeIntegrationType_;

        /************************************************/
        /***************  BOUNDARY LAYER  ***************/
        /************************************************/

        //Time-step related to the min of all the relaxation times
        double dtRelax_;

        //Number of time-step of duration dtRelax_ during the boundary layer resolution 
        int NRelax_;

        //Courant number for the time boundary layer resolution
        double CourantBL_;

        /************************************************/
        /*****************  CONVECTION  *****************/
        /************************************************/

        //Convective Time-step
        double TimeStep_;

        //Discrete time-steps at face face_id
        MatrixXd TimeMatrix_;

        //Courant number for the convection step
        double CourantConv_;

        //Time of the Simulation
        double SimulationTime_;
        double TimeElapsed_;
        int print_freq_;
        double CPUTime_;

        //Numerical scheme for the conservative flux
        string SchemeTypeCons_;
        
        //Discretization for the non-conservative part of the flux
        string SchemeTypeNCons_;

        //constructor:
        Solver_Isen(\
                bool FractionalStep,\
                string SourceTermType,\
                string TimeIntegrationType,\
                int NRelax, double CourantBL,\
                string SchemeTypeCons, string SchemeTypeNCons,\
                double CourantConv,\
                double SimulationTime, int print_freq,\
                string FileOutputFormat,\
                int Nfaces\
                );

        //constructor by copy:
        Solver_Isen( Solver_Isen& solver); 

        //constructor by default:
        Solver_Isen();

        //Methods

        /************************************************/
        /***************  BOUNDARY LAYER  ***************/
        /************************************************/

        /**********************/
        /* BEREUX-SAINSAULIEU */
        /**********************/

        //Resolution of the time boundary layer related to the relaxation processes
        void BoundaryLayerUpdate(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
                double& dtRelax, double& dtRelax_estim, string VariableType);

        //BN Resolution of the time boundary layer related to the relaxation processes
        void BN_BoundaryLayerUpdate(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
                double& dtRelax, double& dtRelax_estim, string VariableType);


        //Update dtRelax_ using the relaxation cofactors and the Jacobian eigenvalues
        void BoundaryLayerTimeStepUpdate(Sol_Isen& sol, Mesh& mesh);

        //Update the solution from t to t+ NRelax_*dtRelax_
        void BoundaryLayer(Sol_Isen& sol, Mesh& mesh);
        void BoundaryLayerTest(Sol_Isen& sol, Mesh& mesh, string Filename);

        /**********************/
        /*** FRACTIONAL-STEP **/
        /**********************/

        //Resolution of the time boundary layer related to the relaxation processes
        void BoundaryLayerUpdateFracStep(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
                double& dtRelax, double& dtRelax_estim, string VariableType);

        //Update the pure relaxation solution from t to t + *dtConv_
        void BoundaryLayerFracStep(Sol_Isen& sol, Mesh& mesh);

        /************************************************/
        /*****************  CONVECTION  *****************/
        /************************************************/

        //Update of the conservative flux
        void ConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh);
        //Returns a NcellExt*Nvariables matrix containing the non-conservative flux
        void NConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh);

        //Updates the conservative vectors using conservative and non-conservative fluxes
        void ConsVarUpdate(Sol_Isen& sol, Mesh& mesh);

        //Updates the conservative vectors using the projective method of P. Lafitte and T. Rey
        void ConsVarUpdateLafitteRey(Sol_Isen& sol, Mesh& mesh);

        //Updates the cell-colocated SourceTerm vector after the boundary layer resolution
        void SourceTermsUpdate(Sol_Isen& sol, Mesh& mesh);

        //Updates the TimeStep_ based the Jacobian eigenvalues
        void TimeStepUpdate(Sol_Isen& sol, Mesh& mesh);

        //Perform the Simulation of SimulationTime
        void Simulation(Sol_Isen& sol, Mesh& mesh,\
                string FileOutputFormat, string FileName\
                );

        //Save the solution 
        void Save_Sol(string FileOutputFormat, string FileName,\
                Sol_Isen& sol, Mesh& mesh, int ite);

        void Save_Sol_Local(string FileName,\
                Sol_Isen& sol, Mesh& mesh, double x_cell, double time);
};

//External Functions

/*************************************************/
/********** BOUNDARY LAYER RESOLUTION ************/
/*************************************************/

/*%%%%%%%%%%  EqRelax Variables %%%%%%%%%%*/

//Returns the solution of the ODE dt(V) = B_eq * V, with B_eq the linearized source term at equilibrium
Vector5d LinearizedSourceTermEqODE(\
        Vector5d& W_state_avr, Vector5d& W_state_ini,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time\
        );

//Returns the CFL constraint related to the eigenvalues of the LinearizedSourceTermEq
Vector2d LocalCourant_LSTEq(\
        Vector5d& W_state0, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauMin, double etaP, double etaU,\
        double SpaceStep,\
        double Damp,\
        int face_id\
        );

//Returns the solution of the system of ODEs shown by Bereux and Sainsaulieu
void BoundaryLayerResolutionLocal(\
        Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin, double NRelax,\
        double SpaceStep\
        );

void BoundaryLayerResolutionLocalTest(\
        string VariableType,\
        Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
        Matrix5d& LinJac_L, Matrix5d& LinJac_R,\
        Matrix5d& LinSouTerm, Vector5d& STerm,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin,\
        double SpaceStep\
        );

void BoundaryLayerResolutionLocalTestNewton(\
        string VariableType,\
        Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
        Matrix5d& LinJac_L, Matrix5d& LinJac_R,\
        Matrix5d& LinSouTermInv, Vector5d& STerm,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin,\
        double SpaceStep\
        );

//Returns the solution of the system of ODEs with the Source Term, Newton Method
void SourceTermResolutionLocalNewton(\
        Vector5d& H_state,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double CourantBL, double tauMin, double NRelax\
        );

//During the time evolution of the Bereux-Sainsaulieu ODE, save the solution H_state, 
//in the filename
void SaveBereuxSainsaulieu(Vector5d& H_state, string filename, double time);

/*%%%%%%%%%%  Conservative Variables %%%%%%%%%%*/

//Returns the solution of the system of ODEs shown by Bereux and Sainsaulieu
void BoundaryLayerResolutionLocalConsVar(\
        Vector5d& H_state,Vector5d& U_state_L, Vector5d& U_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin, double NRelax,\
        double SpaceStep\
        );

//Returns the non-conservative flux à la Ambroso Galié Chalons
Vector5d NConsFluxLoc(Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        );

/********** exact solution functions ************/

//Returns the solution of the boundary layer dynamics when A(jacobian) = c0*Id
Vector5d BoundaryLayerSolSingleCharact(\
        Vector5d& W_state_avr, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time,\
        double c0, double dx\
        );

//Returns the solution of the boundary layer dynamics when A(jacobian) = c0*Id + Perth
//Perth is related to: tau_pert*dx(P) + p_perth*dx(U)
Vector5d BoundaryLayerSolSingleCharact_Pert(\
        Vector5d& W_state_avr, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time,\
        double c0, double rho_pert, double p_pert, double dx\
        );

/*****************************************************************************/
/************** HIGH ORDER IMPLICIT TIME INTEGRATION METHOD ******************/
/*****************************************************************************/

//Fourth order time integration
void TimeIntegration(Sol_Isen& sol, double SimulationTime, double dtRelax,\
        string FileNameInput, double CourantBL\
        );

void TimeIntegrationLoc(Vector5d& U_state_ini,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Matrix5d& TauMat,\
        Vector5d& U_state_ref, Vector5d& U_state_eq,\
        double dtRelax, double tauMin);

//Rosenbrock of order 4 method
void RosenBrockFourthOrder(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini, Vector5d& JacVector_ini, Vector5d& One,\
        double& dtRelax_ini, double& dtRelax_end, double& dtRelax_estim,\
        int& TimeStepIteCounter\
        );

//Euler Implicit of order 1 method
void ImplicitEuler(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini, Vector5d& JacVector_ini, Vector5d& One,\
        double& dtRelax_ini, double& dtRelax_end, double& dtRelax_estim,\
        int& TimeStepIteCounter\
        );

//Explicit Runge-Kutta of order 4 method
void RungeKuttaFourthOrder(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini,\
        double& dtRelax_ini, double& dtRelax_end, double& dtRelax_estim,\
        int& TimeStepIteCounter\
        );

//Update the source term corresponding to a linear spring
void LinearSpring(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq, Vector5d& STerm,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv);

//Returns the diagonal vector of the source term jacobian related to a linear spring
Vector5d LinearSpringGradient(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq,\
        Matrix5d& EigenVectorBasisInv);

//Update the source term corresponding to a non-linear cubic spring
void NonLinearSpring(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq, Vector5d& STerm,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv);

//Derive the real source term of the isentropic BN system
void IsenBNSourceTerm(double tauU, double tauP, double pref, double rhoref,\
        Vector5d& U_state, Vector5d& STerm,\
        ThermoLaw& Therm
        );

//Derive the real source term gradient of the isentropic BN system
Matrix5d IsenBNSourceTermGradient(double tauU, double tauP, double pref, double rhoref,\
        Vector5d& U_state,\
        ThermoLaw& Therm
        );

//Returns the diagonal vector of the source term jacobian related to a non-linear cubic spring
Vector5d NonLinearSpringGradient(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq,\
        Matrix5d& EigenVectorBasisInv);

//Returns the exact solution at time t of the non-linear cubic relaxation
void NonLinearSpringExactSolution(Vector5d& U_state_exact, Vector5d& U_state_ini,\
        Vector5d& U_state_eq, Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& TimeRefVector, double time);

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

//Returns the conservative flux according to the scheme choice SchemeTypeCons_
Vector5d ConsVarFluxUpdateLoc(\
        Vector5d& W_state_L, Vector5d&  W_state_R,\
        ThermoLaw& Therm,\
        Matrix5d& JacConvFrozen,\
        double EigenvaluesFrozen,\
        string SchemeTypeCons\
        );

// Returns the spectral radius used as diffusion coefficient in the Rusanov scheme
double SpectralRadiusRusanov(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        );

//Returns the local timestep based on the Jacobian eigenvalues
double LocalCourant_Conv(Vector5d& W_state_L, Vector5d& W_state_R,\
                ThermoLaw& Therm,\
                double SpaceStep,\
                double Courant\
                );
//Returns the convergence curves, orders and efficiency curves of the simulation
void Convergence_Curve(\
        double Length, int NGhostCells, \
        string LeftBCType, string RightBCType, \
        ThermoLaw& therm,\
        string SolType,\
        string VariableType,\
        Vector5d& InitL, Vector5d& InitR,\
        int Nbr_Areas, double x0,\
        double SimulationTime, \
        bool CFL_ramp, int CFL_ramp_range, \
        double CourantBL, double CourantConv,\
        bool FractionalStep,\
        string SourceTermType,\
        string TimeIntegrationType,\
        int NRelax,\
        double pRef, double mRef,\
        Vector3d& tauRelax,\
        string SchemeTypeCons, string SchemeTypeNCons,\
        Array<int, Eigen::Dynamic, 1>& CellsTab, \
        string FileOutputFormat, string CV_curve, int print_freq);
#endif
