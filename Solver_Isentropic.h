#ifndef SOLVER_ISENTROPIC
#define SOLVER_ISENTROPIC
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
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
using Eigen::Vector3d; //For the function Compute_k_Star_State
using Eigen::MatrixXd;
using Eigen::SparseLU; //For acoustic implicitation

//Defining a type for sparse matrix used for acoustic implicitation
typedef Eigen::SparseMatrix<double> SpMat; 

//Defining a type for matrix double 4 by 4
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;

//Defining a type for triplets used for acoustic implicitation
typedef Eigen::Triplet<double> T;

//: public Sol_Isen //Heritage: the Solver class heritates of the methods of the Sol 
	                 //class for the sake of simplicity
class Solver_Isen
{
    private:

    public:

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

        //Courant number for the convection step
        double CourantConv_;


        //Time of the Simulation
        double SimulationTime_;
        int print_freq_;

        //Numerical scheme for the conservative flux
        string SchemeTypeCons_;
        
        //Discretization for the non-conservative part of the flux
        string SchemeTypeNCons_;

        //constructor:
        Solver_Isen( double dtRelax,\
                int NRelax, double CourantBL,\
                string SchemeTypeCons, string SchemeTypeNCons,\
                double TimeStep, double CourantConv,\
                double SimulationTime, int print_freq\
                );

        //constructor by copy:
        Solver_Isen( Solver_Isen& solver); 

        //constructor by default:
        Solver_Isen();

        //Methods

        /************************************************/
        /***************  BOUNDARY LAYER  ***************/
        /************************************************/

        //Resolution of the time boundary layer related to the relaxation processes
        void BoundaryLayerUpdate(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
                double dtRelax, int NRelax\
                );

        //Update dtRelax_ using the relaxation cofactors and the Jacobian eigenvalues
        void BoundaryLayerTimeStepUpdate(Sol_Isen& sol, Mesh& mesh);

        //Update the solution from t to t+ NRelax_*dtRelax_
        void BoundaryLayer(Sol_Isen& sol, Mesh& mesh);

        /************************************************/
        /*****************  CONVECTION  *****************/
        /************************************************/

        //Update of the conservative flux
        void ConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh);
        //Returns a NcellExt*Nvariables matrix containing the non-conservative flux
        void NConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh);

        //Updates the conservative vectors using conservative and non-conservative fluxes
        void ConsVarUpdate(Sol_Isen& sol, Mesh& mesh);

        //Updates the TimeStep_ based the Jacobian eigenvalues
        void TimeStepUpdate(Sol_Isen& sol, Mesh& mesh);

        //Perform the Simulation of SimulationTime
        void Simulation(Sol_Isen& sol, Mesh& mesh,\
                string FileOutputFormat, string FileName\
                );

        //Save the solution 
        void Save_Sol(string FileOutputFormat, string FileName,\
                Sol_Isen& sol, Mesh& mesh, int ite);
};

//External Functions

/*************************************************/
/********** BOUNDARY LAYER RESOLUTION ************/
/*************************************************/

//Returns the solution of the ODE dt(V) = B_eq * V, with B_eq the linearized source term at equilibrium
Vector5d LinearizedSourceTermEqODE(\
        Vector5d& W_state_avr, Vector5d& W_state_ini,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time\
        );

//Returns the CFL constraint related to the eigenvalues of the LinearizedSourceTermEq
double LocalCourant_LSTEq(Vector5d& W_state_avr,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauMin, double etaP, double etaU,\
        double SpaceStep,\
        double Damp\
        );

//Returns the solution of the system of ODEs shown by Bereux and Sainsaulieu
void BoundaryLayerResolutionLocal(\
        Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin, double NRelax,\
        double SpaceStep\
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

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

//Returns the conservative flux according to the scheme choice SchemeTypeCons_
Vector5d ConsVarFluxUpdateLoc(\
        Vector5d& W_state_L, Vector5d&  W_state_R,\
        ThermoLaw& Therm,\
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
#endif
