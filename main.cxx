#include <iostream>
#include "Mesh.h"
#include "Sol.h"
#include "ThermoLaw.h"
//#include "Solver.h"
//#include "Stability.h"
#include <Eigen/SparseLU>
#include <limits> //file reading
#include <vector> //storing the different strings of Mach numbers

using namespace std;
using Eigen::MatrixXd;
using Eigen::Array;
using Eigen::VectorXd;

int main(int argc, char* argv[])
{

	cout.precision(COUT_PRECISION);
	cout<<std::scientific;

    //Mesh parameters
    double Length   = ONE;
    int Ncells      = 5;
    int NGhostcells = 2;

    Mesh  mesh_try(Length, Ncells, NGhostcells);

    //Thermodynamics parameters
    string ThermoLawType1 = "PG";
    string ThermoLawType2 = "SG";
    double Gamma1         = 1.4;
    double Gamma2         = 7.5;
    double PiSG1          = ZERO;
    double PiSG2          = 3e8;

    ThermoLaw therm_try(\
            ThermoLawType1, ThermoLawType2,\
            Gamma1, Gamma2,\
            PiSG1,  PiSG2\
            );

    //Solution parameters
    Vector7d InitL, InitR;
    InitL<< 8e-1,
            1e3,
            0.0,
            1e5,
            1e0,
            0.0,
            1e4;
    InitR<< 2e-1,
            9e2,
            0.0,
            8e4,
            9e-1,
            0.0,
            9e3;
    string SolType     = "Poulet";
    string LeftBCType  = "Transparent";
    string RightBCType = "Transparent";

    double x_0 = 5e-1;

    Sol sol_try(mesh_try,\
		therm_try,\
        InitL, InitR,\
		SolType,\
		LeftBCType, RightBCType,\
        x_0\
		);


	/***** SIMULATION END *****/

    cout << "Hello World!\n"<<endl;

    return 0;

}


