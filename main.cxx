#include <iostream>
#include "Mesh.h"
#include "Sol.h"
#include "ThermoLaw.h"
#include "Solver.h"
#include "Stability.h"
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

	/***** SIMULATION END *****/

    cout << "Hello World!\n"<<endl;

    return 0;

}


