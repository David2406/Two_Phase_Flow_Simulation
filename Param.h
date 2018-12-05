#ifndef PARAM
#define PARAM
#include <math.h>
#include <complex> //For complex numbers manipulations
#include <iomanip> //For setprecision in the ofstreams

using namespace std;

//Index Parameters

extern int NoIndex;

//Cout Precision integer

extern int COUT_PRECISION;

//Numbers

extern double ZERO;
extern double ONE_OVER_TWO;
extern double ONE;
extern double TWO;
extern double THREE;
extern double FOUR;
extern double EIGHT;
extern double TEN;
extern double TWENTY;
extern double HUNDRED;
extern double THOUSAND;

extern double Pi;

//Complex Numbers

extern complex<double> I;

//Numerical Parameters

extern double eps_check_cst;
extern double epsZero;
extern double epsDicho;
extern double Big;
extern double Inf;

//Physical Parameters

//Involved in eps0sq Update: max(Mach_ref, min(MachMax_,1))
extern double Mach_ref;

//Stiffened Gas reference pressure
extern double PiSG;    


#endif
