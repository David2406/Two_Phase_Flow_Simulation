#include <iostream>
#include "Mesh.h"
#include "Sol.h"
#include "Sol_Isentropic.h"
#include "Solver_Isentropic.h"
#include "ThermoLaw.h"
#include "Stability.h"
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

    /* FILE READING PARAMETERS */
    string file_name="BN_Simulation.input";
    char input_comment='%';
    char character;
    string key_word;
    string key_str_info;
    char delimiter('=');
    char end_of_line('\n');

    /* SIMULATION VARIABLES */

    string Simulation_Name;
    //Solution choice
    string SolType;

    //Thermodynamics choice
    string ThermoLawType1; 
    string ThermoLawType2; 
    double Gamma1;         
    double Gamma2;         
    double PiSG1;          
    double PiSG2;          
    double kBAR1;          
    double kBAR2; 

    //Time Relaxation Parameters
    int NbrTauRelax;
    Vector3d tauRelax; 
    double pRef;
    double mRef;

    bool SourceTerm;
    int NRelax;
    
    //Boundary Conditions 
    string LeftBCType;
    string RightBCType;

    //Initial Riemann Problem
    double alpha1_L, alpha1_R;
    double p1_L, p1_R;
    double u1_L, u1_R;
    double p2_L, p2_R;
    double u2_L, u2_R;

    //Solution parameters
    Vector5d InitL, InitR;

    //Scheme choice
    string SchemeTypeCons;
    string SchemeTypeNCons;
    
    //CFL choice
    bool CFL_ramp;
    int CFL_ramp_range;

    double CourantBL;
    double CourantConv;
    double TimeStep = Big;
    
    //Mesh
    double Length;
    double TranslationLength;
    int NGhostCells;
    int nbNcells_cases;
    Array<int, Eigen::Dynamic, 1> CellsTab;
    
    //Double Riemann Problem
    int Nbr_Areas;
    double x0;

    //We print every print_freq iteration
    int print_freq;
    string FileOutputFormat;

    ifstream fichier( file_name.c_str(), ios::in);  // on ouvre le fichier en lecture

    if(fichier)  // si l'ouverture a réussi
    { 
        int nblines(0);
        string line;
        while( getline(fichier, line) ) {

            nblines++;

        }     
        fichier.clear();
        //Going to octet zero from the begining
        fichier.seekg(0, ios::beg);

        cout<<endl;
        cout<<"%%%%%%%%%%%%%%%%"<<endl;
        cout<<"File Parameters: "<<endl;
        cout<<"%%%%%%%%%%%%%%%%"<<endl<<endl;

        for(int i=1; i<=nblines; i++){

            fichier.get(character);
            if(character==input_comment){

                fichier.ignore(numeric_limits<int>::max(), '\n');
            }
            else{

                //the first get has moved the cursor of 1
                //octet, then we have to replace it
                fichier.seekg(-1, ios::cur);
                getline(fichier,key_word,delimiter);

                //Analysis of the different key words:
                if(key_word=="CASE_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    Simulation_Name+=key_str_info;
                    //cout<<"Simulation_Name= "<<Simulation_Name<<endl;

                }
                else if(key_word=="SOLUTION_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    SolType=key_str_info;
                    cout<<"Sol Type: "<<SolType<<endl;

                }
                else if(key_word=="THERMODYNAMICS1"){

                    getline(fichier,key_str_info, ' ');
                    ThermoLawType1=key_str_info;
                    //cout<<"ThermoLawType: "<<ThermoLawType<<endl;
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            ' ');
                    Gamma1=stod(key_str_info);
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            ' ');
                    PiSG1=stod(key_str_info);
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            end_of_line);
                    kBAR1=stod(key_str_info);

                    Simulation_Name+="_"+ThermoLawType1;
                    cout<<"Thermodynamics1: "<<ThermoLawType1<<", "<<"Gamma1= "<<Gamma1<<\
                        ", PiSG1=  "<<PiSG1<<", kBAR1=  "<<kBAR1<<endl;
                }
                else if(key_word=="THERMODYNAMICS2"){

                    getline(fichier,key_str_info, ' ');
                    ThermoLawType2=key_str_info;
                    //cout<<"ThermoLawType: "<<ThermoLawType<<endl;
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            ' ');
                    Gamma2=stod(key_str_info);
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            ' ');
                    PiSG2=stod(key_str_info);
                    fichier.ignore(numeric_limits<int>::max(),
                            delimiter);
                    getline(fichier,key_str_info,
                            end_of_line);
                    kBAR2=stod(key_str_info);

                    cout<<"Thermodynamics2: "<<ThermoLawType2<<", "<<"Gamma2= "<<Gamma2<<\
                        ", PiSG2=  "<<PiSG2<<", kBAR2=  "<<kBAR2<<endl;
                }
                else if(key_word=="NBR_RELAXATION"){
                    getline(fichier,key_str_info, end_of_line);
                    NbrTauRelax=stoi(key_str_info);
                }
                else if(key_word=="RELAXATION_SCALES"){

                    for(int i_tau=0;i_tau<NbrTauRelax-1;i_tau++){
                        getline(fichier,key_str_info, ' ');
                        tauRelax(i_tau) = stod(key_str_info); 
                    }
                    //The last relaxation time
                    getline(fichier,key_str_info, end_of_line);
                    tauRelax(NbrTauRelax-1) = stod(key_str_info); 

                    cout<<"tauRelax: "<<tauRelax.transpose()<<endl;
                }
                else if(key_word=="PREF"){
                    getline(fichier,key_str_info, end_of_line);
                    pRef = stod(key_str_info);
                }
                else if(key_word=="MREF"){
                    getline(fichier,key_str_info, end_of_line);
                    mRef = stod(key_str_info);
                    cout<<"pRef: "<<pRef<<", mRef: "<<mRef<<endl;
                }
                else if(key_word=="SOURCE_TERM"){
                    getline(fichier,key_str_info, end_of_line);
                    if(key_str_info=="true"){
                        SourceTerm = true;
                    }
                    else{
                        SourceTerm = false;
                    }
                }
                else if(key_word=="NRELAX"){
                    getline(fichier,key_str_info, end_of_line);
                    NRelax = stoi(key_str_info);
                    cout<<"SourceTerm: "<<SourceTerm<<", NRelax: "<<NRelax<<endl;
                }
                else if(key_word=="LEFT_BC_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    LeftBCType=key_str_info;
                }
                else if(key_word=="RIGHT_BC_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    RightBCType=key_str_info;
                    cout<<"Left BC Type: "<<LeftBCType<<", "<<"Right BC Type: "<<RightBCType<<endl;
                    cout<<endl;
                }
                else if(key_word=="ALPHA1"){

                    getline(fichier,key_str_info, ' ');
                    alpha1_L = stod(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    alpha1_R = stod(key_str_info);
                }
                else if(key_word=="P1"){

                    getline(fichier,key_str_info, ' ');
                    p1_L = stod(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    p1_R = stod(key_str_info);
                }
                else if(key_word=="U1"){

                    getline(fichier,key_str_info, ' ');
                    u1_L = stod(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    u1_R = stod(key_str_info);
                }
                else if(key_word=="P2"){

                    getline(fichier,key_str_info, ' ');
                    p2_L = stod(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    p2_R = stod(key_str_info);
                }
                else if(key_word=="U2"){

                    getline(fichier,key_str_info, ' ');
                    u2_L = stod(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    u2_R = stod(key_str_info);

                    InitL<< alpha1_L,
                        p1_L,
                        u1_L,
                        p2_L,
                        u2_L;

                    InitR<< alpha1_R,
                        p1_R,
                        u1_R,    
                        p2_R,    
                        u2_R;

                    ThermoLaw therm_print(\
                            ThermoLawType1, ThermoLawType2,\
                            Gamma1, Gamma2,\
                            PiSG1,  PiSG2,\
                            kBAR1,  kBAR2\
                            );

                    Vector5d V_EqRelax_L = NConsVarToEqRelaxLoc(InitL,therm_print);
                    Vector5d V_EqRelax_R = NConsVarToEqRelaxLoc(InitR,therm_print);

                    cout<<"VARIABLES  |      LEFT      |      RIGHT     |"<<endl;
                    cout<<" alpha1    |"<<alpha1_L<<"|"<<alpha1_R<<"|"<<endl;
                    cout<<"    p1     |"<<p1_L<<"|"<<p1_R<<"|"<<endl;
                    cout<<"    u1     |"<<u1_L<<"|"<<u1_R<<"|"<<endl;
                    cout<<"    p2     |"<<p2_L<<"|"<<p2_R<<"|"<<endl;
                    cout<<"    u2     |"<<u2_L<<"|"<<u2_R<<"|"<<endl;
                    cout<<endl<<endl;

                    cout<<"VARIABLES  |      LEFT      |      RIGHT     |"<<endl;
                    cout<<" alpha1    |"<<V_EqRelax_L(0)<<"|"<<V_EqRelax_R(0)<<"|"<<endl;
                    cout<<"    U      |"<<V_EqRelax_L(1)<<"|"<<V_EqRelax_R(1)<<"|"<<endl;
                    cout<<"    P      |"<<V_EqRelax_L(2)<<"|"<<V_EqRelax_R(2)<<"|"<<endl;
                    cout<<"    du     |"<<V_EqRelax_L(3)<<"|"<<V_EqRelax_R(3)<<"|"<<endl;
                    cout<<"    dp     |"<<V_EqRelax_L(4)<<"|"<<V_EqRelax_R(4)<<"|"<<endl;
                    cout<<endl;

                }
                else if(key_word=="SCHEME_TYPE_CONS"){
                    getline(fichier,key_str_info, end_of_line);
                    SchemeTypeCons = key_str_info;
                }
                else if(key_word=="SCHEME_TYPE_NCONS"){
                    getline(fichier,key_str_info, end_of_line);
                    SchemeTypeNCons = key_str_info;

                    cout<<"Scheme type, CONS: "<<SchemeTypeCons<<", NCONS: "<<SchemeTypeNCons<<endl;
                }
                else if(key_word=="CFL_RAMP"){

                    getline(fichier,key_str_info, end_of_line);
                    if(key_str_info=="true"){
                        CFL_ramp=true;
                    }
                    else{
                        CFL_ramp=false;
                    }
                }
                else if(key_word=="CFL_RAMP_RANGE"){

                    getline(fichier,key_str_info, end_of_line);
                    if(CFL_ramp==true){

                        CFL_ramp_range=stoi(key_str_info);
                        //cout<<"CFL_ramp_range: "<<CFL_ramp_range<<endl;
                    }
                }
                else if(key_word=="CFL_CONV"){

                    getline(fichier,key_str_info, end_of_line);
                    CourantConv=stod(key_str_info);

                }
                else if(key_word=="CFL_BL"){

                    getline(fichier,key_str_info, end_of_line);
                    CourantBL=stod(key_str_info);
                    cout<<"CFL, CONV: "<<CourantConv<<", BL: "<<CourantBL<<endl;

                }
                else if(key_word=="LENGTH"){

                    getline(fichier,key_str_info, end_of_line);
                    Length=stod(key_str_info);

                }
                else if(key_word=="TRANSLATION_LENGTH"){

                    getline(fichier,key_str_info, end_of_line);
                    TranslationLength=stod(key_str_info);
                    cout<<"Length: "<<Length<<", TranslationLength: "<<TranslationLength<<endl;

                }
                else if(key_word=="NBR_AREAS"){

                    getline(fichier,key_str_info, ' ');
                    Nbr_Areas = stoi(key_str_info);
                    getline(fichier,key_str_info, end_of_line);
                    x0        = stod(key_str_info);
                    cout<<"Nbr_Areas: "<<Nbr_Areas<<", "<<"x0: "<<x0<<endl;

                }
                else if(key_word=="NGHOSTCELLS"){

                    getline(fichier,key_str_info, end_of_line);
                    NGhostCells=stoi(key_str_info);

                }
                else if(key_word=="NCELLS_LIST_NUMBER"){

                    getline(fichier,key_str_info, end_of_line);
                    nbNcells_cases=stoi(key_str_info);

                }
                else if(key_word=="NCELLS_LIST"){

                    CellsTab.resize(nbNcells_cases);
                    for (int j=0; j<nbNcells_cases;j++){

                        if(j<nbNcells_cases-1){

                            getline(fichier,key_str_info, ' ');
                        }
                        else{
                            getline(fichier,key_str_info, end_of_line);
                        }
                        CellsTab(j)=stod(key_str_info);
                    }
                }
                else if(key_word=="PRINT_FREQ"){

                    getline(fichier,key_str_info, end_of_line);
                    print_freq=stoi(key_str_info);
                }
                else if(key_word=="FILE_OUTPUT_FORMAT"){

                    getline(fichier,key_str_info, end_of_line);
                    FileOutputFormat=key_str_info;
                }
            }
        }

        fichier.close();  // on ferme le fichier
    }
    else{cerr << "Impossible d'ouvrir le fichier !" << endl;}

	/***** SIMULATION *****/
    //Mesh parameters
    int Ncells = CellsTab(0);
    Mesh  mesh_try(Length, Ncells, NGhostCells);

    //EOS parameters
    ThermoLaw therm_try(\
            ThermoLawType1, ThermoLawType2,\
            Gamma1, Gamma2,\
            PiSG1,  PiSG2,\
            kBAR1,  kBAR2\
            );

    //Solution parameters
    if(SolType=="Pure Contact"){

        //Single contact discontinuity test case
        InitR = WstateContactResolution(InitL, alpha1_R,\
                therm_try, epsDicho);

        alpha1_R = InitR(0);
            p1_R = InitR(1);
            u1_R = InitR(2);
            p2_R = InitR(3);
            u2_R = InitR(4);

        Vector5d V_EqRelax_L = NConsVarToEqRelaxLoc(InitL,therm_try);
        Vector5d V_EqRelax_R = NConsVarToEqRelaxLoc(InitR,therm_try);

        cout<<"Inside main, Single Contact Detected, New Right State: "<<endl<<endl;

        cout<<"VARIABLES  |      LEFT      |      RIGHT     |"<<endl;
        cout<<" alpha1    |"<<alpha1_L<<"|"<<alpha1_R<<"|"<<endl;
        cout<<"    p1     |"<<p1_L<<"|"<<p1_R<<"|"<<endl;
        cout<<"    u1     |"<<u1_L<<"|"<<u1_R<<"|"<<endl;
        cout<<"    p2     |"<<p2_L<<"|"<<p2_R<<"|"<<endl;
        cout<<"    u2     |"<<u2_L<<"|"<<u2_R<<"|"<<endl;
        cout<<endl<<endl;

        cout<<"VARIABLES  |      LEFT      |      RIGHT     |"<<endl;
        cout<<" alpha1    |"<<V_EqRelax_L(0)<<"|"<<V_EqRelax_R(0)<<"|"<<endl;
        cout<<"    U      |"<<V_EqRelax_L(1)<<"|"<<V_EqRelax_R(1)<<"|"<<endl;
        cout<<"    P      |"<<V_EqRelax_L(2)<<"|"<<V_EqRelax_R(2)<<"|"<<endl;
        cout<<"    du     |"<<V_EqRelax_L(3)<<"|"<<V_EqRelax_R(3)<<"|"<<endl;
        cout<<"    dp     |"<<V_EqRelax_L(4)<<"|"<<V_EqRelax_R(4)<<"|"<<endl;
        cout<<endl;
    }

    //Solver parameters
    Sol_Isen sol_try(mesh_try,\
		therm_try,\
        pRef, mRef,\
        tauRelax,\
        InitL, InitR,\
		SolType,\
		LeftBCType, RightBCType,\
        x0\
		);

    double dtRelax = Big;
    double SimulationTime = 5.e-2;

    Solver_Isen solver_try(dtRelax, NRelax, CourantBL,\
            SchemeTypeCons, SchemeTypeNCons,\
            TimeStep, CourantConv,\
            SimulationTime, print_freq\
            );
    solver_try.Simulation(sol_try, mesh_try, FileOutputFormat, Simulation_Name);

    // */

    /**********************************************************************************/
    /*****************  ANALYTICAL SOLUTION RHS OF THE H-DYNAMICS ****************/
    /**********************************************************************************/

    /*
    Vector5d W_state = sol_try.NConsVar_.row(0).transpose();

    Matrix5d LSTEq = LinearizedSourceTermsEq(W_state,\
            sol_try.SolTherm_,\
            pRef, mRef,\
            etaP, etaU\
            );

    cout<<"Inside main LinearizedSourceTermEq"<<endl;
    cout<<LSTEq<<endl<<endl;

    Matrix5d LJEqRel = LinearizedJacobianEqRelax(W_state,\
            sol_try.SolTherm_\
            );

    cout<<"Inside main LinearizedJacobianEqRelax"<<endl;
    cout<<LJEqRel<<endl<<endl;

    //Plotting the relaxation dynamics:
    Vector5d LST_EqODE;
    double alpha1, U, P, du, dp;

    string filename = "ODE_LinearizedSourceTermEq.dat";
    ofstream file((filename).c_str(), ios::out);
    if(file){

        file<<"Time"<<" "<<"alpha1"<<" "<<"U"<<" "<<"P"<<" "<<"du"<<" "<<"dp"<<endl;

        for (double t=ZERO; t<=30*tauMin; t+=tauMin/HUNDRED){

            LST_EqODE = LinearizedSourceTermEqODE(\
                    W_state_avr, W_state_ini,\
                    sol_try.SolTherm_, pRef, mRef,\
                    tauMin, etaP, etaU, t\
                    );

            alpha1  = LST_EqODE(0);
            U       = LST_EqODE(1);
            P       = LST_EqODE(2);
            du      = LST_EqODE(3);
            dp      = LST_EqODE(4);

            file<<std::setprecision(COUT_PRECISION)<<std::scientific\
              <<t<<" "<<alpha1<<" "<<U<<" "<<P<<" "<<du<<" "<<dp<<endl;
             }
    }
    else{cerr << "Linearized Source Term Equilibrium ODE, Can't open file !" << endl;}

    */

    /**********************************************************************************/
    /*****************  BOUNDARY LAYER RESOLUTION ****************/
    /**********************************************************************************/

    /*
    Vector5d W_state_L   = InitL;
    Vector5d W_state_R   = InitR;
    W_state_avr = NConsVarAveraged(\
        W_state_L, W_state_R,\
        ONE_OVER_TWO);

    Vector5d H_state = NConsVarToEqRelax(W_state_avr,\
        sol_try.SolTherm_\
        );

    BoundaryLayerResolutionLocal(\
            H_state, W_state_L, W_state_R,\
            sol_try.SolTherm_, pRef, mRef,\
            etaP, etaU,\
            dtRelax, tauMin, NRelax,\
            SpaceStep);
    */



    /**********************************************************************************/
    /*****************  STABILITY CONVECTION SOURCE BN: U-P RELAXATION ****************/
    /**********************************************************************************/

    /*

    //Equilibrium definition

    double alpha0_eq = 0.8;
    double alpha1_eq = 0.8;
    double alpha2_eq = 0.2;

    double p0_eq = 2e6; //(Pa)
    double s1_eq = 2.44e3; //(J.kg-1.K-1)
    double s2_eq = 6.32e3; //(J.kg-1.K-1)

    double rho1_eq = 8.46e2; //(kg.m-3)
    double rho2_eq = 1.20e1; //(kg.m-3)

    double u0_eq   = 1.e0;   //(m.s-1)
    //double c1_eq   = 1.50e3; //(m.s-1)
    double c1_eq   = 3.50e2; //(m.s-1)
    double c2_eq   = 3.00e2; //(m.s-1)

    //Thermodynamics parameters
    ThermoLawType1 = "SG";
    ThermoLawType2 = "SG";
    Gamma1         = 1.4;
    //Gamma2         = 7.5;
    Gamma2         = 1.4;
    PiSG1          = ZERO;
    //PiSG2          = 3e8;
    PiSG2          = ZERO;

    ThermoLaw therm_try_UP_relax(\
            ThermoLawType1, ThermoLawType2,\
            Gamma1, Gamma2,\
            PiSG1,  PiSG2,\
            kBAR1,  kBAR2\
            );

    //Space-step and Fourier mode parameters
    double dx = 1e-3; //(m)
    double k  = ONE/(TWO*dx); //(m-1)

    //Time-step and relaxation time parameters
    double dt_tauPbar = ONE_OVER_TWO;
    double dt_Big     = dx/max(c1_eq, c2_eq); 
    double tauPbar    = 1.e-12; //(s)
    double tauPbar_tauU = 1.e-3;

    double K_integrator = 3.;

    Row4cd Fourier_R1_UP_relax = Fourier_R1_BN_UP(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            k, dx);

    Row4cd Fourier_W1_UP_relax = Fourier_W1_BN_UP(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            k, dx);

    Row4cd Fourier_R2_UP_relax = Fourier_R2_BN_UP(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            k, dx);

    Row4cd Fourier_W2_UP_relax = Fourier_W2_BN_UP(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            k, dx);

    cout<<"Inside main, Fourier_R1_UP_relax: "<<endl;
    cout<<Fourier_R1_UP_relax<<endl<<endl;
    cout<<"Inside main, Fourier_W1_UP_relax: "<<endl;
    cout<<Fourier_W1_UP_relax<<endl<<endl;
    cout<<"Inside main, Fourier_R2_UP_relax: "<<endl;
    cout<<Fourier_R2_UP_relax<<endl<<endl;
    cout<<"Inside main, Fourier_W2_UP_relax: "<<endl;
    cout<<Fourier_W2_UP_relax<<endl<<endl;

    Mat4cd Fourier_Matrix_BN_UP_relax = Fourier_Matrix_BN_UP(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            k, dx);

    cout<<"Inside main, Fourier_Matrix_BN_UP_relax: "<<endl;
    cout<<Fourier_Matrix_BN_UP_relax<<endl<<endl;

    int dk = 25;

    Fourier_Matrix_BN_UP_EigenValues(\
            alpha0_eq, u0_eq, p0_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            dt_tauPbar, tauPbar, tauPbar_tauU,\
            dt_Big, K_integrator,\
            dk, dx);

    Matrix_BN_UP_EigenValues(\
            u0_eq, p0_eq,\
            alpha1_eq, alpha2_eq,\
            rho1_eq, rho2_eq, s1_eq, s2_eq,\
            c1_eq, c2_eq,\
            therm_try_UP_relax,\
            tauPbar_tauU);


    string filename = "Wood_Velocity.dat";

    rho1_eq = 1.20e1; //(kg.m-3)
    rho2_eq = 8.46e2; //(kg.m-3)

    c1_eq   = 3.50e2; //(m.s-1)
    c2_eq   = 1.00e3; //(m.s-1)

    double C1 = rho1_eq*pow(c1_eq,TWO);
    double C2 = rho2_eq*pow(c2_eq,TWO);
    double c1sq = pow(c1_eq,TWO);
    double c2sq = pow(c2_eq,TWO);

    double M      = ZERO;
    double Tau    = ZERO;
    double Y      = ZERO;
    double Cwood  = ZERO; 
    double Ctilde = ZERO; 
    double Cstar  = ZERO; 

    //Plotting the Wood curve:

    ofstream file((filename).c_str(), ios::out);
    if(file){

        file<<"Alpha_1"<<" "<<"Cwood"<<" "<<"Ctilde"<<" "<<"Cstar"<<endl;

        for (double alpha1=ZERO; alpha1<=ONE; alpha1+=0.005){

            M      = alpha1*rho1_eq + (ONE-alpha1)*rho2_eq;
            Y      = (alpha1*rho1_eq)/M;
            Tau    = (ONE-alpha1)/rho1_eq + alpha1/rho2_eq;
            Cwood  = sqrt( (ONE/M) * (C1*C2 / ((ONE-alpha1)*C1 + alpha1*C2)) );  
            Ctilde = sqrt( c1sq*c2sq / ((ONE-alpha1)*c1sq + alpha1*c2sq));
            Cstar  = sqrt(Tau*((ONE-Y)*C1 + Y*C2));

            file<<std::setprecision(COUT_PRECISION)<<std::scientific\
              <<alpha1<<" "<<Cwood<<" "<<Ctilde<<" "<<Cstar<<endl;
             }
    }
    else{cerr << "Wood_Velocity, Can't open file !" << endl;}

    //Plotting the modified characteristic curves when linear source term:

    filename = "Characteristic_Source_Term.dat";
    
    double X_L(ZERO), X_R(ZERO);
    double u_eq = 1.;
    double u_L = 0.5;
    double u_R = 0.75;
    double tau = 0.1;
    double Time = 0.5;

    ofstream file2((filename).c_str(), ios::out);
    if(file){

        file2<<"t"<<" "<<"X_L"<<" "<<"X_R"<<endl;

        for (double t=ZERO; t<=Time; t+=0.01){
            
            X_L = u_eq*t + (u_L*tau) * (ONE -  exp(-(t/tau)));
            X_R = u_eq*t + (u_R*tau) * (ONE -  exp(-(t/tau)));

            file2<<std::setprecision(COUT_PRECISION)<<std::scientific\
              <<t<<" "<<X_L<<" "<<X_R<<endl;
             }
    }
    else{cerr << "Characteristic_Source_Term, Can't open file !" << endl;}


    // */

	/***** SIMULATION END *****/

    cout << "Hello World!\n"<<endl;

    return 0;

}


