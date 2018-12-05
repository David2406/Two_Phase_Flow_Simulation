#include <iostream>
#include "Mesh.h"
#include "Sol_Isentropic.h"
#include "Solver_Isentropic.h"
#include "ThermoLaw.h"
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

    /* READING PARAMETER FILE */
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
    string VariableType;

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

    string SourceTermType;
    bool FractionalStep;
    string TimeIntegrationType;
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
    double SimulationTime;
    
    //Mesh
    double Length;
    double TranslationLength;
    int NGhostCells;
    int nbNcells_cases;
    Array<int, Eigen::Dynamic, 1> CellsTab;
    
    //Riemann Problem parameters
    int Nbr_Areas;
    double x0;

    //We print every print_freq iteration
    int print_freq;
    string FileOutputFormat;

    // file is opened in reading mode
    ifstream fichier( file_name.c_str(), ios::in);  

    if(fichier)  
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

                }
                else if(key_word=="SOLUTION_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    SolType=key_str_info;
                    cout<<"Sol Type: "<<SolType<<endl;

                }
                else if(key_word=="VARIABLE_TYPE"){

                    getline(fichier,key_str_info,
                            end_of_line);
                    VariableType=key_str_info;
                    cout<<"Variable Type: "<<VariableType<<endl;

                }
                else if(key_word=="THERMODYNAMICS1"){

                    getline(fichier,key_str_info, ' ');
                    ThermoLawType1=key_str_info;
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

                    Simulation_Name+="_"+ThermoLawType1+"_";
                    cout<<"Thermodynamics1: "<<ThermoLawType1<<", "<<"Gamma1= "<<Gamma1<<\
                        ", PiSG1=  "<<PiSG1<<", kBAR1=  "<<kBAR1<<endl;
                }
                else if(key_word=="THERMODYNAMICS2"){

                    getline(fichier,key_str_info, ' ');
                    ThermoLawType2=key_str_info;
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
                else if(key_word=="SOURCE_TERM_TYPE"){
                    getline(fichier,key_str_info, end_of_line);
                    SourceTermType=key_str_info;
                }
                else if(key_word=="FRACTIONAL_STEP"){
                    getline(fichier,key_str_info, end_of_line);
                    if(key_str_info=="true"){
                        FractionalStep = true;
                        Simulation_Name +="FracStep";
                    }
                    else{
                        FractionalStep = false;
                    }
                }
                else if(key_word=="TIME_INTEGRATION_TYPE"){
                    getline(fichier,key_str_info, end_of_line);
                    TimeIntegrationType=key_str_info;
                }
                else if(key_word=="NRELAX"){
                    getline(fichier,key_str_info, end_of_line);
                    NRelax = stoi(key_str_info);
                    if(FractionalStep==true){

                        cout<<"FractionalStep: "<<FractionalStep<<", SourceTermType: "<<SourceTermType<<\
                            ", NO BS RELAXATION"<<endl;
                    }
                    else{
                        cout<<"FractionalStep: "<<FractionalStep<<", SourceTermType: "<<SourceTermType<<\
                            ", NRelax: "<<NRelax<<", Time-integration: "<<TimeIntegrationType<<endl;
                    }
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
                else if(key_word=="SIMULATION_TIME"){

                    getline(fichier,key_str_info, end_of_line);
                    SimulationTime=stod(key_str_info);
                    cout<<"Simulation Time: "<<SimulationTime<<endl;

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

        fichier.close();  
    }
    else{cerr << "Impossible to open the file !" << endl;}

    cout<<endl;

    /**********************/
	/***** SIMULATION *****/
    /**********************/

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
        if(fabs(u1_L-u2_L)> epsZero){

        InitR = WstateContactResolution(InitL, alpha1_R,\
                therm_try, epsDicho);
        }

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
    else if(SolType=="Shock Contact"){

        double rho1_L     = Density_EOS(1, therm_try, p1_L, ZERO);
        double rho2_L     = Density_EOS(2, therm_try, p2_L, ZERO);

        double sigma1Shock = TWO*u1_L - Sound_Speed_EOS(1, therm_try, rho1_L, p1_L);
        double sigma2Shock = TWO*u2_L - Sound_Speed_EOS(2, therm_try, rho2_L, p2_L);

        cout<<"Inside main, Shock-Contact detected, sigma1Shock = "<<sigma1Shock\
            <<", sigma2Shock = "<<sigma2Shock<<endl<<endl;

        InitL = WstateRHJumpCondition(\
                InitL,\
                sigma1Shock, sigma2Shock,\
                therm_try, epsDicho);

        alpha1_L = InitL(0);
            p1_L = InitL(1);
            u1_L = InitL(2);
            p2_L = InitL(3);
            u2_L = InitL(4);

        cout<<"New Left State: "<<endl<<endl;

        cout<<"VARIABLES  |      LEFT      |      RIGHT     |"<<endl;
        cout<<" alpha1    |"<<alpha1_L<<"|"<<alpha1_R<<"|"<<endl;
        cout<<"    p1     |"<<p1_L<<"|"<<p1_R<<"|"<<endl;
        cout<<"    u1     |"<<u1_L<<"|"<<u1_R<<"|"<<endl;
        cout<<"    p2     |"<<p2_L<<"|"<<p2_R<<"|"<<endl;
        cout<<"    u2     |"<<u2_L<<"|"<<u2_R<<"|"<<endl;
        cout<<endl<<endl;

    }

    /* Convergence_Curve creates the mandatory objects for the simulation: */
    // -a Mesh              object the computational domain is 1D and structured,
    // -a ThermoLaw         object storing the functions and parameters of the fluid Equation of State,
    // -a Sol_Isentropic    object which shelters the shape of the computed solution at any time,
    // -a Solver_Isentropic object which updates the solution for each time step using a given numerical scheme 

    Convergence_Curve(\
        Length, NGhostCells, \
        LeftBCType, RightBCType, \
        therm_try,\
        SolType,\
        VariableType,\
        InitL, InitR,\
        Nbr_Areas, x0,\
        SimulationTime, \
        CFL_ramp, CFL_ramp_range, \
        CourantBL, CourantConv,\
        FractionalStep,\
        SourceTermType,\
        TimeIntegrationType,\
        NRelax,\
        pRef, mRef,\
        tauRelax,\
        SchemeTypeCons, SchemeTypeNCons,\
        CellsTab, \
        FileOutputFormat, Simulation_Name, print_freq);

    /***************************/
	/****** SIMULATION END *****/
    /***************************/

    cout << "Hello World!\n"<<endl;

    return 0;

}


