#include "Solver_Isentropic.h"

//constructor:
Solver_Isen::Solver_Isen( double dtRelax,\
        int NRelax, double CourantBL,\
        string SchemeTypeCons, string SchemeTypeNCons,\
        double TimeStep, double CourantConv,\
        double SimulationTime, int print_freq\
        ){

    dtRelax_         = dtRelax;
    NRelax_          = NRelax;
    CourantBL_       = CourantBL;
    CourantConv_     = CourantConv;
    SchemeTypeCons_  = SchemeTypeCons;
    SchemeTypeNCons_ = SchemeTypeNCons;
    TimeStep_        = TimeStep;
    SimulationTime_  = SimulationTime;
    print_freq_      = print_freq;
}

Solver_Isen::Solver_Isen( Solver_Isen& solver){

    dtRelax_         = solver.dtRelax_;
    NRelax_          = solver.NRelax_;
    CourantBL_       = solver.CourantBL_;
    CourantConv_     = solver.CourantConv_;
    SchemeTypeCons_  = solver.SchemeTypeCons_;
    SchemeTypeNCons_ = solver.SchemeTypeNCons_;
    TimeStep_        = solver.TimeStep_;
    SimulationTime_  = solver.SimulationTime_;
    print_freq_      = solver.print_freq_;
}

Solver_Isen::Solver_Isen(){

    dtRelax_   = ZERO;
    NRelax_    = 0;
    CourantBL_ = ZERO;
}

/************************************************/
/***************  BOUNDARY LAYER  ***************/
/************************************************/

void Solver_Isen::BoundaryLayerUpdate(\
        Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double dtRelax, int NRelax\
        ){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R, W_state_avr;
    Vector5d H_state_avr;

    //Function

    //mesh inputs
    double SpaceStep = mesh.Get_SpaceStep();

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    double tauMin = sol.etaRelax_(0);
    double etaU   = sol.etaRelax_(1);
    double etaP   = sol.etaRelax_(2);

    //Initialization of H_state_avr:
    H_state_avr<<ZERO,
                 ZERO,
                 ZERO,
                 ZERO,
                 ZERO;

    int Nfaces      = mesh.Get_Nfaces();
    double NcellExt = mesh.Get_NcellExt();

    for(int face_id = 0; face_id < Nfaces; face_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            L = mesh.FaceIndex_(face_id,1);
            R = mesh.FaceIndex_(face_id,2);

            W_state_L = sol.NConsVar_.row(L).transpose();
            W_state_R = sol.NConsVar_.row(R).transpose();
        }
        //Boundary Layer Resolution performed on the faces of the DUAL mesh
        else{

            L = mesh.FaceIndexDual_(face_id,1);
            R = mesh.FaceIndexDual_(face_id,2);

            W_state_L = sol.NConsVarDual_.row(L).transpose();
            W_state_R = sol.NConsVarDual_.row(R).transpose();
        }

        double weight_L = ONE_OVER_TWO;
        W_state_avr = NConsVarAveraged(W_state_L, W_state_R, weight_L);
        H_state_avr = NConsVarToEqRelaxLoc(W_state_avr, sol.SolTherm_);

        //Resolution of the time-boundary layer

        BoundaryLayerResolutionLocal(\
                H_state_avr, W_state_L, W_state_R,\
                sol.SolTherm_, sol.pRef_, sol.mRef_,\
                etaP, etaU,\
                dtRelax, tauMin, NRelax,\
                SpaceStep\
                );

        //Update of the staggered solution

        //PRIMAL mesh --> DUAL update
        if (MeshTag =="Primal"){

            sol.NConsVarDual_.row(face_id) = EqRelaxToNConsVar(H_state_avr,\
                    sol.SolTherm_).transpose();
        }
        //DUAL mesh --> PRIMAL update + Imposed boundary conditions
        else{

            sol.NConsVar_.row(face_id+1) = EqRelaxToNConsVar(H_state_avr,\
                    sol.SolTherm_).transpose();

            if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                sol.NConsVar_.row(0)          = sol.NConsVar_.row(1);
                sol.NConsVar_.row(NcellExt-1) = sol.NConsVar_.row(NcellExt-2);

            }
        }
    }

}

void Solver_Isen::BoundaryLayerTimeStepUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R, W_state_avr;

    //Function

    //Always initialize dtRelax_ to a big value
    dtRelax_ = Big;

    //mesh inputs
    double SpaceStep = mesh.Get_SpaceStep();

    //relaxation inputs
    double tauMin = sol.etaRelax_(0);
    double etaU   = sol.etaRelax_(1);
    double etaP   = sol.etaRelax_(2);

    int Nfaces      = mesh.Get_Nfaces();

    for(int face_id = 0; face_id < Nfaces; face_id++){

        L = mesh.FaceIndex_(face_id,1);
        R = mesh.FaceIndex_(face_id,2);

        W_state_L = sol.NConsVar_.row(L).transpose();
        W_state_R = sol.NConsVar_.row(R).transpose();

        double weight_L = ONE_OVER_TWO;
        W_state_avr = NConsVarAveraged(W_state_L, W_state_R, weight_L);

        //local timestep provided by the relaxation process
        //and the jacobian eigenvalues
        double dtRelaxConvLoc = LocalCourant_LSTEq(W_state_avr,\
                sol.SolTherm_, sol.pRef_, sol.mRef_,\
                tauMin, etaP, etaU,\
                SpaceStep,\
                CourantBL_\
                );

        dtRelax_ = min(dtRelax_, dtRelaxConvLoc);
    }

}

void Solver_Isen::BoundaryLayer(Sol_Isen& sol, Mesh& mesh){

    cout<<"Inside BoundaryLayer, dt initial = "<<dtRelax_<<endl;

    //Update of the timestep in order to update the time boundary layer
    BoundaryLayerTimeStepUpdate(sol, mesh);

    cout<<"Inside BoundaryLayer, dt updated = "<<dtRelax_<<endl;

    string MeshTag = "Primal";
    int NRelaxHalf = floor(NRelax_/2);

    cout<<"Inside BoundaryLayer, NConsVar_     = "<<endl;
    cout<<sol.NConsVar_<<endl<<endl;
    cout<<"Inside BoundaryLayer, NConsVarDual_ = "<<endl;
    cout<<sol.NConsVarDual_<<endl<<endl;

    //Update starting from PRIMAL face
    BoundaryLayerUpdate(sol, mesh, MeshTag,\
            dtRelax_, NRelaxHalf\
            );

    cout<<"Inside BoundaryLayer PRIMAL UPDATE, NConsVar_     = "<<endl;
    cout<<sol.NConsVar_<<endl<<endl;
    cout<<"Inside BoundaryLayer PRIMAL UPDATE, NConsVarDual_ = "<<endl;
    cout<<sol.NConsVarDual_<<endl<<endl;

    MeshTag = "Dual";

    //Update starting from DUAL face
    BoundaryLayerUpdate(sol, mesh, MeshTag,\
            dtRelax_, NRelaxHalf\
            );

    cout<<"Inside BoundaryLayer DUAL UPDATE, NConsVar_     = "<<endl;
    cout<<sol.NConsVar_<<endl<<endl;
    cout<<"Inside BoundaryLayer DUAL UPDATE, NConsVarDual_ = "<<endl;
    cout<<sol.NConsVarDual_<<endl<<endl;
}

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

void Solver_Isen::ConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R, W_state_avr;
    double alpha1_L, alpha1_R;

    //Function

    int Nfaces = mesh.Get_Nfaces();    

    for(int face_id = 0; face_id < Nfaces; face_id++){

            L = mesh.FaceIndex_(face_id,1);
            R = mesh.FaceIndex_(face_id,2);

            W_state_L = sol.NConsVar_.row(L).transpose();
            W_state_R = sol.NConsVar_.row(R).transpose();

            //Storing the face quantities needed for the non-conservative 
            //discrete terms
            alpha1_L = W_state_L(0);
            alpha1_R = W_state_R(0);

            sol.NConsVarFace_(face_id) = ONE_OVER_TWO*(alpha1_L + alpha1_R);

            //Building the conservative flux
            sol.ConsVarFlux_.row(face_id) = ConsVarFluxUpdateLoc(W_state_L, W_state_R,\
                    sol.SolTherm_,\
                    SchemeTypeCons_\
                    ).transpose();
        }

}

void Solver_Isen::NConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int Lf,Rf;
    Vector5d W_state, NConsTerm_state;
    double alpha1_state, alpha1_Lf, alpha1_Rf;
    double uI_state, pI_state;

    //Function

    int NcellExt = mesh.Get_NcellExt();

    for(int cell_id = 0; cell_id < NcellExt; cell_id++){

        W_state   = sol.NConsVar_.row(cell_id).transpose();
        //uI = u2
        uI_state  = W_state(4); 
        //pI = p1
        pI_state  = W_state(1);

        Lf        = mesh.CellIndex_(cell_id,1);
        Rf        = mesh.CellIndex_(cell_id,2);

        //Left boundary condition
        if(cell_id==0){

            alpha1_Rf    = sol.NConsVarFace_(Rf);
            alpha1_state = W_state(0);
            alpha1_Lf    = alpha1_state;

        }
        //Right boundary condition
        else if(cell_id==NcellExt-1){

            alpha1_Lf    = sol.NConsVarFace_(Lf);
            alpha1_state = W_state(0);
            alpha1_Rf    = alpha1_state;
        }
        //normal cell
        else{

            alpha1_Lf = sol.NConsVarFace_(Lf);
            alpha1_Rf = sol.NConsVarFace_(Rf);
        }

        NConsTerm_state<<uI_state*(alpha1_Rf - alpha1_Lf),
            ZERO,
            -pI_state*(alpha1_Rf - alpha1_Lf),
            ZERO,
            pI_state*(alpha1_Rf - alpha1_Lf);
        sol.NConsVarFlux_.row(cell_id) = NConsTerm_state.transpose();
    }

}

void Solver_Isen::ConsVarUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int Lf,Rf;

    //Function

    int Ncells       = mesh.Get_Ncells();
    int NcellExt     = mesh.Get_NcellExt();
    double SpaceStep = mesh.Get_SpaceStep();

    string  LeftBCType   =sol.LeftBCType_;
    string  RightBCType  =sol.RightBCType_;

    for(int cell_id = 1; cell_id <= Ncells; cell_id++){

        Lf = mesh.CellIndex_(cell_id,1);
        Rf = mesh.CellIndex_(cell_id,2);

        sol.ConsVar_.row(cell_id) +=-(TimeStep_/SpaceStep)*(
            sol.ConsVarFlux_.row(Rf) - sol.ConsVarFlux_.row(Lf) +\
            sol.NConsVarFlux_.row(cell_id) );
    }

    //Treatment of the boundary conditions
    if(LeftBCType=="transparent"){
        sol.ConsVar_.row(0) = sol.ConsVar_.row(1);

    }
    if(RightBCType=="transparent"){
        sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

    }
    else{
        sol.ConsVar_.row(0)          = MatrixXd::Zero(1,5);
        sol.ConsVar_.row(NcellExt-1) = MatrixXd::Zero(1,5);
    }

}

void Solver_Isen::TimeStepUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R;

    //Function

    //Always initialize TimeStep to a big value
    TimeStep_ = Big;

    //mesh inputs
    double SpaceStep = mesh.Get_SpaceStep();
    int Nfaces       = mesh.Get_Nfaces();

    for(int face_id = 0; face_id < Nfaces; face_id++){

        L = mesh.FaceIndex_(face_id,1);
        R = mesh.FaceIndex_(face_id,2);

        W_state_L = sol.NConsVar_.row(L).transpose();
        W_state_R = sol.NConsVar_.row(R).transpose();

        //local timestep provided by the jacobian eigenvalues
        double dtConvLoc = LocalCourant_Conv(W_state_L, W_state_R,\
                sol.SolTherm_,\
                SpaceStep,\
                CourantConv_\
                );
        TimeStep_ = min(TimeStep_, dtConvLoc);
    }

}

void Solver_Isen::Simulation(Sol_Isen& sol, Mesh& mesh,\
                string FileOutputFormat, string FileName\
                ){


    //Local variables
    double time(ZERO);
    int ite(0);
   
    //Saving solution
    Save_Sol(FileOutputFormat, FileName,\
                sol, mesh, ite);

    while(time < SimulationTime_){

        //Update of the discrete timestep
        TimeStepUpdate(sol, mesh);

        time+=TimeStep_;
        ite++;

        //cout<<"Inside Simulation: TimeStep_ = "<<TimeStep_<<endl;

        //Update of the conservative flux
        ConsVarFluxUpdate(sol, mesh);

        //Update of the non-conservative flux
        NConsVarFluxUpdate(sol, mesh);

        //Update of the conservative variables
        ConsVarUpdate(sol, mesh);

        //Update NConsVar_ using the matrix ConsVar_
        sol.ConsVarToNConsVar();

        //Update NConsVarEqRelax_ using the matrix NConsVar_
        sol.NConsVarToEqRelax();

        //Updates SoundSpeed_  of Sol using new values
        sol.SoundSpeed_Update();

        //Updates MachMax_ , MachMin_ using the new Sol values
        sol.Mach_Update();
        
        //Updating the exact solution
        sol.SolExact_Update(mesh, time);

        //Saving solution
        Save_Sol(FileOutputFormat, FileName,\
                sol, mesh, ite);

        if(ite%print_freq_==0){

            //Saving results at last iteration
            cout<<"Inside Simulation: "<<time/SimulationTime_<<"  completed"<<endl; 
        }
    }

}

void Solver_Isen::Save_Sol(string FileOutputFormat, string FileName,\
        Sol_Isen& sol, Mesh& mesh, int ite){

    int NcellExt= mesh.Get_NcellExt();
    int Ncells= mesh.Get_Ncells();

    char numstrite[20]; // enough to hold all numbers up to 64-bits
    char numstrNcells[20]; // enough to hold all numbers up to 64-bits
    sprintf(numstrite, "%d", ite);
    sprintf(numstrNcells, "%d",Ncells);

    double x,y,z;
    string FileName_ite;
    string FileName_Face_ite;

    string Folder_location="./Output/";

    //Get the right precision for all the prints
    std::setprecision(COUT_PRECISION);

    if(FileOutputFormat==".vtk"){
        FileName_ite=Folder_location+FileName+"_Cell"+numstrNcells+"_ite"+numstrite+".vtk";
    }

    else if(FileOutputFormat==".dat"){
        FileName_ite=Folder_location+FileName+"_Cell"+numstrNcells+"_ite"+numstrite+".dat";
    }

    ofstream file(FileName_ite.c_str(), ios::out);  //file flux declaration and file opening

    if(file){  // If the opening is successful

        if(FileOutputFormat==".vtk"){

            file<<"# vtk DataFile Version 2.0"<<endl;
            file<<"Simulation Data"<<endl;
            file<<"ASCII"<<endl;
            file<<"DATASET STRUCTURED_GRID"<<endl;
            file<<"DIMENSIONS 1 1 1"<<endl;
            file<<"POINTS "<<NcellExt<<"  double"<<endl;

            for (int i=0; i<NcellExt;i++){

                x=mesh.CellCoordsTab_(i,1);
                y=mesh.CellCoordsTab_(i,2);
                z=mesh.CellCoordsTab_(i,3);
                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<x<<" "<<y<<" "<<z<<endl;
            }
            file<<endl;

            file<<"POINT_DATA "<<NcellExt<<endl;

            /*********************************/
            /***** COMPUTED SOLUTION  ********/
            /*********************************/

            file<<"SCALARS "<<"Alpha1"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVar_(i,0)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P1"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVar_(i,1)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U1"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVar_(i,2)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P2"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVar_(i,3)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U2"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVar_(i,4)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVarEqRelax_(i,1)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVarEqRelax_(i,2)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"dU"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVarEqRelax_(i,3)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"dP"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.NConsVarEqRelax_(i,4)<<endl;

            }
            file<<endl;

            /*********************************/
            /******* EXACT SOLUTION  *********/
            /*********************************/

            file<<"SCALARS "<<"Alpha1_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExact_(i,0)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P1_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExact_(i,1)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U1_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExact_(i,2)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P2_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExact_(i,3)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U2_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExact_(i,4)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"U_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExactEqRelax_(i,1)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"P_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExactEqRelax_(i,2)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"dU_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExactEqRelax_(i,3)<<endl;

            }
            file<<endl;
            
            file<<"SCALARS "<<"dP_ex"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.SolExactEqRelax_(i,4)<<endl;

            }
            file<<endl;
        }

        else if(FileOutputFormat==".dat"){

            //Writting the columns names:

            file<<"x"<<" "<<"y"<<" "<<"z"<<" "\
                <<"Alpha1"<<" "<<"P1"<<" "<<"U1"<<" "<<"P2"<<" "<<"U2"<<" "\
                <<"U"<<" "<<"P"<<" "<<"dU"<<" "<<"dP"<<" "\
                <<"Alpha1_ex"<<" "<<"P1_ex"<<" "<<"U1_ex"<<" "<<"P2_ex"<<" "<<"U2_ex"<<" "\
                <<"U_ex"<<" "<<"P_ex"<<" "<<"dU_ex"<<" "<<"dP_ex"<<endl;

            //Writing values
            for (int i=0; i<NcellExt;i++){

                x=mesh.CellCoordsTab_(i,1);
                y=mesh.CellCoordsTab_(i,2);
                z=mesh.CellCoordsTab_(i,3);
                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<x<<" "<<y<<" "<<z<<" ";

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<\
                    sol.NConsVar_(i,0)<<" "<<sol.NConsVar_(i,1)<<" "\
                    <<sol.NConsVar_(i,2)<<" "<<sol.NConsVar_(i,3)<<" "\
                    <<sol.NConsVar_(i,4)<<" "<<sol.NConsVarEqRelax_(i,1)<<" "\
                    <<sol.NConsVarEqRelax_(i,2)<<" "<<sol.NConsVarEqRelax_(i,3)<<" "\
                    <<sol.NConsVarEqRelax_(i,4)<<" "\
                    <<sol.SolExact_(i,0)<<" "<<sol.SolExact_(i,1)<<" "\
                    <<sol.SolExact_(i,2)<<" "<<sol.SolExact_(i,3)<<" "\
                    <<sol.SolExact_(i,4)<<" "<<sol.SolExactEqRelax_(i,1)<<" "\
                    <<sol.SolExactEqRelax_(i,2)<<" "<<sol.SolExactEqRelax_(i,3)<<" "\
                    <<sol.SolExactEqRelax_(i,4)<<endl;
            }

        }

        file.close();  // the file is closed

    }
    else  {cerr << "Error Save_Sol, Can't open file !" << endl;}

}

//External Functions

/********** BOUNDARY LAYER RESOLUTION ************/

Vector5d LinearizedSourceTermEqODE(\
        Vector5d& W_state_avr, Vector5d& W_state_ini,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time\
        ){

    //Output
    Vector5d LinSouTerEq;

    //Local variables
    double alpha1_avr, alpha2_avr, rho1_avr, rho2_avr, p1_avr, p2_avr;
    double m1_avr, m2_avr, m_avr, c1_avr, c2_avr, C1_avr, C2_avr;

    double alpha1_ini; 
    double U_ini, P_ini, du_ini, dp_ini;

    //Deriving the equilibrium/relaxation variables for the initial condition
    Vector5d EqRelax_ini = NConsVarToEqRelaxLoc(W_state_ini, Therm);

    alpha1_ini = EqRelax_ini(0);
    U_ini      = EqRelax_ini(1);
    P_ini      = EqRelax_ini(2);
    du_ini     = EqRelax_ini(3);
    dp_ini     = EqRelax_ini(4);

    //Getting the other non-conservative quantities for the averaged state
    alpha1_avr = W_state_avr(0);
    alpha2_avr = ONE - alpha1_avr;
    p1_avr     = W_state_avr(1);
    p2_avr     = W_state_avr(3);

    rho1_avr = Density_EOS(1, Therm, p1_avr, ZERO);
    rho2_avr = Density_EOS(2, Therm, p2_avr, ZERO);
    m1_avr = alpha1_avr*rho1_avr;
    m2_avr = alpha2_avr*rho2_avr;
    m_avr  = m1_avr + m2_avr;
    c1_avr = Sound_Speed_EOS(1, Therm, rho1_avr, p1_avr);
    c2_avr = Sound_Speed_EOS(2, Therm, rho2_avr, p2_avr);

    C1_avr = rho1_avr*pow(c1_avr,TWO);
    C2_avr = rho2_avr*pow(c2_avr,TWO);

    //Deriving the co-factors related to the pressure and velocity relaxations
    double lambda_du = etaU*Ku_Coeff(mRef, alpha1_avr)*( m_avr/(m1_avr*m2_avr)  );
    double lambda_dp = etaP*Kp_Coeff(pRef, alpha1_avr)*( C1_avr/alpha1_avr + C2_avr/alpha2_avr  );

    //cout<<"Inside LinearizedSourceTermEqODE: lambda_du = "<<lambda_du<<", lambda_dp = "<<lambda_dp<<endl;

    //Time-integrators
    double TIdu_decr = exp(-lambda_du*( time/tauRelax  ) );
    double TIdp_decr = exp(-lambda_dp*( time/tauRelax  ) );
    double TIdp_plat = ONE - TIdp_decr;

    //Cofactors contributions
    double Cofac_avr = ( alpha1_avr*alpha2_avr  )/( alpha2_avr*C1_avr + alpha1_avr*C2_avr  );

    //Filling the output vector
    LinSouTerEq<<alpha1_ini - TIdp_plat*Cofac_avr*dp_ini,
                 U_ini,
                 P_ini      -  TIdp_plat*Cofac_avr*(C2_avr - C1_avr)*dp_ini,
                 TIdu_decr*du_ini,
                 TIdp_decr*dp_ini;

    return LinSouTerEq;
}

double LocalCourant_LSTEq(Vector5d& W_state_avr,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauMin, double etaP, double etaU,\
        double SpaceStep,\
        double Damp\
        ){

    //Local variables
    double alpha1_avr, alpha2_avr, rho1_avr, rho2_avr, p1_avr, p2_avr;
    double m1_avr, m2_avr, m_avr, c1_avr, c2_avr, C1_avr, C2_avr;

    double u1_avr, u2_avr;

    //Getting the other non-conservative quantities for the averaged state
    alpha1_avr = W_state_avr(0);
    alpha2_avr = ONE - alpha1_avr;
    p1_avr     = W_state_avr(1);
    p2_avr     = W_state_avr(3);

    u1_avr = W_state_avr(2);
    u2_avr = W_state_avr(4);

    rho1_avr = Density_EOS(1, Therm, p1_avr, ZERO);
    rho2_avr = Density_EOS(2, Therm, p2_avr, ZERO);
    m1_avr = alpha1_avr*rho1_avr;
    m2_avr = alpha2_avr*rho2_avr;
    m_avr  = m1_avr + m2_avr;
    c1_avr = Sound_Speed_EOS(1, Therm, rho1_avr, p1_avr);
    c2_avr = Sound_Speed_EOS(2, Therm, rho2_avr, p2_avr);

    C1_avr = rho1_avr*pow(c1_avr,TWO);
    C2_avr = rho2_avr*pow(c2_avr,TWO);

    //Deriving the co-factors related to the pressure and velocity relaxations
    double lambda_du = etaU*Ku_Coeff(mRef, alpha1_avr)*( m_avr/(m1_avr*m2_avr)  );
    double lambda_dp = etaP*Kp_Coeff(pRef, alpha1_avr)*( C1_avr/alpha1_avr + C2_avr/alpha2_avr  );

    double dtRelax = min( ONE/lambda_du, ONE/lambda_dp  )*tauMin;
    double dtConv  = (SpaceStep/TWO)*(ONE/(max(fabs(u2_avr)+c2_avr, fabs(u1_avr)+c1_avr)));

    /*
       cout<<"Inside LocalCourant_LSTEq: lambda_du = "<<lambda_du<<", lambda_dp = "<<lambda_dp<<endl<<endl;
     */

    return Damp*min(dtRelax, dtConv);

}

void BoundaryLayerResolutionLocal(\
        Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin, double NRelax,\
        double SpaceStep\
        ){

    //Local variables
    double time;
    Vector5d LSTEqODE_L, LSTEqODE_R;

    //Function
    int nrows = H_state.rows();
    Matrix5d Id = MatrixXd::Identity(nrows, nrows);

    //Averaged state
    Vector5d W_state_avr = NConsVarAveraged(\
        W_state_L, W_state_R,\
        ONE_OVER_TWO);

    //Linearized source term 
    Matrix5d LinSouTermEq = LinearizedSourceTermsEq(W_state_avr,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );

    //Constant Jacobian
    Matrix5d LinJacEqRelax = LinearizedJacobianEqRelax(W_state_avr,\
            Therm\
            );

    //Initialization
    Vector5d H_ini = H_state;

        time = ZERO;

        for (int ite=1;ite<=NRelax;ite++){

            time = ite*dtRelax;

            //Non homogeneous term of the ODE
            LSTEqODE_L = LinearizedSourceTermEqODE(\
                    W_state_avr, W_state_L,\
                    Therm, pRef, mRef,\
                    tauMin, etaP, etaU, time\
                    );

            //Non homogeneous term of the ODE
            LSTEqODE_R = LinearizedSourceTermEqODE(\
                    W_state_avr, W_state_R,\
                    Therm, pRef, mRef,\
                    tauMin, etaP, etaU, time\
                    );

            //Dynamics of H
            H_state =(Id + (dtRelax/tauMin)*LinSouTermEq )*H_ini\
                     -(dtRelax/SpaceStep)*LinJacEqRelax*( LSTEqODE_R - LSTEqODE_L  );

            //Update of H_ini
            H_ini = H_state;

        }

}

/*
   void BoundaryLayerResolution(\
   Vector5d& H_state,Vector5d& W_state_L, Vector5d& W_state_R,\
   ThermoLaw& Therm, double pRef, double mRef,\
   double etaP, double etaU,\
   double dtRelax, double tauMin, double NRelax,\
   double SpaceStep\
   ){

//Local variables
double time;
Vector5d LSTEqODE_L, LSTEqODE_R;

//Function
int nrows = H_state.rows();
Matrix5d Id = MatrixXd::Identity(nrows, nrows);

//Averaged state
Vector5d W_state_avr = NConsVarAveraged(\
W_state_L, W_state_R,\
ONE_OVER_TWO);

//Linearized source term 
Matrix5d LinSouTermEq = LinearizedSourceTermsEq(W_state_avr,\
Therm,\
pRef, mRef,\
etaP, etaU\
);

//Constant Jacobian
//Matrix5d LinJacEqRelax = LinearizedJacobianEqRelax(W_state_avr,\
Therm\
);

double c0       = 3.e2;
double rho_pert = 1.e3;
double p_pert   = 5.e6;
Matrix5d LinJacEqRelax = c0*Id;
LinJacEqRelax(3,2) = ONE/rho_pert;
LinJacEqRelax(4,1) = p_pert;

cout<<"Inside ODE H: Linearized Jacobian = "<<endl;
cout<<LinJacEqRelax<<endl<<endl;

//Initialization
Vector5d H_ini = H_state;

string filename     = "ODE_H.dat";
string filename_err = "ODE_H_Pert_Error.dat";
ofstream file((filename).c_str(), ios::out);
ofstream file_err((filename_err).c_str(), ios::app);

//printed variables
double alpha1, U, P, du, dp;

//Exact solution
Vector5d SolSingleCharact;
double alpha1_ex, U_ex, P_ex, du_ex, dp_ex;

//Error Curves Quantitites
double alpha1_err(ZERO), U_err(ZERO), P_err(ZERO), du_err(ZERO), dp_err(ZERO);
double norm_alpha1_ex(ZERO), norm_U_ex(ZERO), norm_P_ex(ZERO), norm_du_ex(ZERO), norm_dp_ex(ZERO);

if(file){

file<<"Time"<<" "<<"alpha1"<<" "<<"U"<<" "<<"P"<<" "<<"du"<<" "<<"dp"\
<<" "<<" "<<"alpha1_ex"<<" "<<"U_ex"<<" "<<"P_ex"<<" "<<"du_ex"<<" "<<"dp_ex"<<endl;

//Save H_state at initial condition
alpha1  = H_state(0);
U       = H_state(1);
P       = H_state(2);
du      = H_state(3);
dp      = H_state(4);
time = ZERO;

file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<time<<" "<<alpha1<<" "<<U<<" "<<P<<" "<<du<<" "<<dp\
        <<" "<<alpha1<<" "<<U<<" "<<P<<" "<<du<<" "<<dp<<endl;

while (time <= NRelax*tauMin){

    time += dtRelax;

    //Non homogeneous term of the ODE
    LSTEqODE_L = LinearizedSourceTermEqODE(\
            W_state_avr, W_state_L,\
            Therm, pRef, mRef,\
            tauMin, etaP, etaU, time\
            );

    //Non homogeneous term of the ODE
    LSTEqODE_R = LinearizedSourceTermEqODE(\
            W_state_avr, W_state_R,\
            Therm, pRef, mRef,\
            tauMin, etaP, etaU, time\
            );

    //Dynamics of H
    H_state =(Id + (dtRelax/tauMin)*LinSouTermEq )*H_ini\
             -(dtRelax/SpaceStep)*LinJacEqRelax*( LSTEqODE_R - LSTEqODE_L  );

    //Update of H_ini
    H_ini = H_state;

    //Save H_state
    alpha1  = H_state(0);
    U       = H_state(1);
    P       = H_state(2);
    du      = H_state(3);
    dp      = H_state(4);

    //Save Exact H (when known)

    SolSingleCharact = BoundaryLayerSolSingleCharact_Pert(\
            W_state_avr, W_state_L, W_state_R,\
            Therm, pRef, mRef,\
            tauMin, etaP, etaU, time,\
            c0, rho_pert, p_pert, SpaceStep\
            );

    alpha1_ex  = SolSingleCharact(0);
    U_ex       = SolSingleCharact(1);
    P_ex       = SolSingleCharact(2);
    du_ex      = SolSingleCharact(3);
    dp_ex      = SolSingleCharact(4);



    alpha1_ex  = alpha1 ;
    U_ex       = U      ;
    P_ex       = P      ;
    du_ex      = du     ;
    dp_ex      = dp     ;


    //Updating the error curves


    alpha1_err  += fabs(alpha1-alpha1_ex);
    U_err       += fabs(U-U_ex);
    P_err       += fabs(P-P_ex);
    du_err      += fabs(du-du_ex);
    dp_err      += fabs(dp-dp_ex);

    norm_alpha1_ex  += fabs(alpha1_ex);
    norm_U_ex       += fabs(U_ex);
    norm_P_ex       += fabs(P_ex);
    norm_du_ex      += fabs(du_ex);
    norm_dp_ex      += fabs(dp_ex);

    file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<time<<" "<<alpha1<<" "<<U<<" "<<P<<" "<<du<<" "<<dp\
        <<" "<<alpha1_ex<<" "<<U_ex<<" "<<P_ex<<" "<<du_ex<<" "<<dp_ex<<endl;

}


alpha1_err  /= norm_alpha1_ex;
U_err       /= norm_U_ex;
P_err       /= norm_P_ex;
du_err      /= norm_du_ex;
dp_err      /= norm_dp_ex;


}

else{cerr << "ODE H, Can't open file !" << endl;}

file_err<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<dtRelax<<" "<<alpha1_err<<" "<<U_err<<" "<<P_err<<" "<<du_err<<" "<<dp_err<<endl;

}

*/


/********** exact solution functions ************/

Vector5d BoundaryLayerSolSingleCharact(\
        Vector5d& W_state_avr, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time,\
        double c0, double dx\
        ){

    //Output
    Vector5d SolSingleCharact;

    //Local variables
    double alpha1_avr, alpha2_avr, rho1_avr, rho2_avr, p1_avr, p2_avr;
    double m1_avr, m2_avr, m_avr, c1_avr, c2_avr, C1_avr, C2_avr;

    double alpha1_ini; 
    double U_ini, P_ini, du_ini, dp_ini;

    double dalpha1_ini, dU_ini, dP_ini, ddu_ini, ddp_ini;

    //Deriving the equilibrium/relaxation variables for the initial condition
    Vector5d EqRelax_ini   = NConsVarToEqRelaxLoc(W_state_avr, Therm);
    Vector5d EqRelax_ini_L = NConsVarToEqRelaxLoc(W_state_L, Therm);
    Vector5d EqRelax_ini_R = NConsVarToEqRelaxLoc(W_state_R, Therm);

    alpha1_ini = EqRelax_ini(0);
    U_ini      = EqRelax_ini(1);
    P_ini      = EqRelax_ini(2);
    du_ini     = EqRelax_ini(3);
    dp_ini     = EqRelax_ini(4);

    dalpha1_ini = EqRelax_ini_R(0) - EqRelax_ini_L(0);
    dU_ini      = EqRelax_ini_R(1) - EqRelax_ini_L(1);
    dP_ini      = EqRelax_ini_R(2) - EqRelax_ini_L(2);
    ddu_ini     = EqRelax_ini_R(3) - EqRelax_ini_L(3);
    ddp_ini     = EqRelax_ini_R(4) - EqRelax_ini_L(4);

    //Getting the other non-conservative quantities for the averaged state
    alpha1_avr = W_state_avr(0);
    alpha2_avr = ONE - alpha1_avr;
    p1_avr     = W_state_avr(1);
    p2_avr     = W_state_avr(3);

    rho1_avr = Density_EOS(1, Therm, p1_avr, ZERO);
    rho2_avr = Density_EOS(2, Therm, p2_avr, ZERO);
    m1_avr = alpha1_avr*rho1_avr;
    m2_avr = alpha2_avr*rho2_avr;
    m_avr  = m1_avr + m2_avr;
    c1_avr = Sound_Speed_EOS(1, Therm, rho1_avr, p1_avr);
    c2_avr = Sound_Speed_EOS(2, Therm, rho2_avr, p2_avr);

    C1_avr = rho1_avr*pow(c1_avr,TWO);
    C2_avr = rho2_avr*pow(c2_avr,TWO);

    //Deriving the co-factors related to the pressure and velocity relaxations
    double lambda_du = etaU*Ku_Coeff(mRef, alpha1_avr)*( m_avr/(m1_avr*m2_avr)  );
    double lambda_dp = etaP*Kp_Coeff(pRef, alpha1_avr)*( C1_avr/alpha1_avr + C2_avr/alpha2_avr  );

    //cout<<"Inside LinearizedSourceTermEqODE: lambda_du = "<<lambda_du<<", lambda_dp = "<<lambda_dp<<endl;

    //Time-integrators
    double TIdu_decr = exp(-lambda_du*( time/tauRelax  ) );
    double TIdp_decr = exp(-lambda_dp*( time/tauRelax  ) );
    double TIdp_plat = ONE - TIdp_decr;
    //double TIdp_lin  = -time + ( tauRelax/lambda_dp  )*TIdp_plat;

    //Cofactors contributions
    double Cofac_avr = ( alpha1_avr*alpha2_avr  )/( alpha2_avr*C1_avr + alpha1_avr*C2_avr  );

    //Filling the output vector
    SolSingleCharact<<alpha1_ini - TIdp_plat*Cofac_avr*dp_ini  - (c0*time/dx)*(dalpha1_ini -TIdp_plat*Cofac_avr*ddp_ini),
                 U_ini - (c0*time/dx)*(dU_ini),
                 P_ini - TIdp_plat*Cofac_avr*(C2_avr - C1_avr)*dp_ini - (c0*time/dx)*(dP_ini - TIdp_plat*Cofac_avr*(C2_avr - C1_avr)*ddp_ini),
                 TIdu_decr*( du_ini- (c0*time/dx)*(ddu_ini) ),
                 TIdp_decr*( dp_ini- (c0*time/dx)*(ddp_ini) );

    return SolSingleCharact;
}

Vector5d BoundaryLayerSolSingleCharact_Pert(\
        Vector5d& W_state_avr, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauRelax, double etaP, double etaU, double time,\
        double c0, double rho_pert, double p_pert, double dx\
        ){

    //Output
    Vector5d SolSingleCharact, SolKconv, SolKrelax, Sol; 

    //Local variables
    double alpha1_avr, alpha2_avr, rho1_avr, rho2_avr, p1_avr, p2_avr;
    double m1_avr, m2_avr, m_avr, c1_avr, c2_avr, C1_avr, C2_avr;

    double dU_ini, dP_ini, ddp_ini;

    //Deriving the equilibrium/relaxation variables for the initial condition
    Vector5d EqRelax_ini_L = NConsVarToEqRelaxLoc(W_state_L, Therm);
    Vector5d EqRelax_ini_R = NConsVarToEqRelaxLoc(W_state_R, Therm);

    dU_ini      = EqRelax_ini_R(1) - EqRelax_ini_L(1);
    dP_ini      = EqRelax_ini_R(2) - EqRelax_ini_L(2);
    ddp_ini     = EqRelax_ini_R(4) - EqRelax_ini_L(4);

    //Getting the other non-conservative quantities for the averaged state
    alpha1_avr = W_state_avr(0);
    alpha2_avr = ONE - alpha1_avr;
    p1_avr     = W_state_avr(1);
    p2_avr     = W_state_avr(3);

    rho1_avr = Density_EOS(1, Therm, p1_avr, ZERO);
    rho2_avr = Density_EOS(2, Therm, p2_avr, ZERO);
    m1_avr = alpha1_avr*rho1_avr;
    m2_avr = alpha2_avr*rho2_avr;
    m_avr  = m1_avr + m2_avr;
    c1_avr = Sound_Speed_EOS(1, Therm, rho1_avr, p1_avr);
    c2_avr = Sound_Speed_EOS(2, Therm, rho2_avr, p2_avr);

    C1_avr = rho1_avr*pow(c1_avr,TWO);
    C2_avr = rho2_avr*pow(c2_avr,TWO);

    //Deriving the co-factors related to the pressure and velocity relaxations
    double lambda_du = etaU*Ku_Coeff(mRef, alpha1_avr)*( m_avr/(m1_avr*m2_avr)  );
    double lambda_dp = etaP*Kp_Coeff(pRef, alpha1_avr)*( C1_avr/alpha1_avr + C2_avr/alpha2_avr  );

    //cout<<"Inside LinearizedSourceTermEqODE: lambda_du = "<<lambda_du<<", lambda_dp = "<<lambda_dp<<endl;

    //Time-integrators
    double TIdu_decr = exp(-lambda_du*( time/tauRelax  ) );
    double TIdp_decr = exp(-lambda_dp*( time/tauRelax  ) );
    double TIdp_plat = ONE - TIdp_decr;
    double TIdu_plat = ONE - TIdu_decr;

    //Cofactors contributions
    double Cofac_avr = ( alpha1_avr*alpha2_avr  )/( alpha2_avr*C1_avr + alpha1_avr*C2_avr  );
    //pertubation contributions
    double time_pert          =  (ONE/(lambda_du-lambda_dp))*( (lambda_dp/lambda_du)*TIdu_plat - TIdp_plat  );
    double time_pert_conv     =  ( TIdp_plat/lambda_dp - time/tauRelax  );
    double Cofac_pert         = (tauRelax/dx)*Cofac_avr*(C2_avr - C1_avr)*time_pert/(rho_pert);
    double Cofac_pert_conv_alpha = (tauRelax/dx)*p_pert*time_pert_conv*Cofac_avr;

    double Cofac_pert_conv_P     = (tauRelax/dx)*p_pert*time_pert_conv*Cofac_avr*(C2_avr - C1_avr);
    double Cofac_pert_conv_du    = (tauRelax/dx)/(rho_pert)*(TIdu_plat/lambda_du);
    double Cofac_pert_conv_dp    = (tauRelax/dx)*p_pert*(TIdp_plat/lambda_dp);

    //Solution of the c0_characteristic
    SolSingleCharact = BoundaryLayerSolSingleCharact(\
        W_state_avr, W_state_L, W_state_R,\
        Therm, pRef, mRef,\
        tauRelax, etaP, etaU, time,\
        c0, dx\
        ); 

    //Delta contribution Kconv matrix
    SolKconv<< - Cofac_pert_conv_alpha*dU_ini,
              ZERO,
              - Cofac_pert_conv_P*dU_ini,
              - Cofac_pert_conv_du*dP_ini,
              - Cofac_pert_conv_dp*dU_ini;

    //Delta contribution Krelax matrix
    SolKrelax<< ZERO,
                ZERO,
                ZERO,
                - Cofac_pert*ddp_ini,
                ZERO;

    //Final solution
    Sol = SolSingleCharact + SolKconv + SolKrelax;

    return Sol;

}

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

Vector5d ConsVarFluxUpdateLoc(\
        Vector5d& W_state_L, Vector5d&  W_state_R,\
        ThermoLaw& Therm,\
        string SchemeTypeCons\
        ){

    //Local Variables
    Vector5d ConsFlux, Flux_L, Flux_R;
    Vector5d U_state_L, U_state_R;
    
    //Function
    if(SchemeTypeCons=="Rusanov"){

        Flux_L    = ConsVarFlux(W_state_L, Therm); 
        Flux_R    = ConsVarFlux(W_state_R, Therm); 
        U_state_L = NConsVarToConsVar(W_state_L, Therm);
        U_state_R = NConsVarToConsVar(W_state_R, Therm);

        double lambda = SpectralRadiusRusanov(\
                W_state_L, W_state_R,\
                Therm);

        ConsFlux = ONE_OVER_TWO*(Flux_L + Flux_R) -\
                   (ONE_OVER_TWO*lambda)*(U_state_R - U_state_L);

    }
    else{

        ConsFlux<<ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO;

    }

    return ConsFlux;

}

double SpectralRadiusRusanov(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        ){

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2, u1, u2;
    double m, m1, m2, Y1, Y2;

    double c1, c2, C1, C2, c_star, c_w;
    double U;

    //Function

    //Left cell

    alpha1 = W_state_L(0);
    alpha2 = ONE - alpha1;
    p1     = W_state_L(1);
    p2     = W_state_L(3);
    u1     = W_state_L(2);
    u2     = W_state_L(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    C1 = rho1*c1*c1;
    C2 = rho2*c2*c2;

    double lambda1_cL = max(fabs(u1 + c1), fabs(u1 - c1));
    double lambda2_cL = max(fabs(u2 + c2), fabs(u2 - c2));

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    Y1 = m1/m;
    Y2 = m2/m;

    U      = Y1*u1 + Y2*u2;
    c_star = sqrt(alpha2*c1*c1 + alpha1*c2*c2);
    c_w    = sqrt( (ONE/m)*( C1*C2/(alpha1*C2 + alpha2*C1)  ) );

    double lambda_cwL    = max(fabs(U + c_w), fabs(U - c_w));
    double lambda_cstarL = max(fabs(U + c_star), fabs(U - c_star));

    //Right cell

    alpha1 = W_state_R(0);
    alpha2 = ONE - alpha1;
    p1     = W_state_R(1);
    p2     = W_state_R(3);
    u1     = W_state_R(2);
    u2     = W_state_R(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    C1 = rho1*c1*c1;
    C2 = rho2*c2*c2;

    double lambda1_cR = max(fabs(u1 + c1), fabs(u1 - c1));
    double lambda2_cR = max(fabs(u2 + c2), fabs(u2 - c2));

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    Y1 = m1/m;
    Y2 = m2/m;

    U      = Y1*u1 + Y2*u2;
    c_star = sqrt(alpha2*c1*c1 + alpha1*c2*c2);
    c_w    = sqrt( (ONE/m)*( C1*C2/(alpha1*C2 + alpha2*C1)  ) );

    double lambda_cwR    = max(fabs(U + c_w), fabs(U - c_w));
    double lambda_cstarR = max(fabs(U + c_star), fabs(U - c_star));

    //Max of the eigenvalues

    double lambda_ckL = max(lambda1_cL, lambda2_cL);
    lambda_ckL = max(lambda_ckL, lambda_cwL);
    lambda_ckL = max(lambda_ckL, lambda_cstarL);
    double lambda_ckR = max(lambda1_cR, lambda2_cR);
    lambda_ckR = max(lambda_ckR, lambda_cwR);
    lambda_ckR = max(lambda_ckR, lambda_cstarR);

    double lambda = max(lambda_ckL, lambda_ckR);

    return lambda;

}

double LocalCourant_Conv(Vector5d& W_state_L, Vector5d& W_state_R,\
                ThermoLaw& Therm,\
                double SpaceStep,\
                double Courant\
                ){

    //Local variables
    double rho1, rho2, p1, p2, u1, u2;
    double c1, c2;

    //Function

    //Left cell

    p1     = W_state_L(1);
    p2     = W_state_L(3);
    u1     = W_state_L(2);
    u2     = W_state_L(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    double lambda1_cL = max(fabs(u1 + c1), fabs(u1 - c1));
    double lambda2_cL = max(fabs(u2 + c2), fabs(u2 - c2));

    //Right cell

    p1     = W_state_R(1);
    p2     = W_state_R(3);
    u1     = W_state_R(2);
    u2     = W_state_R(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    double lambda1_cR = max(fabs(u1 + c1), fabs(u1 - c1));
    double lambda2_cR = max(fabs(u2 + c2), fabs(u2 - c2));

    //Max of the eigenvalues

    double lambda_ckL = max(lambda1_cL, lambda2_cL);
    double lambda_ckR = max(lambda1_cR, lambda2_cR);
    double lambda = max(lambda_ckL, lambda_ckR);

    return Courant*(SpaceStep/lambda);

}
