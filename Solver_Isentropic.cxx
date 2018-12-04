#include "Solver_Isentropic.h"

//constructor:
Solver_Isen::Solver_Isen(\
        bool FractionalStep,\
        string SourceTermType,\
        string TimeIntegrationType,\
        int NRelax, double CourantBL,\
        string SchemeTypeCons, string SchemeTypeNCons,\
        double CourantConv,\
        double SimulationTime, int print_freq,\
        string FileOutputFormat,\
        int Nfaces\
        ){

    FractionalStep_        = FractionalStep;
    SourceTermType_        = SourceTermType;
    TimeIntegrationType_   = TimeIntegrationType;
    dtRelax_               = Big;
    NRelax_                = NRelax;
    CourantBL_             = CourantBL;
    CourantConv_           = CourantConv;
    SchemeTypeCons_        = SchemeTypeCons;
    SchemeTypeNCons_       = SchemeTypeNCons;
    TimeStep_              = Big;
    SimulationTime_        = SimulationTime;
    TimeElapsed_           = ZERO;
    print_freq_            = print_freq;
    FileOutputFormat_      = FileOutputFormat; 
    CPUTime_               = Big;

    //Allocating the size of the TimeMatrix_
    //dt, dt_BL, dt_Conv
    TimeMatrix_ = MatrixXd::Zero(Nfaces,3);
}

Solver_Isen::Solver_Isen( Solver_Isen& solver){

    FractionalStep_        = solver.FractionalStep_;
    SourceTermType_        = solver.SourceTermType_;
    TimeIntegrationType_   = solver.TimeIntegrationType_;
    dtRelax_               = solver.dtRelax_;
    NRelax_                = solver.NRelax_;
    CourantBL_             = solver.CourantBL_;
    CourantConv_           = solver.CourantConv_;
    SchemeTypeCons_        = solver.SchemeTypeCons_;
    SchemeTypeNCons_       = solver.SchemeTypeNCons_;
    TimeStep_              = solver.TimeStep_;
    TimeElapsed_           = solver.TimeElapsed_;
    SimulationTime_        = solver.SimulationTime_;
    print_freq_            = solver.print_freq_;
    FileOutputFormat_      = solver.FileOutputFormat_;
    CPUTime_               = solver.CPUTime_;

    TimeMatrix_ = solver.TimeMatrix_;
}

Solver_Isen::Solver_Isen(){

    dtRelax_   = ZERO;
    NRelax_    = 0;
    CourantBL_ = ZERO;
}

/************************************************/
/***************  BOUNDARY LAYER  ***************/
/************************************************/

/**********************/
/* BEREUX-SAINSAULIEU */
/**********************/

void Solver_Isen::BoundaryLayerUpdate(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double& dtRelax, double& dtRelax_estim, string VariableType){

    //Local variables
    int L,R;
    double x_face;
    Vector5d X_state_L, X_state_R;
    Vector5d H_state_avr, H_state_copy;
    //Matrix5d LinJac_L, LinJac_R;
    Matrix5d LinJac_R;
    double weight_L = ONE_OVER_TWO;

    double dtRelax_face_ini = ZERO;

    //Function

    //mesh inputs
    double SpaceStep = mesh.Get_SpaceStep();

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    //double tauMin = sol.etaRelax_(0);
    //double etaU   = sol.etaRelax_(1);
    //double etaP   = sol.etaRelax_(2);

    //Cubic relaxation 
    Vector5d STerm, STerm_L, STerm_R;
    Vector5d JacVector, JacVector_L, JacVector_R;
    Vector5d DiracMass;

    int Nfaces      = mesh.Get_Nfaces();
    double NcellExt = mesh.Get_NcellExt();

    LinJac_R  = sol.JacConvFrozen_;
    DiracMass = sol.DiracMassFrozen_;

    //Reference state
    Vector5d U_state_ref = NConsVarToConsVarLoc(sol.W_ref_, sol.SolTherm_);

    //Time matrix
    Matrix5d TimeMat = sol.TimeMat_;

    //FIXME
    //Relax_mu *= -ONE;

    //Vector5d ZeroTab = VectorXd::Zero(5);
    //FIXME
    //Linear relaxation
    //Matrix5d RelaxDiagMat = DiagonalMatrix<double, 5>(Relax_mu); 

    //Cubic relaxation
    Vector5d One = VectorXd::Constant(5, ONE);

    //Eigenvector base matrix
    Matrix5d RmatEigenVectors   = sol.EigenVectorBasis_;
   
    //Eigenvector base inverse matrix
    Matrix5d RmatEigenVectorsInv = sol.EigenVectorBasisInv_;

    //Equilibrium state
    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    //Time-step fixed
    dtRelax_face_ini   = dtRelax;

    for(int face_id = 0; face_id < Nfaces; face_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            L = mesh.FaceIndex_(face_id,1);
            R = mesh.FaceIndex_(face_id,2);

            x_face = mesh.FaceIndex_(face_id,3);

            X_state_L = sol.ConsVar_.row(L).transpose();
            X_state_R = sol.ConsVar_.row(R).transpose();
        }
        //Boundary Layer Resolution performed on the faces of the DUAL mesh
        else{

            L = mesh.FaceIndexDual_(face_id,1);
            R = mesh.FaceIndexDual_(face_id,2);

            x_face = mesh.FaceIndexDual_(face_id,3);

            X_state_L = sol.ConsVarDual_.row(L).transpose();
            X_state_R = sol.ConsVarDual_.row(R).transpose();
        }

        if(VariableType=="ConsVar"){

            H_state_avr = NConsVarAveraged(X_state_L, X_state_R, weight_L);

            H_state_copy = H_state_avr;

            //No-source terms active
            if(SourceTermType_=="none"){

                if(abs(x_face-sol.x_0_-sol.uI_Frozen_*TimeElapsed_)<\
                        ONE_OVER_TWO*SpaceStep){

                    //Gathering of the solution
                    H_state_avr = H_state_copy\
                                  -(dtRelax/SpaceStep)*LinJac_R*(X_state_R - X_state_L)\
                                  +(dtRelax/SpaceStep)*DiracMass;
                }
                else{

                    //Gathering of the solution
                    H_state_avr = H_state_copy\
                                  -(dtRelax/SpaceStep)*LinJac_R*(X_state_R - X_state_L);
                }
            }
            //Source terms are active
            else{

                //Source term initialization L
                if(sol.SolType_=="Linear Convection Relaxation BN"){

                    //Source term U_L
                    LinearSpring(TimeMat, X_state_L, U_state_eq, STerm_L,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector_L = LinearSpringGradient(TimeMat, X_state_L, U_state_eq,\
                            RmatEigenVectorsInv);

                    //Source term U_R
                    LinearSpring(TimeMat, X_state_R, U_state_eq, STerm_R,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector_R = LinearSpringGradient(TimeMat, X_state_R, U_state_eq,\
                            RmatEigenVectorsInv);

                    //Source term H
                    LinearSpring(TimeMat, H_state_copy, U_state_eq, STerm,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector = LinearSpringGradient(TimeMat, H_state_copy, U_state_eq,\
                            RmatEigenVectorsInv);
                }
                else if(sol.SolType_=="Linear Convection Relaxation Cubic BN"){

                    //Source term U_L
                    NonLinearSpring(TimeMat, X_state_L, U_state_eq, STerm_L,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector_L = NonLinearSpringGradient(TimeMat, X_state_L, U_state_eq,\
                            RmatEigenVectorsInv);

                    //Source term U_R
                    NonLinearSpring(TimeMat, X_state_R, U_state_eq, STerm_R,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector_R = NonLinearSpringGradient(TimeMat, X_state_R, U_state_eq,\
                            RmatEigenVectorsInv);

                    //Source term H
                    NonLinearSpring(TimeMat, H_state_copy, U_state_eq, STerm,\
                            RmatEigenVectors, RmatEigenVectorsInv);

                    JacVector = NonLinearSpringGradient(TimeMat, H_state_copy, U_state_eq,\
                            RmatEigenVectorsInv);
                }

                if(TimeIntegrationType_=="Rosenbrock4"){

                    //Source term U_L
                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_L, U_state_ref, U_state_eq,\
                            STerm_L, JacVector_L, One,\
                            dtRelax_face_ini\
                            );

                    //Source term U_R
                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_R, U_state_ref, U_state_eq,\
                            STerm_R, JacVector_R, One,\
                            dtRelax_face_ini\
                            );
                    //Source term H
                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_copy, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
                else if(TimeIntegrationType_=="ImplicitEuler1"){

                    //Source term U_L
                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_L, U_state_ref, U_state_eq,\
                            STerm_L, JacVector_L, One,\
                            dtRelax_face_ini\
                            );

                    //Source term U_R
                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_R, U_state_ref, U_state_eq,\
                            STerm_R, JacVector_R, One,\
                            dtRelax_face_ini\
                            );

                    //Source term H
                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_copy, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
                else if(TimeIntegrationType_=="ExplicitRungeKutta4"){

                    //Source term U_L
                    RungeKuttaFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_L, U_state_ref, U_state_eq,\
                            STerm_L,\
                            dtRelax_face_ini\
                            );

                    //Source term U_R
                    RungeKuttaFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            X_state_R, U_state_ref, U_state_eq,\
                            STerm_R,\
                            dtRelax_face_ini\
                            );

                    //Source term H
                    RungeKuttaFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_copy, U_state_ref, U_state_eq,\
                            STerm,\
                            dtRelax_face_ini\
                            );
                }

                if(abs(x_face-sol.x_0_-sol.uI_Frozen_*TimeElapsed_)<\
                        ONE_OVER_TWO*SpaceStep){

                    //Gathering of the solution
                    H_state_avr = H_state_copy\
                                  -(dtRelax/SpaceStep)*LinJac_R*(X_state_R - X_state_L)\
                                  +(dtRelax/SpaceStep)*DiracMass;
                }
                else{

                    //Gathering of the solution
                    H_state_avr = H_state_copy\
                                  -(dtRelax/SpaceStep)*LinJac_R*(X_state_R - X_state_L);
                }
            }


            //Update of the staggered solution

            //PRIMAL mesh --> DUAL update
            if (MeshTag =="Primal"){

                sol.ConsVarDual_.row(face_id) = H_state_avr.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVarDual_.row(0)          = sol.ConsVarDual_.row(1);
                    sol.ConsVarDual_.row(NcellExt-1) = sol.ConsVarDual_.row(NcellExt-2);

                }

            }
            //DUAL mesh --> PRIMAL update + Imposed boundary conditions
            else{

                sol.ConsVar_.row(face_id+1) = H_state_avr.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVar_.row(0)          = sol.ConsVar_.row(1);
                    sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

                }
            }
        }

    }

}

void Solver_Isen::BN_BoundaryLayerUpdate(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double& dtRelax, double& dtRelax_estim, string VariableType){

    //Local variables
    int L,R;
    double alpha1_mass_n;
    Vector5d X_state_L, X_state_R;
    Vector5d H_state_avr, H_state_copy;
    Matrix5d LinJac_R;
    double weight_L = ONE_OVER_TWO;

    //double dtRelax_face_ini = ZERO;

    //Function

    Matrix5d Id;
    Id = Eigen::Matrix<double, 5, 5>::Identity();

    //mesh inputs
    double SpaceStep = mesh.Get_SpaceStep();

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    //double tauMin = sol.etaRelax_(0);
    //double etaU   = sol.etaRelax_(1);
    //double etaP   = sol.etaRelax_(2);
    double tauP     = sol.tauRelax_(0);
    double tauU     = sol.tauRelax_(1);
    //double pref     = sol.pRef_;
    //double rhoref   = sol.mRef_;

    //Non-conservative products discretization
    Vector4d WaveSpeeds;
    Vector6d Star_UP; 

    //Cubic relaxation 
    Vector5d STerm, STerm_L, STerm_R;
    Matrix5d JacMatrix, JacMatrix_L, JacMatrix_R;

    Vector5d W_state_L, W_state_R, NConsFlux, Flux_L, Flux_R;

    int Nfaces      = mesh.Get_Nfaces();
    double NcellExt = mesh.Get_NcellExt();

    //Time matrix
    //Matrix5d TimeMat = sol.TimeMat_;

    //Cubic relaxation
    //Vector5d One = VectorXd::Constant(5, ONE);

    for(int face_id = 0; face_id < Nfaces; face_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            L = mesh.FaceIndex_(face_id,1);
            R = mesh.FaceIndex_(face_id,2);

            X_state_L = sol.ConsVar_.row(L).transpose();
            X_state_R = sol.ConsVar_.row(R).transpose();
        }
        //Boundary Layer Resolution performed on the faces of the DUAL mesh
        else{

            L = mesh.FaceIndexDual_(face_id,1);
            R = mesh.FaceIndexDual_(face_id,2);

            X_state_L = sol.ConsVarDual_.row(L).transpose();
            X_state_R = sol.ConsVarDual_.row(R).transpose();
        }

        if(VariableType=="ConsVar"){

            W_state_L = ConsVarToNConsVarLoc(X_state_L, sol.SolTherm_);
            W_state_R = ConsVarToNConsVarLoc(X_state_R, sol.SolTherm_);

            WaveSpeeds =  WaveSpeedEstimate(W_state_L, W_state_R, sol.SolTherm_);
            Star_UP =  U_P_star_states(W_state_L, W_state_R, WaveSpeeds, sol.SolTherm_);

            NConsFlux = NConsFluxLoc(W_state_L, W_state_R,\
                    WaveSpeeds, Star_UP, SchemeTypeNCons_, sol.SolTherm_);
        }

        H_state_avr = NConsVarAveraged(X_state_L, X_state_R, weight_L);

        H_state_copy = H_state_avr;

        //No-source terms active
        if(SourceTermType_=="none"){

            Flux_L    = ConsVarFlux(W_state_L, sol.SolTherm_); 
            Flux_R    = ConsVarFlux(W_state_R, sol.SolTherm_); 

        }
        else{

            /*
            //Source term U_L
            IsenBNSourceTerm(tauU, tauP, pref, rhoref,\
                    X_state_L, STerm_L,\
                    sol.SolTherm_);

            JacMatrix_L = IsenBNSourceTermGradient(tauU, tauP, pref, rhoref,\
                    X_state_L,\
                    sol.SolTherm_\
                    );

            //Source term U_R
            IsenBNSourceTerm(tauU, tauP, pref, rhoref,\
                    X_state_L, STerm_R,\
                    sol.SolTherm_);

            JacMatrix_R = IsenBNSourceTermGradient(tauU, tauP, pref, rhoref,\
                    X_state_R,\
                    sol.SolTherm_\
                    );

            //Source term H
            IsenBNSourceTerm(tauU, tauP, pref, rhoref,\
                    H_state_copy, STerm,\
                    sol.SolTherm_);

            JacMatrix = IsenBNSourceTermGradient(tauU, tauP, pref, rhoref,\
                    H_state_copy,\
                    sol.SolTherm_\
                    );
                    */

            if(TimeIntegrationType_=="Rosenbrock4"){

                BN_RosenBrockFourthOrder(\
                        tauP, tauU,\
                        sol.SolTherm_,\
                        X_state_L,\
                        dtRelax\
                        );

                BN_RosenBrockFourthOrder(\
                        tauP, tauU,\
                        sol.SolTherm_,\
                        X_state_R,\
                        dtRelax\
                        );

                BN_RosenBrockFourthOrder(\
                        tauP, tauU,\
                        sol.SolTherm_,\
                        H_state_copy,\
                        dtRelax\
                        );

            }
            else if(TimeIntegrationType_=="ImplicitEuler1"){

                BN_ImplicitVelocityPressure(\
                        X_state_L, sol.SolTherm_,\
                        tauU, tauP, dtRelax\
                        );
                BN_ImplicitVelocityPressure(\
                        X_state_R, sol.SolTherm_,\
                        tauU, tauP, dtRelax\
                        );
                BN_ImplicitVelocityPressure(\
                        H_state_copy, sol.SolTherm_,\
                        tauU, tauP, dtRelax\
                        );

                /*
                   BN_ImplicitEuler(\
                   X_state_L,\
                   STerm_L, JacMatrix_L, Id,\
                   dtRelax\
                   );

                   BN_ImplicitEuler(\
                   X_state_R,\
                   STerm_R, JacMatrix_R, Id,\
                   dtRelax\
                   );

                   BN_ImplicitEuler(\
                   H_state_copy,\
                   STerm, JacMatrix, Id,\
                   dtRelax\
                   );
                 */
            }
            else if(TimeIntegrationType_=="ExplicitRungeKutta4"){
            }

            W_state_L = ConsVarToNConsVarLoc(X_state_L, sol.SolTherm_);
            W_state_R = ConsVarToNConsVarLoc(X_state_R, sol.SolTherm_);
            Flux_L    = ConsVarFlux(W_state_L, sol.SolTherm_); 
            Flux_R    = ConsVarFlux(W_state_R, sol.SolTherm_); 

            //FIXME
            alpha1_mass_n = NConsFlux(0);
            NConsFlux(0)  = ZERO;
        }


            //Gathering of the solution
            H_state_avr = H_state_copy -\
                          (dtRelax/SpaceStep)*(Flux_R - Flux_L)+\
                          (dtRelax/SpaceStep)*NConsFlux;

            //FIXME
            //Special treatment of the alpha1 equation
            if(SourceTermType_!="none"){

                BN_Pressure_DichotomyBS(\
                        H_state_avr,\
                        tauP, dtRelax,\
                        SpaceStep,\
                        sol.SolTherm_,\
                        alpha1_mass_n,\
                        epsDicho);
            }
                    

            //Update of the staggered solution

            //PRIMAL mesh --> DUAL update
            if (MeshTag =="Primal"){

                sol.ConsVarDual_.row(face_id) = H_state_avr.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVarDual_.row(0)          = sol.ConsVarDual_.row(1);
                    sol.ConsVarDual_.row(NcellExt-1) = sol.ConsVarDual_.row(NcellExt-2);

                }

            }
            //DUAL mesh --> PRIMAL update + Imposed boundary conditions
            else{

                sol.ConsVar_.row(face_id+1) = H_state_avr.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVar_.row(0)          = sol.ConsVar_.row(1);
                    sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

                }
            }
        }

}

Vector5d NConsFluxLoc(Vector5d& W_state_L, Vector5d& W_state_R,\
        Vector4d& WaveSpeeds, Vector6d& Star_UP, string SchemeTypeNCons,\
        ThermoLaw& Therm){

    Vector5d NConsFlux;
    double alpha1_L = W_state_L(0);
    double alpha2_L = ONE - alpha1_L;
    double alpha1_R = W_state_R(0);
    double alpha2_R = ONE - alpha1_R;

    if(SchemeTypeNCons=="Default"){

        //Local variables
        //double u1_L     = W_state_L(2);
        double u2_L     = W_state_L(4);
        double p1_L     = W_state_L(1);
        //double p2_L     = W_state_L(3);

        //double u1_R     = W_state_R(2);
        double u2_R     = W_state_R(4);
        double p1_R     = W_state_R(1);
        //double p2_R     = W_state_R(3);

        double uI = ONE_OVER_TWO*(u2_L + u2_R); 
        double pI = ONE_OVER_TWO*(p1_L + p1_R); 

        //    NConsFlux<<-uI*(alpha1_R - alpha1_L),
        //               ZERO,
        //               -(alpha2_R*p2_R - alpha2_L*p2_L),
        //               ZERO,
        //                (alpha2_R*p2_R - alpha2_L*p2_L);

        NConsFlux<<-uI*(alpha1_R - alpha1_L),
            ZERO,
            +pI*(alpha1_R - alpha1_L),
            ZERO,
            -pI*(alpha1_R - alpha1_L);

        /* Test Non-conservative product of Ambroso, Chalons, Galié */
        /*
        //Construction de la constante de relaxation
        double rho1_L = Density_EOS(1, Therm, p1_L, ZERO); 
        double rho1_R = Density_EOS(1, Therm, p1_R, ZERO); 
        double c1_L   = Sound_Speed_EOS(1, Therm, rho1_L, p1_L);
        double c1_R   = Sound_Speed_EOS(1, Therm, rho1_R, p1_R);
        double a1 = 1.1*max(rho1_L*c1_L, rho1_R*c1_R);

        //Construction des quantités de la discontinuité de contact
        //On suppose U_L donné on construit U+
        double m = alpha1_L*rho1_L*(u1_L - u2_L);
        double cofac = (m/(alpha1_R*a1) - ONE)/(m/(alpha1_L*a1) - ONE);

        if (cofac < ZERO){

        cout<<"Inside NConsFlux BS: cofac < 0, rho1_p can not be computed"<<endl;
        exit(EXIT_FAILURE);
        }

        double rho1_p = sqrt(cofac)*rho1_L;
        double u1_p   = m/(alpha1_R*rho1_p) + u2_L;

        double a2p2_jump = m*(u1_L - u1_p) - p1_L*(alpha1_R - alpha1_L)\
        -alpha1_R*a1*a1*(ONE/rho1_L - ONE/rho1_p);

        NConsFlux<<-uI*(alpha1_R - alpha1_L),
        ZERO,
        -a2p2_jump,
        ZERO,
        +a2p2_jump;

         */
    }
    //HLLAC non-conservative products discretization
    else{

        double u2_star   = Star_UP(3);
        double p2_star_L = Star_UP(4);
        double p2_star_R = Star_UP(5);

        NConsFlux<<-u2_star*(alpha1_R - alpha1_L),
            ZERO,
            -(alpha2_R*p2_star_R - alpha2_L*p2_star_L),
            ZERO,
            +(alpha2_R*p2_star_R - alpha2_L*p2_star_L);
    }

    return NConsFlux;
}

void Solver_Isen::BoundaryLayerTimeStepUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R, W_state_avr;
    double dtRelaxLoc, dtConvLoc, dtRelaxConvLoc;

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

    Vector5d W_state0 = ONE_OVER_TWO*(sol.InitL_ + sol.InitR_);

    for(int face_id = 0; face_id < Nfaces; face_id++){

        L = mesh.FaceIndex_(face_id,1);
        R = mesh.FaceIndex_(face_id,2);

        W_state_L = sol.NConsVar_.row(L).transpose();
        W_state_R = sol.NConsVar_.row(R).transpose();

        //local timestep provided by the relaxation process
        //and the jacobian eigenvalues
        Vector2d TimeLoc = LocalCourant_LSTEq(\
                W_state0, W_state_L, W_state_R,\
                sol.SolTherm_, sol.pRef_, sol.mRef_,\
                tauMin, etaP, etaU,\
                SpaceStep,\
                CourantBL_,\
                face_id\
                );

        dtRelaxLoc = TimeLoc(0);
        dtConvLoc  = TimeLoc(1);

        dtRelaxConvLoc = min(dtRelaxLoc, dtConvLoc);

        dtRelax_ = min(dtRelax_, dtRelaxConvLoc);

        //Updating the time-steps
        TimeMatrix_(face_id, 0) = dtRelaxConvLoc;
        TimeMatrix_(face_id, 1) = dtRelaxLoc;
        TimeMatrix_(face_id, 2) = dtConvLoc;
        
    }

}

void Solver_Isen::BoundaryLayer(Sol_Isen& sol, Mesh& mesh){

    //Local Variables
    string MeshTag;
    string VariableType = sol.VariableType_;

    //Function

    double TimeStep       = dtRelax_/TWO;
    double TimeStep_estim = Big;

    for(int ite = 1; ite <=NRelax_; ite++){

        MeshTag = "Primal";

        //Update starting from PRIMAL face

        if(sol.SolType_=="BN Relaxation"){

            BN_BoundaryLayerUpdate(sol, mesh, MeshTag,\
                    TimeStep, TimeStep_estim, VariableType);
        }

        else{

            BoundaryLayerUpdate(sol, mesh, MeshTag,\
                    TimeStep, TimeStep_estim, VariableType);
        }

        MeshTag = "Dual";

        //Update starting from DUAL face

        if(sol.SolType_=="BN Relaxation"){

            BN_BoundaryLayerUpdate(sol, mesh, MeshTag,\
                    TimeStep, TimeStep_estim, VariableType);
        }

        else{

            BoundaryLayerUpdate(sol, mesh, MeshTag,\
                    TimeStep, TimeStep_estim, VariableType);
        }

        //Pauline Lafitte & Thomas Rey Idea
        //Storage of the solution at time NRelax-1 and NRelax
        //in order to reconstruct the time-gradient
        /*
        if(ite==NRelax_-1){

            sol.ConsVarOld_ = sol.ConsVar_;

        }
        */

    }

}

void Solver_Isen::BoundaryLayerTest(Sol_Isen& sol, Mesh& mesh, string Filename){

    //Local Variables
    string MeshTag;
    string VariableType = sol.VariableType_;

    //Function

    string FileOutputFormat = FileOutputFormat_;
    string FileName  = Filename;
    //string FileName2 = "BS_Extended_CFL005_adapt.dat";
    char numstrCFL[20]; // enough to hold all numbers up to 64-bits
    sprintf(numstrCFL, "%f", CourantBL_);
    string FileName2 = Filename+"_CFL"+numstrCFL+".dat";

    
    Save_Sol(FileOutputFormat, FileName,\
            sol, mesh, 0);

    double x_cell = 0.5;
    double TimeStep       = Big;
    double TimeStep_estim = Big;

    string Folder_location="./Output/";
    string FileName_ite;
    FileName_ite=Folder_location+FileName2;

    ofstream file(FileName_ite.c_str(), ios::out);  //file flux declaration and file opening

    file<<"Time"<<" "\
        <<"Alpha1"<<" "<<"P1"<<" "<<"U1"<<" "<<"P2"<<" "<<"U2"<<" "\
        <<"Alpha1_ex"<<" "<<"P1_ex"<<" "<<"U1_ex"<<" "<<"P2_ex"<<" "<<"U2_ex"<<endl;
    
    file.close();

    Save_Sol_Local(FileName2,\
            sol, mesh, x_cell, TimeElapsed_);

    int ite = 0;

    //FIXME
    //for (int ite = 1; ite <= NRelax_;ite++){
    //for (int ite = 0; ite <= 2;ite++){

    while (TimeElapsed_ < SimulationTime_){

    //Update of the timestep in order to update the time boundary layer
    BoundaryLayerTimeStepUpdate(sol, mesh);

    TimeStep = dtRelax_/TWO;

    MeshTag = "Primal";

    //Update starting from PRIMAL face

    BoundaryLayerUpdate(sol, mesh, MeshTag,\
            TimeStep, TimeStep_estim, VariableType);

    //FIXME
    TimeElapsed_ += TimeStep;
    
    MeshTag = "Dual";

    //Update starting from DUAL face

    BoundaryLayerUpdate(sol, mesh, MeshTag,\
    TimeStep, TimeStep_estim, VariableType);

    //FIXME
    TimeElapsed_ += TimeStep;

    //TimeElapsed_ += dtRelax_;
    //TimeElapsed_ += Big;
    ite++;

    if(ite%300==0){

        //Saving results at last iteration
        cout<<"Inside BoundaryLayer: "<<TimeElapsed_/SimulationTime_<<"  completed"<<endl; 
    }

    //Exact solution update
    sol.SolExact_Update(mesh, TimeElapsed_);

    //Save computed and exact solutions
   
    
       Save_Sol(FileOutputFormat, FileName,\
       sol, mesh, ite);
    
            
    //Save locally the time evolution of the solution

    Save_Sol_Local(FileName2,\
            sol, mesh, x_cell, TimeElapsed_);
     
    }

    //Save computed and exact solutions
   
    
    Save_Sol(FileOutputFormat, FileName,\
            sol, mesh, ite);
    sol.Compute_L1_err(mesh);
    

}


/**********************/
/*** FRACTIONAL-STEP **/
/**********************/

void Solver_Isen::BoundaryLayerFracStep(Sol_Isen& sol, Mesh& mesh){

    //Local Variables
    string MeshTag;
    string VariableType = sol.VariableType_;

    //Function

    double TimeStep       = TimeStep_;
    MeshTag = "Primal";

    //Pure time-relaxation update in a cell-loop
    /*
       BoundaryLayerUpdateFracStep(sol, mesh, MeshTag,\
       TimeStep, VariableType);
     */

    //Pure velocity time-relaxation update in a cell-loop
    BoundaryLayerUpdateFracStepVelocity(sol, mesh, MeshTag,\
            TimeStep, VariableType);

    //Pure pressure time-relaxation update in a cell-loop
    BoundaryLayerUpdateFracStepPressure(sol, mesh, MeshTag,\
            TimeStep, VariableType);
            

}

void Solver_Isen::BoundaryLayerUpdateFracStep(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double& dtRelax, string VariableType){

    //Local variables
    Vector5d W_state_cell;
    Vector5d H_state_cell;

    double dtRelax_face_ini = ZERO;
    Matrix5d Id;
    Id = Eigen::Matrix<double, 5, 5>::Identity();

    //Function

    //mesh inputs

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    //double tauMin = sol.etaRelax_(0);
    //double etaU   = sol.etaRelax_(1);
    //double etaP   = sol.etaRelax_(2);

    //relaxation inputs
    double tauP     = sol.tauRelax_(0);
    double tauU     = sol.tauRelax_(1);
    double pref     = sol.pRef_;
    double rhoref   = sol.mRef_;

    //Cubic relaxation 
    Vector5d STerm, STerm_L, STerm_R;
    Vector5d JacVector;
    Matrix5d JacMatrix, JacMatrix_L, JacMatrix_R;

    double NcellExt = mesh.Get_NcellExt();

    //Reference state
    Vector5d U_state_ref = NConsVarToConsVarLoc(sol.W_ref_, sol.SolTherm_);

    //Equilibrium state
    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    //Time matrix
    Matrix5d TimeMat = sol.TimeMat_;

    //Cubic relaxation
    Vector5d One = VectorXd::Constant(5, ONE);

    //Eigenvector base matrix
    Matrix5d RmatEigenVectors   = sol.EigenVectorBasis_;
   
    //Eigenvector base inverse matrix
    Matrix5d RmatEigenVectorsInv = sol.EigenVectorBasisInv_;

    //Time-step fixed
    dtRelax_face_ini   = dtRelax;

    //Loop on the real-cells
    for(int cell_id = 1; cell_id < NcellExt-1; cell_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            W_state_cell = sol.NConsVar_.row(cell_id).transpose();
        }
        else{

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);

        }

        if(VariableType=="ConsVar"){

            H_state_cell = NConsVarToConsVarLoc(W_state_cell, sol.SolTherm_);

                if(sol.SolType_=="Convection Relaxation Cubic BN"){

                NonLinearSpring(TimeMat, H_state_cell, U_state_eq, STerm,\
                        RmatEigenVectors, RmatEigenVectorsInv);

                JacVector = NonLinearSpringGradient(TimeMat, H_state_cell, U_state_eq,\
                        RmatEigenVectorsInv);

                if (TimeIntegrationType_=="Rosenbrock4"){

                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
                else if (TimeIntegrationType_=="ImplicitEuler1"){

                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
            }
            else if(sol.SolType_=="BN Relaxation"){

                IsenBNSourceTerm(tauU, tauP, pref, rhoref,\
                        H_state_cell, STerm,\
                        sol.SolTherm_);

                JacMatrix = IsenBNSourceTermGradient(tauU, tauP, pref, rhoref,\
                        H_state_cell,\
                        sol.SolTherm_\
                        );


                if(TimeIntegrationType_=="Rosenbrock4"){
                }

                else if(TimeIntegrationType_=="ImplicitEuler1"){

                    BN_ImplicitEuler(\
                            H_state_cell,\
                            STerm, JacMatrix, Id,\
                            dtRelax_face_ini\
                            );
                }
            }

            //Update of the staggered solution

            //PRIMAL mesh --> PRIMAL mesh update
            //Because we are in a fractional step strategy
            if (MeshTag =="Primal"){

                sol.ConsVar_.row(cell_id) = H_state_cell.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVar_.row(0)          = sol.ConsVar_.row(1);
                    sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

                }

            }
            //DUAL mesh --> PRIMAL update + Imposed boundary conditions
            else{

                cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
                exit(EXIT_FAILURE);

            }

        }
    }

    //Update of the variable arrays
    if(VariableType=="EqRelaxVar"){

        if(MeshTag=="Primal"){

            sol.EqRelaxToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(VariableType=="ConsVar"){

        if(MeshTag=="Primal"){

            sol.ConsVarToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }

}

void Solver_Isen::BoundaryLayerUpdateFracStepVelocity(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double& dtRelax, string VariableType){

    //Local variables
    Vector5d W_state_cell;
    Vector5d H_state_cell;

    double dtRelax_face_ini = ZERO;
    Matrix5d Id;
    Id = Eigen::Matrix<double, 5, 5>::Identity();

    //Function

    //mesh inputs

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    //double tauMin = sol.etaRelax_(0);
    //double etaU   = sol.etaRelax_(1);
    //double etaP   = sol.etaRelax_(2);

    //relaxation inputs
    double tauU     = sol.tauRelax_(1);
    //double pref     = sol.pRef_;
    //double rhoref   = sol.mRef_;

    //Cubic relaxation 
    Vector5d STerm, STerm_L, STerm_R;
    Vector5d JacVector;
    Matrix5d JacMatrix, JacMatrix_L, JacMatrix_R;

    double NcellExt = mesh.Get_NcellExt();

    //Reference state
    Vector5d U_state_ref = NConsVarToConsVarLoc(sol.W_ref_, sol.SolTherm_);

    //Equilibrium state
    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    //Time matrix
    Matrix5d TimeMat = sol.TimeMat_;

    //Cubic relaxation
    Vector5d One = VectorXd::Constant(5, ONE);

    //Eigenvector base matrix
    Matrix5d RmatEigenVectors   = sol.EigenVectorBasis_;
   
    //Eigenvector base inverse matrix
    Matrix5d RmatEigenVectorsInv = sol.EigenVectorBasisInv_;

    //Time-step fixed
    dtRelax_face_ini   = dtRelax;

    //Loop on the real-cells
    for(int cell_id = 1; cell_id < NcellExt-1; cell_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            W_state_cell = sol.NConsVar_.row(cell_id).transpose();
        }
        else{

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);

        }

        if(VariableType=="ConsVar"){

            H_state_cell = NConsVarToConsVarLoc(W_state_cell, sol.SolTherm_);

                if(sol.SolType_=="Convection Relaxation Cubic BN"){

                NonLinearSpring(TimeMat, H_state_cell, U_state_eq, STerm,\
                        RmatEigenVectors, RmatEigenVectorsInv);

                JacVector = NonLinearSpringGradient(TimeMat, H_state_cell, U_state_eq,\
                        RmatEigenVectorsInv);

                if (TimeIntegrationType_=="Rosenbrock4"){

                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
                else if (TimeIntegrationType_=="ImplicitEuler1"){

                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
            }
            else if(sol.SolType_=="BN Relaxation"){


                if(TimeIntegrationType_=="Rosenbrock4"){
                }

                else if(TimeIntegrationType_=="ImplicitEuler1"){

                    BN_ImplicitVelocity(\
                            H_state_cell,\
                            tauU, dtRelax_face_ini\
                            );
                }
            }

            //Update of the staggered solution

            //PRIMAL mesh --> PRIMAL mesh update
            //Because we are in a fractional step strategy
            if (MeshTag =="Primal"){

                sol.ConsVar_.row(cell_id) = H_state_cell.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVar_.row(0)          = sol.ConsVar_.row(1);
                    sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

                }

            }
            //DUAL mesh --> PRIMAL update + Imposed boundary conditions
            else{

                cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
                exit(EXIT_FAILURE);

            }

        }
    }

    //Update of the variable arrays
    if(VariableType=="EqRelaxVar"){

        if(MeshTag=="Primal"){

            sol.EqRelaxToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(VariableType=="ConsVar"){

        if(MeshTag=="Primal"){

            sol.ConsVarToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }

}

void Solver_Isen::BoundaryLayerUpdateFracStepPressure(Sol_Isen& sol, Mesh& mesh, string MeshTag,\
        double& dtRelax, string VariableType){

    //Local variables
    Vector5d W_state_cell;
    Vector5d H_state_cell;

    double dtRelax_face_ini = ZERO;
    Matrix5d Id;
    Id = Eigen::Matrix<double, 5, 5>::Identity();

    //Function

    //mesh inputs

    //sol inputs
    string LeftBCType  = sol.LeftBCType_;
    string RightBCType = sol.RightBCType_;

    //relaxation inputs
    //double tauMin = sol.etaRelax_(0);
    //double etaU   = sol.etaRelax_(1);
    //double etaP   = sol.etaRelax_(2);

    //relaxation inputs
    double tauP     = sol.tauRelax_(0);

    //Cubic relaxation 
    Vector5d STerm, STerm_L, STerm_R;
    Vector5d JacVector;
    Matrix5d JacMatrix, JacMatrix_L, JacMatrix_R;

    double NcellExt = mesh.Get_NcellExt();

    //Reference state
    Vector5d U_state_ref = NConsVarToConsVarLoc(sol.W_ref_, sol.SolTherm_);

    //Equilibrium state
    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    //Time matrix
    Matrix5d TimeMat = sol.TimeMat_;

    //Cubic relaxation
    Vector5d One = VectorXd::Constant(5, ONE);

    //Eigenvector base matrix
    Matrix5d RmatEigenVectors   = sol.EigenVectorBasis_;
   
    //Eigenvector base inverse matrix
    Matrix5d RmatEigenVectorsInv = sol.EigenVectorBasisInv_;

    //Time-step fixed
    dtRelax_face_ini   = dtRelax;

    //Loop on the real-cells
    for(int cell_id = 1; cell_id < NcellExt-1; cell_id++){

        //Boundary Layer Resolution performed on the faces of the PRIMAL mesh
        if (MeshTag =="Primal"){

            W_state_cell = sol.NConsVar_.row(cell_id).transpose();
        }
        else{

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);

        }

        if(VariableType=="ConsVar"){

            H_state_cell = NConsVarToConsVarLoc(W_state_cell, sol.SolTherm_);

                if(sol.SolType_=="Convection Relaxation Cubic BN"){

                NonLinearSpring(TimeMat, H_state_cell, U_state_eq, STerm,\
                        RmatEigenVectors, RmatEigenVectorsInv);

                JacVector = NonLinearSpringGradient(TimeMat, H_state_cell, U_state_eq,\
                        RmatEigenVectorsInv);

                if (TimeIntegrationType_=="Rosenbrock4"){

                    RosenBrockFourthOrder(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
                else if (TimeIntegrationType_=="ImplicitEuler1"){

                    ImplicitEuler(\
                            TimeMat,\
                            RmatEigenVectors, RmatEigenVectorsInv,\
                            H_state_cell, U_state_ref, U_state_eq,\
                            STerm, JacVector, One,\
                            dtRelax_face_ini\
                            );
                }
            }
            else if(sol.SolType_=="BN Relaxation"){


                if(TimeIntegrationType_=="Rosenbrock4"){
                }

                else if(TimeIntegrationType_=="ImplicitEuler1"){


                    BN_ImplicitPressure(\
                            H_state_cell, sol.SolTherm_,\
                            tauP, dtRelax_face_ini\
                            );

                    /*
                    BN_Pressure_Dichotomy(\
                            H_state_cell,\
                            tauP, dtRelax_face_ini,\
                            sol.SolTherm_,\
                            epsDicho);
                            */
                }
            }

            //Update of the staggered solution

            //PRIMAL mesh --> PRIMAL mesh update
            //Because we are in a fractional step strategy
            if (MeshTag =="Primal"){

                sol.ConsVar_.row(cell_id) = H_state_cell.transpose();

                if (LeftBCType=="transparent" && RightBCType=="transparent" ){

                    sol.ConsVar_.row(0)          = sol.ConsVar_.row(1);
                    sol.ConsVar_.row(NcellExt-1) = sol.ConsVar_.row(NcellExt-2);

                }

            }
            //DUAL mesh --> PRIMAL update + Imposed boundary conditions
            else{

                cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
                exit(EXIT_FAILURE);

            }

        }
    }

    //Update of the variable arrays
    if(VariableType=="EqRelaxVar"){

        if(MeshTag=="Primal"){

            sol.EqRelaxToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(VariableType=="ConsVar"){

        if(MeshTag=="Primal"){

            sol.ConsVarToNConsVar("Primal");
        }
        else if(MeshTag=="Dual"){

            cout<<"Inside BoundaryLayerUpdateFracStep: Dual mesh used!"<<endl;
            exit(EXIT_FAILURE);
        }
    }

}

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

void Solver_Isen::ConsVarFluxUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables
    int L,R;
    Vector5d W_state_L, W_state_R, W_state_avr;
    Vector5d ConsFlux     = VectorXd::Zero(5);
    Vector5d NConsVarFluxR = VectorXd::Zero(5);
    Vector5d NConsVarFluxL = VectorXd::Zero(5);
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
            ConsFlux      = VectorXd::Zero(5);
            NConsVarFluxR = VectorXd::Zero(5);
            NConsVarFluxL = VectorXd::Zero(5);

            ConsVarFluxUpdateLoc(\
                    ConsFlux,\
                    NConsVarFluxR, NConsVarFluxL,\
                    W_state_L, W_state_R,\
                    sol.SolTherm_,\
                    sol.JacConvFrozen_,\
                    sol.EigenvaluesFrozen_,\
                    SchemeTypeCons_\
                    );

            sol.ConsVarFlux_.row(face_id)   = ConsFlux.transpose();
            sol.NConsVarFluxR_.row(face_id) = NConsVarFluxR.transpose();
            sol.NConsVarFluxL_.row(face_id) = NConsVarFluxL.transpose();

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

        if(SchemeTypeCons_=="Rusanov"){

            sol.ConsVar_.row(cell_id) +=-(TimeStep_/SpaceStep)*(
                    sol.ConsVarFlux_.row(Rf) - sol.ConsVarFlux_.row(Lf) +\
                    sol.NConsVarFlux_.row(cell_id) ) +\
                TimeStep_*(sol.SourceTerms_.row(cell_id));
        }
        else if(SchemeTypeCons_=="HLLAC"){

            sol.ConsVar_.row(cell_id) +=-(TimeStep_/SpaceStep)*(
                    sol.ConsVarFlux_.row(Rf) - sol.ConsVarFlux_.row(Lf) +\
                    sol.NConsVarFluxL_.row(Rf) + sol.NConsVarFluxR_.row(Lf) ) +\
                TimeStep_*(sol.SourceTerms_.row(cell_id));
        }
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

void Solver_Isen::ConsVarUpdateLafitteRey(Sol_Isen& sol, Mesh& mesh){

    //Local variables

    //Function

    int Ncells       = mesh.Get_Ncells();
    int NcellExt     = mesh.Get_NcellExt();

    string  LeftBCType   =sol.LeftBCType_;
    string  RightBCType  =sol.RightBCType_;

    for(int cell_id = 1; cell_id <= Ncells; cell_id++){

        sol.ConsVar_.row(cell_id) += (TimeStep_-NRelax_*dtRelax_)*(\
                sol.ConsVar_.row(cell_id) - sol.ConsVarOld_.row(cell_id) )/dtRelax_;
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

void Solver_Isen::SourceTermsUpdate(Sol_Isen& sol, Mesh& mesh){

    //Local variables

    //Function

    int Ncells       = mesh.Get_Ncells();
    int NcellExt     = mesh.Get_NcellExt();

    string  LeftBCType   =sol.LeftBCType_;
    string  RightBCType  =sol.RightBCType_;

    Vector5d U_state;
    Vector5d STerm;

    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    for(int cell_id = 1; cell_id <= Ncells; cell_id++){

        U_state = sol.ConsVar_.row(cell_id).transpose();

        NonLinearSpring(sol.TimeMat_, U_state, U_state_eq, STerm,\
                sol.EigenVectorBasis_, sol.EigenVectorBasisInv_);

        sol.SourceTerms_.row(cell_id) = STerm.transpose(); 

    }

    //Treatment of the boundary conditions
    if(LeftBCType=="transparent"){
        sol.SourceTerms_.row(0) = sol.SourceTerms_.row(1);

    }
    if(RightBCType=="transparent"){
        sol.SourceTerms_.row(NcellExt-1) = sol.SourceTerms_.row(NcellExt-2);

    }
    else{
        sol.SourceTerms_.row(0)          = MatrixXd::Zero(1,5);
        sol.SourceTerms_.row(NcellExt-1) = MatrixXd::Zero(1,5);
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
    int ite(0);

    //CPU time variables
    clock_t t;

    //Saving solution
    Save_Sol(FileOutputFormat, FileName,\
            sol, mesh, ite);

    //Measuring the wall clock time for the simulation
    t = clock();

    //Definition of the relaxation time step
    double TauMin = sol.etaRelax_(0);
    dtRelax_ = CourantBL_*TauMin;

    //double TimeCounter = ONE;
    //double TimeSlope   = HUNDRED;

    while(TimeElapsed_ < SimulationTime_){
    //while(TimeCounter < THREE){

        //The fractional step strategy is triggered
        if(FractionalStep_==true){

            //Update of the discrete timestep for the convection part
            TimeStepUpdate(sol, mesh);

            //Time modulation
            /*
            if(TimeCounter < TimeSlope){

                cout<<"old TimeStep = "<<TimeStep_<<endl;
                TimeStep_   *= (TimeCounter/TimeSlope);
                TimeCounter +=ONE;
                cout<<"TimeCounter  = "<<TimeCounter<<endl;
                cout<<"new TimeStep = "<<TimeStep_<<endl;
            }
            */

            //Classical Update of the time step
            TimeElapsed_+=TimeStep_;

            if(SourceTermType_!="none"){

                //Boundary Layer Resolution part during: TimeStep_ = convection time scale
                BoundaryLayerFracStep(sol, mesh);

            }

            //Update of the conservative flux
            ConsVarFluxUpdate(sol, mesh);

            //Update of the non-conservative flux
            if(SchemeTypeCons_=="Rusanov"){

                NConsVarFluxUpdate(sol, mesh);

            }

            //Update of the conservative variables
            ConsVarUpdate(sol, mesh);

        }

        //The Bereux-sainsaulieu strategy is triggered
        else{

            //Lafitte Rey Technique
            TimeStepUpdate(sol, mesh);

            //Classical Update of the time step
            TimeElapsed_+=TimeStep_;

            //cout<<"TimeStep = "<<TimeStep_<<endl;

            //FIXME
            //WARNING after BoundaryLayer the
            dtRelax_ = TimeStep_;

            //Boundary Layer Resolution part during: NRelax_*dtRelax_
            BoundaryLayer(sol, mesh);

            //FIXME
            //Update of the conservative variables
            //ConsVarUpdateLafitteRey(sol, mesh);

            //Update of the source terms
            //FIXME
            //SourceTermsUpdate(sol, mesh);

            //Update of the discrete timestep for the convection part
            //FIXME
            //TimeStepUpdate(sol, mesh);

        }


        /*
        if(TimeElapsed_ > SimulationTime_){
            TimeStep_ = TimeElapsed_ - SimulationTime_;
            TimeElapsed_      = SimulationTime_;
        }
        */

        ite++;

       /* 
        cout<<"Inside Simulation: dtRelax = "<<dtRelax_<<", TimeStep = "<<TimeStep_<<\
        ", TimeElapsed = "<<TimeElapsed_<<endl;
        */
        

        //Update NConsVar_ using the matrix ConsVar_
        sol.ConsVarToNConsVar("Primal");

        //Update NConsVarEqRelax_ using the matrix NConsVar_
        sol.NConsVarToEqRelax("Primal");

        //Updates SoundSpeed_  of Sol using new values
        sol.SoundSpeed_Update();

        //Updates MachMax_ , MachMin_ using the new Sol values
        sol.Mach_Update();

        //Updating the exact solution
        sol.SolExact_Update(mesh, TimeElapsed_);

        if(ite%print_freq_==0){

            //Saving solution
            Save_Sol(FileOutputFormat, FileName,\
                    sol, mesh, ite);
        }
        if(ite%500==0){

            //Saving results at last iteration
            cout<<"Inside Simulation: "<<TimeElapsed_/SimulationTime_<<"  completed"<<endl; 
        }
    }

    //Measuring the CPU time
    t        = clock()-t;
    CPUTime_ = ((double) t)/ ((double)CLOCKS_PER_SEC);

    //Saving solution
    Save_Sol(FileOutputFormat, FileName,\
            sol, mesh, ite);

    //Computing the L1 errors and L1 exact norms at final time
    sol.Compute_L1_err(mesh);
}

void Solver_Isen::Save_Sol(string FileOutputFormat, string FileName,\
        Sol_Isen& sol, Mesh& mesh, int ite){

    int NcellExt = mesh.Get_NcellExt();
    int Ncells   = mesh.Get_Ncells();
    int Nfaces   = mesh.Get_Nfaces();

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
        FileName_Face_ite=Folder_location+FileName+"_Cell"+numstrNcells+"_FACE_ite"+numstrite+".vtk";
    }

    else if(FileOutputFormat==".dat"){
        FileName_ite=Folder_location+FileName+"_Cell"+numstrNcells+"_ite"+numstrite+".dat";
        FileName_Face_ite=Folder_location+FileName+"_Cell"+numstrNcells+"_FACE_ite"+numstrite+".dat";
    }

    ofstream file(FileName_ite.c_str(), ios::out);  //file flux declaration and file opening
    ofstream file_Face(FileName_Face_ite.c_str(), ios::out);  //file flux declaration and file opening

    if(file){  // If the opening is successful

        if(FileOutputFormat==".vtk"){

            file<<"# vtk DataFile Version 2.0"<<endl;
            file<<"Simulation Data"<<endl;
            file<<"ASCII"<<endl;
            file<<"DATASET STRUCTURED_GRID"<<endl;
            file<<"DIMENSIONS 1 1 1"<<endl;
            file<<"POINTS "<<NcellExt<<"  double"<<endl;


            /********* CELL-LOCATED OUTPUTS **********/

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

            file<<"SCALARS "<<"M1"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.ConsVar_(i,1)<<endl;

            }
            file<<endl;


            file<<"SCALARS "<<"M2"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.ConsVar_(i,3)<<endl;

            }
            file<<endl;

            file<<"SCALARS "<<"M"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.ConsVar_(i,1) + \
                    sol.ConsVar_(i,3)<<endl;

            }
            file<<endl;

            file<<"SCALARS "<<"MU"<<" double 1"<<endl;
            file<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<NcellExt;i++){

                file<<std::setprecision(COUT_PRECISION)<<std::scientific<<sol.ConsVar_(i,2) + \
                    sol.ConsVar_(i,4)<<endl;

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

            /********* FACE-LOCATED OUTPUTS **********/

            file_Face<<"# vtk DataFile Version 2.0"<<endl;
            file_Face<<"Simulation Data"<<endl;
            file_Face<<"ASCII"<<endl;
            file_Face<<"DATASET STRUCTURED_GRID"<<endl;
            file_Face<<"DIMENSIONS 1 1 1"<<endl;
            file_Face<<"POINTS "<<Nfaces<<"  double"<<endl;

            for (int i=0; i<Nfaces;i++){

                x=mesh.FaceIndex_(i,3);
                y=mesh.FaceIndex_(i,4);
                z=mesh.FaceIndex_(i,5);
                file_Face<<std::setprecision(COUT_PRECISION)<<std::scientific<<x<<" "<<y<<" "<<z<<endl;
            }

            file_Face<<endl;

            file_Face<<"POINT_DATA "<<Nfaces<<endl;

            /*********************************/
            /***** COMPUTED TIME STEPS  ******/
            /*********************************/

            file_Face<<"SCALARS "<<"Dt_Loc"<<" double 1"<<endl;
            file_Face<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<Nfaces;i++){

                file_Face<<std::setprecision(COUT_PRECISION)<<std::scientific<<TimeMatrix_(i,0)<<endl;

            }
            file_Face<<endl;

            file_Face<<"SCALARS "<<"DtRelax_Loc"<<" double 1"<<endl;
            file_Face<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<Nfaces;i++){

                file_Face<<std::setprecision(COUT_PRECISION)<<std::scientific<<TimeMatrix_(i,1)<<endl;

            }
            file_Face<<endl;

            file_Face<<"SCALARS "<<"DtConv_Loc"<<" double 1"<<endl;
            file_Face<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<Nfaces;i++){

                file_Face<<std::setprecision(COUT_PRECISION)<<std::scientific<<TimeMatrix_(i,0)<<endl;

            }
            file_Face<<endl;

            file_Face<<"SCALARS "<<"Dt_Min"<<" double 1"<<endl;
            file_Face<<"LOOKUP_TABLE default"<<endl;
            for (int i=0; i<Nfaces;i++){

                file_Face<<std::setprecision(COUT_PRECISION)<<std::scientific<<dtRelax_<<endl;

            }
            file_Face<<endl;
        }

        else if(FileOutputFormat==".dat"){

            //Writting the columns names:

            file<<"x"<<" "<<"y"<<" "<<"z"<<" "\
                <<"Alpha1"<<" "<<"P1"<<" "<<"U1"<<" "<<"P2"<<" "<<"U2"<<" "\
                <<"U"<<" "<<"P"<<" "<<"dU"<<" "<<"dP"<<" "\
                <<"M1"<<" "<<"M"<<" "<<"MU"<<" "\
                <<"Alpha1_ex"<<" "<<"P1_ex"<<" "<<"U1_ex"<<" "<<"P2_ex"<<" "\
                <<"U2_ex"<<" "\
                <<"U_ex"<<" "<<"P_ex"<<" "<<"dU_ex"<<" "<<"dP_ex"<<endl;

            //Writing values
            for (int i=0; i<NcellExt;i++){

                x=mesh.CellCoordsTab_(i,1);
                y=mesh.CellCoordsTab_(i,2);
                z=mesh.CellCoordsTab_(i,3);
                file<<std::setprecision(COUT_PRECISION)<<std::scientific\
                    <<x<<" "<<y<<" "<<z<<" "\
                    <<sol.NConsVar_(i,0)<<" "<<sol.NConsVar_(i,1)<<" "\
                    <<sol.NConsVar_(i,2)<<" "<<sol.NConsVar_(i,3)<<" "\
                    <<sol.NConsVar_(i,4)<<" "<<sol.NConsVarEqRelax_(i,1)<<" "\
                    <<sol.NConsVarEqRelax_(i,2)<<" "<<sol.NConsVarEqRelax_(i,3)<<" "\
                    <<sol.NConsVarEqRelax_(i,4)<<" "\
                    <<sol.ConsVar_(i,1)<<" "<<sol.ConsVar_(i,1) + sol.ConsVar_(i,3)<<" "\
                    <<sol.ConsVar_(i,2) + sol.ConsVar_(i,4)<<" "\
                    <<sol.SolExact_(i,0)<<" "<<sol.SolExact_(i,1)<<" "\
                    <<sol.SolExact_(i,2)<<" "<<sol.SolExact_(i,3)<<" "\
                    <<sol.SolExact_(i,4)<<" "<<sol.SolExactEqRelax_(i,1)<<" "\
                    <<sol.SolExactEqRelax_(i,2)<<" "<<sol.SolExactEqRelax_(i,3)<<" "\
                    <<sol.SolExactEqRelax_(i,4)<<endl;
            }

        }

        file.close();  // the file is closed
        file_Face.close();  // the file is closed

    }
    else  {cerr << "Error Save_Sol, Can't open file !" << endl;}

}

void Solver_Isen::Save_Sol_Local(string FileName,\
        Sol_Isen& sol, Mesh& mesh, double x_cell, double time){

    double SpaceStep = mesh.Get_SpaceStep();

    string Folder_location="./Output/";

    //Get the right precision for all the prints
    std::setprecision(COUT_PRECISION);

    string FileName_ite;
    FileName_ite=Folder_location+FileName;

    ofstream file(FileName_ite.c_str(), ios::app);  //file flux declaration and file opening

    if(file){  // If the opening is successful

        //Writing values
        int i = floor(x_cell/SpaceStep);
        file<<std::setprecision(COUT_PRECISION)<<std::scientific\
            <<time<<" "\
            <<sol.NConsVar_(i,0)<<" "<<sol.NConsVar_(i,1)<<" "\
            <<sol.NConsVar_(i,2)<<" "<<sol.NConsVar_(i,3)<<" "\
            <<sol.NConsVar_(i,4)<<" "<<sol.SolExact_(i,0)<<" "\
            <<sol.SolExact_(i,1)<<" "<<sol.SolExact_(i,2)<<" "\
            <<sol.SolExact_(i,3)<<" "<<sol.SolExact_(i,4)<<endl;

        file.close();  // the file is closed

    }
    else  {cerr << "Error Save_Sol_local, Can't open file !" << endl;}

}

//External Functions

/********** BOUNDARY LAYER RESOLUTION ************/

/*%%%%%%%%%%  EqRelax Variables %%%%%%%%%%*/

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

Vector2d  LocalCourant_LSTEq(Vector5d& W_state0, Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double tauMin, double etaP, double etaU,\
        double SpaceStep,\
        double Damp,\
        int face_id\
        ){

    /*
    //Local variables
    Vector2d TimeVector;
    double alpha1_avr, alpha2_avr, rho1_avr, rho2_avr, p1_avr, p2_avr;
    double m1_avr, m2_avr, m_avr, c1_avr, c2_avr, C1_avr, C2_avr;

    double u1_avr, u2_avr;

    Vector2d TimeLoc;

    //Getting the other non-conservative quantities for the averaged state

    //FIXME CUBIC RELAXATION
    double weight_L = ONE_OVER_TWO;
    Vector5d W_state_avr = NConsVarAveraged(W_state_L, W_state_R, weight_L);

    Vector5d U_state_L   = NConsVarToConsVarLoc(W_state_L, Therm);
    Vector5d U_state_R   = NConsVarToConsVarLoc(W_state_R, Therm);
    Vector5d U_state_avr = NConsVarAveraged(U_state_L, U_state_R, weight_L);

    //Reference state
    Vector5d W_ref;
    W_ref<<0.5,
           1.e5,
           ONE,
           3.e6,
           ONE;

    Vector5d U_state_ref = NConsVarToConsVarLoc(W_ref, Therm);

    //Equilibrium state
    Vector5d W_eq;
        W_eq<<0.5,
              3.e5,
              6.,
              4.e6,
              -0.5;

    Vector5d U_eq = NConsVarToConsVarLoc(W_eq, Therm);

    //Time relaxation exponential operator
    double mu_0 = ONE;
    double mu_1 = ONE;
    double mu_2 = ONE;
    double mu_3 = ONE;
    double mu_4 = ONE;

    Vector5d Relax_mu;
    Relax_mu <<mu_0,
               mu_1,
               mu_2,
               mu_3,
               mu_4;

    string VariableType = "ConsVar";

    //Eigenvector base inverse matrix
    Matrix5d RmatEigenVectorsInv = IsentropicBN_EigenvectorsInv(\
            VariableType,\
            W_state0,\
            Therm\
            );

    Vector5d a_ref = RmatEigenVectorsInv*(U_state_ref);
    Vector5d a_n   = RmatEigenVectorsInv*(U_state_avr);
    Vector5d a_inf = RmatEigenVectorsInv*(U_eq);

    //double CofacAlpha1 = ZERO;
    double a_min = ZERO;
    double a_max = ONE;

    //Time step constraint to garantee the alpha1 in [0, 1]

    //Case very close to the equilibrium state
    if ( fabs(a_n(0) - a_inf(0)) < epsZero){

        //CofacAlpha1 = Big;

    }

    //Condition to be in [0,..
    else if ( a_n(0) > a_inf(0)){

        //cout<<" a_n > a_inf --> [0,.."<<endl;
        //CofacAlpha1 = ((a_min - a_n(0))*pow(a_ref(0), TWO))/pow(a_inf(0) - a_n(0), THREE);
        //cout<<"CofacAlpha1 = "<<CofacAlpha1<<endl;

    }
    //Condition to be in ..,1]
    else{

        //cout<<" a_n < a_inf --> ..,1]"<<endl;
        //CofacAlpha1 = ((a_max - a_n(0))*pow(a_ref(0), TWO))/pow(a_inf(0) - a_n(0), THREE); 
        //cout<<"CofacAlpha1 = "<<CofacAlpha1<<endl;

    }

    //Function
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
 
    //FIXME
    double mu_max = ONE;
    //double mu_max = ONE/CofacAlpha1;
    lambda_du = mu_max;
    lambda_dp = mu_max;
    double dtRelax = min( ONE/lambda_du, ONE/lambda_dp  )*tauMin;

    double dtConv  = (SpaceStep)*(ONE/(max(max(fabs(u2_avr+c2_avr), fabs(u1_avr+c1_avr)), max(fabs(u2_avr-c2_avr), fabs(u1_avr-c1_avr)))));

    TimeLoc(0) = Damp*dtRelax;
    TimeLoc(1) = Damp*dtConv;

    */

    Vector2d TimeLoc;
    TimeLoc(0) =Big;
    TimeLoc(1) =Big;
    return TimeLoc;

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

    Vector5d V_state_L   = NConsVarToEqRelaxLoc(W_state_L, Therm); 
    Vector5d V_state_R   = NConsVarToEqRelaxLoc(W_state_R, Therm); 

    //Linearized source term 
    Matrix5d LinSouTermEq   = LinearizedSourceTermsEq(W_state_avr,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );

    Matrix5d LinSouTermOutEq = LinearizedSourceTermsOutOfEq(W_state_avr,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );

    Matrix5d LinSouTerm = LinSouTermEq + LinSouTermOutEq;
    //Matrix5d LinSouTerm = LinSouTermEq;

    //Constant Jacobian
    Matrix5d LinJacEqRelax = LinearizedJacobianEqRelax(W_state_avr,\
            Therm\
            );

    //FIXME
    LinJacEqRelax = MatrixXd::Zero(5,5);

    /*
    VectorXcd eivals = LinSouTerm.eigenvalues();
    cout << "The eigenvalues are:" << endl << eivals << endl;
    MatrixXcd Test = LinSouTerm.cast <complex<double>> ();
    ComplexEigenSolver<MatrixXcd> ces(Test);
    cout << "The eigenvectors are:" 
         << endl << ces.eigenvectors() << endl;

    */

    //Initialization
    Vector5d H_ini = H_state;
    Vector5d RHS_L = V_state_L;
    Vector5d RHS_R = V_state_R;

    string filename = "BS_Trajectories.dat";
    ofstream file((filename).c_str(), ios::out);
    file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<"Time"<<" "<<"Alpha1"<<" "<<"U"<<" "<<"P"<<" "<<"du"<<" "<<"dp"<<endl;

    time = ZERO;
    SaveBereuxSainsaulieu(H_ini, filename, time);

    for (int ite=1;ite<=NRelax;ite++){

        time = ite*dtRelax;

        /*
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
         */

        //Non homogeneous term of the left ODE
        RHS_L = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_L;

        //Non homogeneous term of the right ODE
        RHS_R = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_R;

        //Dynamics of H
        H_state =(Id + (dtRelax/tauMin)*LinSouTerm )*H_ini\
                 -(dtRelax/SpaceStep)*LinJacEqRelax*( RHS_R - RHS_L  );

        //Update of H_ini
        H_ini = H_state;
        SaveBereuxSainsaulieu(H_ini, filename, time);
    }

    file.close();
}

void BoundaryLayerResolutionLocalTest(\
        string VariableType,\
        Vector5d& H_state, Vector5d& W_state_L, Vector5d& W_state_R,\
        Matrix5d& LinJac_L, Matrix5d& LinJac_R,\
        Matrix5d& LinSouTerm, Vector5d& STerm,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin,\
        double SpaceStep\
        ){

    //Local variables
    Vector5d V_state_L, V_state_R;

    //Function
    //int nrows = H_state.rows();
    //Matrix5d Id = MatrixXd::Identity(nrows, nrows);

    //Averaged state

    /*
    Vector5d W_state_avr = NConsVarAveraged(\
            W_state_L, W_state_R,\
            ONE_OVER_TWO);
            */

    if(VariableType=="EqRelaxVar"){

        V_state_L   = NConsVarToEqRelaxLoc(W_state_L, Therm); 
        V_state_R   = NConsVarToEqRelaxLoc(W_state_R, Therm); 

    }
    else if(VariableType=="ConsVar"){

        V_state_L   = NConsVarToConsVarLoc(W_state_L, Therm); 
        V_state_R   = NConsVarToConsVarLoc(W_state_R, Therm); 

    }

    //Gradient of the source term
    /*
    Matrix5d LinSouTerm = LinearizedSourceTerm(W_state_avr,\
            VariableType,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );
            */

    //Source Term
    /*
    Vector5d STerm = SourceTerm(W_state_avr,\
            VariableType,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );
            */

    //FIXME
    //LinSouTerm = MatrixXd::Zero(nrows,nrows);
    //STerm      = VectorXd::Zero(nrows);

    //Jacobian matrix
    //FIXME
    /*
       Matrix5d LinJac = LinearizedJacobian(W_state_avr,\
       VariableType,\
       Therm\
       );

     */

    //Initialization
    //Vector5d H_ini = H_state;
    Vector5d RHS_L = -(V_state_R - V_state_L)/TWO;
    Vector5d RHS_R = (V_state_R - V_state_L)/TWO;

    //Non homogeneous term of the left ODE
    //RHS_L = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_L + (dtRelax/tauMin)*STerm;

    //Non homogeneous term of the right ODE
    //RHS_R = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_R + (dtRelax/tauMin)*STerm;

    //Dynamics of H
    /*
       H_state =(Id + (dtRelax/tauMin)*LinSouTerm )*H_ini-\
       (dtRelax/SpaceStep)*(LinJac_R*RHS_R - LinJac_L*RHS_L) +\
       (dtRelax/tauMin)*STerm;
     */

    /*
    cout<<"H_state BEFORE = "<<endl;
    cout<<H_state<<endl<<endl;
    */
    //cout<<"U_R - U_L = "<<RHS_R - RHS_L<<endl<<endl;

    H_state = H_state + LinSouTerm*(dtRelax*STerm) -\
              LinSouTerm*((dtRelax/SpaceStep)*LinJac_R)*LinSouTerm*(RHS_R - RHS_L);

    /*
    H_state = H_state + LinSouTerm*(dtRelax*STerm) -\
              LinSouTerm*((dtRelax/SpaceStep)*LinJac_R)*(RHS_R - RHS_L);
    */
    /*
    cout<<"H_state AFTER = "<<endl;
    cout<<H_state<<endl<<endl;
    */
}

void SourceTermResolutionLocalNewton(\
        Vector5d& H_state,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double CourantBL, double tauMin, double NRelax\
        ){

    //Local variables
    double time, dtRelax;
    Vector5d LSTEqODE_L, LSTEqODE_R;
    Matrix5d LinSouTerm, LinSouTermEq, LinSouTermOutEq;
    Vector5d RHS;

    //Function
    int nrows = H_state.rows();
    Matrix5d Id = MatrixXd::Identity(nrows, nrows);

    //Initialization
    Vector5d H_updated = VectorXd::Zero(nrows);
    Vector5d H_current = H_state;
    Vector5d W_current = EqRelaxToNConsVarLoc(H_current, Therm);

    string filename = "Newton_Source_Term_Trajectories.dat";
    ofstream file((filename).c_str(), ios::out);
    file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<"Time"<<" "<<"Alpha1"<<" "<<"U"<<" "<<"P"<<" "<<"du"<<" "<<"dp"<<endl;

    time = ZERO;
    SaveBereuxSainsaulieu(H_current, filename, time);

    for (int ite=1;ite<=NRelax;ite++){

        //FIXME
        /*
        dtRelax = LocalCourant_LSTEq(W_current, W_current, W_current,\
                Therm, pRef, mRef,\
                tauMin, etaP, etaU,\
                Big,\
                CourantBL, ite);
                */
        dtRelax = ZERO;

        time += dtRelax;

        //Linearized source term 
        LinSouTermEq   = LinearizedSourceTermsEq(W_current,\
                Therm,\
                pRef, mRef,\
                etaP, etaU\
                );

        LinSouTermOutEq = LinearizedSourceTermsOutOfEq(W_current,\
                Therm,\
                pRef, mRef,\
                etaP, etaU\
                );
        
        string VariableType = "EqRelaxVar";

        RHS = SourceTerm(H_current,\
                VariableType,\
                Therm,\
                pRef, mRef,\
                etaP, etaU\
                );
        //RHS = VectorXd::Zero(5);

        LinSouTerm = LinSouTermEq + LinSouTermOutEq;

        //Dynamics of H
        H_updated = (Id + (dtRelax/tauMin)*LinSouTerm )*H_current +\
                    (dtRelax/tauMin)*RHS;

        //Update of H_ini
        H_current = H_updated;
        W_current = EqRelaxToNConsVarLoc(H_current, Therm);

        SaveBereuxSainsaulieu(H_current, filename, time);
    }

    file.close();

}

void SaveBereuxSainsaulieu(Vector5d& H_state, string filename, double time){

//Local variables
double alpha1, U, P, du, dp;

//Function
ofstream file((filename).c_str(), ios::app);

if(file){

//Save H_state at initial condition
alpha1  = H_state(0);
U       = H_state(1);
P       = H_state(2);
du      = H_state(3);
dp      = H_state(4);

file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<time<<" "<<alpha1<<" "<<U<<" "<<P<<" "<<du<<" "<<dp<<endl;

}

else{cerr << "SaveBereuxSainsaulieu, Can't open file !" << endl;}

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

/*%%%%%%%%%%  Conservative Variables %%%%%%%%%%*/

void BoundaryLayerResolutionLocalConsVar(\
        Vector5d& H_state,Vector5d& U_state_L, Vector5d& U_state_R,\
        ThermoLaw& Therm, double pRef, double mRef,\
        double etaP, double etaU,\
        double dtRelax, double tauMin, double NRelax,\
        double SpaceStep\
        ){

    //Local variables

    //Function
    int nrows = H_state.rows();
    Matrix5d Id = MatrixXd::Identity(nrows, nrows);

    //Averaged state
    Vector5d U_state_avr = NConsVarAveraged(\
        U_state_L, U_state_R,\
        ONE_OVER_TWO);

    //Linearized source term
    Vector5d W_state_avr = ConsVarToNConsVarLoc(U_state_avr, Therm);

    Matrix5d LinSouTerm = LinearizedSourceTermsCons(W_state_avr,\
            Therm,\
            pRef, mRef,\
            etaP, etaU\
            );

    //Constant Jacobian
    Matrix5d LinJacCons = LinearizedJacobianCons(W_state_avr,\
            Therm\
            );

    //Initialization
    Vector5d H_ini = H_state;
    Vector5d RHS_L = U_state_L;
    Vector5d RHS_R = U_state_R;

        for (int ite=1;ite<=NRelax;ite++){

            //Non homogeneous term of the left ODE
            RHS_L = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_L;

            //Non homogeneous term of the right ODE
            RHS_R = (Id + (dtRelax/tauMin)*LinSouTerm)*RHS_R;

            //Dynamics of H
            H_state =(Id + (dtRelax/tauMin)*LinSouTerm )*H_ini\
                     -(dtRelax/SpaceStep)*LinJacCons*( RHS_R - RHS_L  );

            //Update of H_ini
            H_ini = H_state;

        }

}
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

/*****************************************************************************/
/************** HIGH ORDER IMPLICIT TIME INTEGRATION METHOD ******************/
/*****************************************************************************/

void TimeIntegration(Sol_Isen& sol, double SimulationTime, double dtRelax,\
        string FileNameInput, double CourantBL\
        ){

    double tauMin = sol.etaRelax_(0);

    cout<<"tauMin = "<<tauMin<<endl;

    Matrix5d TimeMat_;
    Vector5d W_eq_;
    Vector5d W_ref_;
    Matrix5d EigenVectorBasis_;
    Matrix5d EigenVectorBasisInv_;

    Matrix5d EigenVectorBasis   = sol.EigenVectorBasis_;

    //Eigenvector base inverse matrix
    Matrix5d EigenVectorBasisInv = sol.EigenVectorBasisInv_;

    //Reference state
    Vector5d U_state_ref = NConsVarToConsVarLoc(sol.W_ref_, sol.SolTherm_);

    //Equilibrium state
    Vector5d U_state_eq = NConsVarToConsVarLoc(sol.W_eq_, sol.SolTherm_);

    //Initial state
    Vector5d W_state_ini = sol.InitL_;
    Vector5d U_state_ini = NConsVarToConsVarLoc(sol.InitL_, sol.SolTherm_);

    //Exact state
    Vector5d W_state_exact = W_state_ini;
    Vector5d U_state_exact_ini = NConsVarToConsVarLoc(sol.InitL_, sol.SolTherm_);
    Vector5d U_state_exact = U_state_ini;

    //Time matrix
    double mu_0 = ONE;
    double mu_1 = ONE;
    double mu_2 = ONE;
    double mu_3 = ONE;
    double mu_4 = ONE;

    Vector5d Relax_mu;
    Relax_mu <<mu_0,
             mu_1,
             mu_2,
             mu_3,
             mu_4;

    Matrix5d TauMat = sol.TimeMat_;

    Vector5d TimeRefVector = Relax_mu.cwiseQuotient(U_state_ref.cwiseProduct(U_state_ref));
    TimeRefVector *=(ONE/tauMin);

    double TimeElapsed = ZERO;

    //Source term initialization
    Vector5d STerm_ini = VectorXd::Zero(5);
    Vector5d One       = VectorXd::Constant(5, ONE);
    NonLinearSpring(TauMat, U_state_ini, U_state_eq, STerm_ini,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Initialization of the gradient vector for the jacobian
    Vector5d JacVector_ini = NonLinearSpringGradient(\
            TauMat, U_state_ini, U_state_eq,\
            EigenVectorBasisInv);

    //File output initialization

    char numstrCFL[20]; // enough to hold all numbers up to 64-bits
    sprintf(numstrCFL, "%f", CourantBL);
    string FileName    = "TimeIntegration_Ros4_"+FileNameInput+"_CFL"+numstrCFL+".dat";
    string FileNameErr = "TimeIntegration_Ros4_"+FileNameInput+"Error.dat";

    string Folder_location="./Output/";
    FileName    = Folder_location+FileName;
    FileNameErr = Folder_location+FileNameErr;

    ofstream file(FileName.c_str(), ios::out);  //file flux declaration and file opening
    ofstream fileErr(FileNameErr.c_str(), ios::app);  //file flux declaration and file opening

    //Initialization of the time-steps
    double dtRelax_ini   = ZERO;

    //Initialization of the iteration counter
    int ite = 0;
    int TimeStepIteCounter = 0;

    file<<"Time"<<" "<<"Iteration"<<" "<<"SubIteration"<<" "\
        <<"dtRelax_ini"<<" "\
        <<"Alpha1"<<" "<<"P1"<<" "<<"U1"<<" "<<"P2"<<" "<<"U2"<<" "\
        <<"Alpha1_ex"<<" "<<"P1_ex"<<" "<<"U1_ex"<<" "<<"P2_ex"<<" "<<"U2_ex"<<endl;
    //FIXME
    
    file<<std::setprecision(COUT_PRECISION)<<std::scientific\
        <<TimeElapsed<<" "<<ite<<" "<<TimeStepIteCounter<<" "\
        <<dtRelax_ini<<" "\
        <<W_state_ini(0)<<" "<<W_state_ini(1)<<" "<<W_state_ini(2)<<" "\
        <<W_state_ini(3)<<" "<<W_state_ini(4)<<" "\
        <<W_state_exact(0)<<" "<<W_state_exact(1)<<" "<<W_state_exact(2)<<" "\
        <<W_state_exact(3)<<" "<<W_state_exact(4)<<endl;
    
    /*
       file<<std::setprecision(COUT_PRECISION)<<std::scientific\
       <<TimeElapsed<<" "<<ite<<" "<<TimeStepIteCounter<<" "\
       <<dtRelax_ini<<" "\
       <<U_state_ini(0)<<" "<<U_state_ini(1)<<" "<<U_state_ini(2)<<" "\
       <<U_state_ini(3)<<" "<<U_state_ini(4)<<" "\
       <<U_state_exact(0)<<" "<<U_state_exact(1)<<" "<<U_state_exact(2)<<" "\
       <<U_state_exact(3)<<" "<<U_state_exact(4)<<endl;
     */

    //Error Initialization

    Vector5d Error   = VectorXd::Zero(5);
    Vector5d NormRef = VectorXd::Zero(5);

    //Once the initial solution saved, update of time-step as
    //a first guess

    dtRelax_ini = dtRelax;
    //cout<<"dtRelax = "<<dtRelax<<endl;

    if(file){  // If the opening is successful

        //FIXME
        //while(ite < 1){
        while(TimeElapsed < SimulationTime){

            ite++;
            TimeStepIteCounter = 0;
            
            
               RosenBrockFourthOrder(\
               TauMat,\
               EigenVectorBasis, EigenVectorBasisInv,\
               U_state_ini, U_state_ref, U_state_eq,\
               STerm_ini, JacVector_ini, One,\
               dtRelax_ini\
               );
             
               /**
               
               ImplicitEuler(\
                       TauMat,\
                       EigenVectorBasis, EigenVectorBasisInv,\
                       U_state_ini, U_state_ref, U_state_eq,\
                       STerm_ini, JacVector_ini, One,\
                       dtRelax_ini\
                       );

                       */

               JacVector_ini = NonLinearSpringGradient(\
                       TauMat, U_state_ini, U_state_eq,\
                       EigenVectorBasisInv);

               /*
                  RungeKuttaFourthOrder(\
                  TauMat,\
                  EigenVectorBasis, EigenVectorBasisInv,\
                  U_state_ini, U_state_ref, U_state_eq,\
                  STerm_ini,\
                  dtRelax_ini\
                  );
                */

            //Update of STerm_ini
            NonLinearSpring(TauMat, U_state_ini, U_state_eq, STerm_ini,\
                    EigenVectorBasis, EigenVectorBasisInv);

            TimeElapsed += dtRelax_ini;

            //Error update
            //FIXME
            
            //Time-step for the next iteration
            //AFTER the error update
            //dtRelax_ini = dtRelax_estim;
            
            NonLinearSpringExactSolution(U_state_exact, U_state_exact_ini,\
                    U_state_eq, EigenVectorBasis, EigenVectorBasisInv,\
                    TimeRefVector, TimeElapsed);

            W_state_ini   = ConsVarToNConsVarLoc(U_state_ini, sol.SolTherm_);
            W_state_exact = ConsVarToNConsVarLoc(U_state_exact, sol.SolTherm_);

            for (int k =0; k<5;k++){
                Error(k)   += dtRelax_ini*fabs(W_state_ini(k) - W_state_exact(k));
                NormRef(k) += dtRelax_ini*fabs(sol.W_ref_(k));
            }
            
            /*
            for (int k =0; k<5;k++){
                Error(k)   += dtRelax_ini*fabs(U_state_ini(k) - U_state_exact(k));
                NormRef(k) += dtRelax_ini*fabs(U_state_ref(k));
            }
            */

            //FIXME
            
            file<<std::setprecision(COUT_PRECISION)<<std::scientific\
                <<TimeElapsed<<" "<<ite<<" "<<TimeStepIteCounter<<" "\
                <<dtRelax_ini<<" "\
                <<W_state_ini(0)<<" "<<W_state_ini(1)<<" "<<W_state_ini(2)<<" "\
                <<W_state_ini(3)<<" "<<W_state_ini(4)<<" "\
                <<W_state_exact(0)<<" "<<W_state_exact(1)<<" "<<W_state_exact(2)<<" "\
                <<W_state_exact(3)<<" "<<W_state_exact(4)<<endl;
            
            /*
            file<<std::setprecision(COUT_PRECISION)<<std::scientific\
                <<TimeElapsed<<" "<<ite<<" "<<TimeStepIteCounter<<" "\
                <<dtRelax_ini<<" "\
                <<U_state_ini(0)<<" "<<U_state_ini(1)<<" "<<U_state_ini(2)<<" "\
                <<U_state_ini(3)<<" "<<U_state_ini(4)<<" "\
                <<U_state_exact(0)<<" "<<U_state_exact(1)<<" "<<U_state_exact(2)<<" "\
                <<U_state_exact(3)<<" "<<U_state_exact(4)<<endl;
                */
        }

        file.close();

        //Normalizing error
        //Saving error

        /*
        for (int k =0; k<5;k++){
            Error(k)   += fabs(W_state_ini(k) - W_state_exact(k));
            NormRef(k) += fabs(W_ref(k));
        }
        */

        for (int k =0; k<5;k++){
            Error(k)   /= NormRef(k);
        }

        fileErr<<std::setprecision(COUT_PRECISION)<<std::scientific\
            <<dtRelax<<" "<<Error(0)<<" "<<Error(1)<<" "<<Error(2)<<" "\
            <<Error(3)<<" "<<Error(4)<<endl;
    }
    else  {cerr << "Error TimeIntegration, Can't open file !" << endl;}

}


void TimeIntegrationLoc(Vector5d& U_state_ini,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Matrix5d& TauMat,\
        Vector5d& U_state_ref, Vector5d& U_state_eq,\
        double dtRelax, double tauMin){

    double TimeElapsed = ZERO;

    //Source term initialization
    Vector5d STerm_ini = VectorXd::Zero(5);
    NonLinearSpring(TauMat, U_state_ini, U_state_eq, STerm_ini,\
            EigenVectorBasis, EigenVectorBasisInv);
    Vector5d One       = VectorXd::Constant(5, ONE);

    //Initialization of the gradient vector for the jacobian
    Vector5d JacVector_ini = NonLinearSpringGradient(\
            TauMat, U_state_ini, U_state_eq,\
            EigenVectorBasisInv);

    //Initialization of the time-steps
    double dtRelax_ini   = dtRelax;

    //Initialization of the iteration counter
    int ite = 0;

        while(TimeElapsed < dtRelax){

            ite++;
            RosenBrockFourthOrder(\
                    TauMat,\
                    EigenVectorBasis, EigenVectorBasisInv,\
                    U_state_ini, U_state_ref, U_state_eq,\
                    STerm_ini, JacVector_ini, One,\
                    dtRelax_ini\
                    );
            
            TimeElapsed += dtRelax_ini;
        }
}

//Explicit Runge-Kutta order 4 method
void RungeKuttaFourthOrder(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini,\
        double& dtRelax_ini\
        ){

    //EXPLICIT-RK COEFFICIENTS FOR ORDER 4 METHOD

    const double A21 = ONE_OVER_TWO;
    const double A32 = ONE_OVER_TWO;
    const double A43 = ONE;

    //ORDER4
    const double B1 = 1./6., B2 = 1./3., B3 = 1./3., B4 = 1./6.;
    //ORDER1
    //const double B1 = 1.;

    //Saving initial state
    Vector5d U_state_sav = U_state_ini;

    //local variables
    Vector5d k1, k2, k3, k4;

    //Deriving k1 vector:
    k1 = STerm_ini;

    //Updating the state using the new slope k1
    //U2
    U_state_ini = U_state_sav + dtRelax_ini*A21*k1;

    //Deriving k2 vector:
    NonLinearSpring(TauMat, U_state_ini, U_eq, k2,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Updating the state using the new slopes k2
    //U3
    U_state_ini = U_state_sav + dtRelax_ini*A32*k2;

    //Deriving k3 vector:
    NonLinearSpring(TauMat, U_state_ini, U_eq, k3,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Updating the state using the new slopes k3
    //U4
    U_state_ini = U_state_sav + dtRelax_ini*A43*k3;

    //Deriving k4 vector:
    NonLinearSpring(TauMat, U_state_ini, U_eq, k4,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Updating U_state_ini
    //ORDER4
    U_state_ini  = U_state_sav + dtRelax_ini*(B1*k1 + B2*k2 + B3*k3 + B4*k4);
    //ORDER1
    //U_state_ini  = U_state_sav + dtRelax_ini*B1*k1;

    return;
}

//Euler Implicit of order 1 method
void ImplicitEuler(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini, Vector5d& JacVector_ini, Vector5d& One,\
        double& dtRelax_ini\
        ){

    //local variables
    Vector5d A_coords_inv;
    Matrix5d DistanceDiagMat; 
    Vector5d g1;

    //Computing [Id/(Gam*dtRelax) - Jac_ini]^-1
    A_coords_inv = ((One/(dtRelax_ini) + JacVector_ini).array()).pow(-ONE);
    DistanceDiagMat = EigenVectorBasis*DiagonalMatrix<double, 5>(A_coords_inv)*EigenVectorBasisInv; 

    //Deriving g1 vector:
    g1 = DistanceDiagMat*STerm_ini;

    U_state_ini  = U_state_ini + g1;

    return;

}

//Euler Implicit of order 1 method for Baer-Nunziato 
void BN_ImplicitEuler(\
        Vector5d& U_state_ini,\
        Vector5d& STerm_ini, Matrix5d& JacMatrix_ini, Matrix5d& Id,\
        double& dtRelax_ini\
        ){

    //local variables
    Matrix5d A;
    Vector5d g1;

    //Computing [Id/(Gam*dtRelax) - Jac_ini]^-1
    A = (ONE/(dtRelax_ini))*Id - JacMatrix_ini;

    //Defining the solver based on a QR decomposition
    Eigen::ColPivHouseholderQR<Matrix5d> QRsolver(A);
    //Deriving g1 vector:
    g1 = QRsolver.solve(STerm_ini);

    U_state_ini  = U_state_ini + g1;

    return;

}

//Implicit time integration of order 1 method for the Velocity Baer-Nunziato 
void BN_ImplicitVelocity(\
        Vector5d& U_state_ini,\
        double tauU, double& dtRelax_ini\
        ){

    //Local variables
    double m1    = U_state_ini(1);
    double m1u1  = U_state_ini(2);
    double m2    = U_state_ini(3);
    double m2u2  = U_state_ini(4);
    double m     = m1 + m2;
    double Y1    = m1/m;
    double Y2    = m2/m;
    double tadim = dtRelax_ini/tauU;

    double TimeDenom = ONE/(ONE + tadim);
    double w1        = Y1*tadim*TimeDenom; 
    double w2        = Y2*tadim*TimeDenom; 

    double m1u1_up = (TimeDenom + w1)*m1u1 + w1*m2u2;
    double m2u2_up = w2*m1u1 + (TimeDenom + w2)*m2u2;

    U_state_ini(2) = m1u1_up;
    U_state_ini(4) = m2u2_up;

    return;

}

//Implicit time integration of order 1 method for the Pressure Baer-Nunziato 
void BN_ImplicitPressure(\
        Vector5d& U_state_ini, ThermoLaw& Therm,\
        double tauP, double& dtRelax_ini\
        ){

    //Local variables
    Vector5d W_state_ini = ConsVarToNConsVarLoc(U_state_ini, Therm); 
    double alpha1 = W_state_ini(0);
    double alpha2 = ONE - alpha1;
    double p1     = W_state_ini(1);
    double p2     = W_state_ini(3);
    double dp     = p2 - p1;
    double P      = p2 + p1;
    double rho1   = Density_EOS(1, Therm, p1, ZERO);
    double rho2   = Density_EOS(2, Therm, p2, ZERO);
    double c1     = Sound_Speed_EOS(1, Therm, rho1, p1);
    double c2     = Sound_Speed_EOS(1, Therm, rho2, p2);
    double C1     = rho1*c1*c1;
    double C2     = rho2*c2*c2;
    double tadim  = dtRelax_ini/tauP;
    double noneq  = (ONE - TWO*alpha1)*(dp/P);
    double Pt     = (alpha1*C2 + alpha2*C1) - (alpha1*C2 - alpha2*C1)*(dp/P);

    // Gallouet et al. alpha1 update
    /*
       double alpha1_up = (alpha1/alpha2)*exp(-tadim*(dp/P))/(ONE +(alpha1/alpha2)*\
       exp(-tadim*(dp/P))); 
     */

    //BS alpha1 update
    double alpha1_up = alpha1 - tadim*(alpha1*alpha2)*(dp/P)/\
                       ( ONE + tadim*(noneq + Pt/P)); 

    if( fabs(alpha1_up - ONE) < 0.0001*epsZero || alpha1_up < 0.0001*epsZero){

        cout<<"Inside BN_ImplicitPressure: alpha1 vanishing phase"<<endl;
        exit(EXIT_FAILURE);
    }
    /*
       double alpha2_up = ONE - alpha1_up;

    double alpha = exp( -tadim*( (alpha1*C2 + alpha2*C1)/P ) );
    double beta  = exp(  tadim*(alpha2*C1/p1 - alpha1*C2/p2)*(dp/P) );
    double beta  = pow(alpha1_up/alpha1, - C1/p1)*\
                   pow(alpha2_up/alpha2, - C2/p2);

    double p1_up = ONE_OVER_TWO*( -(dp*alpha) + sqrt( (dp*alpha*dp*alpha) \
                + FOUR*beta*p1*p2 ) );

    double p2_up = (p1*p2)*beta/p1_up;
    */

    //W_state_ini update
    //FIXME
    //W_state_ini(0) = alpha1_up;
    U_state_ini(0) = alpha1_up;
    //W_state_ini(1) = p1_up;
    //W_state_ini(3) = p2_up;

    //U_state_ini = NConsVarToConsVarLoc(W_state_ini, Therm);

    return;

}

//Simultaneous implicit time integration of order 1 for the Velocity-Pressure Baer-Nunziato
void BN_ImplicitVelocityPressure(\
        Vector5d& U_state_ini, ThermoLaw& Therm,\
        double tauU, double tauP, double& dtRelax_ini\
        ){

    //Pressure relaxation
    //Vector5d W_state_ini = ConsVarToNConsVarLoc(U_state_ini, Therm); 
    //double alpha1 = W_state_ini(0);
    //double alpha2 = ONE - alpha1;
    //double p1     = W_state_ini(1);
    //double p2     = W_state_ini(3);
    //double dp     = p2 - p1;
    //double P      = p2 + p1;
    //double rho1   = Density_EOS(1, Therm, p1, ZERO);
    //double rho2   = Density_EOS(2, Therm, p2, ZERO);
    //double c1     = Sound_Speed_EOS(1, Therm, rho1, p1);
    //double c2     = Sound_Speed_EOS(1, Therm, rho2, p2);
    //double C1     = rho1*c1*c1;
    //double C2     = rho2*c2*c2;

    double tadim  = dtRelax_ini/tauP;
    //double noneq  = (ONE - TWO*alpha1)*(dp/P);
    //double Pt     = (alpha1*C2 + alpha2*C1) - (alpha1*C2 - alpha2*C1)*(dp/P);
    //noneq = ZERO;

    //Velocity relaxation
    double m1   = U_state_ini(1);
    double m1u1 = U_state_ini(2);
    double m2   = U_state_ini(3);
    double m2u2 = U_state_ini(4);
    double m    = m1 + m2;
    double Y1   = m1/m;
    double Y2   = m2/m;

    double TimeDenom = ONE/(ONE + tadim);
    double w1        = Y1*tadim*TimeDenom; 
    double w2        = Y2*tadim*TimeDenom; 

    double m1u1_up = (TimeDenom + w1)*m1u1 + w1*m2u2;
    double m2u2_up = w2*m1u1 + (TimeDenom + w2)*m2u2;

    //m1u1 m2u2 updates
    U_state_ini(2) = m1u1_up;
    U_state_ini(4) = m2u2_up;

    //FIXME
    //Pressure relaxation Bereux-Sainsaulieu
    
    /* 
       double alpha1_up = alpha1 - tadim*(alpha1*alpha2)*(dp/P)/\
       ( ONE + tadim*(noneq + Pt/P)); 
     */
                       
    /*
    BN_Pressure_Dichotomy(\
            U_state_ini,\
            tauP, dtRelax_ini,\
            Therm,\
            epsDicho);
            */
            

    //U_state_ini update
    //U_state_ini(0) = alpha1_up;

    return;

}

double  AlphaRelaxFunction(\
        double alpha1,\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        ThermoLaw& Therm\
        ){

    //local variables

    //State at time t^n
    double alpha1_n = U_state_ini(0);
    double alpha2_n = ONE - alpha1_n;
    double m1_n     = U_state_ini(1);
    double m2_n     = U_state_ini(3);
    double p1_n     = Pressure_EOS(1, Therm, m1_n/alpha1_n, ZERO);
    double p2_n     = Pressure_EOS(2, Therm, m2_n/alpha2_n, ZERO);
    double cofac    = (dtRelax_ini/tauP)*(alpha1_n*alpha2_n/(p2_n + p1_n)); 

    //New pressure
    double p1  = Pressure_EOS(1, Therm, m1_n/alpha1, ZERO);
    double p2  = Pressure_EOS(2, Therm, m2_n/(ONE - alpha1), ZERO);
    double dp  = p2 - p1;

    return alpha1 + cofac*dp - alpha1_n;
    
}

double  AlphaRelaxFunctionBS(\
        double& alpha1,\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        double SpaceStep,\
        ThermoLaw& Therm,\
        double alpha1_mass_n\
        ){

    //local variables

    //State at time t^n
    double alpha1_n = U_state_ini(0);
    double alpha2_n = ONE - alpha1_n;
    double m1_n     = U_state_ini(1);
    double m2_n     = U_state_ini(3);
    double p1_n     = Pressure_EOS(1, Therm, m1_n/alpha1_n, ZERO);
    double p2_n     = Pressure_EOS(2, Therm, m2_n/alpha2_n, ZERO);
    double cofac    = (dtRelax_ini/tauP)*(alpha1_n*alpha2_n/(p2_n + p1_n)); 

    //New pressure
    double p1  = Pressure_EOS(1, Therm, m1_n/alpha1, ZERO);
    double p2  = Pressure_EOS(2, Therm, m2_n/(ONE - alpha1), ZERO);
    double dp  = p2 - p1;

    //cout<<"Alpha1 = "<<alpha1<<endl;
    /*
       cout<<"Inside AlphaRelaxFunctionBS: dp            = "<<dp<<endl;
       cout<<"Inside AlphaRelaxFunctionBS: frozen alpha1 = "<<alpha1_n + (dtRelax_ini/SpaceStep)*alpha1_mass_n<<endl;
       cout<<"Inside AlphaRelaxFunctionBS: cofac*dp      = "<<cofac*dp<<endl;
       cout<<"Inside AlphaRelaxFunctionBS: all           = "<<alpha1 + cofac*dp - alpha1_n - (dtRelax_ini/SpaceStep)*alpha1_mass_n<<endl;
     */

    return alpha1 + cofac*dp - alpha1_n - (dtRelax_ini/SpaceStep)*alpha1_mass_n;
    
}

double AlphaDichotomy(\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        ThermoLaw& Therm,\
        double& alpha1_inf, double& alpha1_sup,\
        double& eps){

    //Local variables
    double alpha1_star = ZERO;
    double alpha1_mid  = (alpha1_inf + alpha1_sup)/TWO;

    double Finf = AlphaRelaxFunction(\
            alpha1_inf,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            Therm\
            );

    double Fmid = AlphaRelaxFunction(\
            alpha1_mid,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            Therm\
            );

    //cout<<"Inside AlphaDichotomy: Fmid = "<<Fmid<<endl;

    double Fsup = AlphaRelaxFunction(\
            alpha1_sup,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            Therm\
            );
           
    //Function

    if(Finf*Fsup>= ZERO){

	cout<<"Alert AlphaDichotomy: Finf Fsup of same sign !"<<endl;
	exit(EXIT_FAILURE);

    }
    else if(fabs(Fmid)<=eps){

        return alpha1_mid;

    }
    else{

        if(Finf*Fmid<ZERO){

            alpha1_star = AlphaDichotomy(\
                    U_state_ini,\
                    tauP, dtRelax_ini,\
                    Therm,\
                    alpha1_inf, alpha1_mid,\
                    eps);
        }
        else {

            alpha1_star = AlphaDichotomy(\
                    U_state_ini,\
                    tauP, dtRelax_ini,\
                    Therm,\
                    alpha1_mid, alpha1_sup,\
                    eps);
        }

    }

    return alpha1_star;
}

double AlphaDichotomyBS(\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        double SpaceStep,\
        ThermoLaw& Therm,\
        double alpha1_inf, double alpha1_sup,\
        double alpha1_mass_n,\
        double eps){

    //Local variables
    double alpha1_star = ZERO;
    double Finf        = ZERO;
    double Fmid        = ZERO;
    double Fsup        = ZERO;
    double alpha1_mid  = (alpha1_inf + alpha1_sup)/TWO;

    /*
    cout<<"Dichotomy BS: alpha1_inf = "<<alpha1_inf<<", alpha1_mid = "<<alpha1_mid<<", alpha1_sup = "<<alpha1_sup<<endl;
*/

    //cout<<"Finf: "<<endl;
    Finf = AlphaRelaxFunctionBS(\
            alpha1_inf,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            SpaceStep,\
            Therm,\
            alpha1_mass_n\
            );

    //cout<<"Fmid: "<<endl;
    Fmid = AlphaRelaxFunctionBS(\
            alpha1_mid,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            SpaceStep,\
            Therm,\
            alpha1_mass_n\
            );

    //cout<<"Fsup: "<<endl;
    Fsup = AlphaRelaxFunctionBS(\
            alpha1_sup,\
            U_state_ini,\
            tauP, dtRelax_ini,\
            SpaceStep,\
            Therm,\
            alpha1_mass_n\
            );
    //cout<<"Inside AlphaDichotomy: Fmid = "<<Fmid<<endl;

           
    //Function

    if(fabs(Fmid)<=eps){

        return alpha1_mid;

    }
    else if(Finf*Fsup>= ZERO){

        cout<<"Alert AlphaDichotomy BS: Finf Fsup of same sign !"<<endl;
        exit(EXIT_FAILURE);

    }
    else{

        if(Finf*Fmid<ZERO){

            alpha1_star = AlphaDichotomyBS(\
                    U_state_ini,\
                    tauP, dtRelax_ini,\
                    SpaceStep,\
                    Therm,\
                    alpha1_inf, alpha1_mid,\
                    alpha1_mass_n,\
                    eps);
        }
        else {

            alpha1_star = AlphaDichotomyBS(\
                    U_state_ini,\
                    tauP, dtRelax_ini,\
                    SpaceStep,\
                    Therm,\
                    alpha1_mid, alpha1_sup,\
                    alpha1_mass_n,\
                    eps);
        }
    }

    return alpha1_star;
}

void BN_Pressure_Dichotomy(\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        ThermoLaw& Therm,\
        double& eps){

    double alpha1_inf_ini = eps;
    double alpha1_sup_ini = ONE - eps;

    //Dichotomy new volume fraction
    double alpha1_ini = AlphaDichotomy(\
            U_state_ini,\
            tauP, dtRelax_ini,\
            Therm,\
            alpha1_inf_ini, alpha1_sup_ini,\
            eps);

    //cout<<"Inside Pressure_Dichotomy: alpha1_new = "<<alpha1_ini<<endl;

    //Update
    U_state_ini(0) = alpha1_ini;
}

void BN_Pressure_DichotomyBS(\
        Vector5d& U_state_ini,\
        double tauP, double dtRelax_ini,\
        double SpaceStep,\
        ThermoLaw& Therm,\
        double alpha1_mass_n,\
        double& eps){

    double alpha1_inf_ini = eps;
    double alpha1_sup_ini = ONE - eps;

    //Dichotomy new volume fraction
    double alpha1_ini = AlphaDichotomyBS(\
            U_state_ini,\
            tauP, dtRelax_ini,\
            SpaceStep,\
            Therm,\
            alpha1_inf_ini, alpha1_sup_ini,\
            alpha1_mass_n,\
            eps);

    //cout<<"Inside Pressure_Dichotomy: alpha1_new = "<<alpha1_ini<<endl;

    //Update
    U_state_ini(0) = alpha1_ini;
}

//Rosenbrock of order 4 method
void RosenBrockFourthOrder(\
        Matrix5d& TauMat,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& U_state_ini, Vector5d& U_state_ref, Vector5d& U_eq,\
        Vector5d& STerm_ini, Vector5d& JacVector_ini, Vector5d& One,\
        double& dtRelax_ini\
        ){

    //SHAMPINE COEFFICIENTS FOR THE ROSENBROCK ORDER 4 METHOD
    //FIXME
    const double GAM = ONE_OVER_TWO;
    //const double GAM = ONE;

    //Rmk: the fourth line of coefficients is not provided because
    //at the fourth sub-step, f(y) used is f(y) of the third sub-step
    const double A21 = TWO;
    const double A31 = 48./25., A32 = 6./25.;

    //FIXME
    //B1 = ONE
    const double B1 = 19./9., B2 = ONE_OVER_TWO, B3 = 25./108., B4 = 125./108.;

    const double C21 = -8.;
    const double C31 = 372./25., C32 = 12./5.;
    const double C41 = -112./125., C42 = -54./125., C43 = -2./5.;

    //Saving initial state
    Vector5d U_state_sav = U_state_ini;
    Vector5d STerm_sav   = STerm_ini;

    //local variables
    Vector5d A_coords_inv;
    Matrix5d DistanceDiagMat; 
    Vector5d g1, g2, g3, g4, g5;

    //Computing [Id/(Gam*dtRelax) - Jac_ini]^-1
    A_coords_inv = ((One/(GAM*dtRelax_ini) + JacVector_ini).array()).pow(-ONE);
    DistanceDiagMat = EigenVectorBasis*DiagonalMatrix<double, 5>(A_coords_inv)*EigenVectorBasisInv; 

    //Deriving g1 vector:
    g1 = DistanceDiagMat*STerm_sav;

    //Updating the state using the new slope g1
    U_state_ini = U_state_sav + A21*g1;

    //Updating STerm 
    NonLinearSpring(TauMat, U_state_ini, U_eq, STerm_ini,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Deriving g2 vector:
    g2 = DistanceDiagMat*(STerm_ini + C21*g1/dtRelax_ini);

    //Updating the state using the new slopes g1, g2
    U_state_ini = U_state_sav + A31*g1 + A32*g2;

    //Updating STerm 
    NonLinearSpring(TauMat, U_state_ini, U_eq, STerm_ini,\
            EigenVectorBasis, EigenVectorBasisInv);

    //Deriving g3 vector:
    g3 = DistanceDiagMat*(STerm_ini + C31*g1/dtRelax_ini + C32*g2/dtRelax_ini);

    //No updates of U for the 4-th sub-step

    //Deriving g4 vector:
    g4 = DistanceDiagMat*\
         (STerm_ini + C41*g1/dtRelax_ini + C42*g2/dtRelax_ini + C43*g3/dtRelax_ini);

    //Updating U_state_ini
    //ORDER4
    U_state_ini  = U_state_sav + B1*g1 + B2*g2 + B3*g3 + B4*g4;

}

//BN Rosenbrock of order 4 method
void BN_RosenBrockFourthOrder(\
        double tauP, double tauU,\
        ThermoLaw& Therm,\
        Vector5d& U_state_ini,\
        double& dtRelax_ini\
        ){

    //SHAMPINE COEFFICIENTS FOR THE ROSENBROCK ORDER 4 METHOD
    //FIXME
    const double GAM = ONE_OVER_TWO;
    //const double GAM = ONE;

    //Rmk: the fourth line of coefficients is not provided because
    //at the fourth sub-step, f(y) used is f(y) of the third sub-step
    const double A21 = TWO;
    const double A31 = 48./25., A32 = 6./25.;

    //FIXME
    //B1 = ONE
    const double B1 = 19./9., B2 = ONE_OVER_TWO, B3 = 25./108., B4 = 125./108.;

    const double C21 = -8.;
    const double C31 = 372./25., C32 = 12./5.;
    const double C41 = -112./125., C42 = -54./125., C43 = -2./5.;

    //Saving initial state
    Vector5d U_state_sav = U_state_ini;

    //Initializing source terms
    Vector5d STerm_ini = VectorXd::Zero(5);
    BN_SourceTerm(tauU, tauP,\
            U_state_ini, STerm_ini,\
            Therm\
            );
    Vector5d STerm_sav   = STerm_ini;

    //local variables
    VectorXd g1 = VectorXd::Zero(5);
    VectorXd g2 = VectorXd::Zero(5);
    VectorXd g3 = VectorXd::Zero(5);
    VectorXd g4 = VectorXd::Zero(5);
    VectorXd g5 = VectorXd::Zero(5);

    Vector5d W_state_ini = ConsVarToNConsVarLoc(U_state_ini, Therm); 
    double alpha1 = W_state_ini(0);
    double alpha2 = ONE - alpha1;
    double p1     = W_state_ini(1);
    double p2     = W_state_ini(3);
    double dp     = p2 - p1;
    double P      = p2 + p1;
    double rho1   = Density_EOS(1, Therm, p1, ZERO);
    double rho2   = Density_EOS(2, Therm, p2, ZERO);
    double c1     = Sound_Speed_EOS(1, Therm, rho1, p1);
    double c2     = Sound_Speed_EOS(1, Therm, rho2, p2);
    double C1     = rho1*c1*c1;
    double C2     = rho2*c2*c2;

    double noneq  = (ONE - TWO*alpha1)*(dp/P);
    double Pt     = (alpha1*C2 + alpha2*C1) - (alpha1*C2 - alpha2*C1)*(dp/P);
    //noneq = ZERO;

    double alpha_relax = -(ONE/tauP)*(noneq + Pt/P); 
    double P_relax = (GAM*dtRelax_ini)/(ONE  - (GAM*dtRelax_ini*alpha_relax));
    double U_relax = (GAM*dtRelax_ini)/(ONE + GAM*dtRelax_ini/tauU);

    //Function

    //Deriving g1 vector:
    g1(0) = P_relax*STerm_sav(0); 
    g1(2) = U_relax*STerm_sav(2);
    g1(4) = U_relax*STerm_sav(4);

    //Updating the state using the new slope g1
    U_state_ini = U_state_sav + A21*g1;

    //Updating STerm 
    BN_SourceTerm(tauU, tauP,\
            U_state_ini, STerm_ini,\
            Therm\
            );

    //Deriving g2 vector:
    g2(0) = P_relax*(STerm_ini(0) + C21*g1(0)/dtRelax_ini); 
    g2(2) = U_relax*(STerm_ini(2) + C21*g1(2)/dtRelax_ini);
    g2(4) = U_relax*(STerm_ini(4) + C21*g1(4)/dtRelax_ini);

    //Updating the state using the new slopes g1, g2
    U_state_ini = U_state_sav + A31*g1 + A32*g2;

    //Updating STerm 
    BN_SourceTerm(tauU, tauP,\
            U_state_ini, STerm_ini,\
            Therm\
            );

    //Deriving g3 vector:
    g3(0) = P_relax*(STerm_ini(0) + C31*g1(0)/dtRelax_ini +\
            C32*g2(0)/dtRelax_ini); 
    g3(2) = U_relax*(STerm_ini(2) + C31*g1(2)/dtRelax_ini +\
            C32*g2(2)/dtRelax_ini);
    g3(4) = U_relax*(STerm_ini(4) + C31*g1(4)/dtRelax_ini +\
            C32*g2(4)/dtRelax_ini);

    //No updates of U for the 4-th sub-step

    //Deriving g4 vector:
    g4(0) = P_relax*(STerm_ini(0) + C41*g1(0)/dtRelax_ini +\
            C42*g2(0)/dtRelax_ini + C43*g3(0)/dtRelax_ini); 
    g4(2) = P_relax*(STerm_ini(2) + C41*g1(2)/dtRelax_ini +\
            C42*g2(2)/dtRelax_ini + C43*g3(2)/dtRelax_ini); 
    g4(4) = P_relax*(STerm_ini(4) + C41*g1(4)/dtRelax_ini +\
            C42*g2(4)/dtRelax_ini + C43*g3(4)/dtRelax_ini); 

    //Updating U_state_ini
    //ORDER4
    U_state_ini  = U_state_sav + B1*g1 + B2*g2 + B3*g3 + B4*g4;

}


void LinearSpring(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq, Vector5d& STerm,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv){

    Vector5d A_coords = EigenVectorBasisInv*(U_state - U_eq);

    STerm = -EigenVectorBasis*TauMat*A_coords;

}

Vector5d LinearSpringGradient(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq,\
        Matrix5d& EigenVectorBasisInv){

    Vector5d One;
    One<<ONE,
         ONE,
         ONE,
         ONE,
         ONE;
    Vector5d Grad   = TauMat*One;

    return Grad;

}

//Update the source term corresponding to a non-linear cubic spring
void NonLinearSpring(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq, Vector5d& STerm,\
        Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv){

    Vector5d A_coords = (EigenVectorBasisInv*(U_state - U_eq)).array().pow(THREE);

    STerm = -EigenVectorBasis*TauMat*A_coords;

}

//Returns the diagonal vector of the source term jacobian related to a non-linear cubic spring
Vector5d NonLinearSpringGradient(Matrix5d& TauMat, Vector5d& U_state, Vector5d& U_eq,\
        Matrix5d& EigenVectorBasisInv){

    Vector5d Coords = THREE*((EigenVectorBasisInv*(U_state - U_eq)).array().pow(TWO));
    Vector5d Grad   = TauMat*Coords;

    return Grad;

}

void NonLinearSpringExactSolution(Vector5d& U_state_exact, Vector5d& U_state_ini,\
        Vector5d& U_state_eq, Matrix5d& EigenVectorBasis, Matrix5d& EigenVectorBasisInv,\
        Vector5d& TimeRefVector, double time){

    Vector5d A_ini  = EigenVectorBasisInv*U_state_ini;
    Vector5d A_eq   = EigenVectorBasisInv*U_state_eq;
    Vector5d A_sol;

    int nVariables = 5;
    double exact_ak;

    for (int k = 0; k< nVariables; k++){

        exact_ak = ONE/((A_ini(k) - A_eq(k))*(A_ini(k) - A_eq(k))) + TWO*TimeRefVector(k)*time;

        if(A_ini(k) <= A_eq(k)){

            A_sol(k) = A_eq(k) - pow(exact_ak, -ONE_OVER_TWO);
        }
        else{

            A_sol(k) = A_eq(k) + pow(exact_ak, -ONE_OVER_TWO);
        }

    }

    U_state_exact = EigenVectorBasis*A_sol;

}

//Derive the real source term of the isentropic BN system
void IsenBNSourceTerm(double tauU, double tauP, double pref, double rhoref,\
        Vector5d& U_state, Vector5d& STerm,\
        ThermoLaw& Therm
        ){

    Vector5d W_state = ConsVarToNConsVarLoc(U_state, Therm);

    double alpha1 = W_state(0);
    double p1     = W_state(1);
    double u1     = W_state(2);
    double p2     = W_state(3);
    double u2     = W_state(4);

    double alpha2 = ONE - alpha1;

    //Source term: Alpha 1 
    STerm(0) = - (alpha1*alpha2/tauP)*(p2 - p1)/pref; 
    STerm(1) = ZERO; 
    STerm(2) =  (alpha1*alpha2/tauU)*(u2 - u1)*rhoref; 
    STerm(3) = ZERO; 
    STerm(4) = - (alpha1*alpha2/tauU)*(u2 - u1)*rhoref; 
}

void BN_SourceTerm(double tauU, double tauP,\
        Vector5d& U_state, Vector5d& STerm,\
        ThermoLaw& Therm\
        ){

    Vector5d W_state = ConsVarToNConsVarLoc(U_state, Therm);

    double alpha1 = W_state(0);
    double alpha2 = ONE - alpha1;
    double m1     = U_state(1);
    double p1     = W_state(1);
    double u1     = W_state(2);
    double m2     = U_state(3);
    double p2     = W_state(3);
    double u2     = W_state(4);

    double pref   = p1 + p2;
    double rhoref = m1*m2/(m1+m2);

    //Source term: Alpha 1 
    STerm(0) = - (alpha1*alpha2/tauP)*(p2 - p1)/pref; 
    STerm(1) = ZERO; 
    STerm(2) =  (alpha1*alpha2/tauU)*(u2 - u1)*rhoref; 
    STerm(3) = ZERO; 
    STerm(4) = - (alpha1*alpha2/tauU)*(u2 - u1)*rhoref; 
}

Matrix5d IsenBNSourceTermGradient(double tauU, double tauP, double pref, double rhoref,\
        Vector5d& U_state,\
        ThermoLaw& Therm
        ){

    double m1 = U_state(1);
    double m2 = U_state(3);

    Vector5d W_state = ConsVarToNConsVarLoc(U_state, Therm);

    double alpha1 = W_state(0);
    double p1     = W_state(1);
    double u1     = W_state(2);
    double p2     = W_state(3);
    double u2     = W_state(4);
    double du     = u2 - u1;
    double dp     = p2 - p1;

    double alpha2 = ONE - alpha1;
    double rho1   = Density_EOS(1, Therm, p1, ZERO);
    double rho2   = Density_EOS(2, Therm, p2, ZERO);
    double c1     = Sound_Speed_EOS(1, Therm, rho1, p1);
    double c2     = Sound_Speed_EOS(1, Therm, rho2, p2);
    double C1     = rho1*c1*c1;
    double C2     = rho2*c2*c2;
    double Prel   = alpha1*C2 + alpha2*C1;

    double K  = alpha1*alpha2;
    double Kt = (ONE - TWO*alpha1);

    Matrix5d BN_Source_Gradient = MatrixXd::Zero(5,5);

    //Row Alpha1
    RowVector5d RowAlpha1;
    RowAlpha1<<-(Kt/tauP)*(dp/pref) -(ONE/tauP)*(Prel/pref),\
        (K/(pref*tauP))*(C1/m1), ZERO, -(K/(pref*tauP))*(C2/m2), ZERO;

    //Row U
    RowVector5d RowU1;
    RowU1<< (Kt/tauU)*rhoref*du,\
        (K/tauU)*(rhoref*u1/m1),-(K/tauU)*(rhoref/m1),\
       -(K/tauU)*(rhoref*u2/m2), (K/tauU)*(rhoref/m2);

    BN_Source_Gradient.row(0) = RowAlpha1;
    BN_Source_Gradient.row(2) = RowU1;
    BN_Source_Gradient.row(4) = -RowU1;

    return BN_Source_Gradient;

}

/************************************************/
/*****************  CONVECTION  *****************/
/************************************************/

void ConsVarFluxUpdateLoc(\
        Vector5d& ConsFlux,\
        Vector5d& NConsVarFluxR, Vector5d& NConsVarFluxL,\
        Vector5d& W_state_L, Vector5d&  W_state_R,\
        ThermoLaw& Therm,\
        Matrix5d& JacConvFrozen,\
        double EigenvaluesFrozen,\
        string SchemeTypeCons\
        ){

    //Local Variables
    Vector5d Flux_L, Flux_R;
    Vector5d U_state_L = NConsVarToConsVarLoc(W_state_L, Therm);
    Vector5d U_state_R = NConsVarToConsVarLoc(W_state_R, Therm);
    
    //Function
    if(SchemeTypeCons=="Rusanov"){

        Flux_L    = ConsVarFlux(W_state_L, Therm); 
        Flux_R    = ConsVarFlux(W_state_R, Therm); 

        double lambda = SpectralRadiusRusanov(\
                W_state_L, W_state_R,\
                Therm);

        ConsFlux = ONE_OVER_TWO*(Flux_L + Flux_R) -\
                   (ONE_OVER_TWO*lambda)*(U_state_R - U_state_L);

    }
    else if(SchemeTypeCons=="Frozen Rusanov"){

        ConsFlux = ONE_OVER_TWO*JacConvFrozen*(U_state_L + U_state_R) -\
                   (ONE_OVER_TWO*EigenvaluesFrozen)*(U_state_R - U_state_L);
    }
    else if(SchemeTypeCons=="HLLAC"){

        Vector4d WaveSpeeds =  WaveSpeedEstimate(W_state_L, W_state_R, Therm);
        Vector6d Star_UP =  U_P_star_states(W_state_L, W_state_R, WaveSpeeds, Therm);
        ConsFlux = HLLAC_Flux(W_state_L, W_state_R, WaveSpeeds, Star_UP, Therm);
        NConsVarFluxR = NConsVarFluxUpdateLocR(W_state_L, W_state_R, WaveSpeeds, Star_UP, Therm);
        NConsVarFluxL = NConsVarFluxUpdateLocL(W_state_L, W_state_R, WaveSpeeds, Star_UP, Therm);

    }
    else{

        ConsFlux<<ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO;

    }

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

Vector4d WaveSpeedEstimate(Vector5d& W_state_L, Vector5d& W_state_R,\
                ThermoLaw& Therm){

    //Local variables
    double p1_L   = W_state_L(1); 
    double p1_R   = W_state_R(1); 
    double rho1_L = Density_EOS(1, Therm, p1_L, ZERO); 
    double rho1_R = Density_EOS(1, Therm, p1_R, ZERO); 
    double c1_L   = Sound_Speed_EOS(1, Therm, rho1_L, p1_L); 
    double c1_R   = Sound_Speed_EOS(1, Therm, rho1_R, p1_R); 
    double u1_L   = W_state_L(2); 
    double u1_R   = W_state_R(2); 

    double p2_L   = W_state_L(3); 
    double p2_R   = W_state_R(3); 
    double rho2_L = Density_EOS(2, Therm, p2_L, ZERO); 
    double rho2_R = Density_EOS(2, Therm, p2_R, ZERO); 
    double c2_L   = Sound_Speed_EOS(2, Therm, rho2_L, p2_L); 
    double c2_R   = Sound_Speed_EOS(2, Therm, rho2_R, p2_R); 
    double u2_L   = W_state_L(4); 
    double u2_R   = W_state_R(4);

    //Function
    double u1_avr = (rho1_R/(rho1_R + rho1_L))*u1_R + (rho1_L/(rho1_R + rho1_L))*u1_L;
    double c1_avr = (rho1_R/(rho1_R + rho1_L))*c1_R + (rho1_L/(rho1_R + rho1_L))*c1_L;

    double u2_avr = (rho2_R/(rho2_R + rho2_L))*u2_R + (rho2_L/(rho2_R + rho2_L))*u2_L;
    double c2_avr = (rho2_R/(rho2_R + rho2_L))*c2_R + (rho2_L/(rho2_R + rho2_L))*c2_L;

    double lambda1_L     = u1_L - c1_L; 
    double lambda1_L_avr = u1_avr - c1_avr; 
    double lambda1_R     = u1_R + c1_R; 
    double lambda1_R_avr = u1_avr + c1_avr; 

    double lambda2_L     = u2_L - c2_L; 
    double lambda2_L_avr = u2_avr - c2_avr; 
    double lambda2_R     = u2_R + c2_R; 
    double lambda2_R_avr = u2_avr + c2_avr;

    Vector4d WaveSpeeds;
    WaveSpeeds(0) = min(lambda1_L, lambda1_L_avr);
    WaveSpeeds(1) = max(lambda1_R, lambda1_R_avr);
    WaveSpeeds(2) = min(lambda2_L, lambda2_L_avr);
    WaveSpeeds(3) = max(lambda2_R, lambda2_R_avr);

    return WaveSpeeds;

}

/*
Vector5d U_P_star_states(Vector5d& W_state_L, Vector5d& W_state_R,\
                         Vector4d& WaveSpeeds,\
                         ThermoLaw& Therm){

    //Local variables
    double alpha1_L  = W_state_L(0);
    double alpha1_R  = W_state_R(0);
    double p1_L      = W_state_L(1); 
    double p1_R      = W_state_R(1); 
    double rho1_L    = Density_EOS(1, Therm, p1_L, ZERO); 
    double rho1_R    = Density_EOS(1, Therm, p1_R, ZERO); 
    double u1_L      = W_state_L(2); 
    double u1_R      = W_state_R(2); 

    double alpha2_L  = ONE - alpha1_L;
    double alpha2_R  = ONE - alpha1_R;
    double p2_L      = W_state_L(3); 
    double p2_R      = W_state_R(3); 
    double rho2_L    = Density_EOS(2, Therm, p2_L, ZERO); 
    double rho2_R    = Density_EOS(2, Therm, p2_R, ZERO); 
    double u2_L      = W_state_L(4); 
    double u2_R      = W_state_R(4);

    double S1_L = WaveSpeeds(0);
    double S1_R = WaveSpeeds(1);
    double S2_L = WaveSpeeds(2);
    double S2_R = WaveSpeeds(3);

    //Function
    double mq1_R =  rho1_R*(S1_R - u1_R);
    double mq1_L = -rho1_L*(S1_L - u1_L);
    double mq1   = mq1_R + mq1_L;

    double u1_star = (mq1_R*u1_R + mq1_L*u1_L - (p1_R - p1_L) )/mq1;
    double p1_star = p1_R + mq1_R*(u1_star - u1_R);

    double mq2_R =  alpha2_R*rho2_R*(S2_R - u2_R);
    double mq2_L = -alpha2_L*rho2_L*(S2_L - u2_L);
    double mq2   = mq2_R + mq2_L;
    double u2_star = (mq2_R*u2_R + mq2_L*u2_L - (alpha2_R*p2_R -alpha2_L*p2_L) \
                                              - (alpha1_R - alpha1_L)*p1_star)/mq2;

    double p2_star_L = p2_L - (mq2_L/alpha2_L)*(u2_star - u2_L);
    double p2_star_R = p2_R + (mq2_R/alpha2_R)*(u2_star - u2_R);

    Vector5d Star_State;
    Star_State(0) = u1_star;
    Star_State(1) = p1_star;
    Star_State(2) = u2_star;
    Star_State(3) = p2_star_L;
    Star_State(4) = p2_star_R;

    return Star_State;
}
*/

Vector6d U_P_star_states(Vector5d& W_state_L, Vector5d& W_state_R,\
                         Vector4d& WaveSpeeds,\
                         ThermoLaw& Therm){

    //Local variables
    double alpha1_L  = W_state_L(0);
    double alpha1_R  = W_state_R(0);
    double p1_L      = W_state_L(1); 
    double p1_R      = W_state_R(1); 
    double rho1_L    = Density_EOS(1, Therm, p1_L, ZERO); 
    double rho1_R    = Density_EOS(1, Therm, p1_R, ZERO); 
    double u1_L      = W_state_L(2); 
    double u1_R      = W_state_R(2); 

    double alpha2_L  = ONE - alpha1_L;
    double alpha2_R  = ONE - alpha1_R;
    double p2_L      = W_state_L(3); 
    double p2_R      = W_state_R(3); 
    double rho2_L    = Density_EOS(2, Therm, p2_L, ZERO); 
    double rho2_R    = Density_EOS(2, Therm, p2_R, ZERO); 
    double u2_L      = W_state_L(4); 
    double u2_R      = W_state_R(4);

    double S1_L = WaveSpeeds(0);
    double S1_R = WaveSpeeds(1);
    double S2_L = WaveSpeeds(2);
    double S2_R = WaveSpeeds(3);

    //Function
    double dAlpha = alpha1_R - alpha1_L;
    double q1_R   =  rho1_R*(S1_R - u1_R);
    double q1_L   = -rho1_L*(S1_L - u1_L);

    double q2_R =  rho2_R*(S2_R - u2_R);
    double q2_L = -rho2_L*(S2_L - u2_L);

    Matrix2d Mat;
    Mat(0,0) = (alpha1_R/q1_R) + (alpha1_L/q1_L);
    Mat(0,1) = -dAlpha;
    Mat(1,0) = +dAlpha;
    Mat(1,1) = (alpha2_R*q2_R) + (alpha2_L*q2_L);

    Vector2d RHS;
    RHS(0) = alpha1_L*(u1_L + p1_L/q1_L) - alpha1_R*(u1_R - p1_R/q1_R);
    RHS(1) = alpha2_L*(u2_L*q2_L + p2_L) - alpha2_R*(-u2_R*q2_R + p2_R);

    //Defining the solver based on a QR decomposition
    Eigen::ColPivHouseholderQR<Matrix2d> QRsolver(Mat);

    //Deriving  vector:  [p1_star, u2_star]^T
    Vector2d StarState;
    StarState = QRsolver.solve(RHS);

    double p1_star = StarState(0);
    double u2_star = StarState(1);
    double u1_star_L = u1_L - (p1_star - p1_L)/q1_L;
    double u1_star_R = u1_R + (p1_star - p1_R)/q1_R;
    double p2_star_L = p2_L - q2_L*(u2_star - u2_L);
    double p2_star_R = p2_R + q2_R*(u2_star - u2_R);

    Vector6d Star_State;
    Star_State(0) = u1_star_L;
    Star_State(1) = u1_star_R;
    Star_State(2) = p1_star;
    Star_State(3) = u2_star;
    Star_State(4) = p2_star_L;
    Star_State(5) = p2_star_R;

    return Star_State;
}

Vector5d HLLAC_Flux(Vector5d& W_state_L, Vector5d& W_state_R,\
                    Vector4d& WaveSpeeds, Vector6d& Star_UP,\
                         ThermoLaw& Therm){

    //Local variables
    double alpha1_L  = W_state_L(0);
    double alpha1_R  = W_state_R(0);
    double p1_L      = W_state_L(1); 
    double p1_R      = W_state_R(1); 
    double rho1_L    = Density_EOS(1, Therm, p1_L, ZERO); 
    double rho1_R    = Density_EOS(1, Therm, p1_R, ZERO); 
    double m1_L      = alpha1_L*rho1_L;
    double m1_R      = alpha1_R*rho1_R;
    double u1_L      = W_state_L(2); 
    double u1_R      = W_state_R(2); 

    double alpha2_L  = ONE - alpha1_L;
    double alpha2_R  = ONE - alpha1_R;
    double p2_L      = W_state_L(3); 
    double p2_R      = W_state_R(3); 
    double rho2_L    = Density_EOS(2, Therm, p2_L, ZERO); 
    double rho2_R    = Density_EOS(2, Therm, p2_R, ZERO); 
    double m2_L      = alpha2_L*rho2_L;
    double m2_R      = alpha2_R*rho2_R;
    double u2_L      = W_state_L(4); 
    double u2_R      = W_state_R(4);

    double S1_L = WaveSpeeds(0);
    double S1_R = WaveSpeeds(1);
    double S2_L = WaveSpeeds(2);
    double S2_R = WaveSpeeds(3);

    double u1_star_L = Star_UP(0);
    double u1_star_R = Star_UP(1);
    double u2_star   = Star_UP(3);

    //FIXME
    double p1_star     = Star_UP(2);
    double p2_star_L   = Star_UP(4);
    double p2_star_R   = Star_UP(5);
    double rho1_star   = Density_EOS(1, Therm, p1_star, ZERO); 
    double rho2_star_L = Density_EOS(2, Therm, p2_star_L, ZERO); 
    double rho2_star_R = Density_EOS(2, Therm, p2_star_R, ZERO); 

    double mk   = ZERO;
    double mkuk = ZERO;

    //Function
    Vector5d Flux   = VectorXd::Zero(5);
    Vector5d Flux_L = ConsVarFlux(W_state_L, Therm); 
    Vector5d Flux_R = ConsVarFlux(W_state_R, Therm); 

    //Subsonic flow test
    if( u2_star >= S1_R || u2_star <= S1_L){

        cout<<"Inside HLLAC Flux: Supersonic configuration!"<<endl;
        //exit(EXIT_FAILURE);
    }

    //Phase1
    if (S1_L > ZERO){

        Flux(1) = Flux_L(1);
        Flux(2) = Flux_L(2);
    }
    else if (ZERO < u2_star){

        //FIXME
        //mk = m1_L*(u1_L - S1_L)/(u1_star_L - S1_L);
        mk = alpha1_L*rho1_star;
        mkuk = mk*u1_star_L;
        Flux(1) = Flux_L(1) + S1_L*(mk - m1_L);
        Flux(2) = Flux_L(2) + S1_L*(mkuk - m1_L*u1_L);
    }
    else if (ZERO < S1_R){

        //FIXME
        //mk = m1_R*(u1_R - S1_R)/(u1_star_R - S1_R);
        mk = alpha1_R*rho1_star;
        mkuk = mk*u1_star_R;
        Flux(1) = Flux_R(1) + S1_R*(mk - m1_R);
        Flux(2) = Flux_R(2) + S1_R*(mkuk - m1_R*u1_R);
    }
    else{

        Flux(1) = Flux_R(1);
        Flux(2) = Flux_R(2);
    }
    //Phase2
    if (S2_L > ZERO){

        Flux(3) = Flux_L(3);
        Flux(4) = Flux_L(4);
    }
    else if (ZERO < u2_star){

        //FIXME
        //mk = m2_L*(u2_L - S2_L)/(u2_star - S2_L);
        mk   = alpha2_L*rho2_star_L;
        mkuk = mk*u2_star;
        Flux(3) = Flux_L(3) + S2_L*(mk - m2_L);
        Flux(4) = Flux_L(4) + S2_L*(mkuk - m2_L*u2_L);
    }
    else if (ZERO < S2_R){

        //FIXME
        //mk = m2_R*(u2_R - S2_R)/(u2_star - S2_R);
        mk   = alpha2_R*rho2_star_R; 
        mkuk = mk*u2_star;
        Flux(3) = Flux_R(3) + S2_R*(mk - m2_R);
        Flux(4) = Flux_R(4) + S2_R*(mkuk - m2_R*u2_R);
    }
    else{

        Flux(3) = Flux_R(3);
        Flux(4) = Flux_R(4);
    }

    return Flux;

}

Vector5d NConsVarFluxUpdateLocR(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        Vector4d& WaveSpeeds, Vector6d& Star_UP,\
        ThermoLaw& Therm\
        ){

    //Local variables
    double alpha1_L  = W_state_L(0);
    double alpha1_R  = W_state_R(0);

    double alpha2_L  = ONE - alpha1_L;
    double alpha2_R  = ONE - alpha1_R;

    double u2_star   = Star_UP(3);
    double p2_star_L = Star_UP(4);
    double p2_star_R = Star_UP(5);

    //Function
    Vector5d NConsVarFluxR = VectorXd::Zero(5);

    if( u2_star > ZERO){

        NConsVarFluxR(0) = u2_star*(alpha1_R - alpha1_L);
        NConsVarFluxR(2) = (alpha2_R*p2_star_R - alpha2_L*p2_star_L);
        NConsVarFluxR(4) = -(alpha2_R*p2_star_R - alpha2_L*p2_star_L);

    }

    return NConsVarFluxR;

}

Vector5d NConsVarFluxUpdateLocL(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        Vector4d& WaveSpeeds, Vector6d& Star_UP,\
        ThermoLaw& Therm\
        ){

    //Local variables
    double alpha1_L  = W_state_L(0);
    double alpha1_R  = W_state_R(0);

    double alpha2_L  = ONE - alpha1_L;
    double alpha2_R  = ONE - alpha1_R;

    double u2_star   = Star_UP(3);
    double p2_star_L = Star_UP(4);
    double p2_star_R = Star_UP(5);

    //Function
    Vector5d NConsVarFluxL = VectorXd::Zero(5);

    if( u2_star <= ZERO){

        NConsVarFluxL(0) = u2_star*(alpha1_R - alpha1_L);
        NConsVarFluxL(2) = (alpha2_R*p2_star_R - alpha2_L*p2_star_L);
        NConsVarFluxL(4) = -(alpha2_R*p2_star_R - alpha2_L*p2_star_L);

    }

    return NConsVarFluxL;

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

void Convergence_Curve(\
        double Length, int NGhostCells, \
        string LeftBCType, string RightBCType, \
        ThermoLaw& therm,\
        string SolType,\
        string VariableType,\
        Vector5d& InitL, Vector5d& InitR,\
        int Nbr_Areas, double x0,\
        double SimulationTime, \
        bool CFL_ramp, int CFL_ramp_range,\
        double CourantBL, double CourantConv,\
        bool FractionalStep,\
        string SourceTermType,\
        string TimeIntegrationType,\
        int NRelax,\
        double pRef, double mRef,\
        Vector3d& tauRelax,\
        string SchemeTypeCons, string SchemeTypeNCons,\
        Array<int, Eigen::Dynamic, 1>& CellsTab, \
        string FileOutputFormat, string CV_curve, int print_freq){

    cout<<"Convergence Curve: Case "<<CV_curve<<endl;
    cout<<"Simulation Time: "<<SimulationTime<<endl;

    int Npoints=CellsTab.size();
    int NbVariables = 5;

    double Order;
    VectorXd X(Npoints);
    //5 columns: err_alpha1, err_p1, err_u1, err_p2, err_u2
    MatrixXd ErrMatrix(Npoints,5);
    double Xa, Xb, Ya, Yb;

    Vector5d normL1_Var_exact_ref;

    //Get the right precision for all the prints
    std::setprecision(COUT_PRECISION);

    //Plotting the convergence curve
    ofstream file((CV_curve+"CV_curve.dat").c_str(), ios::app);  //file flux declaration and file opening

    //Calculating the order of convergence
    ofstream file_Order((CV_curve+"CV_order.dat").c_str(), ios::app);  //file flux declaration and file opening

    //CPU effiency curve
    ofstream file_CPU((CV_curve+"CPU_efficiency.dat").c_str(), ios::app);  //file flux declaration and file opening
    if(file && file_Order && file_CPU)  // If the opening is successful

    {
	//Writting column names
	file<<"Ncells"<<" "<<"Alpha1_L1_err"<<" "\
		<<"P1_L1_err"<<" "<<"U1_L1_err"<<" "\
		<<"P2_L1_err"<<" "<<"U2_L1_err"<<" "\
		<<"Alpha1_ex_L1_norm"<<" "<<"P1_ex_L1_norm"<<" "\
		<<"U1_ex_L1_norm"<<" "<<"P2_ex_L1_norm"<<" "\
        <<"U2_ex_L1_norm"<<endl;

	file_Order<<"Alpha1_Order"<<" "<<"P1_Order"<<" "\
		<<"U1_Order"<<" "<<"P2_Order"<<" "<<"U2_Order"<<endl;

	file_CPU<<"Ncells"<<" "<<"Alpha1_L1_err"<<" "\
		<<"P1_L1_err"<<" "<<"U1_L1_err"<<" "\
		<<"P2_L1_err"<<" "<<"U2_L1_err"<<" "\
		<<"CPU_time"<<endl;

	for (int i=0; i< Npoints; i++){

	    int Ncells=CellsTab(i);

	    cout<<"Convergence curve: Ncells= "<<Ncells<<endl;
	    //Mesh building
	    Mesh mesh_try(Length, Ncells, NGhostCells);

        int Nfaces = mesh_try.Get_Nfaces();

	    //Test of the Sol class building:
        Sol_Isen sol_try(mesh_try,\
                therm,\
                pRef, mRef,\
                tauRelax,\
                VariableType,\
                InitL, InitR,\
                SolType,\
                LeftBCType, RightBCType,\
                x0\
                );

	    //Solver building
        Solver_Isen solver_try(\
                FractionalStep,\
                SourceTermType,\
                TimeIntegrationType,\
                NRelax, CourantBL,\
                SchemeTypeCons, SchemeTypeNCons,\
                CourantConv,\
                SimulationTime, print_freq,\
                FileOutputFormat,\
                Nfaces\
                );

        //FIXME
        //double CourantBL = solver_try.CourantBL_;
        //double TauMin    = sol_try.etaRelax_(0);
        //TimeIntegration(sol_try, solver_try.SimulationTime_, CourantBL*TauMin, CV_curve, CourantBL);
        //exit(EXIT_FAILURE);

        //solver_try.BoundaryLayerTest(sol_try, mesh_try, CV_curve);

	    string FileName=CV_curve;        
       
        solver_try.Simulation(sol_try, mesh_try,\
                FileOutputFormat, FileName); 
                                
	    //Getting the exact error for the first Ncells as reference in the error calculation

	    X(i)=Ncells;

        for (int nVar=0; nVar < NbVariables; nVar++){

            //Storing the exact L1 norm for the coarser mesh
            if(i==0){
                normL1_Var_exact_ref(nVar) = sol_try.normL1_exact_(nVar);
            }

            //Storing the L1 error
            if(normL1_Var_exact_ref(nVar) < epsZero){

                ErrMatrix(i,nVar) = sol_try.errL1_(nVar);
            }
            else{
                ErrMatrix(i,nVar) = sol_try.errL1_(nVar)/normL1_Var_exact_ref(nVar);
            }
        }

	    //Get the right precision for all the prints
	    
	    file<<std::setprecision(COUT_PRECISION)<<std::scientific<<Ncells<<" "\
        <<ErrMatrix(i,0)<<" "\
        <<ErrMatrix(i,1)<<" "\
        <<ErrMatrix(i,2)<<" "\
        <<ErrMatrix(i,3)<<" "\
        <<ErrMatrix(i,4)<<" "\
        <<sol_try.normL1_exact_(0)<<" "\
        <<sol_try.normL1_exact_(1)<<" "\
        <<sol_try.normL1_exact_(2)<<" "\
        <<sol_try.normL1_exact_(3)<<" "\
        <<sol_try.normL1_exact_(4)<<endl;

	    file_CPU<<std::setprecision(COUT_PRECISION)<<std::scientific<<Ncells<< " "\
        <<ErrMatrix(i,0)<<" "\
        <<ErrMatrix(i,1)<<" "\
        <<ErrMatrix(i,2)<<" "\
        <<ErrMatrix(i,3)<<" "\
        <<ErrMatrix(i,4)<<" "\
		<<solver_try.CPUTime_<<endl;
	}

    //Writting the CV order:
    if(Npoints>=2){

        Xa=log(X(Npoints-2));
        Xb=log(X(Npoints-1));

        for(int nVar=0; nVar< NbVariables ; nVar++){

            Ya=log(ErrMatrix(Npoints-2,nVar));
            Yb=log(ErrMatrix(Npoints-1,nVar));
            Order= abs((Yb-Ya)/(Xb-Xa));
            file_Order<<std::setprecision(COUT_PRECISION)<<std::scientific<<Order<<" ";
        }
    }

    file.close();  // the file is closed
    file_Order.close();  // the file_Order is closed
    file_CPU.close();  // the file_CPU is closed
    }

    else{cerr << "Error Convergence_curve, Can't open file !" << endl;}

    cout<<endl;
}
