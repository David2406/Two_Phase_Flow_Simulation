#include "Sol_Isentropic.h"

//constructor:

Sol_Isen::Sol_Isen(Mesh& M,\
		ThermoLaw& Therm,\
        double pRef, double mRef,\
        Vector3d tauRelax,\
        string VariableType,\
        Vector5d& InitL, Vector5d& InitR,\
		string SolType,\
		string LeftBCType, string RightBCType,\
        double x0\
		): SolTherm_(Therm)
{

    LeftBCType_ =LeftBCType;
    RightBCType_=RightBCType;
    SolType_    =SolType;

    int Ncells    =M.Get_Ncells();
    int NcellExt  =M.Get_NcellExt();
    int Nfaces    =M.Get_Nfaces();
    double Length =M.Get_Length();

    x_0_ = x0;

    //Copying the variable type
    VariableType_ = VariableType;

    //Copying the initial Rieman state (theoretical solution)
    InitL_ = InitL;
    InitR_ = InitR;


    //Copying the relaxation times
    tauRelax_ = tauRelax;

    //Relaxation parameters
    etaRelax_(0) = min(tauRelax_(0), tauRelax_(1));
    etaRelax_(1) = etaRelax_(0)/tauRelax_(0);
    etaRelax_(2) = etaRelax_(0)/tauRelax_(1);

    //Copying the reference pressure and density
    pRef_ = pRef;
    mRef_ = mRef;

    //Initializing over the overall mesh
    NConsVarInit(InitR, InitL, tauRelax, pRef, mRef,\
            NcellExt, Ncells, Length);

    ConsVarInit(InitR, InitL, tauRelax, pRef, mRef,\
            NcellExt, Ncells, Length);

    //Initializing the L1 errors 
    errL1_<<Big,
            Big,
            Big,
            Big,
            Big;

    //Initializing the L1 norms 
    normL1_exact_<<ZERO,
                   ZERO,
                   ZERO,
                   ZERO,
                   ZERO;
    
    //Initializing the ConsVarDual_ to ZERO 
    ConsVarDual_     = MatrixXd::Zero(NcellExt,5);
    ConsVarOld_      = MatrixXd::Zero(NcellExt,5);

    //Initializing the ConsVarDual_ to ZERO 
    NConsVarDual_    = MatrixXd::Zero(NcellExt,5);

    //Initializing the NConsVarEqRelaxDual_ to ZERO 
    NConsVarEqRelaxDual_ = MatrixXd::Zero(NcellExt,5);

    //Initializing the conservative flux terms
    ConsVarFlux_ = MatrixXd::Zero(Nfaces,5);

    //Initializing the non-conservative terms
    NConsVarFlux_  = MatrixXd::Zero(NcellExt,5);
    NConsVarFluxR_ = MatrixXd::Zero(NcellExt,5);
    NConsVarFluxL_ = MatrixXd::Zero(NcellExt,5);

    //Initializing the face values used for the non-conservative terms
    NConsVarFace_ = MatrixXd::Zero(Nfaces,1);

    //Initializing the exact solution
    SolExact_        = NConsVar_;
    SolExactEqRelax_ = NConsVarEqRelax_;

    /////////////////////////////////////////////////////////////////
    //FROZEN PART AND CUBIC RELAXATION//
    Vector5d Init_avr = ONE_OVER_TWO*(InitR_ + InitL_); 

    //Reference state
    //Out of Equilibrium
    /*
       W_ref_<<0.5,
       1.e5,
       ONE,
       3.e6,
       ONE;
     */
    //At Equilibrium
    W_ref_<<0.6,
           1.5e7,
           1.,
           1.5e7,
           1.;

    //Equilibrium state
    
    W_eq_<<0.5,
        3.e5,
        3.,
        3.e5,
        3.;
        
    //W_eq_ = VectorXd::Zero(5);

    //Time matrix
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

    Vector5d U_state_ref = NConsVarToConsVarLoc(W_ref_, Therm);

    Vector5d One;
    One<<ONE,
         ONE,
         ONE,
         ONE,
         ONE;
    //TimeMat_ = TimeRelaxMat(U_state_ref, Relax_mu);
    TimeMat_ = TimeRelaxMat(One, Relax_mu);
    TimeMat_         *= (ONE/etaRelax_(0));

    //Eigenvector base matrix
    EigenVectorBasis_ = IsentropicBN_Eigenvectors(\
            VariableType,\
            Init_avr,\
            Therm\
            );
   
    //Eigenvector base inverse matrix
    EigenVectorBasisInv_ = IsentropicBN_EigenvectorsInv(\
            VariableType,\
            Init_avr,\
            Therm\
            );

    JacConvFrozen_ = LinearizedJacobian(Init_avr,\
        VariableType,\
        Therm\
        );

    //Dirac MassFrozen
    Vector5d Uinit_avr = NConsVarToConsVarLoc(Init_avr, Therm);
    
    double alpha1_L = InitL_(0);
    double alpha1_R = InitR_(0);

    double m1_avr = Uinit_avr(1);
    double u1_avr = Init_avr(2);
    double p1_avr = Init_avr(1);

    double m2_avr = Uinit_avr(3);
    double u2_avr = Init_avr(4);
    double p2_avr = Init_avr(3);

    double m_avr  = m1_avr + m2_avr;
    double uI_frozen = (m1_avr/m_avr)*u1_avr + (m2_avr/m_avr)*u2_avr; 
    double pI_frozen = (m2_avr/m_avr)*p1_avr + (m1_avr/m_avr)*p2_avr;

    //uI_Frozen_=uI_frozen;
    uI_Frozen_=100.;

    DiracMassFrozen_ << -uI_frozen*(alpha1_R - alpha1_L),
                        ZERO,
                        +pI_frozen*(alpha1_R - alpha1_L),
                        ZERO,
                        -pI_frozen*(alpha1_R - alpha1_L);

    Vector5d EigenvaluesFrozen = IsentropicBN_Eigenvalues(\
                Init_avr,\
                Therm\
                );

    EigenvaluesFrozen_ = (EigenvaluesFrozen.array().abs()).maxCoeff();

    //////////////////////////////////////////////////////////////////
}

Sol_Isen::Sol_Isen( Sol_Isen& solution){

    LeftBCType_   = solution.LeftBCType_;
    RightBCType_  = solution.RightBCType_;
    VariableType_ = solution.VariableType_;

    SolTherm_=solution.SolTherm_;

    ConsVar_            = solution.ConsVar_;
    ConsVarDual_        = solution.ConsVarDual_;
    NConsVar_           = solution.NConsVar_;
    NConsVarDual_       = solution.NConsVarDual_;
    NConsVarEqRelax_    = solution.NConsVarEqRelax_;
    NConsVarEqRelaxDual_    = solution.NConsVarEqRelaxDual_;
    SoundSpeed_         = solution.SoundSpeed_;
    Mach_               = solution.Mach_;
    ConsVarFlux_        = solution.ConsVarFlux_;
    NConsVarFlux_       = solution.NConsVarFlux_;
    NConsVarFluxL_      = solution.NConsVarFluxL_;
    NConsVarFluxR_      = solution.NConsVarFluxR_;
    SourceTerms_        = solution.SourceTerms_;
    SourceTermsEqRelax_ = solution.SourceTermsEqRelax_;
    SolExact_           = solution.SolExact_;

}

Sol_Isen::Sol_Isen(){

    LeftBCType_   ="None";
    RightBCType_  ="None";

}

//methods:

void Sol_Isen::NConsVarInit(Vector5d& InitR, Vector5d& InitL,\
            Vector3d tauRelax, double pRef, double mRef,\
            int NcellExt, int Ncells, double Length){

    //local variables

    double tauP = tauRelax(0);
    double tauU = tauRelax(1);
    
    double rho1, u1, p1, c1; 
    double rho2, u2, p2, c2;

    double alpha1, alpha2;
    double m1, m2, m;
    double C1, C2, du, dp, dC;

    double SpaceStep=Length/Ncells;
    double x_0 = (*this).x_0_;
    double Ncols = InitR.rows();

    //Resizing the NConsVar_ array
    NConsVar_.resize(NcellExt, Ncols);
    NConsVarEqRelax_.resize(NcellExt, Ncols);
    SourceTermsEqRelax_.resize(NcellExt, Ncols);
    SoundSpeed_.resize(NcellExt, 2);
    Mach_.resize(NcellExt, 2);

    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

        if(i<=Riemann_index){

            ((*this).NConsVar_).block(i,0,1,Ncols) = InitL.transpose();

            alpha1 = InitL(0);
            p1     = InitL(1);
            u1     = InitL(2);
            alpha2 = ONE - InitL(0);
            p2     = InitL(3);
            u2     = InitL(4);
        }

        else{

            ((*this).NConsVar_).block(i,0,1,Ncols) = InitR.transpose(); 

            alpha1 = InitR(0);
            p1     = InitR(1);
            u1     = InitR(2);
            alpha2 = ONE - InitR(0);
            p2     = InitR(3);
            u2     = InitR(4);

        }

        rho1   = Density_EOS(1,(*this).SolTherm_, p1, ZERO);
        m1     = alpha1*rho1;
        rho2   = Density_EOS(2,(*this).SolTherm_, p2, ZERO);
        m2     = alpha2*rho2;
        m      = m1 + m2;

        c1 = Sound_Speed_EOS(1,(*this).SolTherm_, rho1, p1);
        c2 = Sound_Speed_EOS(2,(*this).SolTherm_, rho2, p2);

        C1 = rho1*pow(c1, TWO);
        C2 = rho2*pow(c2, TWO);
        dC = C2 - C1;

        dp = p2 - p1;
        du = u2 - u1;

        ((*this).SoundSpeed_)(i,0) = c1; 
        ((*this).SoundSpeed_)(i,1) = c2;

        ((*this).Mach_)(i,0) = fabs(u1)/c1; 
        ((*this).Mach_)(i,1) = fabs(u2)/c2; 

        //NConsVarEqRelax Vector Filled
        ((*this).NConsVarEqRelax_)(i,0) = alpha1;
        ((*this).NConsVarEqRelax_)(i,1) = (m1/m)*u1 + (m2/m)*u2;
        ((*this).NConsVarEqRelax_)(i,2) = alpha1*p1 + alpha2*p2;
        ((*this).NConsVarEqRelax_)(i,3) = du;
        ((*this).NConsVarEqRelax_)(i,4) = dp;

        //SourceTermsEqRelax Vector Filled
        ((*this).SourceTermsEqRelax_)(i,0) = -(Kp_Coeff(pRef, alpha1)/tauP)*dp;
        ((*this).SourceTermsEqRelax_)(i,1) = ZERO;
        ((*this).SourceTermsEqRelax_)(i,2) = (Kp_Coeff(pRef, alpha1)/tauP)*(dp - dC)*dp;
        ((*this).SourceTermsEqRelax_)(i,3) = -(Ku_Coeff(mRef, alpha1)/tauU)*(m/(m1*m2))*du;
        ((*this).SourceTermsEqRelax_)(i,4) = -(Kp_Coeff(pRef, alpha1)/tauP)*(C1/alpha1 + C2/alpha2)*dp;
    }
}

void Sol_Isen::ConsVarInit(Vector5d& InitR, Vector5d& InitL,\
            Vector3d tauRelax, double pRef, double mRef,\
            int NcellExt, int Ncells, double Length){

    double SpaceStep=Length/Ncells;
    double x_0 = (*this).x_0_;
    double Ncols = InitR.rows();

    double tauP = tauRelax(0);
    double tauU = tauRelax(1);

    double du, dp;

    //local variables
    double alpha1, rho1, u1, p1; 
    double alpha2, rho2, u2, p2;

    double m1, m2;

    //Resizing the NConsVar_ array
    ConsVar_.resize(NcellExt, Ncols);
    SourceTerms_.resize(NcellExt, Ncols);

    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

        if(i<=Riemann_index){

            alpha1 = InitL(0);
            p1     = InitL(1);
            u1     = InitL(2);
            alpha2 = ONE-InitL(0);
            p2     = InitL(3);
            u2     = InitL(4);

        }

        else{

            alpha1 = InitR(0);
            p1     = InitR(1);
            u1     = InitR(2);
            alpha2 = ONE-InitR(0);
            p2     = InitR(3);
            u2     = InitR(4);

        }

        rho1 = Density_EOS(1,(*this).SolTherm_, p1, ZERO);
        rho2 = Density_EOS(2,(*this).SolTherm_, p2, ZERO);
        m1   = alpha1*rho1;
        m2   = alpha2*rho2;

        dp = p2 - p1;
        du = u2 - u1;

        //Conservative Vector Filled
        ((*this).ConsVar_)(i,0) = alpha1;
        ((*this).ConsVar_)(i,1) = m1;
        ((*this).ConsVar_)(i,2) = m1*u1;
        ((*this).ConsVar_)(i,3) = m2;
        ((*this).ConsVar_)(i,4) = m2*u2;

        //SourceTerms Vector Filled
        //FIXME
        ((*this).SourceTerms_)(i,0) = ZERO;
        ((*this).SourceTerms_)(i,1) = ZERO;
        ((*this).SourceTerms_)(i,2) = ZERO;
        ((*this).SourceTerms_)(i,3) = ZERO;
        ((*this).SourceTerms_)(i,4) = ZERO;
        /*
        ((*this).SourceTerms_)(i,0) = -(Kp_Coeff(pRef, alpha1)/tauP)*dp;
        ((*this).SourceTerms_)(i,1) = ZERO;
        ((*this).SourceTerms_)(i,2) = +(Ku_Coeff(mRef, alpha1)/tauU)*du;
        ((*this).SourceTerms_)(i,3) = ZERO;
        ((*this).SourceTerms_)(i,4) = -(Ku_Coeff(mRef, alpha1)/tauU)*du;
        */
    }
}

//Update methods:

void Sol_Isen::SoundSpeed_Update(){

    //local variables
    int size = NConsVar_.cols();
    VectorXd ZeroTab = VectorXd::Zero(size);

    //rho1, p1
    VectorXd P   = NConsVar_.col(1);
    VectorXd Rho = Density_EOS_Tab(1, (*this).SolTherm_, P, ZeroTab);
    SoundSpeed_.col(0) = Sound_Speed_EOS_Tab(1, (*this).SolTherm_, Rho,P);                        

    //rho2, p2
    P   = NConsVar_.col(3);
    Rho = Density_EOS_Tab(2, (*this).SolTherm_, P, ZeroTab);
    SoundSpeed_.col(1) = Sound_Speed_EOS_Tab(2, (*this).SolTherm_, Rho,P);
}

void Sol_Isen::Mach_Update(){

    //u1
    VectorXd U   = NConsVar_.col(2);
    Mach_.col(0) = U.cwiseQuotient(SoundSpeed_.col(0));                         

    //u2
    U   = NConsVar_.col(4);
    Mach_.col(1) = U.cwiseQuotient(SoundSpeed_.col(1));                        
}

void Sol_Isen::NConsVarToEqRelax(string MeshTag){

    //local variables
    int nrows = NConsVar_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

    //Update using the Primal mesh values
    if(MeshTag == "Primal"){

        VectorXd alpha1 = NConsVar_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd p1   = NConsVar_.col(1);
        VectorXd p2   = NConsVar_.col(3);
        VectorXd u1   = NConsVar_.col(2);
        VectorXd u2   = NConsVar_.col(4);
        VectorXd rho1 = Density_EOS_Tab(1, SolTherm_, p1, ZeroTab);
        VectorXd rho2 = Density_EOS_Tab(2, SolTherm_, p2, ZeroTab);
        VectorXd m1  = alpha1.cwiseProduct(rho1);
        VectorXd m2  = alpha2.cwiseProduct(rho2);
        VectorXd m   = m1 + m2;
        VectorXd Y1  = m1.cwiseQuotient(m); 
        VectorXd Y2  = m2.cwiseQuotient(m); 

        //alpha1
        NConsVarEqRelax_.col(0)  = NConsVar_.col(0);

        //U
        NConsVarEqRelax_.col(1)  = Y1.cwiseProduct(u1) +\
                                   Y2.cwiseProduct(u2);
        //P
        NConsVarEqRelax_.col(2)  = alpha1.cwiseProduct(p1)+\
                                   alpha2.cwiseProduct(p2);

        //du
        NConsVarEqRelax_.col(3) = u2 - u1;

        //dp
        NConsVarEqRelax_.col(4) = p2 - p1;

    }
    //Update using the Dual mesh values
    else{

        VectorXd alpha1 = NConsVarDual_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd p1   = NConsVarDual_.col(1);
        VectorXd p2   = NConsVarDual_.col(3);
        VectorXd u1   = NConsVarDual_.col(2);
        VectorXd u2   = NConsVarDual_.col(4);
        VectorXd rho1 = Density_EOS_Tab(1, SolTherm_, p1, ZeroTab);
        VectorXd rho2 = Density_EOS_Tab(2, SolTherm_, p2, ZeroTab);
        VectorXd m1  = alpha1.cwiseProduct(rho1);
        VectorXd m2  = alpha2.cwiseProduct(rho2);
        VectorXd m   = m1 + m2;
        VectorXd Y1  = m1.cwiseQuotient(m); 
        VectorXd Y2  = m2.cwiseQuotient(m); 

        //alpha1
        NConsVarEqRelaxDual_.col(0)  = NConsVarDual_.col(0);

        //U
        NConsVarEqRelaxDual_.col(1)  = Y1.cwiseProduct(u1) +\
                                   Y2.cwiseProduct(u2);
        //P
        NConsVarEqRelaxDual_.col(2)  = alpha1.cwiseProduct(p1)+\
                                   alpha2.cwiseProduct(p2);

        //du
        NConsVarEqRelaxDual_.col(3) = u2 - u1;

        //dp
        NConsVarEqRelaxDual_.col(4) = p2 - p1;
    }

}

void Sol_Isen::ConsVarToNConsVar(string MeshTag){

    //local variables
    int nrows = NConsVar_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

    //Update using the Primal mesh values
    if(MeshTag == "Primal"){

        VectorXd alpha1 = ConsVar_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd m1     = ConsVar_.col(1);
        VectorXd m2     = ConsVar_.col(3);
        VectorXd m1u1   = ConsVar_.col(2);
        VectorXd m2u2   = ConsVar_.col(4);

        VectorXd rho1 = m1.cwiseQuotient(alpha1);
        VectorXd rho2 = m2.cwiseQuotient(alpha2);
        VectorXd p1 = Pressure_EOS_Tab(1, SolTherm_, rho1, ZeroTab);
        VectorXd p2 = Pressure_EOS_Tab(2, SolTherm_, rho2, ZeroTab);
        VectorXd u1 = m1u1.cwiseQuotient(m1);
        VectorXd u2 = m2u2.cwiseQuotient(m2);

        //alpha1
        NConsVar_.col(0)  = alpha1;

        //p1
        NConsVar_.col(1)  = p1;

        //u1
        NConsVar_.col(2)  = u1;

        //p2
        NConsVar_.col(3)  = p2;

        //u2
        NConsVar_.col(4)  = u2;
    }
    //Update using the Dual mesh values
    else{

        VectorXd alpha1 = ConsVarDual_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd m1     = ConsVarDual_.col(1);
        VectorXd m2     = ConsVarDual_.col(3);
        VectorXd m1u1   = ConsVarDual_.col(2);
        VectorXd m2u2   = ConsVarDual_.col(4);

        VectorXd rho1 = m1.cwiseQuotient(alpha1);
        VectorXd rho2 = m2.cwiseQuotient(alpha2);
        VectorXd p1 = Pressure_EOS_Tab(1, SolTherm_, rho1, ZeroTab);
        VectorXd p2 = Pressure_EOS_Tab(2, SolTherm_, rho2, ZeroTab);
        VectorXd u1 = m1u1.cwiseQuotient(m1);
        VectorXd u2 = m2u2.cwiseQuotient(m2);

        //alpha1
        NConsVarDual_.col(0)  = alpha1;

        //p1
        NConsVarDual_.col(1)  = p1;

        //u1
        NConsVarDual_.col(2)  = u1;

        //p2
        NConsVarDual_.col(3)  = p2;

        //u2
        NConsVarDual_.col(4)  = u2;
    }

}

void Sol_Isen::EqRelaxToNConsVar(string MeshTag){

    //local variables
    int nrows = NConsVarEqRelax_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

    //Update using the Primal mesh values
    if(MeshTag == "Primal"){

        VectorXd alpha1 = NConsVarEqRelax_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd U      = NConsVarEqRelax_.col(1);
        VectorXd P      = NConsVarEqRelax_.col(2);
        VectorXd du     = NConsVarEqRelax_.col(3);
        VectorXd dp     = NConsVarEqRelax_.col(4);

        VectorXd p1     = P - alpha2.cwiseProduct(dp);
        VectorXd p2     = P + alpha1.cwiseProduct(dp);

        VectorXd rho1 = Density_EOS_Tab(1, SolTherm_, p1, ZeroTab);
        VectorXd rho2 = Density_EOS_Tab(2, SolTherm_, p2, ZeroTab);
        VectorXd m1  = alpha1.cwiseProduct(rho1);
        VectorXd m2  = alpha2.cwiseProduct(rho2);
        VectorXd m   = m1 + m2;
        VectorXd Y1  = m1.cwiseQuotient(m); 
        VectorXd Y2  = m2.cwiseQuotient(m); 

        //alpha1
        NConsVar_.col(0)  = alpha1;

        //p1
        NConsVar_.col(1)  = p1;

        //u1
        NConsVar_.col(2)  = U - Y2.cwiseProduct(du);

        //p2
        NConsVar_.col(3)  = p2;

        //dp
        NConsVar_.col(4)  = U + Y1.cwiseProduct(du);
    }
    //Update using the Dual mesh values
    else{

        VectorXd alpha1 = NConsVarEqRelaxDual_.col(0);
        VectorXd alpha2 = VectorXd::Constant(nrows, ONE) - alpha1;
        VectorXd U      = NConsVarEqRelaxDual_.col(1);
        VectorXd P      = NConsVarEqRelaxDual_.col(2);
        VectorXd du     = NConsVarEqRelaxDual_.col(3);
        VectorXd dp     = NConsVarEqRelaxDual_.col(4);

        VectorXd p1     = P - alpha2.cwiseProduct(dp);
        VectorXd p2     = P + alpha1.cwiseProduct(dp);

        VectorXd rho1 = Density_EOS_Tab(1, SolTherm_, p1, ZeroTab);
        VectorXd rho2 = Density_EOS_Tab(2, SolTherm_, p2, ZeroTab);
        VectorXd m1  = alpha1.cwiseProduct(rho1);
        VectorXd m2  = alpha2.cwiseProduct(rho2);
        VectorXd m   = m1 + m2;
        VectorXd Y1  = m1.cwiseQuotient(m); 
        VectorXd Y2  = m2.cwiseQuotient(m); 

        //alpha1
        NConsVarDual_.col(0)  = alpha1;

        //p1
        NConsVarDual_.col(1)  = p1;

        //u1
        NConsVarDual_.col(2)  = U - Y2.cwiseProduct(du);

        //p2
        NConsVarDual_.col(3)  = p2;

        //dp
        NConsVarDual_.col(4)  = U + Y1.cwiseProduct(du);
    }

}

void Sol_Isen::SolExact_Update(Mesh& mesh, double time){

    //Local variables
    double xcell;

    //Function
    int NcellExt = mesh.Get_NcellExt();

    Vector5d SolExactEqRelax_L = NConsVarToEqRelaxLoc(InitL_, SolTherm_);
    Vector5d SolExactEqRelax_R = NConsVarToEqRelaxLoc(InitR_, SolTherm_);

    double u2_L = InitL_(4);

    if(SolType_== "Pure Contact"){

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);

            if(xcell<=x_0_+time*u2_L){

                SolExact_.row(cell) = InitL_.transpose();
                SolExactEqRelax_.row(cell) = SolExactEqRelax_L.transpose();

            }
            else{

                SolExact_.row(cell) = InitR_.transpose();
                SolExactEqRelax_.row(cell) = SolExactEqRelax_R.transpose();
            }


        }

    }
    else if(SolType_== "Shock Contact" || SolType_=="BN Relaxation"){

       double sigma1Shock =  -839.323232317123598;
       double sigma2Shock =  -528.970548257695896;

       double p1_star     = 5.e5;
       double u1_star     = 5.e0;
       double p2_star     = 3.e5;
       double u2_star     = 4.e0;

       Vector5d W_state1_star, W_state2_star;

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);

            if(xcell<=x_0_+time*sigma1Shock){

                SolExact_.row(cell) = InitL_.transpose();
                SolExactEqRelax_.row(cell) = SolExactEqRelax_L.transpose();

            }
            //after the shock wave of phase 1 
            else if(xcell<=x_0_+time*sigma2Shock){

                W_state1_star << InitL_(0),
                              p1_star,
                              u1_star,
                              InitL_(3),
                              InitL_(4);

                SolExact_.row(cell) = W_state1_star.transpose();
                SolExactEqRelax_.row(cell) = NConsVarToEqRelaxLoc(W_state1_star,\
                                              SolTherm_).transpose();

            }
            //after the shock wave of phase 2 
            else if(xcell<=x_0_+time*u2_star){

                W_state2_star << InitL_(0),
                              p1_star,
                              u1_star,
                              p2_star,
                              u2_star;

                SolExact_.row(cell)        = W_state2_star.transpose();
                SolExactEqRelax_.row(cell) = NConsVarToEqRelaxLoc(W_state2_star,\
                                              SolTherm_).transpose();
            }
            //After the contact discontinuity wave
            else{

                SolExact_.row(cell) = InitR_.transpose();
                SolExactEqRelax_.row(cell) = SolExactEqRelax_R.transpose();
            }


        }

    }
    else if(SolType_== "Linear BN"){

        double weight_L      = ONE_OVER_TWO;
        Vector5d W_state_avr = NConsVarAveraged(InitL_, InitR_, weight_L);
        Vector5d U_update, W_update, V_update;

        //Used to compute the exact solution
        double aCoordsIni;
        Vector5d PureConvVector;

        //Eigenvector base matrix
        Matrix5d RmatEigenVectors = EigenVectorBasis_;

        //EigenvectorInv base matrix
        Matrix5d RmatEigenVectorsInv = EigenVectorBasisInv_;

        //Coords of InitR_ - InitL_ in the eigenvector base

        Vector5d U_state_L = NConsVarToConsVarLoc(InitL_, SolTherm_);
        Vector5d U_state_R = NConsVarToConsVarLoc(InitR_, SolTherm_);
        Vector5d DeltaCoords = RmatEigenVectorsInv*(U_state_R - U_state_L);
        Vector5d InitLCoords = RmatEigenVectorsInv*(U_state_L);

        //Heaviside function according to the eigenvalues
        Vector5d H_Eigen;
        double H_uI_Frozen;

        //Ordering of the eigenvalues
        Vector5d EigenValues = IsentropicBN_Eigenvalues(\
                W_state_avr,\
                SolTherm_\
                );

        //Dirac mass vector
        Vector5d Dirac_Jump = RmatEigenVectorsInv*DiracMassFrozen_;
        for(int nVar=0;nVar<5;nVar++){

            Dirac_Jump(nVar) /= (EigenValues(nVar) - uI_Frozen_);

        }

        //Result
        Vector5d Result;
        Result<<ZERO,
                ZERO,
                ZERO,
                ZERO,
                ZERO;

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);
            H_Eigen<<Heaviside(xcell-x_0_, time, EigenValues(0)),
                     Heaviside(xcell-x_0_, time, EigenValues(1)),
                     Heaviside(xcell-x_0_, time, EigenValues(2)),
                     Heaviside(xcell-x_0_, time, EigenValues(3)),
                     Heaviside(xcell-x_0_, time, EigenValues(4));

            H_uI_Frozen = Heaviside(xcell-x_0_, time, uI_Frozen_); 

            PureConvVector = InitLCoords + DeltaCoords.cwiseProduct(H_Eigen);

            for(int nVar=0;nVar<5;nVar++){

                aCoordsIni   = PureConvVector(nVar);

                Result(nVar) = aCoordsIni -\
                        Dirac_Jump(nVar)*Heaviside(xcell-x_0_, time, EigenValues(nVar));
            }

            U_update = RmatEigenVectors*Result + RmatEigenVectors*Dirac_Jump*\
                                               H_uI_Frozen;

            //Update of the solution
            W_update = ConsVarToNConsVarLoc(U_update, SolTherm_);

            V_update = NConsVarToEqRelaxLoc(W_update, SolTherm_);

            SolExact_.row(cell) = W_update.transpose();
            SolExactEqRelax_.row(cell) = V_update.transpose(); 

        }

    }

    else if(SolType_== "Linear Convection Relaxation BN"){

        double weight_L      = ONE_OVER_TWO;
        Vector5d W_state_avr = NConsVarAveraged(InitL_, InitR_, weight_L);
        Vector5d U_update, W_update, V_update;

        //Used to compute the exact solution
        double aInf;
        double aCoord_L, aCoord_R, aRelaxL, aRelaxR;
        double aRelaxDelta, aDirac_Jump;
        double lambda, H_eigen_loc, Delta_H;

        //Eigenvector base matrix
        Matrix5d RmatEigenVectors = EigenVectorBasis_;

        //EigenvectorInv base matrix
        Matrix5d RmatEigenVectorsInv = EigenVectorBasisInv_;
         
        //Exponential time operator
        double Tk;
        Vector5d One;
        One<<ONE,
             ONE,
             ONE,
             ONE,
             ONE;
        Vector5d tauvec = TimeMat_*One;

        //Coords of InitR_ - InitL_ in the eigenvector base

        Vector5d U_state_L = NConsVarToConsVarLoc(InitL_, SolTherm_);
        Vector5d U_state_R = NConsVarToConsVarLoc(InitR_, SolTherm_);

        Vector5d InitLCoords = RmatEigenVectorsInv*(U_state_L);
        Vector5d InitRCoords = RmatEigenVectorsInv*(U_state_R);
        
        Vector5d U_state_inf = NConsVarToConsVarLoc(W_eq_, SolTherm_);
        Vector5d U_infCoords = RmatEigenVectorsInv*U_state_inf;
        Vector5d InitLCoordsRelax,InitRCoordsRelax;

        //Heaviside function according to the eigenvalues
        Vector5d H_Eigen;
        double H_uI_Frozen;

        //Ordering of the eigenvalues
        Vector5d EigenValues = IsentropicBN_Eigenvalues(\
                W_state_avr,\
                SolTherm_\
                );

        //Dirac mass vector
        double dirac_weight;
        Vector5d Dirac_Jump = RmatEigenVectorsInv*DiracMassFrozen_;
        for(int nVar=0;nVar<5;nVar++){

            Dirac_Jump(nVar) /= (EigenValues(nVar) - uI_Frozen_);

        }

        //Result
        Vector5d Result;
        Result<<ZERO,
                ZERO,
                ZERO,
                ZERO,
                ZERO;

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);
            H_Eigen<<Heaviside(xcell-x_0_, time, EigenValues(0)),
                     Heaviside(xcell-x_0_, time, EigenValues(1)),
                     Heaviside(xcell-x_0_, time, EigenValues(2)),
                     Heaviside(xcell-x_0_, time, EigenValues(3)),
                     Heaviside(xcell-x_0_, time, EigenValues(4));

            H_uI_Frozen = Heaviside(xcell-x_0_, time, uI_Frozen_); 

            /*
            cout<<"xcell = "<<xcell<<endl;
            cout<<"H_Eigen_uI_frozen = "<<H_uI_Frozen<<endl;
            */

            for(int nVar=0;nVar<5;nVar++){

                aCoord_L    = InitLCoords(nVar);
                aCoord_R    = InitRCoords(nVar);
                aInf        = U_infCoords(nVar);
                dirac_weight = Dirac_Jump(nVar);
                Tk = tauvec(nVar);
                lambda      = EigenValues(nVar);
                H_eigen_loc = Heaviside(xcell-x_0_, time, lambda);
                Delta_H = H_eigen_loc - H_uI_Frozen;

                InitLCoordsRelax(nVar) = (ONE-exp(-Tk*time))*aInf +\
                                         exp(-Tk*time)*aCoord_L;  
                aRelaxL = InitLCoordsRelax(nVar);

                InitRCoordsRelax(nVar) = (ONE-exp(-Tk*time))*aInf +\
                                         exp(-Tk*time)*aCoord_R;
                aRelaxR = InitRCoordsRelax(nVar);

                aRelaxDelta = aRelaxL + (aRelaxR - aRelaxL)*H_eigen_loc; 

                aDirac_Jump  = AkRelaxDiracTest(xcell - x_0_, time,\
                        uI_Frozen_, lambda,\
                        dirac_weight, Tk); 

                Result(nVar) = -aDirac_Jump*Delta_H;
                Result(nVar) += aRelaxDelta;

            }

            U_update = RmatEigenVectors*Result;

            //Update of the solution
            W_update = ConsVarToNConsVarLoc(U_update, SolTherm_);

            V_update = NConsVarToEqRelaxLoc(W_update, SolTherm_);

            SolExact_.row(cell) = W_update.transpose();
            SolExactEqRelax_.row(cell) = V_update.transpose(); 

        }

    }
    else if(SolType_== "Linear Convection Relaxation Cubic BN"){

        double weight_L      = ONE_OVER_TWO;
        Vector5d W_state_avr = NConsVarAveraged(InitL_, InitR_, weight_L);
        Vector5d U_update, W_update, V_update;

        //Used to compute the exact solution
        double aCoordsIni, aInf;
        Vector5d PureConvVector;

        //Eigenvector base matrix
        Matrix5d RmatEigenVectors = EigenVectorBasis_;

        //EigenvectorInv base matrix
        Matrix5d RmatEigenVectorsInv = EigenVectorBasisInv_;

        //Reference state
        Vector5d U_state_ref = NConsVarToConsVarLoc(W_ref_, SolTherm_);
         
        Vector5d aRef = U_state_ref;

        //Exponential time operator
        double tau      = etaRelax_(0);

        //Coords of InitR_ - InitL_ in the eigenvector base

        Vector5d U_state_L = NConsVarToConsVarLoc(InitL_, SolTherm_);
        Vector5d U_state_R = NConsVarToConsVarLoc(InitR_, SolTherm_);
        Vector5d DeltaCoords = RmatEigenVectorsInv*(U_state_R - U_state_L);
        Vector5d InitLCoords = RmatEigenVectorsInv*(U_state_L);
        
        Vector5d U_state_inf = NConsVarToConsVarLoc(W_eq_, SolTherm_);

        Vector5d U_infCoords = RmatEigenVectorsInv*U_state_inf;

        //Heaviside function according to the eigenvalues
        Vector5d H_Eigen;

        //Ordering of the eigenvalues
        Vector5d EigenValues = IsentropicBN_Eigenvalues(\
                W_state_avr,\
                SolTherm_\
                );

        //Result
        Vector5d Result;
        Result<<ZERO,
                ZERO,
                ZERO,
                ZERO,
                ZERO;

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);
            H_Eigen<<Heaviside(xcell-x_0_, time, EigenValues(0)),
                     Heaviside(xcell-x_0_, time, EigenValues(1)),
                     Heaviside(xcell-x_0_, time, EigenValues(2)),
                     Heaviside(xcell-x_0_, time, EigenValues(3)),
                     Heaviside(xcell-x_0_, time, EigenValues(4));

            PureConvVector = InitLCoords + DeltaCoords.cwiseProduct(H_Eigen);

            for(int nVar=0;nVar<5;nVar++){

                aCoordsIni   = PureConvVector(nVar);
                aInf         = U_infCoords(nVar);

                if(aCoordsIni < aInf){

                Result(nVar) = aInf - ONE/sqrt( ONE/pow(aCoordsIni - aInf, TWO) + TWO*(time/tau)*(ONE/pow(aRef(nVar), TWO)));

                }
                else{

                Result(nVar) = aInf + ONE/sqrt( ONE/pow(aCoordsIni - aInf, TWO) + TWO*(time/tau)*(ONE/pow(aRef(nVar), TWO)));

                }

            }

            U_update = RmatEigenVectors*Result;

            //Update of the solution
            W_update = ConsVarToNConsVarLoc(U_update, SolTherm_);

            V_update = NConsVarToEqRelaxLoc(W_update, SolTherm_);

            SolExact_.row(cell) = W_update.transpose();
            SolExactEqRelax_.row(cell) = V_update.transpose(); 

        }

    }

}

void Sol_Isen::Compute_L1_err(Mesh& mesh){

    //Local variable
    int NbVariables = 5;

    //Function
    double SpaceStep = mesh.Get_SpaceStep();
    int Ncells       = mesh.Get_Ncells();

    VectorXd Var_disc(Ncells), Var_exact(Ncells), Delta_Var(Ncells);
    VectorXd Var_exact_abs(Ncells);

    for (int nVar = 0; nVar< NbVariables; nVar++){

        //Extracting the variable without the fictitious cells
        Var_disc = (NConsVar_).block(1,nVar,Ncells,1);

        //Extracting the exact variable without the fictitious cells
        Var_exact     = (SolExact_).block(1,nVar,Ncells,1);
        Var_exact_abs = Var_exact.array().abs();

        Delta_Var = (Var_disc-Var_exact).array().abs();

        //norm L1 of Var exact
        normL1_exact_(nVar) = SpaceStep*Var_exact_abs.sum();

        //error in norm L1 of the variable
        errL1_(nVar) = (SpaceStep*(Delta_Var.sum()));
    }

}

//Checking methods

void Sol_Isen::Check_Conservativity(int nVar, Mesh& mesh){

    //Local variable
    int NcellExt         = mesh.Get_NcellExt();
    double SpaceStep     = mesh.Get_SpaceStep();
    double IntegratedVar = ZERO;

    //Function

    for(int icell =0; icell< NcellExt; icell++){

        IntegratedVar += ConsVar_(icell,nVar);

    }

    IntegratedVar *= SpaceStep;

    cout<<IntegratedVar<<endl;

}

//External functions:

Vector5d NConsVarAveraged(\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        double weight_L){

    return weight_L*W_state_L + (ONE-weight_L)*W_state_R;

}


double Kp_Coeff(double p0, double alpha1){

    return alpha1*(ONE - alpha1)/p0;

}

double Ku_Coeff(double m0, double alpha1){

    return alpha1*(ONE - alpha1)*m0;

}
double Da1_Kp_Coeff(double p0, double alpha1){

    return (ONE - TWO*alpha1)/p0;

}

double Da1_Ku_Coeff(double m0, double alpha1){

    return (ONE - TWO*alpha1)*m0;

}

/************************************************/
/*************  CHANGE OF VARIABLE  *************/
/************************************************/

Vector5d NConsVarToConsVarLoc(Vector5d& W_state,\
        ThermoLaw& Therm\
        ){

    Vector5d U_state;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2, u1, u2;
    double m1, m2;

    //Function

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);
    u1     = W_state(2);
    u2     = W_state(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;

    U_state(0) = alpha1;
    U_state(1) = m1;
    U_state(2) = m1*u1;
    U_state(3) = m2;
    U_state(4) = m2*u2;

    return U_state;

}

Vector5d ConsVarToNConsVarLoc(Vector5d& U_state,\
        ThermoLaw& Therm\
        ){

    Vector5d W_state;

    //Local variables
    double alpha1, alpha2, m1, m2, p1, p2, m1u1, m2u2;

    //Function

    alpha1   = U_state(0);
    alpha2   = ONE - alpha1;
    m1       = U_state(1);
    m2       = U_state(3);
    m1u1     = U_state(2);
    m2u2     = U_state(4);

    p1 = Pressure_EOS(1, Therm, m1/alpha1, ZERO);
    p2 = Pressure_EOS(2, Therm, m2/alpha2, ZERO);

    W_state(0) = alpha1;
    W_state(1) = p1;
    W_state(2) = m1u1/m1;
    W_state(3) = p2;
    W_state(4) = m2u2/m2;

    return W_state;

}

Vector5d NConsVarToEqRelaxLoc(Vector5d& W_state,\
        ThermoLaw& Therm\
        ){

    Vector5d V_state;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2, u1, u2;
    double m, m1, m2, Y1, Y2;

    //Function

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);
    u1     = W_state(2);
    u2     = W_state(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    Y1 = m1/m;
    Y2 = m2/m;

    V_state(0) = alpha1;
    V_state(1) = Y1*u1 + Y2*u2;
    V_state(2) = alpha1*p1 + alpha2*p2;
    V_state(3) = u2 - u1;
    V_state(4) = p2 - p1;

    return V_state;

}

Vector5d EqRelaxToNConsVarLoc(Vector5d& V_state,\
        ThermoLaw& Therm\
        ){

    Vector5d W_state;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2, u1, u2;
    double U, P, du, dp;
    double m, m1, m2, Y1, Y2;

    //Function

    alpha1 = V_state(0);
    U  = V_state(1);
    P  = V_state(2);
    du = V_state(3); 
    dp = V_state(4); 

    alpha2 = ONE - alpha1;
    p1     = P - alpha2*dp;
    p2     = P + alpha1*dp;

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    Y1 = m1/m;
    Y2 = m2/m;

    u1     = U - Y2*du;
    u2     = U + Y1*du;

    W_state(0) = alpha1;
    W_state(1) = p1;
    W_state(2) = u1;
    W_state(3) = p2;
    W_state(4) = u2;

    return W_state;

}

/************************************************/
/**********  BOUNDARY LAYER RESOLUTION  *********/
/************************************************/

/*%%%%%%%%%%  All Variables Functions %%%%%%%%%%*/

Vector5d SourceTerm(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        ){

    //Local variables
    Vector5d SourceTerm;
    double alpha1, rho1, p1, alpha2, rho2, p2;
    double m1, m2, m, c1, c2, C1, C2;
    double u1, u2;

    double Coeff_du, Coeff_dp; 

    //Function
    alpha1 = W_state(0);
    alpha2 = ONE - W_state(0);

    p1   = W_state(1);
    rho1 = Density_EOS(1, Therm, p1, ZERO);
    c1   = Sound_Speed_EOS(1, Therm, rho1, p1);
    p2   = W_state(3);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c2   = Sound_Speed_EOS(2, Therm, rho2, p2);

    u1   = W_state(2);
    u2   = W_state(4);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    C1 = rho1*c1*c1; 
    C2 = rho2*c2*c2; 

    double du     = u2 - u1;
    double dp     = p2 - p1;

    Coeff_du = etaU*Ku_Coeff(mRef, alpha1); 
    Coeff_dp = etaP*Kp_Coeff(pRef, alpha1);

    if(VariableType=="EqRelaxVar"){

        SourceTerm<< - Coeff_dp*dp,
            ZERO,
            Coeff_dp*(dp - (C2-C1))*dp,
            - Coeff_du*(m/(m1*m2))*du,
            - Coeff_dp*(C1/alpha1 + C2/alpha2)*dp;
    }
    else if(VariableType=="ConsVar"){

        SourceTerm<< - Coeff_dp*dp,
            ZERO,
            Coeff_du*du,
            ZERO,
            -Coeff_du*du;
    }

    return SourceTerm;

}

Matrix5d LinearizedSourceTerm(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        ){

    //Local variables
    Matrix5d SourceTermGradient;

    //Function
    if(VariableType=="EqRelaxVar"){

        Matrix5d LinSourceTermEq = LinearizedSourceTermsEq(W_state,\
                Therm,\
                pRef, mRef,\
                etaP, etaU);

        Matrix5d LinSourceTermOutOfEq = LinearizedSourceTermsOutOfEq(W_state,\
                Therm,\
                pRef, mRef,\
                etaP, etaU);

        SourceTermGradient = LinSourceTermEq + LinSourceTermOutOfEq;

    }
    else if(VariableType=="ConsVar"){

        SourceTermGradient = LinearizedSourceTermsCons(W_state,\
                Therm,\
                pRef, mRef,\
                etaP, etaU\
                );

    }

    return SourceTermGradient;

}

Matrix5d LinearizedJacobian(Vector5d& W_state,\
        string VariableType,\
        ThermoLaw& Therm\
        ){

    //Local variables
    Matrix5d JacobianMatrix;

    //Function
    if(VariableType=="EqRelaxVar"){

        JacobianMatrix = LinearizedJacobianEqRelax(W_state,\
                Therm);

    }
    else if(VariableType=="ConsVar"){

        JacobianMatrix = LinearizedJacobianCons(W_state,\
                Therm\
                );
    }

    return JacobianMatrix;

}


/*%%%%%%%%%%  EqRelax Variables %%%%%%%%%%*/

Matrix5d LinearizedSourceTermsEq(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        ){

    //
    Matrix5d LinSouTerEq;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2;
    double m, m1, m2, c1, c2, C1, C2;

    double Coeff_du, Coeff_dp;

    //Function

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    C1 = rho1*pow(c1,TWO);
    C2 = rho2*pow(c2,TWO);

    Coeff_du = etaU*Ku_Coeff(mRef, alpha1); 
    Coeff_dp = etaP*Kp_Coeff(pRef, alpha1); 

    LinSouTerEq<<ZERO, ZERO, ZERO, ZERO, -Coeff_dp,
                 ZERO, ZERO, ZERO, ZERO, ZERO,
                 ZERO, ZERO, ZERO, ZERO, -Coeff_dp*(C2 - C1),
                 ZERO, ZERO, ZERO, -Coeff_du*(m/(m1*m2)), ZERO,
                 ZERO, ZERO, ZERO, ZERO, -Coeff_dp*(C1/alpha1 + C2/alpha2);

    return LinSouTerEq;

}

Matrix5d LinearizedSourceTermsOutOfEq(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        ){

    //
    Matrix5d LinSouTerOutEq;
    Matrix5d LinSouTerOutEqAlpha1Var;
    Matrix5d LinSouTerOutEqThermoVar;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2;
    double m, m1, m2, c1, c2, C1, C2;
    double C1p, C2p;

    double Coeff_du, Coeff_dp;
    double Da1_Coeff_du, Da1_Coeff_dp;
    double du, dp;

    //Function

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);
    dp     = p2 - p1;
    du     = W_state(4) - W_state(2);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    C1 = rho1*pow(c1,TWO);
    C2 = rho2*pow(c2,TWO);
    double dC = C2 - C1;

    C1p = Dp_Rho_Sound_Speed_Squared_EOS(1, Therm, rho1, p1);
    C2p = Dp_Rho_Sound_Speed_Squared_EOS(2, Therm, rho2, p2);

    //Computing the thermodynamical variations
    double Peq     = (alpha2*C1 + alpha1*C2)/(alpha1*alpha2);
    double dC_da1  = dp*(C2p - C1p); 
    double dC_dP   = (C2p - C1p); 
    double dC_ddp  = (alpha1*C2p +  alpha2*C1p);

    double dm_inv_da1  = ONE/(alpha2*m2) - ONE/(alpha1*m1) -dp*(ONE/(C1*m1) + ONE/(C2*m2)); 
    double dm_inv_dP   =-(ONE/(C1*m1) + ONE/(C2*m2));
    double dm_inv_ddp  =(alpha2/(C1*m1) - alpha1/(C2*m2));

    double dCa_da1 = C2/(alpha2*alpha2) - C1/(alpha1*alpha1) +dp*(C1p/alpha1 + C2p/alpha2);
    double dCa_dP  = (C1p/alpha1 + C2p/alpha2);
    double dCa_ddp = (-(alpha2/alpha1)*C1p + (alpha1/alpha2)*C2p);

    Coeff_du     = etaU*Ku_Coeff(mRef, alpha1); 
    Coeff_dp     = etaP*Kp_Coeff(pRef, alpha1); 
    Da1_Coeff_du = etaU*Da1_Ku_Coeff(mRef, alpha1); 
    Da1_Coeff_dp = etaP*Da1_Kp_Coeff(pRef, alpha1); 

    LinSouTerOutEqThermoVar<<ZERO, ZERO, ZERO, ZERO, ZERO,
        ZERO, ZERO, ZERO, ZERO, ZERO,
        -Coeff_dp*dp*dC_da1, ZERO, -Coeff_dp*dp*dC_dP, ZERO, -Coeff_dp*dp*dC_ddp,
        -Coeff_du*du*dm_inv_da1, ZERO, -Coeff_du*du*dm_inv_dP, ZERO, -Coeff_du*du*dm_inv_ddp,
        -Coeff_dp*dp*dCa_da1, ZERO, -Coeff_dp*dp*dCa_dP, ZERO, -Coeff_dp*dp*dCa_ddp;

    LinSouTerOutEqAlpha1Var<<-Da1_Coeff_dp*dp, ZERO, ZERO, ZERO, ZERO,
        ZERO, ZERO, ZERO, ZERO, ZERO,
        Da1_Coeff_dp*dp*(dp - dC), ZERO, ZERO, ZERO, Da1_Coeff_dp*dp*TWO,
        -Da1_Coeff_du*du*(m/(m1*m2)), ZERO, ZERO, ZERO, ZERO,
        -Da1_Coeff_dp*dp*Peq, ZERO, ZERO, ZERO, ZERO;

    LinSouTerOutEq = LinSouTerOutEqThermoVar + LinSouTerOutEqAlpha1Var;

    return LinSouTerOutEq;
}

Matrix5d LinearizedJacobianEqRelax(Vector5d& W_state,\
        ThermoLaw& Therm\
        ){

    Matrix5d LinJacEqRel;

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2;
    double m, m1, m2, Y1, Y2, c1, c2, C1, C2;

    double U, du, dp;

    //local variables matrix-subfunctions

    double W_alpha, W_P, W_dp;

    double w_Ualpha1, w_Palpha1, w_dualpha1, w_dpalpha1;
    double w_UP, w_PP, w_dpP;
    double w_Pdu, w_dudu;
    double w_Udp, w_Pdp, w_dpdp;

    //Function

    //Deriving the equilibrium/relaxation variables
    Vector5d EqRelax = NConsVarToEqRelaxLoc(W_state, Therm);

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    m  = m1 + m2;
    Y1 = m1/m;
    Y2 = m2/m;
    C1 = rho1*pow(c1,TWO);
    C2 = rho2*pow(c2,TWO);

    U  = EqRelax(1);
    du = EqRelax(3);
    dp = EqRelax(4);

    //Filling each column of the jacobian
    Vector5d Col_alpha1, Col_U, Col_P, Col_du, Col_dp;

    double mclinsq = alpha1*C1 + alpha2*C2;
    W_P     = ( (m1*m2) / (m*m) )*( (C2 - C1)/(C1*C2)  );
    W_alpha = W_P*dp + ( ONE/(alpha1*alpha2) )*( (m1*m2) / (m*m) );
    W_dp    = - ( (alpha1*alpha2)/m )*( (rho1*Y2/C1) + (rho2*Y1/C2)  );

    //col alpha1
    w_Ualpha1  = Y1*Y2*((Y2/C1 + Y1/C2)*dp + (Y2/alpha1 - Y1/alpha2)  );
    w_Palpha1  = - ( alpha1*dp + C1 - mclinsq*W_alpha );
    w_dualpha1 = - ( ONE/m2 + ONE/rho1 - ONE/rho2  );
    w_dpalpha1 = dp + (C1/alpha1) + (C2 - C1)*W_alpha;

    Col_alpha1<<U + Y1*du,
                w_Ualpha1*du*du,
                w_Palpha1*du,
                W_alpha*du*du + w_dualpha1*dp,
                w_dpalpha1*du;

    //col U
    Col_U<<ZERO,
           U,
           alpha1*C1 + alpha2*C2,
           du,
           C2 - C1;

    //col P
    w_UP  = m*Y1*Y2*( Y2/C1 + Y1/C2  );
    w_PP  = -alpha1*Y2 + alpha2*Y1 + mclinsq*W_P;
    w_dpP = ONE + (C2 - C1)*W_P;

    Col_P<<ZERO,
           (ONE/m)*( ONE + w_UP*du*du  ),
           U + w_PP*du,
           W_P*du*du + ONE/rho2 - ONE/rho1,
           w_dpP*du;

    //col du
    w_Pdu  = - ( alpha1*Y2*C1 - alpha2*Y1*C2  );
    w_dudu = TWO*Y1 - ONE;

    Col_du<<ZERO,
           TWO*Y1*Y2*du,
           w_Pdu,
           U + w_dudu*du,
           Y2*C1 + Y1*C2;

    //col dp
    w_Udp  = - Y1*Y2*( alpha2*Y2/C1 - alpha1*Y1/C2  );
    w_Pdp  = alpha1*alpha2 + mclinsq*W_dp;
    w_dpdp = -alpha2*Y2 + alpha1*Y1 + (C2 - C1)*W_dp;

    Col_dp<<ZERO,
           w_Udp*du*du,
           w_Pdp*du,
           alpha2/rho1 + alpha1/rho2 + W_dp*du*du,
           U + w_dpdp*du;

    LinJacEqRel.col(0) = Col_alpha1;
    LinJacEqRel.col(1) = Col_U;
    LinJacEqRel.col(2) = Col_P;
    LinJacEqRel.col(3) = Col_du;
    LinJacEqRel.col(4) = Col_dp;

    return LinJacEqRel;

}

/*%%%%%%%%%%  Conservative Variables %%%%%%%%%%*/

Matrix5d LinearizedSourceTermsCons(Vector5d& W_state,\
        ThermoLaw& Therm,\
        double pRef, double mRef,\
        double etaP, double etaU\
        ){

    //

    //Local variables
    double alpha1, alpha2, rho1, rho2, p1, p2;
    double m1, m2, c1, c2, C1, C2;
    double u1, u2;

    double Coeff_du, Coeff_dp;
    double lambda_dp;
    double alpha12;
    double du, dp;

    //Function

    int nVar = W_state.rows();
    Matrix5d LinSouTerCons = MatrixXd::Zero(nVar,nVar);

    alpha1  = W_state(0);
    alpha2  = ONE - alpha1;
    alpha12 = ONE - TWO*alpha1;
    p1      = W_state(1);
    p2      = W_state(3);
    dp = p2 - p1;
    u1      = W_state(2);
    u2      = W_state(4);
    du = u2 - u1;

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;
    C1 = rho1*pow(c1,TWO);
    C2 = rho2*pow(c2,TWO);

    Coeff_du = etaU*Ku_Coeff(mRef, alpha1); 
    Coeff_dp = etaP*Kp_Coeff(pRef, alpha1); 

    lambda_dp = Coeff_dp*( C1/alpha1 + C2/alpha2  );

    //Filling the matrix row by row

    //Row alpha1
    Vector5d R_alpha1;
    R_alpha1<<-(etaP/pRef)*alpha12*dp - lambda_dp,
              Coeff_dp*(C1/m1),
              ZERO,
             -Coeff_dp*(C2/m2),
              ZERO;

    //Row m1 u1
    Vector5d R_m1u1;
    R_m1u1<< (etaU*mRef)*alpha12*du,
              Coeff_du*(u1/m1),
             -Coeff_du/m1,
             -Coeff_du*(u2/m2),
              Coeff_du/m2;

    LinSouTerCons.row(0) = R_alpha1.transpose();
    LinSouTerCons.row(2) = R_m1u1.transpose();
    LinSouTerCons.row(4) = -R_m1u1.transpose();
                 
    return LinSouTerCons;

}

Matrix5d LinearizedJacobianCons(Vector5d& W_state,\
        ThermoLaw& Therm\
        ){

    Matrix5d LinJacCons;

    //Local variables
    double rho1, rho2, p1, p2, dp;
    double c1, c2, C1, C2;
    double u1, u2;

    //Function

    p1     = W_state(1);
    p2     = W_state(3);
    dp     = p2 - p1;
    u1     = W_state(2);
    u2     = W_state(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c1 = Sound_Speed_EOS(1, Therm, rho1, p1);
    c2 = Sound_Speed_EOS(2, Therm, rho2, p2);

    C1 = rho1*pow(c1,TWO);
    C2 = rho2*pow(c2,TWO);

    //Filling each column of the jacobian
    Vector5d Col_alpha1, Col_m1, Col_m1u1, Col_m2, Col_m2u2;

    //col alpha1
    Col_alpha1<<u2,
                ZERO,
                -C1,
                ZERO,
                -(dp - C2);

    //col m1
    Col_m1<<ZERO,
            ZERO,
            c1*c1 - u1*u1,
            ZERO,
            ZERO;

    //col m1u1
    Col_m1u1<<ZERO,
              ONE,
              TWO*u1,
              ZERO,
              ZERO;

    //col m2
    Col_m2<<ZERO,
            ZERO,
            ZERO,
            ZERO,
            c2*c2 - u2*u2;

    //col m2u2
    Col_m2u2<<ZERO,
              ZERO,
              ZERO,
              ONE,
              TWO*u2;

    LinJacCons.col(0) = Col_alpha1;
    LinJacCons.col(1) = Col_m1;
    LinJacCons.col(2) = Col_m1u1;
    LinJacCons.col(3) = Col_m2;
    LinJacCons.col(4) = Col_m2u2;

    return LinJacCons;

}

/************************************************/
/**************  CONSERVATIVE FLUX  *************/
/************************************************/

Vector5d ConsVarFlux(Vector5d& W_state, ThermoLaw& Therm){

    //Local variables
    Vector5d Flux;

    double alpha1, alpha2, rho1, rho2, p1, p2, u1, u2;
    double m1, m2;

    //Function

    alpha1 = W_state(0);
    alpha2 = ONE - alpha1;
    p1     = W_state(1);
    p2     = W_state(3);
    u1     = W_state(2);
    u2     = W_state(4);

    rho1 = Density_EOS(1, Therm, p1, ZERO);
    rho2 = Density_EOS(2, Therm, p2, ZERO);

    m1 = alpha1*rho1;
    m2 = alpha2*rho2;

    Flux<<ZERO,
          m1*u1,
          m1*u1*u1 + alpha1*p1,
          m2*u2,
          m2*u2*u2 + alpha2*p2;

    return Flux;
}

/************************************************/
/********  ANALYTICAL SOLUTION FUNCIONS  ********/
/************************************************/

double EntropyJumpFunction(double mL, double alpha1_in, double rho1, ThermoLaw& Therm){

    //Local variables
    double p1, h1;

    //Function
    p1     = Pressure_EOS(1,Therm, rho1, ZERO);
    h1     = Internal_Energy_EOS(1, Therm, rho1, p1) + p1/rho1;

    return (mL*mL/TWO)*(ONE/pow(alpha1_in*rho1, TWO)) + h1;
}

double Rho1UpperBoundEntropyJumpFunction(double mL, double alpha1_in, ThermoLaw& Therm){

    double Gamma1 = Therm.Gamma1_;
    double kBAR1  = Therm.kBAR1_;

    double rhoGamma = (ONE/(Gamma1*kBAR1))*pow(mL/alpha1_in,TWO);

    return pow(rhoGamma, ONE/(Gamma1 + ONE));
}

double EntropyJumpConditionDichotomy(\
        double& rho_inf, double& rho_sup, ThermoLaw& Therm,\
        Vector5d& W_state_L, double alpha1_R,\
        double& eps){

    //Local variables
    double alpha1_L, p1_L, rho1_L, u1_L, u2_L; 
           
    //Function
    alpha1_L = W_state_L(0);
    p1_L     = W_state_L(1);
    u1_L     = W_state_L(2);
    u2_L     = W_state_L(4);

    rho1_L   = Density_EOS(1, Therm, p1_L, ZERO);

    double mL = alpha1_L*rho1_L*(u1_L - u2_L);
    double HL = EntropyJumpFunction(mL, alpha1_L, rho1_L, Therm);

    //Initial density guess
    double rho_star(ZERO);
    double rho_mid  = (rho_inf + rho_sup)/TWO;

    double  Hinf = EntropyJumpFunction(mL, alpha1_R, rho_inf, Therm) - HL;
    double  Hmid = EntropyJumpFunction(mL, alpha1_R, rho_mid, Therm) - HL;
    double  Hsup = EntropyJumpFunction(mL, alpha1_R, rho_sup, Therm) - HL;

    if(Hinf*Hsup>= ZERO){

	cout<<"Alert EntropyJumpConditionDichotomy: H(rho_inf) H(rho_sup) of same sign !"<<endl;
	exit(EXIT_FAILURE);

    }

    else if(fabs(Hmid)<=eps){

        return rho_mid;

    }
    else{

        if(Hinf*Hmid<ZERO){

            rho_star = EntropyJumpConditionDichotomy(\
                    rho_inf, rho_mid, Therm,\
                    W_state_L, alpha1_R,\
                    eps);
        }
        else {

            rho_star = EntropyJumpConditionDichotomy(\
                    rho_mid, rho_sup, Therm,\
                    W_state_L, alpha1_R,\
                    eps);
        }

    }

    return rho_star;

}


double VelocityJumpFunction(\
        double alpha1_L, double rho1_L, double u1_L, double u2_L,\
        double alpha1_R, double rho1_R\
        ){

    return u2_L + (alpha1_L/alpha1_R)*(rho1_L/rho1_R)*(u1_L - u2_L);

}

double PressureJumpFunction(\
        double mL, double alpha1_L, double p1_L, double u1_L, double p2_L,\
        double alpha1_R, double p1_R, double u1_R\
        ){

    double alpha2_L = ONE - alpha1_L;
    double alpha2_R = ONE - alpha1_R;
    double P_L = alpha1_L*p1_L + alpha2_L*p2_L;

    return (ONE/alpha2_R)*(mL*(u1_L - u1_R) + P_L - alpha1_R*p1_R);

}

Vector5d WstateContactResolution(Vector5d W_state_L, double alpha1_R,\
        ThermoLaw& Therm, double eps){

    //Local variables
    Vector5d W_state_R;

    //Function
    double alpha1_L = W_state_L(0);
    double p1_L     = W_state_L(1);
    double rho1_L   = Density_EOS(1, Therm, p1_L, ZERO);
    double u1_L     = W_state_L(2);
    double p2_L     = W_state_L(3);
    double u2_L     = W_state_L(4);
    double mL       = alpha1_L*rho1_L*(u1_L - u2_L);

    //Dichotomy on the subsonic branch
    double rho_sup  = Rho1UpperBoundEntropyJumpFunction(mL, alpha1_R, Therm) + epsZero;

    //Entropy Jump Condition Resolution (rho1_R)
    double rho1_R = EntropyJumpConditionDichotomy(\
            rho_sup, Big, Therm,\
            W_state_L, alpha1_R,\
            epsDicho);

    double p1_R = Pressure_EOS(1, Therm, rho1_R, ZERO);

    //Mass1 Jump Resolution (u1_R)
   double u1_R = VelocityJumpFunction(\
        alpha1_L, rho1_L, u1_L, u2_L,\
        alpha1_R, rho1_R\
        );

    //Mixture momentum Jump Resolution (p2_R)
   double p2_R = PressureJumpFunction(\
        mL, alpha1_L, p1_L, u1_L, p2_L,\
        alpha1_R, p1_R, u1_R\
        ); 

   W_state_R<<alpha1_R,
              p1_R,
              u1_R,
              p2_R,
              u2_L;

   return W_state_R;

}

double TauUpperBoundPressureJumpFunction(int Phase_Id,\
        double Jshock,\
        ThermoLaw& Therm){

    //local variables
    double kBAR, Gamma;

    //Function
    if(Phase_Id == 1){
        kBAR  = Therm.kBAR1_;
        Gamma = Therm.Gamma1_;
    }
    else{
        kBAR  = Therm.kBAR2_;
        Gamma = Therm.Gamma2_;
    }

    double tau_res = pow( ((Jshock*Jshock)/(Gamma*kBAR)), -ONE/(Gamma+ONE) );

    return tau_res;

}

double PressureRHJumpFunction(\
        int Phase_Id,\
        double Jshock,\
        double tau_R,\
        ThermoLaw& Therm){

    //local variables
    double rho_R = ONE/tau_R;

    return Jshock*Jshock*tau_R  + Pressure_EOS(Phase_Id, Therm, rho_R, ZERO);

}

double RHJumpConditionDichotomy(\
        int& Phase_Id,\
        double& tau_inf, double& tau_sup, ThermoLaw& Therm,\
        Vector5d& W_state_R,\
        double& Jshock,\
        double& eps){

    //Local variables
    double p_R, tau_R;
    double PR;
    double Pinf, Pmid, Psup;
           
    //Function

    //Initial density guess
    double tau_star(ZERO);
    double tau_mid  = (tau_inf + tau_sup)/TWO;


    if(Phase_Id == 1){
        p_R     = W_state_R(1);
        tau_R   = ONE/Density_EOS(1, Therm, p_R, ZERO);
        PR      = PressureRHJumpFunction(1, Jshock, tau_R, Therm);
        Pinf    = PressureRHJumpFunction(1, Jshock, tau_inf, Therm) - PR;
        Pmid    = PressureRHJumpFunction(1, Jshock, tau_mid, Therm) - PR;
        Psup    = PressureRHJumpFunction(1, Jshock, tau_sup, Therm) - PR;
    }
    else{
        p_R     = W_state_R(3);
        tau_R   = ONE/Density_EOS(2, Therm, p_R, ZERO);
        PR      = PressureRHJumpFunction(2, Jshock, tau_R, Therm);
        Pinf    = PressureRHJumpFunction(2, Jshock, tau_inf, Therm) - PR;
        Pmid    = PressureRHJumpFunction(2, Jshock, tau_mid, Therm) - PR;
        Psup    = PressureRHJumpFunction(2, Jshock, tau_sup, Therm) - PR;
    }
    /*
    cout<<"Inside RHJumpConditionDichotomy, p_in = "<<p_R<<endl;
    cout<<"Inside RHJumpConditionDichotomy, tau_inf = "<<tau_inf<<", tau_R = "<<tau_R<<", tau_sup = "<<tau_sup<<endl;
    cout<<"Inside RHJumpConditionDichotomy, P_inf   = "<<Pinf<<", P_R = "<<PR<<", P_sup = "<<Psup<<endl;
*/


    if(Pinf*Psup>= ZERO){

	cout<<"Alert RHJumpConditionDichotomy: P(tau_inf) P(tau_sup) of same sign !"<<endl;
	exit(EXIT_FAILURE);

    }

    else if(fabs(Pmid)<=eps){

        return tau_mid;

    }
    else{

        if(Pinf*Pmid<ZERO){

            tau_star = RHJumpConditionDichotomy(\
                    Phase_Id,\
                    tau_inf, tau_mid, Therm,\
                    W_state_R,\
                    Jshock,\
                    eps);
        }
        else {

            tau_star = RHJumpConditionDichotomy(\
                    Phase_Id,\
                    tau_mid, tau_sup, Therm,\
                    W_state_R,\
                    Jshock,\
                    eps);
        }

    }

    return tau_star;
}

double VelocityRHJumpFunction(\
        double Jshock,\
        double tau_L,\
        double tau_R, double u_R\
        ){

    return u_R - Jshock*( tau_R - tau_L );


}

Vector5d WstateRHJumpCondition(\
        Vector5d& W_state_R,\
        double sigma1Shock, double sigma2Shock,\
        ThermoLaw& Therm, double eps){

    //Local variables
    double p_R, p_L, rho_R, tau_R, tau_L, u_R, u_L;
    double Jshock, tau_inf;
    Vector5d W_state_L;

    //Function
    double alpha_L  = W_state_R(0);
    double tau_sup  = Big;

    int Phase_Id = 1;

    p_R      = W_state_R(1);
    u_R      = W_state_R(2);
    rho_R    = Density_EOS(1, Therm, p_R, ZERO);
    tau_R    = ONE/rho_R;
    Jshock   = rho_R*(u_R - sigma1Shock);
    tau_inf  = TauUpperBoundPressureJumpFunction(Phase_Id, Jshock, Therm) + epsZero; 

    tau_L   = RHJumpConditionDichotomy(Phase_Id,\
            tau_inf, tau_sup, Therm,\
            W_state_R,\
            Jshock,\
            eps);

    p_L     = Pressure_EOS(Phase_Id, Therm, ONE/tau_L, ZERO);
    u_L     = VelocityRHJumpFunction(Jshock, tau_L, tau_R, u_R);

    W_state_L(0) = alpha_L;
    W_state_L(1) = p_L;
    W_state_L(2) = u_L;

    Phase_Id = 2;

    p_R      = W_state_R(3);
    u_R      = W_state_R(4);
    rho_R    = Density_EOS(2, Therm, p_R, ZERO);
    tau_R    = ONE/rho_R;
    Jshock   = rho_R*(u_R - sigma2Shock);
    tau_inf  = TauUpperBoundPressureJumpFunction(Phase_Id, Jshock, Therm) + epsZero; 

    cout<<"tau_inf = "<<tau_inf<<", tau_in = "<<tau_R<<", tau_in - tau_inf = "<<tau_R - tau_inf<<endl;

    tau_L   = RHJumpConditionDichotomy(Phase_Id,\
            tau_inf, tau_sup, Therm,\
            W_state_R,\
            Jshock,\
            eps);

    p_L     = Pressure_EOS(Phase_Id, Therm, ONE/tau_L, ZERO);
    u_L     = VelocityRHJumpFunction(Jshock, tau_L, tau_R, u_R);

    W_state_L(3) = p_L;
    W_state_L(4) = u_L;

    return W_state_L;

}


Vector5d IsentropicBN_Eigenvalues(\
                Vector5d& W_state_avr,\
                ThermoLaw& Therm\
                ){

    //Local variables
    Vector5d EigenValues;
    double rho1, p1, rho2, p2;
    double c1, c2;
    double u1, u2;

    //Function
    p1   = W_state_avr(1);
    rho1 = Density_EOS(1, Therm, p1, ZERO);
    c1   = Sound_Speed_EOS(1, Therm, rho1, p1);
    p2   = W_state_avr(3);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c2   = Sound_Speed_EOS(2, Therm, rho2, p2);

    u1   = W_state_avr(2);
    u2   = W_state_avr(4);

    EigenValues<<u2,
                 u1+c1,
                 u1-c1,
                 u2+c2,
                 u2-c2;

    return EigenValues;

}

Matrix5d IsentropicBN_Eigenvectors(\
        string VariableType,\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        ){

    Matrix5d EigenVectorsBase;

    if (VariableType=="ConsVar"){

        EigenVectorsBase = IsentropicBN_EigenvectorsCons(\
                W_state_avr,\
                Therm\
                ); 
    }

    return EigenVectorsBase;
}

Matrix5d IsentropicBN_EigenvectorsCons(\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        ){

    //Local variables
    double rho1, p1, rho2, p2;
    double c1, c2, C1, C2;
    double u1, u2;

    Matrix5d EigenVectorsBase;

    //Function
    p1   = W_state_avr(1);
    rho1 = Density_EOS(1, Therm, p1, ZERO);
    c1   = Sound_Speed_EOS(1, Therm, rho1, p1);
    p2   = W_state_avr(3);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c2   = Sound_Speed_EOS(2, Therm, rho2, p2);

    u1   = W_state_avr(2);
    u2   = W_state_avr(4);

    C1 = rho1*c1*c1; 
    C2 = rho2*c2*c2; 

    double du     = u2 - u1;
    double dp     = p2 - p1;

    double res    = c1*c1 - du*du;
    double Coeff1 = C1/res;
    double Coeff2 = (dp - C2)/(c2*c2);

    EigenVectorsBase<<ONE, ZERO , ZERO , ZERO , ZERO ,
                     Coeff1, ONE ,ONE , ZERO , ZERO ,
                     Coeff1*u2, (u1 + c1), (u1-c1) , ZERO, ZERO,
                     Coeff2, ZERO , ZERO , ONE , ONE ,
                     Coeff2*u2, ZERO , ZERO ,(u2+c2) ,(u2-c2);

    return EigenVectorsBase;

}

Matrix5d IsentropicBN_EigenvectorsInv(\
        string VariableType,\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        ){

    Matrix5d EigenVectorsBaseInv;

    if (VariableType=="ConsVar"){

        EigenVectorsBaseInv = IsentropicBN_EigenvectorsInvCons(\
                W_state_avr,\
                Therm\
                ); 
    }

    return EigenVectorsBaseInv;
}

Matrix5d IsentropicBN_EigenvectorsInvCons(\
        Vector5d& W_state_avr,\
        ThermoLaw& Therm\
        ){

    //Local variables
    Matrix5d EigenVectorsBaseInv;
    double rho1, p1, rho2, p2;
    double c1, c2, C1, C2;
    double u1, u2;
    double du, dp;

    //Function

    //Averaged quantities
    p1   = W_state_avr(1);
    rho1 = Density_EOS(1, Therm, p1, ZERO);
    c1   = Sound_Speed_EOS(1, Therm, rho1, p1);
    p2   = W_state_avr(3);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c2   = Sound_Speed_EOS(2, Therm, rho2, p2);

    u1   = W_state_avr(2);
    u2   = W_state_avr(4);

    C1 = rho1*c1*c1; 
    C2 = rho2*c2*c2; 

    du = u2 - u1;
    dp = p2 - p1;

    //Row alpha1

    RowVector5d R_alpha1;
    R_alpha1 << TWO, ZERO, ZERO, ZERO, ZERO;

    //Row m1
    double Ca1 = -(ONE/c1)*(C1/(c1-du));
    double Cm  = -(u1-c1)/c1;
    double Cmu = ONE/c1;

    RowVector5d R_m1;
    R_m1 << Ca1, Cm, Cmu, ZERO, ZERO;

    //Row m1u1

    Ca1 = -(ONE/c1)*(C1/(c1+du));
    Cm  = (u1+c1)/c1;
    Cmu = -ONE/c1;

    RowVector5d R_m1u1;
    R_m1u1 << Ca1, Cm, Cmu, ZERO, ZERO;

    //Row m2

    Ca1 = -(ONE/(c2*c2))*(dp - C2);
    Cm  = -(u2-c2)/c2;
    Cmu = ONE/c2;

    RowVector5d R_m2;
    R_m2 << Ca1, ZERO, ZERO, Cm, Cmu;

    //Row m2u2

    Ca1 = -(ONE/(c2*c2))*(dp - C2);
    Cm  = (u2 + c2)/c2;
    Cmu = -ONE/c2;

    RowVector5d R_m2u2;
    R_m2u2 << Ca1, ZERO, ZERO, Cm, Cmu;

    //Building of the matrix
    EigenVectorsBaseInv.row(0) = R_alpha1;
    EigenVectorsBaseInv.row(1) = R_m1;
    EigenVectorsBaseInv.row(2) = R_m1u1;
    EigenVectorsBaseInv.row(3) = R_m2;
    EigenVectorsBaseInv.row(4) = R_m2u2;

    return ONE_OVER_TWO*EigenVectorsBaseInv;
}

Vector5d IsentropicBN_EigenvectorsBaseProjection(\
        string VariableType,\
        Vector5d& W_state_avr,\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        ){

    Vector5d CoordsEigenVectorsBase;

    if (VariableType=="ConsVar"){

        CoordsEigenVectorsBase = IsentropicBN_EigenvectorsBaseProjectionCons(\
                W_state_avr,\
                W_state_L, W_state_R,\
                Therm\
                ); 
    }

    return CoordsEigenVectorsBase;
}

Vector5d IsentropicBN_EigenvectorsBaseProjectionCons(\
        Vector5d& W_state_avr,\
        Vector5d& W_state_L, Vector5d& W_state_R,\
        ThermoLaw& Therm\
        ){

    //Local variables
    Vector5d CoordsEigenVectorsBase;
    double rho1, p1, rho2, p2;
    double c1, c2, C1, C2;
    double u1, u2;

    //Function

    //U_R - U_L
    Vector5d U_R = NConsVarToConsVarLoc(W_state_R, Therm);
    Vector5d U_L = NConsVarToConsVarLoc(W_state_L, Therm);
    Vector5d dU = U_R - U_L;

    //Averaged quantities
    p1   = W_state_avr(1);
    rho1 = Density_EOS(1, Therm, p1, ZERO);
    c1   = Sound_Speed_EOS(1, Therm, rho1, p1);
    p2   = W_state_avr(3);
    rho2 = Density_EOS(2, Therm, p2, ZERO);
    c2   = Sound_Speed_EOS(2, Therm, rho2, p2);

    u1   = W_state_avr(2);
    u2   = W_state_avr(4);

    C1 = rho1*c1*c1; 
    C2 = rho2*c2*c2; 

    double du     = u2 - u1;
    double dp     = p2 - p1;

    //Row alpha1, m1

    double Ca1 = -(ONE/c1)*(C1/(c1-du));
    double Cm  = -(u1-c1)/c1;
    double Cmu = ONE/c1;

    CoordsEigenVectorsBase(0) = TWO*dU(0);
    CoordsEigenVectorsBase(1) = Ca1*dU(0) + Cm*dU(1) + Cmu*dU(2);

    //Row m1u1

    Ca1 = -(ONE/c1)*(C1/(c1+du));
    Cm  = (u1+c1)/c1;
    Cmu = -ONE/c1;

    CoordsEigenVectorsBase(2) = Ca1*dU(0) + Cm*dU(1) + Cmu*dU(2);

    //Row m2

    Ca1 = -(ONE/(c2*c2))*(dp - C2);
    Cm  = -(u2-c2)/c2;
    Cmu = ONE/c2;

    CoordsEigenVectorsBase(3) = Ca1*dU(0) + Cm*dU(3) + Cmu*dU(4);

    //Row m2u2

    Ca1 = -(ONE/(c2*c2))*(dp - C2);
    Cm  = (u2 + c2)/c2;
    Cmu = -ONE/c2;

    CoordsEigenVectorsBase(4) = Ca1*dU(0) + Cm*dU(3) + Cmu*dU(4);

    return ONE_OVER_TWO*CoordsEigenVectorsBase;
}

double Heaviside(double x, double time, double lambda){

    //Local variables

    //Function
    if(x <= lambda*time){

        return ZERO;
    }
    else{

        return ONE;
    }

    return Big;
}

double AkRelaxDirac(double x, double t,\
        double uI_Frozen, double lambda_k,\
        double ak, double ak_eq,\
        double dirac_weight, double tauvec){

    /*
    double ak_star = ak + dirac_weight/fabs(lambda_k - uI_Frozen);
    double tcofac = tauvec*(x - uI_Frozen*t)/(lambda_k - uI_Frozen);

    double result = ak_star*exp(-tcofac) + ak_eq*(ONE - exp(-tcofac));
    */

    double result = ZERO;
    double vmax = max(uI_Frozen, lambda_k);
    double vmin = min(uI_Frozen, lambda_k);

    if( x< vmax*t && x> vmin*t){

        double tcofac  = tauvec*(x -lambda_k*t)/(uI_Frozen - lambda_k);
        double tcofac2 = tauvec*(x - uI_Frozen*t)/(lambda_k - uI_Frozen);

        double ak_relax = ak*exp(-tcofac) + ak_eq*(ONE - exp(-tcofac)); 
        double ak_star  = ak_relax + dirac_weight/fabs(lambda_k - uI_Frozen);

        result = ak_star*exp(-tcofac2) + ak_eq*(ONE - exp(-tcofac2));
    }

    return result;

}

double AkRelaxDiracTest(double x, double t,\
        double uI_Frozen, double lambda_k,\
        double dirac_weight, double tauvec){

    double result = ZERO;
    double vmax = max(uI_Frozen, lambda_k);
    double vmin = min(uI_Frozen, lambda_k);
    double tcofac  = ZERO;

    if( x< vmax*t && x> vmin*t){

        tcofac = tauvec*(x - uI_Frozen*t)/(lambda_k - uI_Frozen);
        result = dirac_weight*exp(-tcofac);
    }

    return result;

}


Matrix5d TimeRelaxMat(Vector5d& U_state_ref, Vector5d& TimeWeight){

    Vector5d Result = TimeWeight.cwiseQuotient(U_state_ref.cwiseProduct(U_state_ref));
    Matrix5d RelaxDiagMat = DiagonalMatrix<double, 5>(Result); 

    return RelaxDiagMat;
}

Matrix5d QuadraticDistance(Vector5d& A_coords){

    Matrix5d DistanceDiagMat = DiagonalMatrix<double, 5>(A_coords); 

    return DistanceDiagMat*DistanceDiagMat;
}

Matrix5d QuadraticDistanceInv(Matrix5d& TimeMat, Vector5d& A_coords, double dtRelax){

    //Local variables
    int nrows = A_coords.rows();
    Vector5d One = VectorXd::Constant(nrows, ONE);
    Vector5d A_coords_sq  = A_coords.cwiseProduct(A_coords);

    //Function

    Vector5d A_coords_inv = ((One+THREE*dtRelax*TimeMat*A_coords_sq).array()).pow(-ONE);
    Matrix5d DistanceDiagMat = DiagonalMatrix<double, 5>(A_coords_inv); 

    return DistanceDiagMat;
}

Vector5d CubicVectorDistance(Vector5d& A_coords){

    Vector5d DistanceDiagMat = ((A_coords).array()).pow(THREE); 

    return DistanceDiagMat;
}
