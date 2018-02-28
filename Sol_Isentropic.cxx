#include "Sol_Isentropic.h"

//constructor:

Sol_Isen::Sol_Isen(Mesh& M,\
		ThermoLaw& Therm,\
        double pRef, double mRef,\
        Vector3d tauRelax,\
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
    
    //Initializing the ConsVarDual_ to ZERO 
    NConsVarDual_ = MatrixXd::Zero(NcellExt,5);

    //Initializing the conservative flux terms
    ConsVarFlux_ = MatrixXd::Zero(Nfaces,5);

    //Initializing the non-conservative terms
    NConsVarFlux_ = MatrixXd::Zero(NcellExt,5);

    //Initializing the face values used for the non-conservative terms
    NConsVarFace_ = MatrixXd::Zero(Nfaces,1);

    //Initializing the exact solution
    SolExact_        = NConsVar_;
    SolExactEqRelax_ = NConsVarEqRelax_;

}

Sol_Isen::Sol_Isen( Sol_Isen& solution){

    LeftBCType_   =solution.LeftBCType_;
    RightBCType_  =solution.RightBCType_;

    SolTherm_=solution.SolTherm_;

    ConsVar_            = solution.ConsVar_;
    NConsVar_           = solution.NConsVar_;
    NConsVarDual_       = solution.NConsVarDual_;
    NConsVarEqRelax_    = solution.NConsVarEqRelax_;
    SoundSpeed_         = solution.SoundSpeed_;
    Mach_               = solution.Mach_;
    ConsVarFlux_        = solution.ConsVarFlux_;
    NConsVarFlux_       = solution.NConsVarFlux_;
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
        ((*this).SourceTerms_)(i,0) = -(Kp_Coeff(pRef, alpha1)/tauP)*dp;
        ((*this).SourceTerms_)(i,1) = ZERO;
        ((*this).SourceTerms_)(i,2) = +(Ku_Coeff(mRef, alpha1)/tauU)*du;
        ((*this).SourceTerms_)(i,3) = ZERO;
        ((*this).SourceTerms_)(i,4) = -(Ku_Coeff(mRef, alpha1)/tauU)*du;
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

void Sol_Isen::NConsVarToEqRelax(){

    //local variables
    int nrows = NConsVar_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

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

void Sol_Isen::ConsVarToNConsVar(){

    //local variables
    int nrows = NConsVar_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

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

void Sol_Isen::EqRelaxToNConsVar(){

    //local variables
    int nrows = NConsVarEqRelax_.rows();

    //Building the mass fraction vectors Y1, Y2:

    //Building the mass vectors m1, m2, m
    VectorXd ZeroTab = VectorXd::Zero(nrows);

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
    NConsVarEqRelax_.col(0)  = alpha1;

    //p1
    NConsVarEqRelax_.col(1)  = p1;

    //u1
    NConsVarEqRelax_.col(2)  = U - Y2.cwiseProduct(du);

    //p2
    NConsVarEqRelax_.col(3)  = p2;

    //dp
    NConsVarEqRelax_.col(4)  = U + Y1.cwiseProduct(du);

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

/************************************************/
/*************  CHANGE OF VARIABLE  *************/
/************************************************/

Vector5d NConsVarToConsVar(Vector5d& W_state,\
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

Vector5d EqRelaxToNConsVar(Vector5d& V_state,\
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
    double rho_sup  = Rho1UpperBoundEntropyJumpFunction(mL, alpha1_R, Therm) - epsZero;

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
