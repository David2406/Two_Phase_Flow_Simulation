#include "Sol.h"

//constructor:

Sol::Sol(Mesh& M,\
		ThermoLaw& Therm,\
        Vector7d& InitL, Vector7d& InitR,\
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

    //Initializing over the overall mesh
    NConsVarInit(InitR, InitL,\
            NcellExt, Ncells, Length);
    ConsVarInit(InitR, InitL,\
            NcellExt, Ncells, Length);

    //Initializing the L1 errors 
    errL1_<<Big,
            Big,
            Big,
            Big,
            Big,
            Big,
            Big;
    
    //Initializing the conservative flux terms
    ConsVarFlux_ = MatrixXd::Zero(Nfaces,7);

    //Initializing the non-conservative terms
    NConsVarFlux_ = MatrixXd::Zero(Nfaces,7);

    //Initializing the source terms
    SourceTerms_ = MatrixXd::Zero(NcellExt,7);

    //Initializing the exact solution
    SolExact_ = NConsVar_;

}

Sol::Sol( Sol& solution){

    LeftBCType_   =solution.LeftBCType_;
    RightBCType_  =solution.RightBCType_;

    SolTherm_=solution.SolTherm_;

    ConsVar_      = solution.ConsVar_;
    NConsVar_     = solution.NConsVar_;
    SoundSpeed_   = solution.SoundSpeed_;
    Mach_         = solution.Mach_;
    ConsVarFlux_  = solution.ConsVarFlux_;
    NConsVarFlux_ = solution.NConsVarFlux_;
    SourceTerms_  = solution.SourceTerms_;

}

Sol::Sol(){

    LeftBCType_   ="None";
    RightBCType_  ="None";

}

//methods:

void Sol::NConsVarInit(Vector7d& InitR, Vector7d& InitL,\
        int NcellExt, int Ncells, double Length){

    //local variables
    double rho1, u1, p1, c1; 
    double rho2, u2, p2, c2;

    double SpaceStep=Length/Ncells;
    double x_0 = (*this).x_0_;
    double Ncols = InitR.rows();

    //Resizing the NConsVar_ array
    NConsVar_.resize(NcellExt, Ncols);
    SoundSpeed_.resize(NcellExt, 2);
    Mach_.resize(NcellExt, 2);

    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

        if(i<=Riemann_index){

            ((*this).NConsVar_).block(i,0,1,Ncols) = InitL.transpose();

            rho1   = InitL(1);
            u1     = InitL(2);
            p1     = InitL(3);
            rho2   = InitL(4);
            u2     = InitL(5);
            p2     = InitL(6);
        }

        else{

            ((*this).NConsVar_).block(i,0,1,Ncols) = InitR.transpose(); 

            rho1   = InitR(1);
            u1     = InitR(2);
            p1     = InitR(3);
            rho2   = InitR(4);
            u2     = InitR(5);
            p2     = InitR(6);
        }
            c1 = Sound_Speed_EOS(1,(*this).SolTherm_, rho1, p1);
            c2 = Sound_Speed_EOS(2,(*this).SolTherm_, rho2, p2);

            ((*this).SoundSpeed_)(i,0) = c1; 
            ((*this).SoundSpeed_)(i,1) = c2;

            ((*this).Mach_)(i,0) = fabs(u1)/c1; 
            ((*this).Mach_)(i,1) = fabs(u2)/c2; 
    }
}

void Sol::ConsVarInit(Vector7d& InitR, Vector7d& InitL,\
        int NcellExt, int Ncells, double Length){

    double SpaceStep=Length/Ncells;
    double x_0 = (*this).x_0_;
    double Ncols = InitR.rows();

    //local variables
    double alpha1, rho1, u1 , p1, e_int1; 
    double alpha2, rho2, u2 , p2, e_int2;

    double m1, m2;

    //Resizing the NConsVar_ array
    ConsVar_.resize(NcellExt, Ncols);

    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

        if(i<=Riemann_index){

            alpha1 = InitL(0);
            rho1   = InitL(1);
            u1     = InitL(2);
            p1     = InitL(3);
            e_int1 = Internal_Energy_EOS(\
                    1,\
                    (*this).SolTherm_, \
                    rho1,p1);
            alpha2 = ONE-InitL(0);
            rho2   = InitL(4);
            u2     = InitL(5);
            p2     = InitL(6);
            e_int2 = Internal_Energy_EOS(\
                    2,\
                    (*this).SolTherm_, \
                    rho2,p2);

            m1     = alpha1*rho1;
            m2     = alpha2*rho2;

        }

        else{

            alpha1 = InitR(0);
            rho1   = InitR(1);
            u1     = InitR(2);
            p1     = InitR(3);
            e_int1 = Internal_Energy_EOS(\
                    1,\
                    (*this).SolTherm_, \
                    rho1,p1);
            alpha2 = ONE-InitR(0);
            rho2   = InitR(4);
            u2     = InitR(5);
            p2     = InitR(6);
            e_int2 = Internal_Energy_EOS(\
                    2,\
                    (*this).SolTherm_, \
                    rho2,p2);

            m1     = alpha1*rho1;
            m2     = alpha2*rho2;
        }

        //Conservative Vector Filled
        ((*this).ConsVar_)(i,0) = alpha1;
        ((*this).ConsVar_)(i,1) = m1;
        ((*this).ConsVar_)(i,2) = m1*u1;
        ((*this).ConsVar_)(i,3) = m1*(pow(u1,TWO)/TWO + e_int1);
        ((*this).ConsVar_)(i,4) = m2;
        ((*this).ConsVar_)(i,5) = m2*u2;
        ((*this).ConsVar_)(i,6) = m2*(pow(u2,TWO)/TWO + e_int2);
    }
}

//Update methods:

void Sol::SoundSpeed_Update(){

    //rho1, p1
    VectorXd Rho = NConsVar_.col(1);
    VectorXd P   = NConsVar_.col(3);

    SoundSpeed_.col(0) = Sound_Speed_EOS_Tab(1, (*this).SolTherm_, Rho,P);                         
    //rho2, p2
    Rho = NConsVar_.col(4);
    P   = NConsVar_.col(6);

    SoundSpeed_.col(1) = Sound_Speed_EOS_Tab(2, (*this).SolTherm_, Rho,P);                         
}

void Sol::Mach_Update(){

    //u1
    VectorXd U   = NConsVar_.col(2);
    Mach_.col(0) = U.cwiseQuotient(SoundSpeed_.col(0));                         

    //u2
    U   = NConsVar_.col(5);
    Mach_.col(1) = U.cwiseQuotient(SoundSpeed_.col(1));                         
}

/*
void Sol::ConsVar_Update( Mesh& mesh, double TimeStep, \
        string SchemeType\
        ){

    //Mesh parameters:
    int NcellExt=mesh.Get_NcellExt();
    int Nfaces=mesh.Get_Nfaces();
    int LeftBCface=mesh.Get_LeftBCface();
    int RightBCface=mesh.Get_RightBCface();
    double SpaceStep=mesh.Get_SpaceStep();
    int L,R;

    string LeftBCType=(*this).Get_LeftBCType();
    string RightBCType=(*this).Get_RightBCType();

    //For the JMH splitting conservative variables are classically updated
    if(!(SchemeType1=="GIRARDIN" && SchemeType2=="GIRARDIN")){

        for(int face_id=0; face_id<Nfaces; face_id++){ //Loop on all the faces, including these located at the boundaries


            L=mesh.FaceIndex_(face_id,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(face_id,2); //Index of the cell located to the right

            (*this).ConsVar_.row(L)+=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_.row(face_id);
            //Conservativity is ensured 
            (*this).ConsVar_.row(R)-=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_.row(face_id);

        }

        //Special treatment for boundary faces:

        if(LeftBCType=="transparent"){

            L=mesh.FaceIndex_(LeftBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(LeftBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(L)=(*this).ConsVar_.row(R); 
        }
        if(RightBCType=="transparent"){

            L=mesh.FaceIndex_(RightBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(RightBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(R)=(*this).ConsVar_.row(L); 
        }

    }

    //In the case Girardin's splitting and acoustic step
    else if(SchemeType2=="GIRARDIN" && SplittingStep2==true){

        //Loop on all the faces, including these located at the boundaries
        for(int face_id=0; face_id<Nfaces; face_id++){ 

            L=mesh.FaceIndex_(face_id,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(face_id,2); //Index of the cell located to the right

            //building 1+delta_t/delta_x * (u_faceR- u_faceL)
            (*this).NLGirardinOp_(L)+=(TimeStep/SpaceStep)*ConsVarFlux_(face_id,0);
            (*this).NLGirardinOp_(R)-=(TimeStep/SpaceStep)*ConsVarFlux_(face_id,0);

            //Momentum first update
            (*this).ConsVar_(L,1)+=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_(face_id,1);
            (*this).ConsVar_(R,1)-=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_(face_id,1);

            //Energy first uppdate
            (*this).ConsVar_(L,2)+=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_(face_id,2);
            (*this).ConsVar_(R,2)-=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_(face_id,2);

	    //No Update for mass fraction Y
        }
        //Loop on all cells, including these located on the boundaries
        for (int cell_id=0; cell_id<NcellExt; cell_id++){

            //Mass Update:
            (*this).ConsVar_(cell_id,0)/=NLGirardinOp_(cell_id);
            //Momentum Update:
            (*this).ConsVar_(cell_id,1)/=NLGirardinOp_(cell_id);
            //Energy Update:
            (*this).ConsVar_(cell_id,2)/=NLGirardinOp_(cell_id);
            //Y Update:
            (*this).ConsVar_(cell_id,3)/=NLGirardinOp_(cell_id);

        }

        //Special treatment for boundary faces:

        if(LeftBCType=="transparent"){

            L=mesh.FaceIndex_(LeftBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(LeftBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(L)=(*this).ConsVar_.row(R); 
        }
        if(RightBCType=="transparent"){

            L=mesh.FaceIndex_(RightBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(RightBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(R)=(*this).ConsVar_.row(L); 
        }
    }
    //In the case Girardin's splitting and convective step
    else if(SchemeType1=="GIRARDIN" && SplittingStep1==true){

        //Loop on all cells, including these located on the boundaries
        for (int cell_id=0; cell_id<NcellExt; cell_id++){

            //Mass Update:
            (*this).ConsVar_(cell_id,0)*=NLGirardinOp_(cell_id);
            //Momentum Update:
            (*this).ConsVar_(cell_id,1)*=NLGirardinOp_(cell_id);
            //Energy Update:
            (*this).ConsVar_(cell_id,2)*=NLGirardinOp_(cell_id);
            //Y Update:
            (*this).ConsVar_(cell_id,3)*=NLGirardinOp_(cell_id);


        }
        //Loop on all the faces, including these located at the boundaries
        for(int face_id=0; face_id<Nfaces; face_id++){ 

            L=mesh.FaceIndex_(face_id,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(face_id,2); //Index of the cell located to the right

            (*this).ConsVar_.row(L)+=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_.row(face_id);
            //Conservativity is ensured 
            (*this).ConsVar_.row(R)-=-(TimeStep/SpaceStep)*(*this).ConsVarFlux_.row(face_id);

        }

        //Special treatment for boundary faces:

        if(LeftBCType=="transparent"){

            L=mesh.FaceIndex_(LeftBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(LeftBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(L)=(*this).ConsVar_.row(R); 
        }
        if(RightBCType=="transparent"){

            L=mesh.FaceIndex_(RightBCface,1); //Index of the cell located to the left
            R=mesh.FaceIndex_(RightBCface,2); //Index of the cell located to the right
            //The Fictitious Boundary State is set equal to the Physical State
            (*this).ConsVar_.row(R)=(*this).ConsVar_.row(L); 
        }
    }
}
*/

/*void Sol::NConsVar_Update(){

    // rho update
    (*this).NConsVar_.col(0)=(*this).ConsVar_.col(0);
    VectorXd rho=(*this).ConsVar_.col(0);

    // u update
    (*this).NConsVar_.col(1)=(*this).ConsVar_.col(1).cwiseQuotient((*this).ConsVar_.col(0)); 
    //Total energy by unit mass
    VectorXd e_T=(*this).ConsVar_.col(2).cwiseQuotient((*this).ConsVar_.col(0));
    //Kinetic energy by unit mass
    VectorXd e_k=(ONE_OVER_TWO)*((*this).NConsVar_.col(1)).array().square();
    //Internal energy by unit mass
    VectorXd e=e_T-e_k;

    //pressure update
    VectorXd p=Pressure_EOS_Tab((*this).SolTherm_, rho, e);
    (*this).NConsVar_.col(2)=p;

    // Y update
    (*this).NConsVar_.col(3)=(*this).ConsVar_.col(3).cwiseQuotient((*this).ConsVar_.col(0)); 
}


void Sol::SolExact_Update(Mesh& mesh, double time){

    string SolType=(*this).SolType_;

    double rho0_L=(*this).rho0_L_;
    double rho0_R=(*this).rho0_R_;

    double u0_L=(*this).u0_L_;
    double u0_R=(*this).u0_R_;

    double p0_L=(*this).p0_L_;
    double p0_R=(*this).p0_R_;

    double Y0_L=(*this).Y0_L_;
    double Y0_R=(*this).Y0_R_;

    double c0_L=(*this).c0_L_;
    double c0_R=(*this).c0_R_;

    //Double Riemann Problem
    double rho_L_starstar, rho_R_starstar;
    double p_starstar, u_starstar;
    double c_L_starstar;
    double sigma_L_starstar, sigma_R_starstar;
    double rho0_interm;
    double u0_interm  ;
    double p0_interm  ;
    double Y0_interm  ;
    double c0_interm  ;
    double x_1        ;

    if((*this).Nbr_Areas_> 2){

        rho0_interm = (*this).rho0_interm_;
        u0_interm   = (*this).u0_interm_;
        p0_interm   = (*this).p0_interm_;
        Y0_interm   = (*this).Y0_interm_;
        c0_interm   = (*this).c0_interm_;
        x_1         = (*this).x_1_;

    }

    double xcell;
    double x_0=(*this).x_0_;

    int NcellExt     = mesh.Get_NcellExt();
    double SpaceStep = mesh.Get_SpaceStep();

    //Exact solution states:
    double rho_L_star, rho_R_star;
    double p_star, u_star;
    double c_L_star, c_R_star;
    double sigma_L_star, sigma_R_star; //Shock or Rarefaction wave velocities

    if(SolType!="Pure Contact" && SolType!="Pure Shock" && SolType!="Pure Rarefaction" &&\
	    SolType!="Sinus Impulse" && SolType!="Double Riemann Problem"){

	//In the case of a pure contact discontinuity we don't need to compute p_star

	double p_inf=ZERO;
	double p_sup=HUNDRED*max(p0_L,p0_R);
	p_star=P_Euler_Dichotomy(p_inf, p_sup, SolType, (*this).SolTherm_, rho0_L, rho0_R,\
		u0_L, u0_R, p0_L, p0_R, c0_L, c0_R, epsDicho);
    }

    if (SolType=="Pure Contact"){
	//The theoretical sol is maid of one contact discontinuity only
	//its speed should be u=u0_L=u0_R

	if(fabs(u0_L-u0_R)>epsZero){

	    cout<<"Alert SolExact_Update: Pure Contact discontinuity but u0_L != u0_R"\
		<<endl;
	    exit(EXIT_FAILURE);
	}

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);
	    //cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;

	    if(xcell<=x_0+time*u0_L){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }


	}

    }
    else if (SolType=="Pure Shock"){

	//Shock Speed
	double sigma=c0_L;

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);

	    if(xcell<=x_0+time*sigma){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }

	}

    }
    else if (SolType=="Pure Rarefaction"){

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);

	    if(xcell<=x_0+time*(u0_L-c0_L)){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if (xcell<=x_0+time*(u0_R-c0_R)){

		double X=xcell-x_0;
		Vector3d State_L_Fan=State_k_Fan(rho0_L, u0_L, p0_L, c0_L,\
			(*this).SolTherm_, 1, \
			X, time); 

		(*this).SolExact_(cell,0)=State_L_Fan(0);
		(*this).SolExact_(cell,1)=State_L_Fan(1);
		(*this).SolExact_(cell,2)=State_L_Fan(2);
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }

	}

    }
    else if (SolType=="Sinus Impulse"){

	int    NPI_sin   = (*this).NPI_sin_;
	int    NPI_gauss = (*this).NPI_gauss_;
	double M0        = (*this).M0_;
	double sigma = NPI_gauss * SpaceStep;
	double k_c   = (TWO*Pi) / (NPI_sin*SpaceStep);

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);

	    if(xcell<=x_0+time*u0_L){

		(*this).SolExact_(cell,0)=\
		rho0_L * (ONE + pow(M0, TWO) * sin(k_c * \
			(xcell-(x_0+time*u0_L-k_c*Pi/TWO))) *\
		    (ONE/(sigma*sqrt(TWO*Pi)))*\
		    exp(-ONE_OVER_TWO * pow((xcell - (x_0+time*u0_L))/sigma, TWO)));

		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else{

		(*this).SolExact_(cell,0)=\
		rho0_L * (ONE + pow(M0, TWO) * sin(k_c * \
			(xcell-(x_0+time*u0_L-k_c*Pi/TWO))) *\
		    (ONE/(sigma*sqrt(TWO*Pi)))*\
		    exp(-ONE_OVER_TWO * pow((xcell - (x_0+time*u0_L))/sigma, TWO)));
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_R;
	    }


	}
    }
    else if (SolType=="Double Riemann Problem"){

        string SubSolType = "Left Rarefaction Right Shock";

        //1st Riemann Problem: p_star
        double p_inf=ZERO;
        double p_sup=HUNDRED*max(p0_L,p0_interm);
        p_star=P_Euler_Dichotomy(p_inf, p_sup, SubSolType,\
                (*this).SolTherm_,\
                rho0_L, rho0_interm,\
                u0_L, u0_interm,\
                p0_L, p0_interm,\
                c0_L, c0_interm, epsDicho);

        //u star estimation
        u_star=u0_L-Euler_k_Resolvant(p_star, (*this).SolTherm_, \
                rho0_L, u0_L, p0_L, c0_L);

        //Left quantities
        rho_L_star=Rho_k_star_Resolvant(rho0_L, p0_L, p_star, (*this).SolTherm_);
        c_L_star=C_k_star_Resolvant( c0_L, p0_L, p_star, (*this).SolTherm_);

        sigma_L_star=u_star-c_L_star;

        //Right quantities
        rho_R_star=Rho_k_star_Resolvant(rho0_interm, p0_interm, p_star, (*this).SolTherm_);

        sigma_R_star=Sigma_k_star_Resolvant(u0_interm, c0_interm, p0_interm, p_star,\
                2, (*this).SolTherm_ );

        //2nd Riemann Problem: p_starstar
        p_inf=ZERO;
        p_sup=HUNDRED*max(p0_interm,p0_R);
        p_starstar=P_Euler_Dichotomy(p_inf, p_sup, SubSolType,\
                (*this).SolTherm_,\
                rho0_interm, rho0_R,\
                u0_interm, u0_R,\
                p0_interm, p0_R,\
                c0_interm, c0_R, epsDicho);

        //u starstar estimation
        u_starstar=u0_interm-Euler_k_Resolvant(p_starstar, (*this).SolTherm_, \
                rho0_interm, u0_interm, p0_interm, c0_interm);

        //Left quantities
        rho_L_starstar=Rho_k_star_Resolvant(rho0_interm, p0_interm, p_starstar, (*this).SolTherm_);
        c_L_starstar=C_k_star_Resolvant( c0_interm, p0_interm, p_starstar, (*this).SolTherm_);

        sigma_L_starstar=u_starstar-c_L_starstar;

        //Right quantities
        rho_R_starstar=Rho_k_star_Resolvant(rho0_R, p0_R, p_starstar, (*this).SolTherm_);

        sigma_R_starstar=Sigma_k_star_Resolvant(u0_R, c0_R, p0_R, p_starstar,\
                2, (*this).SolTherm_ );

        for (int cell=0; cell<NcellExt; cell++){

            xcell=mesh.CellCoordsTab_(cell,1);

            if(xcell<=x_0+(u0_L-c0_L)*time){

                (*this).SolExact_(cell,0)=rho0_L;
                (*this).SolExact_(cell,1)=u0_L;
                (*this).SolExact_(cell,2)=p0_L;
                (*this).SolExact_(cell,3)=Y0_L;

            }
            else if(xcell<=x_0+sigma_L_star*time){
                //Inside the Left Rarefaction wave

                double X=xcell-x_0;
                Vector3d State_L_Fan=State_k_Fan(rho0_L, u0_L, p0_L, c0_L,\
                        (*this).SolTherm_, 1, \
                        X, time); 

                (*this).SolExact_(cell,0)=State_L_Fan(0);
                (*this).SolExact_(cell,1)=State_L_Fan(1);
                (*this).SolExact_(cell,2)=State_L_Fan(2);
                (*this).SolExact_(cell,3)=Y0_L;

            }
            else if(xcell<=x_0+u_star*time){
                //Inside the Left star region

                (*this).SolExact_(cell,0)=rho_L_star;
                (*this).SolExact_(cell,1)=u_star;
                (*this).SolExact_(cell,2)=p_star;
                (*this).SolExact_(cell,3)=Y0_L;

            }
            else if(xcell<=x_0+sigma_R_star*time){

                (*this).SolExact_(cell,0)=rho_R_star;
                (*this).SolExact_(cell,1)=u_star;
                (*this).SolExact_(cell,2)=p_star;
                (*this).SolExact_(cell,3)=Y0_interm;

            }
            else if(xcell<=x_1+(u0_interm-c0_interm)*time){

                (*this).SolExact_(cell,0)=rho0_interm;
                (*this).SolExact_(cell,1)=u0_interm;
                (*this).SolExact_(cell,2)=p0_interm;
                (*this).SolExact_(cell,3)=Y0_interm;

            }
            else if(xcell<=x_1+sigma_L_starstar*time){
                //Inside the Left Rarefaction wave

                double X=xcell-x_1;
                Vector3d State_L_Fan=State_k_Fan(\
                        rho0_interm, u0_interm,\
                        p0_interm, c0_interm,\
                        (*this).SolTherm_, 1, \
                        X, time); 

                (*this).SolExact_(cell,0)=State_L_Fan(0);
                (*this).SolExact_(cell,1)=State_L_Fan(1);
                (*this).SolExact_(cell,2)=State_L_Fan(2);
                (*this).SolExact_(cell,3)=Y0_interm;

            }
            else if(xcell<=x_1+u_starstar*time){
                //Inside the Left star region

                (*this).SolExact_(cell,0)=rho_L_starstar;
                (*this).SolExact_(cell,1)=u_starstar;
                (*this).SolExact_(cell,2)=p_starstar;
                (*this).SolExact_(cell,3)=Y0_interm;

            }
            else if(xcell<=x_1+sigma_R_starstar*time){

                (*this).SolExact_(cell,0)=rho_R_starstar;
                (*this).SolExact_(cell,1)=u_starstar;
                (*this).SolExact_(cell,2)=p_starstar;
                (*this).SolExact_(cell,3)=Y0_R;

            }
            else{

                (*this).SolExact_(cell,0)=rho0_R;
                (*this).SolExact_(cell,1)=u0_R;
                (*this).SolExact_(cell,2)=p0_R;
                (*this).SolExact_(cell,3)=Y0_R;
            }

        }
    }

    else if (SolType=="Left Rarefaction Right Rarefaction"){
	//The theoretical sol is maid of 2 rarefaction waves

	//u star estimation
	u_star=u0_L-Euler_k_Resolvant(p_star, (*this).SolTherm_, \
		rho0_L, u0_L, p0_L, c0_L);

	//Left quantities
	rho_L_star=Rho_k_star_Resolvant(rho0_L, p0_L, p_star, (*this).SolTherm_);
	c_L_star=C_k_star_Resolvant( c0_L, p0_L, p_star, (*this).SolTherm_);

	sigma_L_star=u_star-c_L_star;
	
	//Right quantities
	rho_R_star=Rho_k_star_Resolvant(rho0_R, p0_R, p_star, (*this).SolTherm_);
	c_R_star=C_k_star_Resolvant( c0_R, p0_R, p_star, (*this).SolTherm_);

	sigma_R_star=u_star+c_R_star;

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);
	    cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;//

	    if(xcell<=x_0+(u0_L-c0_L)*time){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_L_star*time){
		//Inside the Left Rarefaction wave

		double X=xcell-x_0;
		Vector3d State_L_Fan=State_k_Fan(rho0_L, u0_L, p0_L, c0_L,\
			(*this).SolTherm_, 1, \
			X, time); 

		(*this).SolExact_(cell,0)=State_L_Fan(0);
		(*this).SolExact_(cell,1)=State_L_Fan(1);
		(*this).SolExact_(cell,2)=State_L_Fan(2);
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+u_star*time){
		//Inside the Left star region

		(*this).SolExact_(cell,0)=rho_L_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_R_star*time){

		(*this).SolExact_(cell,0)=rho_R_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else if(xcell<=x_0+(u0_R+c0_R)*time){
		//Inside the Right rarefaction wave 

		double X=xcell-x_0;
		Vector3d State_R_Fan=State_k_Fan(rho0_R, u0_R, p0_R, c0_R,\
			(*this).SolTherm_, 2, \
			X, time); 

		(*this).SolExact_(cell,0)=State_R_Fan(0);
		(*this).SolExact_(cell,1)=State_R_Fan(1);
		(*this).SolExact_(cell,2)=State_R_Fan(2);
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }

	}

    }
    else if (SolType=="Left Shock Right Shock"){
	//The theoretical sol is maid of 2 shock waves

	//u star estimation
	u_star=u0_L-Euler_k_Resolvant(p_star, (*this).SolTherm_, \
		rho0_L, u0_L, p0_L, c0_L);

	//Left quantities
	rho_L_star=Rho_k_star_Resolvant(rho0_L, p0_L, p_star, (*this).SolTherm_);
	
	sigma_L_star=Sigma_k_star_Resolvant(u0_L, c0_L, p0_L, p_star,\
		1, (*this).SolTherm_ );

	//Right quantities
	rho_R_star=Rho_k_star_Resolvant(rho0_R, p0_R, p_star, (*this).SolTherm_);

	sigma_R_star=Sigma_k_star_Resolvant(u0_R, c0_R, p0_R, p_star,\
		2, (*this).SolTherm_ );

	for (int cell=0; cell<NcellExt; cell++){


	    xcell=mesh.CellCoordsTab_(cell,1);
	    cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;//

	    if(xcell<=x_0+sigma_L_star*time){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+u_star*time){

		(*this).SolExact_(cell,0)=rho_L_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_R_star*time){

		(*this).SolExact_(cell,0)=rho_R_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }


	}

    }
    else if (SolType=="Left Rarefaction Right Shock"){
	//The theoretical sol is maid of 2 rarefaction waves

	//u star estimation
	u_star=u0_L-Euler_k_Resolvant(p_star, (*this).SolTherm_, \
		rho0_L, u0_L, p0_L, c0_L);

	//Left quantities
	rho_L_star=Rho_k_star_Resolvant(rho0_L, p0_L, p_star, (*this).SolTherm_);
	c_L_star=C_k_star_Resolvant( c0_L, p0_L, p_star, (*this).SolTherm_);

	sigma_L_star=u_star-c_L_star;
	
	//Right quantities
	rho_R_star=Rho_k_star_Resolvant(rho0_R, p0_R, p_star, (*this).SolTherm_);
	
	sigma_R_star=Sigma_k_star_Resolvant(u0_R, c0_R, p0_R, p_star,\
		2, (*this).SolTherm_ );

	for (int cell=0; cell<NcellExt; cell++){

	    xcell=mesh.CellCoordsTab_(cell,1);
	    cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;//

	    if(xcell<=x_0+(u0_L-c0_L)*time){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_L_star*time){
		//Inside the Left Rarefaction wave

		double X=xcell-x_0;
		Vector3d State_L_Fan=State_k_Fan(rho0_L, u0_L, p0_L, c0_L,\
			(*this).SolTherm_, 1, \
			X, time); 

		(*this).SolExact_(cell,0)=State_L_Fan(0);
		(*this).SolExact_(cell,1)=State_L_Fan(1);
		(*this).SolExact_(cell,2)=State_L_Fan(2);
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+u_star*time){
		//Inside the Left star region

		(*this).SolExact_(cell,0)=rho_L_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_R_star*time){

		(*this).SolExact_(cell,0)=rho_R_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }


	}

    }
    else if (SolType=="Left Shock Right Rarefaction"){
	//The theoretical sol is maid of 2 rarefaction waves

	//u star estimation
	u_star=u0_L-Euler_k_Resolvant(p_star, (*this).SolTherm_, \
		rho0_L, u0_L, p0_L, c0_L);

	//Left quantities
	rho_L_star=Rho_k_star_Resolvant(rho0_L, p0_L, p_star, (*this).SolTherm_);
	
	sigma_L_star=Sigma_k_star_Resolvant(u0_L, c0_L, p0_L, p_star,\
		1, (*this).SolTherm_ );

	//Right quantities
	rho_R_star=Rho_k_star_Resolvant(rho0_R, p0_R, p_star, (*this).SolTherm_);
	c_R_star=C_k_star_Resolvant( c0_R, p0_R, p_star, (*this).SolTherm_);
	
	sigma_R_star=u_star+c_R_star;

	for (int cell=0; cell<NcellExt; cell++){


	    xcell=mesh.CellCoordsTab_(cell,1);
	    cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;//

	    if(xcell<=x_0+sigma_L_star*time){

		(*this).SolExact_(cell,0)=rho0_L;
		(*this).SolExact_(cell,1)=u0_L;
		(*this).SolExact_(cell,2)=p0_L;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+u_star*time){
		//Inside the Left star region

		(*this).SolExact_(cell,0)=rho_L_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_L;

	    }
	    else if(xcell<=x_0+sigma_R_star*time){

		(*this).SolExact_(cell,0)=rho_R_star;
		(*this).SolExact_(cell,1)=u_star;
		(*this).SolExact_(cell,2)=p_star;
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else if(xcell<=x_0+(u0_R+c0_R)*time){
		//Inside the Right rarefaction wave 

		double X=xcell-x_0;
		Vector3d State_R_Fan=State_k_Fan(rho0_R, u0_R, p0_R, c0_R, \
			(*this).SolTherm_, 2, \
			X, time); 

		(*this).SolExact_(cell,0)=State_R_Fan(0);
		(*this).SolExact_(cell,1)=State_R_Fan(1);
		(*this).SolExact_(cell,2)=State_R_Fan(2);
		(*this).SolExact_(cell,3)=Y0_R;

	    }
	    else{

		(*this).SolExact_(cell,0)=rho0_R;
		(*this).SolExact_(cell,1)=u0_R;
		(*this).SolExact_(cell,2)=p0_R;
		(*this).SolExact_(cell,3)=Y0_R;
	    }


	}

    }

}


void Sol::Compute_L1_err(Mesh& mesh){

 double SpaceStep=mesh.Get_SpaceStep();
 int Ncells=mesh.Get_Ncells();

//Rho error:

 VectorXd Rho_disc(Ncells), Rho_exact(Ncells), Delta_Rho(Ncells), Rho_exact_abs(Ncells);
 //Extracting the variable without the fictitious cells
 Rho_disc=((*this).ConsVar_).block(1,0,Ncells,1);

 //Extracting the exact variable without the fictitious cells
 Rho_exact=((*this).SolExact_).block(1,0,Ncells,1);
 Rho_exact_abs=Rho_exact.array().abs();

 Delta_Rho=(Rho_disc-Rho_exact).array().abs();

 //norm L1 of rho exact
 (*this).normL1_rho_exact_=SpaceStep*Rho_exact_abs.sum();

 //error in norm L1 of the variable
 (*this).errL1_rho_=(SpaceStep*(Delta_Rho.sum()));

//U error:
 
 VectorXd U_disc(Ncells), U_exact(Ncells), Delta_U(Ncells), U_exact_abs(Ncells);
 //Extracting the variable without the fictitious cells
 U_disc=((*this).NConsVar_).block(1,1,Ncells,1);

 //Extracting the exact variable without the fictitious cells
 U_exact=((*this).SolExact_).block(1,1,Ncells,1);
 U_exact_abs=U_exact.array().abs();

// Delta_U=(U_disc).array().abs();
Delta_U=(U_disc-U_exact).array().abs(); 

 //norm L1 of u exact
 (*this).normL1_u_exact_=SpaceStep*U_exact_abs.sum();

 //error in norm L1 of the variable
 (*this).errL1_u_=(SpaceStep*(Delta_U.sum()));

//P error:

 VectorXd P_disc(Ncells), P_exact(Ncells), Delta_P(Ncells), P_exact_abs(Ncells);
 //Extracting the variable without the fictitious cells
 P_disc=((*this).NConsVar_).block(1,2,Ncells,1);

 //Extracting the exact variable without the fictitious cells
 P_exact=((*this).SolExact_).block(1,2,Ncells,1);
 P_exact_abs=P_exact.array().abs();

 //Delta_P=(P_disc).array().abs();
 Delta_P=(P_disc-P_exact).array().abs(); 

 //norm L1 of p exact
 (*this).normL1_p_exact_=SpaceStep*P_exact_abs.sum();

 //error in norm L1 of the variable
 (*this).errL1_p_=(SpaceStep*(Delta_P.sum()));

//Y error:
 
 VectorXd Y_disc(Ncells), Y_exact(Ncells), Delta_Y(Ncells), Y_exact_abs(Ncells);
 //Extracting the variable without the fictitious cells
 Y_disc=((*this).NConsVar_).block(1,3,Ncells,1);

 //Extracting the exact variable without the fictitious cells
 Y_exact=((*this).SolExact_).block(1,3,Ncells,1);
 Y_exact_abs=Y_exact.array().abs();

 Delta_Y=(Y_disc-Y_exact).array().abs();

 //norm L1 of Y exact
 (*this).normL1_Y_exact_=SpaceStep*Y_exact_abs.sum();

 //error in norm L1 of the variable
 (*this).errL1_Y_=(SpaceStep*(Delta_Y.sum()));
}

*/

//External functions:

VectorXd Riemann_Field(double val_L, double val_R, double Ncells, double NcellExt, double Length, double x_0){

    VectorXd Field(NcellExt);
    double SpaceStep=Length/Ncells;
    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

	if(i<=Riemann_index){

	    Field(i)=val_L;

	}
	else{

	    Field(i)=val_R;
	}

    }
    return Field;
}

VectorXd Double_Riemann_Field(\
        double val_L,\
        double val_R,\
        double val_interm,\
        double Ncells,\
        double NcellExt,\
        double Length,\
        double x_0, double x_1){

    VectorXd Field(NcellExt);
    double SpaceStep=Length/Ncells;
    int Riemann_index_0=floor(x_0/SpaceStep); 
    int Riemann_index_1=floor(x_1/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

	if(i<=Riemann_index_0){

	    Field(i)=val_L;

	}
    else if(i<=Riemann_index_1){

	    Field(i)=val_interm;

	}
	else{

	    Field(i)=val_R;
	}

    }
    return Field;
}

VectorXd SinusImpulse_Field(\
	double val_L, double M0,\
	int NPI_sin, int NPI_gauss,\
	double Ncells, double NcellExt,\
	double Length, double x_0){

    VectorXd Field(NcellExt);
    double SpaceStep=Length/Ncells;

    double sigma = NPI_gauss * SpaceStep;
    double k_c   = (TWO*Pi) / (NPI_sin*SpaceStep);

    for(int i=0;i<NcellExt;i++){

	    double x_i = (i*ONE-ONE_OVER_TWO) * SpaceStep;
	    Field(i)=val_L * (ONE + pow(M0, ONE) * sin(k_c * \
			(x_i-(x_0-k_c*Pi/TWO))) *\
		    (ONE/(sigma*sqrt(TWO*Pi)))*\
		    exp(-ONE_OVER_TWO * pow((x_i - x_0)/sigma, TWO)));
	}

    return Field;
}

double Euler_k_Resolvant(double& p, ThermoLaw& SolTherm, double& rho_k, double& u_k, double& p_k, double& c_k){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();

    if(ThermoLawType=="PG"){

	if(p>p_k){ //Shock function

	    double A_k=TWO/((Gamma+ONE)*rho_k);
	    double B_k=((Gamma-ONE)/(Gamma+ONE))*p_k;

	    return (p-p_k)*sqrt(A_k/(p+B_k));

	}

	else{  //Rarefaction function

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    return ((TWO*c_k)/(Gamma-ONE))*(pow(p/p_k,z)-ONE);

	}
    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	double alphaPiSG=(TWO*Gamma/(Gamma+1))*PiSG;

	if(p>p_k){ //Shock function

	    double A_k=TWO/((Gamma+ONE)*rho_k);
	    double B_k=((Gamma-ONE)/(Gamma+ONE))*p_k;

	    return (p-p_k)*sqrt(A_k/(p+B_k+alphaPiSG));

	}

	else{  //Rarefaction function

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    return ((TWO*c_k)/(Gamma-ONE))*(pow((p+PiSG)/(p_k+PiSG),z)-ONE);

	}
    }

    return ZERO;

}

double Rho_k_star_Resolvant(double rho0_k, double p0_k, double p_star, ThermoLaw& SolTherm ){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double beta=(Gamma-ONE)/(Gamma+ONE);

    if(ThermoLawType=="PG"){

	if(p_star>p0_k){ //Shock function

	    return rho0_k*((beta+(p_star/p0_k))/(ONE+(beta*p_star/p0_k)));

	}

	else{  //Rarefaction function

	    return rho0_k*pow(p_star/p0_k, ONE/Gamma);

	}
    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	double alphaPiSG=(TWO*Gamma/(Gamma+1))*PiSG;

	if(p_star>p0_k){ //Shock function

	    return rho0_k*((beta+(p_star/p0_k)+(alphaPiSG/p0_k))/(ONE+(beta*p_star/p0_k)+(alphaPiSG/p0_k)));

	}

	else{  //Rarefaction function

	    return rho0_k*pow((p_star+PiSG)/(p0_k+PiSG), ONE/Gamma);

	}
    }

    return ZERO;

}

double C_k_star_Resolvant(double c0_k, double p0_k, double p_star, ThermoLaw& SolTherm ){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double z=(Gamma-ONE)/(TWO*Gamma);

    if(ThermoLawType=="PG"){

	return c0_k*pow(p_star/p0_k,z);

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();

	return c0_k*pow((p_star+PiSG)/(p0_k+PiSG),z);

    }

    return ZERO;

}

double Sigma_k_star_Resolvant(double u0_k, double c0_k, \
	double p0_k, double p_star, \
	int sign, ThermoLaw& SolTherm ){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double A=(Gamma+ONE)/(TWO*Gamma);
    double B=(Gamma-ONE)/(TWO*Gamma);

    if(ThermoLawType=="PG"){

	return u0_k+pow(-ONE,sign)*c0_k*sqrt(A*(p_star/p0_k)+B);

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();

	return u0_k+pow(-ONE,sign)*c0_k*sqrt(A*((p_star+PiSG)/(p0_k+PiSG))+B);

    }

    return ZERO;



}

double Rho_pure_shock(double rho0_L, double u0_L, double p0_L,\
	double sigma, ThermoLaw& SolTherm){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double Beta=(Gamma-ONE)/(Gamma+ONE);

    double tau0_L=ONE/rho0_L;
    double tau_R_shock;

    if(ThermoLawType=="PG"){

	tau_R_shock=tau0_L*((((ONE+Beta)*tau0_L*p0_L)/pow(sigma-u0_L,2))+Beta);
	return ONE/tau_R_shock;

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	tau_R_shock=tau0_L*(((ONE+Beta)*tau0_L*(p0_L+PiSG))/pow(sigma-u0_L,TWO)+Beta);
	return ONE/tau_R_shock;

    }

    return ZERO;
}

double U_pure_shock(double rho0_L, double u0_L, double p0_L,\
	double rho_R_shock, ThermoLaw& SolTherm){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double Beta=(Gamma-ONE)/(Gamma+ONE);

    double tau_R_shock=ONE/rho_R_shock;
    double u_R_shock;
    double tau0_L=ONE/rho0_L;

    if(ThermoLawType=="PG"){

	u_R_shock=u0_L-(tau_R_shock-tau0_L)*\
		  sqrt(((ONE+Beta)*p0_L)/(tau_R_shock-Beta*tau0_L));
	return u_R_shock;

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	u_R_shock=u0_L-(tau_R_shock-tau0_L)*\
		  sqrt(((ONE+Beta)*(p0_L+PiSG))/(tau_R_shock-Beta*tau0_L));
	return u_R_shock;

    }

    return ZERO;
}

double P_pure_shock(double rho0_L, double p0_L,\
	double rho_R_shock,ThermoLaw& SolTherm){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    double Beta=(Gamma-ONE)/(Gamma+ONE);

    double tau_R_shock=ONE/rho_R_shock;
    double p_R_shock;
    double tau0_L=ONE/rho0_L;

    if(ThermoLawType=="PG"){

	p_R_shock=p0_L-((ONE+Beta)*(tau_R_shock-tau0_L)*p0_L)/(tau_R_shock-Beta*tau0_L);
	return p_R_shock;

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	p_R_shock=p0_L-((ONE+Beta)*(tau_R_shock-tau0_L)*(p0_L+PiSG))/(tau_R_shock-Beta*tau0_L);
	return p_R_shock;

    }

    return ZERO;
}

double U_pure_rarefaction(double rho0_L, double u0_L, double p0_L,\
	double rho_R_rarefaction, ThermoLaw& SolTherm){

    if(rho_R_rarefaction> rho0_L){

	cout<<"Inside U_pure_rarefaction: rho_R_rarefaction > rho0_L"<<endl;
	exit(EXIT_FAILURE);
	
    }

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();

    double u_R_rarefaction;

    if(ThermoLawType=="PG"){

	u_R_rarefaction=u0_L-\
			(TWO/(Gamma-ONE))*Sound_Speed_EOS(ONE,SolTherm, rho0_L, p0_L)*\
			(pow(rho_R_rarefaction/rho0_L,(Gamma-ONE)/TWO)-ONE);
	return u_R_rarefaction;

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	u_R_rarefaction=u0_L-\
			(TWO/(Gamma-ONE))*Sound_Speed_EOS(ONE,SolTherm, rho0_L, p0_L)*\
			(pow(rho_R_rarefaction/rho0_L,(Gamma-ONE)/TWO)-ONE);
	return u_R_rarefaction;

    }

    return ZERO;

}
double P_pure_rarefaction(double rho0_L, double p0_L,\
	double rho_R_rarefaction, ThermoLaw& SolTherm){

    if(rho_R_rarefaction> rho0_L){

	cout<<"Inside P_pure_rarefaction: rho_R_rarefaction > rho0_L"<<endl;
	exit(EXIT_FAILURE);
	
    }
    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();

    double p_R_rarefaction;

    if(ThermoLawType=="PG"){

	p_R_rarefaction=pow(rho_R_rarefaction/rho0_L,Gamma)*p0_L;
	return p_R_rarefaction;

    }
    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	p_R_rarefaction=pow(rho_R_rarefaction/rho0_L,Gamma)*(p0_L+PiSG)-PiSG;
	return p_R_rarefaction;

    }

    return ZERO;

}

double Gamma_Update_Impulse(double rho0_L, double u0_L, double p0_L, \
	double M0, double Cv, double T0, string ThermoType){

    double Gamma_Up;

    //Perfect Gas Thermodynamics
    if(ThermoType=="PG"){

	Gamma_Up = (pow(u0_L, TWO) * rho0_L) / (p0_L * pow(M0, TWO));
    }
    //Stiffened Gas Thermodynamics
    else if(ThermoType=="SG"){

	double Beta = pow( u0_L/M0, TWO) / (Cv*T0);
	Gamma_Up = (ONE + sqrt(ONE + FOUR * Beta)) / TWO;
    }
    
    return Gamma_Up;
}

double PiSG_Update_Impulse(double Gamma, double rho0_L, double u0_L, double p0_L,\
	double M0, string ThermoType){

    double PiSG_Up;

    //Perfect Gas Thermodynamics
    if(ThermoType=="PG"){

	PiSG_Up = PiSG;
    }
    //Stiffened Gas Thermodynamics
    else if(ThermoType=="SG"){

	PiSG_Up = pow( u0_L/M0, TWO) * (rho0_L/Gamma) - p0_L;
    }
    
    return PiSG_Up;
}

double Euler_Resolvant(double& p, string& SolType, ThermoLaw& SolTherm, double& rho_L, double& rho_R, double& u_L, double& u_R, double& p_L, double& p_R, double& c_L, double& c_R){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();

    if(ThermoLawType=="PG"){

	if(SolType=="Left Rarefaction Right Rarefaction"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    return (u_R-u_L)+((TWO*c_L)/(Gamma-ONE))*(pow(p/p_L,z)-ONE)+((TWO*c_R)/(Gamma-ONE))*\
		(pow(p/p_R,z)-ONE);


	}
	else if(SolType=="Left Shock Right Shock"){

	    double A_L=TWO/((Gamma+ONE)*rho_L);
	    double A_R=TWO/((Gamma+ONE)*rho_R);

	    double B_L=((Gamma-ONE)/(Gamma+ONE))*p_L;
	    double B_R=((Gamma-ONE)/(Gamma+ONE))*p_R;

	    return (u_R-u_L)+(p-p_L)*sqrt(A_L/(p+B_L))+(p-p_R)*sqrt(A_R/(p+B_R));

	}
	else if(SolType=="Left Rarefaction Right Shock"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    double A_R=TWO/((Gamma+ONE)*rho_R);

	    double B_R=((Gamma-ONE)/(Gamma+ONE))*p_R;

	    return (u_R-u_L)+((TWO*c_L)/(Gamma-ONE))*(pow(p/p_L,z)-ONE)+\
		(p-p_R)*sqrt(A_R/(p+B_R));

	}
	else if(SolType=="Left Shock Right Rarefaction"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    double A_L=TWO/((Gamma+ONE)*rho_L);

	    double B_L=((Gamma-ONE)/(Gamma+ONE))*p_L;

	    return (u_R-u_L)+(p-p_L)*sqrt(A_L/(p+B_L))+\
		((TWO*c_R)/(Gamma-ONE))*(pow(p/p_R,z)-ONE);

	}

    }

    //Stiffened Gas EOS
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();
	double alphaPiSG=(TWO*Gamma/(Gamma+1))*PiSG;

	if(SolType=="Left Rarefaction Right Rarefaction"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    return (u_R-u_L)+((TWO*c_L)/(Gamma-ONE))*(pow((p+PiSG)/(p_L+PiSG),z)-ONE)+\
		((TWO*c_R)/(Gamma-ONE))*\
		(pow((p+PiSG)/(p_R+PiSG),z)-ONE);


	}
	else if(SolType=="Left Shock Right Shock"){

	    double A_L=TWO/((Gamma+ONE)*rho_L);
	    double A_R=TWO/((Gamma+ONE)*rho_R);

	    double B_L=((Gamma-ONE)/(Gamma+ONE))*p_L;
	    double B_R=((Gamma-ONE)/(Gamma+ONE))*p_R;

	    return (u_R-u_L)+(p-p_L)*sqrt(A_L/(p+B_L+alphaPiSG))+\
		(p-p_R)*sqrt(A_R/(p+B_R+alphaPiSG));

	}
	else if(SolType=="Left Rarefaction Right Shock"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    double A_R=TWO/((Gamma+ONE)*rho_R);

	    double B_R=((Gamma-ONE)/(Gamma+ONE))*p_R;

	    return (u_R-u_L)+((TWO*c_L)/(Gamma-ONE))*(pow((p+PiSG)/(p_L+PiSG),z)-ONE)+\
		(p-p_R)*sqrt(A_R/(p+B_R+alphaPiSG));

	}
	else if(SolType=="Left Shock Right Rarefaction"){

	    double z=(Gamma-ONE)/(TWO*Gamma);
	    double A_L=TWO/((Gamma+ONE)*rho_L);

	    double B_L=((Gamma-ONE)/(Gamma+ONE))*p_L;

	    return (u_R-u_L)+(p-p_L)*sqrt(A_L/(p+B_L+alphaPiSG))+\
		((TWO*c_R)/(Gamma-ONE))*(pow((p+PiSG)/(p_R+PiSG),z)-ONE);

	}

    }

    return ZERO;

}

double P_Euler_Dichotomy(double& p_inf, double& p_sup, string& SolType, ThermoLaw& SolTherm, double& rho_L, double& rho_R, double& u_L, double& u_R, double& p_L, double& p_R, double& c_L, double& c_R, double& eps){

    //cout<<"[p_inf; p_sup]= "<<p_inf<<", "<<p_sup<<endl;

    double p_star; 
    double p_mid=(p_inf+p_sup)/TWO;
    double ResolInf, ResolSup, ResolMid;

    ResolInf=Euler_Resolvant(p_inf, SolType, SolTherm, rho_L, rho_R, \
     u_L, u_R, p_L, p_R, c_L, c_R);

    ResolSup=Euler_Resolvant(p_sup, SolType, SolTherm, rho_L, rho_R, \
     u_L, u_R, p_L, p_R, c_L, c_R);

    ResolMid=Euler_Resolvant(p_mid, SolType, SolTherm, rho_L, rho_R, \
	    u_L, u_R, p_L, p_R, c_L, c_R);

    //cout<<"ResolMid= "<<ResolMid<<endl;

    if(ResolInf*ResolSup>= ZERO){

	cout<<"Alert P_Euler_Dichotomy: Resol(p_inf) Resol(p_sup) of same sign !"<<endl;
	exit(EXIT_FAILURE);

    }


    if(fabs(ResolMid)<=eps){

	return p_mid;

    }

    else{

	if(ResolInf*ResolMid<ZERO){

	    p_star=P_Euler_Dichotomy(p_inf, p_mid, SolType, SolTherm, rho_L, rho_R, \
		    u_L, u_R, p_L, p_R, c_L, c_R, eps);

	}
	else {

	    p_star=P_Euler_Dichotomy(p_mid, p_sup, SolType, SolTherm, rho_L, rho_R, \
		    u_L, u_R, p_L, p_R, c_L, c_R, eps);
	}

    }

    return p_star;

}

Vector3d State_k_Fan(double& rho_k, double& u_k, double& p_k, double& c_k, ThermoLaw& SolTherm, int sign, double& x, double& time){

    string ThermoLawType=SolTherm.Get_ThermoLawType();
    double Gamma=SolTherm.Get_Gamma();
    Vector3d NConsVarFan; //In variable rho, u , p

    double A=TWO/(Gamma+ONE);
    double B=TWO/(Gamma-ONE);
    double z=(TWO*Gamma)/(Gamma-ONE);
    double beta=(Gamma-ONE)/(Gamma+ONE);


    if(ThermoLawType=="PG"){

	NConsVarFan(0)=rho_k*pow(A-pow(-ONE,sign)*(beta/c_k)*(u_k-x/time),B);
	NConsVarFan(1)=A*(-pow(-ONE,sign)*c_k+(ONE/B)*u_k+x/time);
	NConsVarFan(2)=p_k*pow(A-pow(-ONE,sign)*(beta/c_k)*(u_k-x/time),z);

    }

    //Stiffened Gas
    else if(ThermoLawType=="SG"){

	double PiSG=SolTherm.Get_PiSG();

	NConsVarFan(0)=rho_k*pow(A-pow(-ONE,sign)*(beta/c_k)*(u_k-x/time),B);
	NConsVarFan(1)=A*(-pow(-ONE,sign)*c_k+(ONE/B)*u_k+x/time);
	NConsVarFan(2)=(p_k+PiSG)*pow(A-pow(-ONE,sign)*(beta/c_k)*(u_k-x/time),z)-PiSG;

    }
    

    return NConsVarFan;

}

