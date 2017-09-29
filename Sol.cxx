#include "Sol.h"

//constructor:

Sol::Sol(Mesh& M,\
		ThermoLaw& Therm,\
        Vector7d& InitL, Vector7d& InitR,\
		string SolType,\
		string LeftBCType, string RightBCType,\
        double x0,\
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

    //Thermodynamical parameters:
//    string ThermoLawType1 =Therm.ThermoLawType1_;
//    double Gamma1         =Therm.Gamma1_;
//    double PiSG1          =Therm.PiSG1_;
//
//    string ThermoLawType2 =Therm.ThermoLawType2_;
//    double Gamma2         =Therm.Gamma2_;
//    double PiSG2          =Therm.PiSG2_;

    //Copying the initial Rieman state (theoretical solution)
    InitL_ = InitL;
    InitR_ = InitR;

//    c0_L_=Sound_Speed_EOS(Therm1, rho0_L, p0_L);
//    c0_R_=Sound_Speed_EOS(Therm1, rho0_R, p0_R);
//    c0_interm_=Sound_Speed_EOS(ONE, Therm, rho0_interm, p0_interm);
//
//    //Initializing the L1 errors 
//    errL1_rho_=Big;
//    errL1_u_=Big;
//    errL1_p_=Big;
//    errL1_Y_=Big;
//
//    //Declaring Initial State Vectors
//    VectorXd Rho0, U0, P0, Y0;
//
//    if(SolType=="Sinus Impulse"){
//
//	//1D Sinus Impulse field: Rho0
//	Rho0=SinusImpulse_Field(rho0_L, M0,\
//		NPI_sin, NPI_gauss,\
//		Ncells, NcellExt,\
//		Length, x_0);
//
//	//1D problem Riemann field: U0
//	U0=Riemann_Field(u0_L, u0_R, Ncells, NcellExt, Length, x_0);
//
//	//1D problem Riemann field: P0
//	P0=Riemann_Field(p0_L, p0_R, Ncells, NcellExt, Length, x_0);
//
//	//1D problem Riemann field: Y0
//	Y0=Riemann_Field(Y0_L, Y0_R, Ncells, NcellExt, Length, x_0);
//    }
//    else if(SolType=="Double Riemann Problem"){
//
//        //1D problem Riemann field: Rho0
//        Rho0=Double_Riemann_Field(\
//                rho0_L, rho0_R, rho0_interm,\
//                Ncells, NcellExt, Length,\
//                x_0, x_1);
//        //1D problem Riemann field: U0
//        U0=Double_Riemann_Field(\
//                u0_L, u0_R, u0_interm,\
//                Ncells, NcellExt, Length,\
//                x_0, x_1);
//        //1D problem Riemann field: P0
//        P0=Double_Riemann_Field(\
//                p0_L, p0_R, p0_interm,\
//                Ncells, NcellExt, Length,\
//                x_0, x_1);
//        //1D problem Riemann field: Y0
//        Y0=Double_Riemann_Field(\
//                Y0_L, Y0_R, Y0_interm,\
//                Ncells, NcellExt, Length,\
//                x_0, x_1);
//    }
//    else{
//	//1D problem Riemann field: Rho0
//	Rho0=Riemann_Field(rho0_L, rho0_R, Ncells, NcellExt, Length, x_0);
//
//	//1D problem Riemann field: U0
//	U0=Riemann_Field(u0_L, u0_R, Ncells, NcellExt, Length, x_0);
//
//	//1D problem Riemann field: P0
//	P0=Riemann_Field(p0_L, p0_R, Ncells, NcellExt, Length, x_0);
//
//	//1D problem Riemann field: Y0
//	Y0=Riemann_Field(Y0_L, Y0_R, Ncells, NcellExt, Length, x_0);
//    }
//
//    VectorXd Qdm0=Rho0.cwiseProduct(U0);
//
//    VectorXd e_intern0=Internal_Energy_EOS_Tab(Therm,  Rho0, P0);
//
//    VectorXd E0=ONE_OVER_TWO*(Qdm0.cwiseProduct(Qdm0)).cwiseQuotient(Rho0)+\
//		Rho0.cwiseProduct(e_intern0);
//
//    ConsVar_.resize(NcellExt, 4);
//    //The first column of ConsVar is filled with Rho0
//    ConsVar_.block(0,0,NcellExt,1)= Rho0;
//    //The second column of ConsVar is filled with Qdm0
//    ConsVar_.block(0,1,NcellExt,1)= Qdm0;
//    //The third column of ConsVar is filled with E0
//    ConsVar_.block(0,2,NcellExt,1)= E0;
//    //The fourth column of ConsVar is filled with RhoY0
//    ConsVar_.block(0,3,NcellExt,1)=Rho0.cwiseProduct(Y0);  
//
//    ConsVarFlux_=MatrixXd::Zero(Nfaces,4);// ConsVarFlux is initialized to 0.
//
//    NConsVar_.resize(NcellExt,4);
//    //The first column of NConsVar is filled with Rho0
//    NConsVar_.block(0,0,NcellExt,1)= Rho0;
//    //The second column of NConsVar is filled with U0
//    NConsVar_.block(0,1,NcellExt,1)= U0;
//    //The third column of NConsVar is filled with p0
//    NConsVar_.block(0,2,NcellExt,1)= P0;  
//    //The fourth column of NConsVar is filled with Y0
//    NConsVar_.block(0,3,NcellExt,1)= Y0;  
//
//    SolExact_.resize(NcellExt,4);
//    //Rho, U, P, Y
//    SolExact_.block(0,0,NcellExt,4)=NConsVar_.block(0,0,NcellExt,4);
//
//    SoundSpeed_=Sound_Speed_EOS_Tab(ONE,Therm,Rho0,P0); 
//
//    Mach_.resize(NcellExt);
//    Mach_=((NConsVar_.col(1)).cwiseQuotient(SoundSpeed_));
//    MachMax_=Mach_.array().abs().maxCoeff();
//    MachMin_=Mach_.array().abs().minCoeff();

}

Sol::Sol( Sol& solution){

    LeftBCType_   =solution.LeftBCType_;
    RightBCType_  =solution.RightBCType_;

    SolTherm_=solution.SolTherm_;

    ConsVar_      =solution.ConsVar_;
    NConsVar_     =solution.NConsVar_;
    SoundSpeed_   =solution.SoundSpeed_;
    ConsVarFlux_  =solution.ConsVarFlux_;

}

Sol::Sol(){

    LeftBCType_   ="None";
    RightBCType_  ="None";

}

//methods:

void Sol::ConsVarInit(Vector7d& InitR, Vector7d& InitL,\
        int NcellExt, int Ncells, double Length){

    double SpaceStep=Length/Ncells;
    double x_0 = (*this).x_0_;

    int Riemann_index=floor(x_0/SpaceStep); 

    for(int i=0;i<NcellExt;i++){

        if(i<=Riemann_index){

            ((*this).ConsVar_).block(i,0,1,Ncols) = InitL; 

        }
        else{

            ((*this).ConsVar_).block(i,0,1,Ncols) = InitR; 
        }

    }
}

//Update methods:

void Sol::SoundSpeed_Update(double eps0sq){

    VectorXd Rhonew=ConsVar_.col(0);
    VectorXd Pnew=NConsVar_.col(2);

    SoundSpeed_=Sound_Speed_EOS_Tab(eps0sq, (*this).SolTherm_, Rhonew,Pnew);                         
}

void Sol::Mach_Update(){

    VectorXd Unew=NConsVar_.col(1);
    Mach_=(Unew.cwiseQuotient(SoundSpeed_));
    MachMax_=Mach_.array().abs().maxCoeff();

    MachMin_=Mach_.array().abs().minCoeff();

}
                                            
void Sol::eps0sq_Update(double Mach){

    double alphaMach=(*this).alphaMach_;
    (*this).eps0sq_=max(pow(Mach_ref, alphaMach), min(pow(Mach, alphaMach),ONE));

}
      
void Sol::eps0sqTab_Update(Mesh& mesh){

    int Nfaces=mesh.Get_Nfaces();
    int L,R;
    double rho_L, rho_R;
    double u_L, u_R, u_star, c_L, c_R, c_max;
    double p_L, p_R;
    double a_Euler;
    double MachLoc;
    double alphaMach=(*this).alphaMach_;
    double Length=mesh.Get_Length();
    double dx=mesh.Get_SpaceStep();

    //eps pour le dectecteur de choc
    double epsDetector = dx/Length;

    double eps0sqTab_max;

    //eps0sqTab is local
    if((*this).eps0_local_==true){

        //Loop on all faces
        for (int face_id=0; face_id < Nfaces; face_id++){

            L=mesh.FaceIndex_(face_id,1);
            R=mesh.FaceIndex_(face_id,2);

            rho_L=(*this).NConsVar_(L,0);
            rho_R=(*this).NConsVar_(R,0);

            u_L=(*this).NConsVar_(L,1);
            u_R=(*this).NConsVar_(R,1);

            p_L=(*this).NConsVar_(L,2);
            p_R=(*this).NConsVar_(R,2);

            c_L=(*this).SoundSpeed_(L);
            c_R=(*this).SoundSpeed_(R);

            a_Euler = (ONE+ONE/TEN)*max(rho_L*c_L, rho_R*c_R);
            u_star  = (u_R + u_L)/TWO -(p_R-p_L)/(TWO*a_Euler); 
            c_max   = max(c_L,c_R);

            MachLoc = abs(u_star)/c_max;
            MachLoc = max(fabs(u_L)/c_L,fabs(u_R)/c_R);

            //eps0sqTab is only based on the Mach number
            if((*this).eps0_type_=="MACH"){

                (*this).eps0sqTab_(face_id)=max(pow(Mach_ref, alphaMach),\
                        min(pow(MachLoc, alphaMach),ONE));
                (*this).ShockFrontTab_(face_id)=ZERO;
                (*this).DeltaUFrontTab_(face_id)=ZERO;
                (*this).MachShockTab_(face_id)=ZERO;
            }
            //eps0sqTab is based on the Mach number and
            //the shock front speed
            if((*this).eps0_type_=="MACH_SHOCK"){

                double shock_speed;

                //Assessing the shock speed front
                //if ((fabs(rho_L-rho_R)/max(rho_R,rho_L))> epsZero){

                //    shock_speed=fabs((rho_R*u_R-rho_L*u_L)/(rho_R-rho_L));
                //}
                if ((fabs(u_L-u_R)/max(fabs(u_R),fabs(u_L)))> epsDetector){

                    shock_speed=fabs((p_R-p_L)/(((rho_R+rho_L)/TWO)*(u_R-u_L)));
                }
                else{ shock_speed=ZERO;}

                (*this).eps0sqTab_(face_id)=max(pow(Mach_ref, alphaMach),\
                        min(pow( max(MachLoc, shock_speed/c_max), alphaMach),ONE));
                (*this).ShockFrontTab_(face_id)=shock_speed;
                (*this).DeltaUFrontTab_(face_id)=fabs(u_L-u_R)/max(fabs(u_R),fabs(u_L));
                (*this).MachShockTab_(face_id)=shock_speed/c_max;
            }

        }

        //FIXME Include a shock detector in eps0sq
        eps0sqTab_max = (*this).eps0sqTab_.array().maxCoeff();
        int nrows=((*this).eps0sqTab_).rows();
        (*this).eps0sqTab_=(eps0sqTab_max)*VectorXd::Constant(nrows,ONE);
        (*this).eps0sq_ = eps0sqTab_max;

    }
    //eps0sqTab is global = eps0sq*Constant_vector
    else{

        int nrows=((*this).eps0sqTab_).rows();
        (*this).eps0sqTab_=((*this).eps0sq_)*VectorXd::Constant(nrows,ONE);
    }

}

void Sol::GammaSplit_Update(double Gamma, double eps0sq){

        (*this).GammaSplit_= eps0sq*(Gamma-ONE)+ONE;

}


void Sol::ConsVar_Update( Mesh& mesh, double TimeStep, \
        string SchemeType1, string SchemeType2,
        bool SplittingStep1, bool SplittingStep2
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


void Sol::NConsVar_Update(){

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
	    /*cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;*/

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
	    /*cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;*/

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
	    /*cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;*/

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
	    /*cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;*/

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
	    /*cout<<"xcell= "<<xcell<<endl;
	    cout<<"u0= "<<u0_L<<endl;
	    cout<<"x_0+t*u0= "<<x_0+time*u0_L<<endl;*/

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

//In case of an isolated contact discontinuity 
//check if the computed velocity is still constant or not

bool Sol::Check_U_Constant(){

    int NcellExt = ((*this).SoundSpeed_).size();

    double u0 = (*this).u0_L_;
    bool is_egal_U = true;

    for (int i=0; i < NcellExt; i++){

	double u = ((*this).NConsVar_)(i,1);
	if (fabs(u-u0)> eps_check_cst * u0){

	    is_egal_U = false;
	    break;
	}

    }

    return is_egal_U;
}

//In case of an isolated contact discontinuity 
//check if the computed pressure is still constant or not

bool Sol::Check_P_Constant(){

    int NcellExt = ((*this).SoundSpeed_).size();

    double p0 = (*this).p0_L_;
    bool is_egal_P = true;

    for (int i=0; i < NcellExt; i++){

	double p = ((*this).NConsVar_)(i,2);
	if (fabs(p-p0)> eps_check_cst * p0){

	    is_egal_P = false;
	    break;
	}

    }

    return is_egal_P;
}
