#include "ThermoLaw.h"

//constructor:

ThermoLaw::ThermoLaw(\
                string ThermoLawType1, string ThermoLawType2,\
                double Gamma1, double Gamma2,\
                double PiSG1 , double PiSG2,\
                double kBAR1 , double kBAR2\
                )
{

    ThermoLawType1_= ThermoLawType1; 
    ThermoLawType2_= ThermoLawType2;
    Gamma1_        = Gamma1;
    Gamma2_        = Gamma2;
    PiSG1_         = PiSG1;
    PiSG2_         = PiSG2;
    kBAR1_         = kBAR1;
    kBAR2_         = kBAR2;

}

//constructor by copy:

ThermoLaw::ThermoLaw(ThermoLaw& Thermo){

    ThermoLawType1_ = Thermo.ThermoLawType1_;
    ThermoLawType2_ = Thermo.ThermoLawType2_;
    Gamma1_         = Thermo.Gamma1_;
    Gamma2_         = Thermo.Gamma2_;
    PiSG1_          = Thermo.PiSG1_;
    PiSG2_          = Thermo.PiSG2_;
    kBAR1_          = Thermo.kBAR1_;
    kBAR2_          = Thermo.kBAR2_;

}

//constructor by default:

ThermoLaw::ThermoLaw(){

    ThermoLawType1_="None";
    ThermoLawType2_="None";
    Gamma1_=ZERO;
    Gamma2_=ZERO;
    PiSG1_=ZERO;
    PiSG2_=ZERO;
    kBAR1_=ZERO;
    kBAR2_=ZERO;

}

//methods

string ThermoLaw::Get_ThermoLawType(){

    return ThermoLawType1_;

}

double ThermoLaw::Get_Gamma(){

    return Gamma1_;

}

double ThermoLaw::Get_PiSG(){

    return PiSG1_;

}

double ThermoLaw::Get_kBAR(){

    return kBAR1_;

}

/********************************/
//External methods:
/********************************/

double  Internal_Energy_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
        double rho, double p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }

    //Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

        return (ONE/(Gamma-ONE))*(p/rho);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

        return (ONE/(Gamma-ONE))*((p+Gamma*PiSG)/rho);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        return (ONE/(Gamma-ONE))*(p/rho);

    }

    return ZERO;

}

VectorXd  Internal_Energy_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;
    VectorXd InternalEnergyTab;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }

    //Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

        InternalEnergyTab=(ONE/(Gamma-ONE))*p.cwiseQuotient(rho);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

        InternalEnergyTab=(ONE/(Gamma-ONE))*(p+(Gamma*PiSG)*VectorXd::Constant(p.rows(),ONE)).cwiseQuotient(rho);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        InternalEnergyTab=(ONE/(Gamma-ONE))*p.cwiseQuotient(rho);

    }

    return InternalEnergyTab;
}

double Entropy_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
        double rho, double p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }

//Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

	return p*pow(rho,-Gamma);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

	return (p+PiSG)*pow(rho, -Gamma);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

	return kBAR;

    }

    return ZERO;
}

VectorXd  Entropy_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;
    VectorXd EntropyTab;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }

//Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

	VectorXd rhoGamma=(rho.array()).pow(-Gamma);
        EntropyTab=p.cwiseProduct(rhoGamma);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

	VectorXd rhoGamma=(rho.array()).pow(-Gamma);
        EntropyTab=(p+(Gamma*PiSG)*\
		VectorXd::Constant(p.rows(),ONE)).cwiseProduct(rhoGamma);
    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        EntropyTab=kBAR*VectorXd::Constant(p.rows(),ONE);

    }

    return EntropyTab;
}

double  Sound_Speed_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
        double rho, double p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    } 

    if (ThermoLawType=="PG"){ //Perfect Gas Law

	return sqrt(Gamma*(p/rho));

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

	return sqrt((Gamma*p+Gamma*PiSG)/rho);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

	return sqrt(Gamma*(p/rho));

    }

    return ZERO;
}

VectorXd  Sound_Speed_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& p){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;
    VectorXd SoundSpeedTab;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }

//Function:
    if (ThermoLawType=="PG"){ //Perfect Gas Law

        SoundSpeedTab=(Gamma*p.cwiseQuotient(rho)).cwiseSqrt();

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

        SoundSpeedTab=((Gamma*p+(Gamma*PiSG)*VectorXd::Constant(p.rows(),ONE)).cwiseQuotient(rho)).cwiseSqrt();

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        SoundSpeedTab=(Gamma*p.cwiseQuotient(rho)).cwiseSqrt();

    }

    return SoundSpeedTab;

}


double Pressure_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
        double rho, double e){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }        //Function: 

        if (ThermoLawType=="PG"){ //Perfect Gas Law

                return (Gamma-ONE)*(rho*e); 
        }
        else if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return (Gamma-ONE)*(rho*e)-Gamma*PiSG;

        }
        else if (ThermoLawType=="BAR"){ //Isentropic EOS

                return kBAR*pow(rho,Gamma); 

        }

        return ZERO;
}

VectorXd Pressure_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& rho, VectorXd& e){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;
    VectorXd PressureTab;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }

        //Function: 

    if (ThermoLawType=="PG"){ //Perfect Gas Law

	PressureTab=(Gamma-ONE)*(rho.cwiseProduct(e));

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

	PressureTab=(Gamma-ONE)*(rho.cwiseProduct(e))-(Gamma*PiSG)*VectorXd::Constant(e.rows(),ONE);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        return kBAR*(rho.array()).pow(Gamma); 

    }

    return PressureTab; 

}

double Density_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double e){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }

    //Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

        return (ONE/(Gamma-ONE))*(p/e);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

        return (ONE/(Gamma-ONE))*((p+Gamma*PiSG)/e);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        return pow(p/kBAR, ONE/Gamma);
    }

    return ZERO;
}

VectorXd Density_EOS_Tab(\
        int PhaseId,\
        ThermoLaw& therm, \
        VectorXd& p, VectorXd& e){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG, kBAR;
    VectorXd DensityTab;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;
        kBAR          = therm.kBAR1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;
        kBAR          = therm.kBAR2_;

    }

    //Function:

    if (ThermoLawType=="PG"){ //Perfect Gas Law

        DensityTab=(ONE/(Gamma-ONE))*p.cwiseQuotient(e);

    }
    else if (ThermoLawType=="SG"){ //Stiffened Gas Law

        DensityTab=(ONE/(Gamma-ONE))*(p+(Gamma*PiSG)*VectorXd::Constant(p.rows(),ONE)).cwiseQuotient(e);

    }
    else if (ThermoLawType=="BAR"){ //Isentropic EOS

        DensityTab=((ONE/kBAR)*p.array()).pow(ONE/Gamma);

    }

    return DensityTab;

}

/*************** CONVECTION SOURCE STABILITY: U-P RELAXATION ************/

double Gruneisen_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s){

    //Local variables:
    string ThermoLawType;
    double Gamma;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }        //Function: 

        if (ThermoLawType=="PG"){ //Perfect Gas Law

                return ONE/(Gamma-ONE); 
        }
        else if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return ONE/(Gamma-ONE);

        }

        return ZERO;

}

double Stiffness_EOS(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }        //Function: 

        if (ThermoLawType=="PG"){ //Perfect Gas Law

                return -ONE/(Gamma-ONE); 
        }
        else if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return -((p+PiSG)/p)*(ONE/(Gamma-ONE));

        }

        return ZERO;

}

double Sound_Speed_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }        //Function: 

    double alpha = ONE/(TWO*Gamma);
    double beta  = (Gamma - ONE)/(TWO*Gamma);

        if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return pow(Gamma,ONE_OVER_TWO)*pow(s,alpha)*pow(p+PiSG,beta); 
        }

        return ZERO;

}

double Rho_Sound_Speed_Squared_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }        //Function: 

        if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return Gamma*(p+PiSG); 
        }

        return ZERO;
}


double Rho_EOS_ps(\
        int PhaseId,\
        ThermoLaw& therm, \
	    double p, double s){

    //Local variables:
    string ThermoLawType;
    double Gamma, PiSG;

    if (PhaseId==1){

        ThermoLawType = therm.ThermoLawType1_;
        Gamma         = therm.Gamma1_;
        PiSG          = therm.PiSG1_;

    }
    else if (PhaseId==2){

        ThermoLawType = therm.ThermoLawType2_;
        Gamma         = therm.Gamma2_;
        PiSG          = therm.PiSG2_;

    }        //Function: 

    double alpha = ONE/(Gamma);

        if (ThermoLawType=="SG"){ //Stiffened Gas Law

                return pow(s,-alpha)*pow(p+PiSG,alpha); 
        }

        return ZERO;
}
 
