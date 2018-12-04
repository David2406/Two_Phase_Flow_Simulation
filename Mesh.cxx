#include "Mesh.h"

//constructor:

Mesh::Mesh( double Length, int Ncells, int NGhostCells)
{
    Length_=Length;
    Ncells_=Ncells;
    int NcellExt=Ncells+NGhostCells;
    NcellExt_=NcellExt;
    int Nfaces=Ncells+1;
    Nfaces_=Nfaces;
    SpaceStep_=Length*ONE/Ncells;

    FaceIndex_.resize(Nfaces,6);
    FaceIndexDual_.resize(Nfaces,6);
    CellCoordsTab_.resize(NcellExt,4);
    CellIndex_.resize(NcellExt,3);

	for (int i=0; i< Nfaces; i++){

	    FaceIndex_(i,0)=i; //Face id is "i"
	    FaceIndex_(i,1)=i; //Left neighbor id of face 'i' is i
	    FaceIndex_(i,2)=i+1; //Right neighbor id of face 'i' is i+1
	    FaceIndex_(i,3)=i*SpaceStep_; //x_face=i*SpaceStep_
	    FaceIndex_(i,4)=ZERO; //y_face=0
	    FaceIndex_(i,5)=ZERO; //z_face=0

	    FaceIndexDual_(i,0)=i; //Face id is "i"
	    FaceIndexDual_(i,1)=i; //Left neighbor id of face 'i' is i
	    FaceIndexDual_(i,2)=i+1; //Right neighbor id of face 'i' is i+1
	    FaceIndexDual_(i,3)=(i + ONE_OVER_TWO)*SpaceStep_; //x_face=(i+1/2)*SpaceStep_
	    FaceIndexDual_(i,4)=ZERO; //y_face=0
	    FaceIndexDual_(i,5)=ZERO; //z_face=0

	}

	for (int i=0; i< NcellExt; i++){

	    //Treatment of the Left ghost cell
	    if (i==0){

		CellIndex_(i,0)=i;
		CellIndex_(i,1)=NoIndex;
		CellIndex_(i,2)=i;

	    }
	    //Treatment of the Right ghost cell
	    else if (i==NcellExt-1){

		CellIndex_(i,0)=i;
		CellIndex_(i,1)=i-1;
		CellIndex_(i,2)=NoIndex;

	    }
	    else{

		CellIndex_(i,0)=i;
		CellIndex_(i,1)=i-1;
		CellIndex_(i,2)=i;
	    }

	    CellCoordsTab_(i,0)=i; //Cell id is "i", Rq: first id is "0" it is
	                         //a ghost cell
	    CellCoordsTab_(i,1)=(i*ONE-ONE_OVER_TWO)*SpaceStep_; // x coords
	    CellCoordsTab_(i,2)=ZERO; // y coords
	    CellCoordsTab_(i,3)=ZERO; // z coords

	}

	LeftBCface_=0; // to treat left Boundary Condition more easily
	RightBCface_=Nfaces_-1;// to treat right Boundary Condition more easily

	/*cout<<"Inside Mesh: CellCoordsTab: "<<endl;
	cout<<endl;
	cout<<CellCoordsTab_<<endl;
	cout<<"Inside Mesh: SpaceStep: "<<SpaceStep_<<endl;*/
}

//methods:

double Mesh::Get_Length(){

    return Length_;
}
double Mesh::Get_SpaceStep(){
    return SpaceStep_;
}
int    Mesh::Get_Ncells(){
    return Ncells_;
}
int    Mesh::Get_NcellExt(){
    return NcellExt_;
}
int    Mesh::Get_Nfaces(){
    return Nfaces_;
}
int    Mesh::Get_LeftBCface(){
    return LeftBCface_;
}
int    Mesh::Get_RightBCface(){
    return RightBCface_;
}


int   Mesh::Get_Left_Neighbor(int face_id){

    return FaceIndex_(face_id,1);

}
int   Mesh::Get_Right_Neighbor(int face_id){

    return FaceIndex_(face_id,2);

}

int Mesh::Get_CoordsTab(int cell_id){

    return cell_id+1;

}

