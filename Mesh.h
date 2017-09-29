#ifndef MESH
#define MESH
#include <iostream>
#include <Eigen/Dense>
#include "Param.h" //in order to have global constant like epsZero or Big

class Mesh
{
    private:
	double  Length_;
	double  SpaceStep_;     //Only one space step for now
	int     Ncells_;        //Number of physical cells
    int     NcellExt_;      //Number of cells including ghost cells
	int     Nfaces_;        //Number of faces
	int     LeftBCface_;    //Id of the left boundary face
	int     RightBCface_;   //Id of the right boundary face
    public:

	//Dynamic matrix, each row stores the cell i
	////and x_cell, y_cell, z_cell
	Eigen::MatrixXd CellCoordsTab_; 
                                
	//Dynamic matrix, each row stores the face id and the left cell neighbor row id 
	//(in CellCoordsTab_), the right cell neighbor row id (in CellCoordsTab_)
	//and the x_face coordinate
        Eigen::MatrixXd FaceIndex_;    
			
	//Dynamic matrix, each row stores 
	//then the left face id , the right face id , for the left 
	//ghost cell the left face id doesn't exist so it takes 
	//the 'NoIndex' value. Same thing for the right ghost cell
        Eigen::MatrixXd CellIndex_;     
                                  
				
	//constructor:
	Mesh(double Length, int Ncells, int NGhostcells); 

	//methods:
	double Get_Length();
	double Get_SpaceStep();
	int    Get_Ncells();
	int    Get_NcellExt();
	int    Get_Nfaces();
	int    Get_LeftBCface();
	int    Get_RightBCface();

	int   Get_Left_Neighbor(int face_id); //return the left neighbor id of the face:
	                                      //face_id
	int   Get_Right_Neighbor(int face_id);//return the right neighbor id of the face:
	                                      //face_id
	int   Get_CoordsTab(int cell_id);//returns the row number in CellCoordsTab_ for a cell whose cell id is 'cell_id'

};


#endif
