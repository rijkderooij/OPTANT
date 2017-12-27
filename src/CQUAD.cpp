#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 
#include "fei_base.hpp"
#include <cmath>

// Constructors:
CQUAD::CQUAD(): Element(Shell){}
CQUAD::~CQUAD(){}

// Build stiffness matrix of the element
int CQUAD::BuildElementMatrix(){     
    // First extract the element properties and interpolation matrix
	Property** prop = Element::Properties();
	int nProp = Element::NumProperties();
	
	// Interpolate the element properties to the center of the triangles, by
	// creating the interpolation matrix P: (actually we store P**T)
	double InterpP[4*nProp];
	if(InterpolateProperties(nProp,InterpP)==1)return 1; 
	
	// Initialize vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle

	// Allocate storage
    int size = 24;
	double storage[576];  // 567, same as tria;
	double* BLm = &storage[0],  *B4m = &storage[27], *B5m = &storage[54], *B6m = &storage[81];
	double* BLb = &storage[108],*B4b = &storage[135],*B5b = &storage[162],*B6b = &storage[189];
	double* Be1 = &storage[216],*Be2 = &storage[243],*Be3 = &storage[270],*Be4 = &storage[297];
	double* Kmm = &storage[324],*Kbb = &storage[405],*Kbm = &storage[486];
	
	// Initialize some matrices and vectors used for each of the four tia elements
	double LambdaT[9]={0,0,0,0,0,0,0,0,0};  // contains ex,ey,ez, i.e. LambdaT == Lambda**T
	double Nodes2d[6]={0,0,0,0,0,0}; // Location of the the three nodes in element coordinates;
	double d1[3] = {0,0,0}, d2[3]={0,0,0}, n[3]={0,0,0}; // Relative nodal position vector and normal vector
	
	
	// Loop through four triangles
	for(std::size_t tria=0; tria<4; tria++){
		// Initialize storage to zero
		for(std::size_t i=0;i<567;i++)storage[i] = 0;
		
		// Get nodes corresponding to this triangle
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute the geometric characteristics of this element
		
		// Compute normal direction of element, which is e_z
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,&LambdaT[6]); // ez=  d1 x d2
		double lN = sqrt(MatrixOperations::Dot(&LambdaT[6],&LambdaT[6],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[6],&LambdaT[6], 3); // ez = 1/norm(ez);
		
		// Compute e_x = emat-(emat'*ez)*ez;
		MatrixOperations::Scale(MatrixOperations::Dot(&LambdaT[6],MaterialDirection_,3),
								0, &LambdaT[6],&LambdaT[0], 3); // ex = dot(emat,ez)*ez;
		MatrixOperations::Add(1,-1,MaterialDirection_, &LambdaT[0],&LambdaT[0],3); //ex = emat-ex;
		lN = sqrt(MatrixOperations::Dot(&LambdaT[0],&LambdaT[0],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[0],&LambdaT[0], 3); // ex = 1/norm(ex);
		
		// Compute ey = cross(ez,ex);
		MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
		
		// Compute nodes2D
		Nodes2d[0] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[1] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[2] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[3] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[4] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[5] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[3],3);
		
		
		// Store coordinate differences in Be1
		Be1[0] = Nodes2d[0]-Nodes2d[2]; Be1[1] = Nodes2d[1]-Nodes2d[3]; // Be1[0] = x12;  Be1[1] = y12 
		Be1[2] = -Be1[0]; 				Be1[3] = -Be1[1]; 				// Be1[2] = x21;  Be1[3] = y21 
		Be1[4] = Nodes2d[2]-Nodes2d[4]; Be1[5] = Nodes2d[3]-Nodes2d[5]; // Be1[4] = x23;  Be1[5] = y23 
		Be1[6] = -Be1[4]; 				Be1[7] = -Be1[5]; 				// Be1[6] = x32;  Be1[7] = y32 
		Be1[8] = Nodes2d[4]-Nodes2d[0]; Be1[9] = Nodes2d[5]-Nodes2d[1]; // Be1[8] = x31;  Be1[9] = y31 
		Be1[10]= -Be1[8]; 				Be1[11]= -Be1[9]; 				// Be1[10]= x13;  Be1[11]= y13
		Be1[12]= (Be1[3]*Be1[10]-Be1[2]*Be1[11])/2; 					// Be1[12]= area
		Be1[13]= Be1[2]*Be1[2]+Be1[3]*Be1[3];							// Be1[13]= L21*L21
		Be1[14]= Be1[6]*Be1[6]+Be1[7]*Be1[7];							// Be1[14]= L32*L32
		Be1[15]= Be1[10]*Be1[10]+Be1[11]*Be1[11];						// Be1[15]= L13*L13
			
		// Calculate membranal B matrices
		CHK_ERR(CTRIA::FelippaTriMembrane(storage,true));
		
		// Calculate bending B matrices
		CTRIA::FelippaTriBending(storage);
		
		// Compute material Am,Bm,Dm matrices in conventional storage
		double *Am = &Be2[0], *Bm = &Be2[9], *Dm = &Be2[18];
		if(MaterialMatrices(Am,Bm,Dm,prop,nProp,InterpP,tria)!=0) return 1; 
		
		// Compute Kmm = BLm**T*Am*BLm*A + (B4m**T*Am*B4m+B5m**T*Am*B5m+B6m**T*Am*B6m)*A/3 * 0.5;
		// NOTE: The 0.5 term is added, to ensure the contribution of each triangle is halved.
		MatrixOperations::dsymm_bact(Be1[12]*0.5,BLm,Am,BLm,Kmm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B4m,Am,B4m,Kmm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B5m,Am,B5m,Kmm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B6m,Am,B6m,Kmm,3,9);
		
		// Compute Kbb = BLb**T*Dm*BLb*A + (B4b**T*Dm*B4b+B5b**T*Dm*B5b+B6b**T*Dm*B6b)*A/3 * 0.5;
		// NOTE: The 0.5 term is added, to ensure the contribution of each triangle is halved.
		MatrixOperations::dsymm_bact(Be1[12]*0.5,BLb,Dm,BLb,Kbb,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B4b,Dm,B4b,Kbb,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B5b,Dm,B5b,Kbb,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B6b,Dm,B6b,Kbb,3,9);
		
		// Compute Kbm = BLb**T*Bm*BLm*A + (B4b**T*Bm*B4m+B5b**T*Bm*B5m+B6b**T*Bm*B6m)*A/3 * 0.5;
		// NOTE: The 0.5 term is added, to ensure the contribution of each triangle is halved.
		MatrixOperations::dsymm_bact(Be1[12]*0.5,BLb,Bm,BLm,Kbm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B4b,Bm,B4m,Kbm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B5b,Bm,B5m,Kbm,3,9);
		MatrixOperations::dsymm_bact(Be1[12]/3*0.5,B6b,Bm,B6m,Kbm,3,9);
		
		// Compute Triangle ElementMatrix Kt, as Kt(mdof,mdof) = Kmm, Kt(bdof,bdof)=Kbb, Kt(bdof,mdof)=Kbm
		// and store immediately in conventional storage. Code is generated by MATLAB.
		// Store K at location of all B matrices:
		double* Kt = &storage[0];
		Kt[0]=Kmm[0];	Kt[1]=Kmm[1];	Kt[18]=Kmm[1];	Kt[2]=Kbm[0];	Kt[36]=Kbm[0];	Kt[3]=Kbm[1];	Kt[54]=Kbm[1];	Kt[4]=Kbm[2];	Kt[72]=Kbm[2];	Kt[5]=Kmm[2];	Kt[90]=Kmm[2];	Kt[6]=Kmm[3];	Kt[108]=Kmm[3];	Kt[7]=Kmm[4];	Kt[126]=Kmm[4];	Kt[8]=Kbm[3];	Kt[144]=Kbm[3];	Kt[9]=Kbm[4];	Kt[162]=Kbm[4];	Kt[10]=Kbm[5];	Kt[180]=Kbm[5];	Kt[11]=Kmm[5];	Kt[198]=Kmm[5];	Kt[12]=Kmm[6];	Kt[216]=Kmm[6];	Kt[13]=Kmm[7];	Kt[234]=Kmm[7];	Kt[14]=Kbm[6];	Kt[252]=Kbm[6];	Kt[15]=Kbm[7];	Kt[270]=Kbm[7];	Kt[16]=Kbm[8];	Kt[288]=Kbm[8];	Kt[17]=Kmm[8];	Kt[306]=Kmm[8];	Kt[19]=Kmm[10];	Kt[20]=Kbm[9];	Kt[37]=Kbm[9];	
		Kt[21]=Kbm[10];	Kt[55]=Kbm[10];	Kt[22]=Kbm[11];	Kt[73]=Kbm[11];	Kt[23]=Kmm[11];	Kt[91]=Kmm[11];	Kt[24]=Kmm[12];	Kt[109]=Kmm[12];	Kt[25]=Kmm[13];	Kt[127]=Kmm[13];	Kt[26]=Kbm[12];	Kt[145]=Kbm[12];	Kt[27]=Kbm[13];	Kt[163]=Kbm[13];	Kt[28]=Kbm[14];	Kt[181]=Kbm[14];	Kt[29]=Kmm[14];	Kt[199]=Kmm[14];	Kt[30]=Kmm[15];	Kt[217]=Kmm[15];	Kt[31]=Kmm[16];	Kt[235]=Kmm[16];	Kt[32]=Kbm[15];	Kt[253]=Kbm[15];	Kt[33]=Kbm[16];	Kt[271]=Kbm[16];	Kt[34]=Kbm[17];	Kt[289]=Kbm[17];	Kt[35]=Kmm[17];	Kt[307]=Kmm[17];	Kt[38]=Kbb[0];	Kt[39]=Kbb[1];	Kt[56]=Kbb[1];	Kt[40]=Kbb[2];	Kt[74]=Kbb[2];	Kt[41]=Kbm[18];	Kt[92]=Kbm[18];	Kt[42]=Kbm[27];	Kt[110]=Kbm[27];	
		Kt[43]=Kbm[36];	Kt[128]=Kbm[36];	Kt[44]=Kbb[3];	Kt[146]=Kbb[3];	Kt[45]=Kbb[4];	Kt[164]=Kbb[4];	Kt[46]=Kbb[5];	Kt[182]=Kbb[5];	Kt[47]=Kbm[45];	Kt[200]=Kbm[45];	Kt[48]=Kbm[54];	Kt[218]=Kbm[54];	Kt[49]=Kbm[63];	Kt[236]=Kbm[63];	Kt[50]=Kbb[6];	Kt[254]=Kbb[6];	Kt[51]=Kbb[7];	Kt[272]=Kbb[7];	Kt[52]=Kbb[8];	Kt[290]=Kbb[8];	Kt[53]=Kbm[72];	Kt[308]=Kbm[72];	Kt[57]=Kbb[10];	Kt[58]=Kbb[11];	Kt[75]=Kbb[11];	Kt[59]=Kbm[19];	Kt[93]=Kbm[19];	Kt[60]=Kbm[28];	Kt[111]=Kbm[28];	Kt[61]=Kbm[37];	Kt[129]=Kbm[37];	Kt[62]=Kbb[12];	Kt[147]=Kbb[12];	Kt[63]=Kbb[13];	Kt[165]=Kbb[13];	Kt[64]=Kbb[14];	Kt[183]=Kbb[14];	Kt[65]=Kbm[46];	Kt[201]=Kbm[46];	
		Kt[66]=Kbm[55];	Kt[219]=Kbm[55];	Kt[67]=Kbm[64];	Kt[237]=Kbm[64];	Kt[68]=Kbb[15];	Kt[255]=Kbb[15];	Kt[69]=Kbb[16];	Kt[273]=Kbb[16];	Kt[70]=Kbb[17];	Kt[291]=Kbb[17];	Kt[71]=Kbm[73];	Kt[309]=Kbm[73];	Kt[76]=Kbb[20];	Kt[77]=Kbm[20];	Kt[94]=Kbm[20];	Kt[78]=Kbm[29];	Kt[112]=Kbm[29];	Kt[79]=Kbm[38];	Kt[130]=Kbm[38];	Kt[80]=Kbb[21];	Kt[148]=Kbb[21];	Kt[81]=Kbb[22];	Kt[166]=Kbb[22];	Kt[82]=Kbb[23];	Kt[184]=Kbb[23];	Kt[83]=Kbm[47];	Kt[202]=Kbm[47];	Kt[84]=Kbm[56];	Kt[220]=Kbm[56];	Kt[85]=Kbm[65];	Kt[238]=Kbm[65];	Kt[86]=Kbb[24];	Kt[256]=Kbb[24];	Kt[87]=Kbb[25];	Kt[274]=Kbb[25];	Kt[88]=Kbb[26];	Kt[292]=Kbb[26];	Kt[89]=Kbm[74];	Kt[310]=Kbm[74];	
		Kt[95]=Kmm[20];	Kt[96]=Kmm[21];	Kt[113]=Kmm[21];	Kt[97]=Kmm[22];	Kt[131]=Kmm[22];	Kt[98]=Kbm[21];	Kt[149]=Kbm[21];	Kt[99]=Kbm[22];	Kt[167]=Kbm[22];	Kt[100]=Kbm[23];	Kt[185]=Kbm[23];	Kt[101]=Kmm[23];	Kt[203]=Kmm[23];	Kt[102]=Kmm[24];	Kt[221]=Kmm[24];	Kt[103]=Kmm[25];	Kt[239]=Kmm[25];	Kt[104]=Kbm[24];	Kt[257]=Kbm[24];	Kt[105]=Kbm[25];	Kt[275]=Kbm[25];	Kt[106]=Kbm[26];	Kt[293]=Kbm[26];	Kt[107]=Kmm[26];	Kt[311]=Kmm[26];	Kt[114]=Kmm[30];	Kt[115]=Kmm[31];	Kt[132]=Kmm[31];	Kt[116]=Kbm[30];	Kt[150]=Kbm[30];	Kt[117]=Kbm[31];	Kt[168]=Kbm[31];	Kt[118]=Kbm[32];	Kt[186]=Kbm[32];	Kt[119]=Kmm[32];	Kt[204]=Kmm[32];	Kt[120]=Kmm[33];	Kt[222]=Kmm[33];	
		Kt[121]=Kmm[34];	Kt[240]=Kmm[34];	Kt[122]=Kbm[33];	Kt[258]=Kbm[33];	Kt[123]=Kbm[34];	Kt[276]=Kbm[34];	Kt[124]=Kbm[35];	Kt[294]=Kbm[35];	Kt[125]=Kmm[35];	Kt[312]=Kmm[35];	Kt[133]=Kmm[40];	Kt[134]=Kbm[39];	Kt[151]=Kbm[39];	Kt[135]=Kbm[40];	Kt[169]=Kbm[40];	Kt[136]=Kbm[41];	Kt[187]=Kbm[41];	Kt[137]=Kmm[41];	Kt[205]=Kmm[41];	Kt[138]=Kmm[42];	Kt[223]=Kmm[42];	Kt[139]=Kmm[43];	Kt[241]=Kmm[43];	Kt[140]=Kbm[42];	Kt[259]=Kbm[42];	Kt[141]=Kbm[43];	Kt[277]=Kbm[43];	Kt[142]=Kbm[44];	Kt[295]=Kbm[44];	Kt[143]=Kmm[44];	Kt[313]=Kmm[44];	Kt[152]=Kbb[30];	Kt[153]=Kbb[31];	Kt[170]=Kbb[31];	Kt[154]=Kbb[32];	Kt[188]=Kbb[32];	Kt[155]=Kbm[48];	Kt[206]=Kbm[48];	
		Kt[156]=Kbm[57];	Kt[224]=Kbm[57];	Kt[157]=Kbm[66];	Kt[242]=Kbm[66];	Kt[158]=Kbb[33];	Kt[260]=Kbb[33];	Kt[159]=Kbb[34];	Kt[278]=Kbb[34];	Kt[160]=Kbb[35];	Kt[296]=Kbb[35];	Kt[161]=Kbm[75];	Kt[314]=Kbm[75];	Kt[171]=Kbb[40];	Kt[172]=Kbb[41];	Kt[189]=Kbb[41];	Kt[173]=Kbm[49];	Kt[207]=Kbm[49];	Kt[174]=Kbm[58];	Kt[225]=Kbm[58];	Kt[175]=Kbm[67];	Kt[243]=Kbm[67];	Kt[176]=Kbb[42];	Kt[261]=Kbb[42];	Kt[177]=Kbb[43];	Kt[279]=Kbb[43];	Kt[178]=Kbb[44];	Kt[297]=Kbb[44];	Kt[179]=Kbm[76];	Kt[315]=Kbm[76];	Kt[190]=Kbb[50];	Kt[191]=Kbm[50];	Kt[208]=Kbm[50];	Kt[192]=Kbm[59];	Kt[226]=Kbm[59];	Kt[193]=Kbm[68];	Kt[244]=Kbm[68];	Kt[194]=Kbb[51];	Kt[262]=Kbb[51];	
		Kt[195]=Kbb[52];	Kt[280]=Kbb[52];	Kt[196]=Kbb[53];	Kt[298]=Kbb[53];	Kt[197]=Kbm[77];	Kt[316]=Kbm[77];	Kt[209]=Kmm[50];	Kt[210]=Kmm[51];	Kt[227]=Kmm[51];	Kt[211]=Kmm[52];	Kt[245]=Kmm[52];	Kt[212]=Kbm[51];	Kt[263]=Kbm[51];	Kt[213]=Kbm[52];	Kt[281]=Kbm[52];	Kt[214]=Kbm[53];	Kt[299]=Kbm[53];	Kt[215]=Kmm[53];	Kt[317]=Kmm[53];	Kt[228]=Kmm[60];	Kt[229]=Kmm[61];	Kt[246]=Kmm[61];	Kt[230]=Kbm[60];	Kt[264]=Kbm[60];	Kt[231]=Kbm[61];	Kt[282]=Kbm[61];	Kt[232]=Kbm[62];	Kt[300]=Kbm[62];	Kt[233]=Kmm[62];	Kt[318]=Kmm[62];	Kt[247]=Kmm[70];	Kt[248]=Kbm[69];	Kt[265]=Kbm[69];	Kt[249]=Kbm[70];	Kt[283]=Kbm[70];	Kt[250]=Kbm[71];	Kt[301]=Kbm[71];	
		Kt[251]=Kmm[71];	Kt[319]=Kmm[71];	Kt[266]=Kbb[60];	Kt[267]=Kbb[61];	Kt[284]=Kbb[61];	Kt[268]=Kbb[62];	Kt[302]=Kbb[62];	Kt[269]=Kbm[78];	Kt[320]=Kbm[78];	Kt[285]=Kbb[70];	Kt[286]=Kbb[71];	Kt[303]=Kbb[71];	Kt[287]=Kbm[79];	Kt[321]=Kbm[79];	Kt[304]=Kbb[80];	Kt[305]=Kbm[80];	Kt[322]=Kbm[80];	Kt[323]=Kmm[80];	

		// Consider shell orientation: Ke:= LambdaT*Ke*LambdaT**T
		MatrixOperations::dsyrk_f(Kt, LambdaT, Kt, 18, 3);
		
		
		// Finally add this stiffness contribution to the elementMatrix in PACKED storage,
		// using the relevant dof of this triangle.
		// Code has been generated by MATLAB.
		double * Kp = &ElementMatrix_[0];
		switch(tria){
			case 0:
				
				Kp[0]+=Kt[0];	Kp[1]+=Kt[1];	Kp[2]+=Kt[2];	Kp[3]+=Kt[3];	Kp[4]+=Kt[4];	Kp[5]+=Kt[5];	Kp[6]+=Kt[6];	Kp[7]+=Kt[7];	Kp[8]+=Kt[8];	Kp[9]+=Kt[9];	Kp[10]+=Kt[10];	Kp[11]+=Kt[11];	Kp[18]+=Kt[12];	Kp[19]+=Kt[13];	Kp[20]+=Kt[14];	Kp[21]+=Kt[15];	Kp[22]+=Kt[16];	Kp[23]+=Kt[17];	Kp[24]+=Kt[19];	Kp[25]+=Kt[20];	Kp[26]+=Kt[21];	Kp[27]+=Kt[22];	Kp[28]+=Kt[23];	Kp[29]+=Kt[24];	Kp[30]+=Kt[25];	Kp[31]+=Kt[26];	Kp[32]+=Kt[27];	Kp[33]+=Kt[28];	Kp[34]+=Kt[29];	Kp[41]+=Kt[30];	Kp[42]+=Kt[31];	Kp[43]+=Kt[32];	Kp[44]+=Kt[33];	Kp[45]+=Kt[34];	Kp[46]+=Kt[35];	Kp[47]+=Kt[38];	Kp[48]+=Kt[39];	Kp[49]+=Kt[40];	Kp[50]+=Kt[41];	Kp[51]+=Kt[42];	
				Kp[52]+=Kt[43];	Kp[53]+=Kt[44];	Kp[54]+=Kt[45];	Kp[55]+=Kt[46];	Kp[56]+=Kt[47];	Kp[63]+=Kt[48];	Kp[64]+=Kt[49];	Kp[65]+=Kt[50];	Kp[66]+=Kt[51];	Kp[67]+=Kt[52];	Kp[68]+=Kt[53];	Kp[69]+=Kt[57];	Kp[70]+=Kt[58];	Kp[71]+=Kt[59];	Kp[72]+=Kt[60];	Kp[73]+=Kt[61];	Kp[74]+=Kt[62];	Kp[75]+=Kt[63];	Kp[76]+=Kt[64];	Kp[77]+=Kt[65];	Kp[84]+=Kt[66];	Kp[85]+=Kt[67];	Kp[86]+=Kt[68];	Kp[87]+=Kt[69];	Kp[88]+=Kt[70];	Kp[89]+=Kt[71];	Kp[90]+=Kt[76];	Kp[91]+=Kt[77];	Kp[92]+=Kt[78];	Kp[93]+=Kt[79];	Kp[94]+=Kt[80];	Kp[95]+=Kt[81];	Kp[96]+=Kt[82];	Kp[97]+=Kt[83];	Kp[104]+=Kt[84];	Kp[105]+=Kt[85];	Kp[106]+=Kt[86];	Kp[107]+=Kt[87];	Kp[108]+=Kt[88];	Kp[109]+=Kt[89];	
				Kp[110]+=Kt[95];	Kp[111]+=Kt[96];	Kp[112]+=Kt[97];	Kp[113]+=Kt[98];	Kp[114]+=Kt[99];	Kp[115]+=Kt[100];	Kp[116]+=Kt[101];	Kp[123]+=Kt[102];	Kp[124]+=Kt[103];	Kp[125]+=Kt[104];	Kp[126]+=Kt[105];	Kp[127]+=Kt[106];	Kp[128]+=Kt[107];	Kp[129]+=Kt[114];	Kp[130]+=Kt[115];	Kp[131]+=Kt[116];	Kp[132]+=Kt[117];	Kp[133]+=Kt[118];	Kp[134]+=Kt[119];	Kp[141]+=Kt[120];	Kp[142]+=Kt[121];	Kp[143]+=Kt[122];	Kp[144]+=Kt[123];	Kp[145]+=Kt[124];	Kp[146]+=Kt[125];	Kp[147]+=Kt[133];	Kp[148]+=Kt[134];	Kp[149]+=Kt[135];	Kp[150]+=Kt[136];	Kp[151]+=Kt[137];	Kp[158]+=Kt[138];	Kp[159]+=Kt[139];	Kp[160]+=Kt[140];	Kp[161]+=Kt[141];	Kp[162]+=Kt[142];	Kp[163]+=Kt[143];	Kp[164]+=Kt[152];	Kp[165]+=Kt[153];	Kp[166]+=Kt[154];	Kp[167]+=Kt[155];	
				Kp[174]+=Kt[156];	Kp[175]+=Kt[157];	Kp[176]+=Kt[158];	Kp[177]+=Kt[159];	Kp[178]+=Kt[160];	Kp[179]+=Kt[161];	Kp[180]+=Kt[171];	Kp[181]+=Kt[172];	Kp[182]+=Kt[173];	Kp[189]+=Kt[174];	Kp[190]+=Kt[175];	Kp[191]+=Kt[176];	Kp[192]+=Kt[177];	Kp[193]+=Kt[178];	Kp[194]+=Kt[179];	Kp[195]+=Kt[190];	Kp[196]+=Kt[191];	Kp[203]+=Kt[192];	Kp[204]+=Kt[193];	Kp[205]+=Kt[194];	Kp[206]+=Kt[195];	Kp[207]+=Kt[196];	Kp[208]+=Kt[197];	Kp[209]+=Kt[209];	Kp[216]+=Kt[210];	Kp[217]+=Kt[211];	Kp[218]+=Kt[212];	Kp[219]+=Kt[213];	Kp[220]+=Kt[214];	Kp[221]+=Kt[215];	Kp[279]+=Kt[228];	Kp[280]+=Kt[229];	Kp[281]+=Kt[230];	Kp[282]+=Kt[231];	Kp[283]+=Kt[232];	Kp[284]+=Kt[233];	Kp[285]+=Kt[247];	Kp[286]+=Kt[248];	Kp[287]+=Kt[249];	Kp[288]+=Kt[250];	
				Kp[289]+=Kt[251];	Kp[290]+=Kt[266];	Kp[291]+=Kt[267];	Kp[292]+=Kt[268];	Kp[293]+=Kt[269];	Kp[294]+=Kt[285];	Kp[295]+=Kt[286];	Kp[296]+=Kt[287];	Kp[297]+=Kt[304];	Kp[298]+=Kt[305];	Kp[299]+=Kt[323];	

				break;
			case 1:
			
				Kp[0]+=Kt[0];	Kp[1]+=Kt[1];	Kp[2]+=Kt[2];	Kp[3]+=Kt[3];	Kp[4]+=Kt[4];	Kp[5]+=Kt[5];	Kp[6]+=Kt[6];	Kp[7]+=Kt[7];	Kp[8]+=Kt[8];	Kp[9]+=Kt[9];	Kp[10]+=Kt[10];	Kp[11]+=Kt[11];	Kp[12]+=Kt[12];	Kp[13]+=Kt[13];	Kp[14]+=Kt[14];	Kp[15]+=Kt[15];	Kp[16]+=Kt[16];	Kp[17]+=Kt[17];	Kp[24]+=Kt[19];	Kp[25]+=Kt[20];	Kp[26]+=Kt[21];	Kp[27]+=Kt[22];	Kp[28]+=Kt[23];	Kp[29]+=Kt[24];	Kp[30]+=Kt[25];	Kp[31]+=Kt[26];	Kp[32]+=Kt[27];	Kp[33]+=Kt[28];	Kp[34]+=Kt[29];	Kp[35]+=Kt[30];	Kp[36]+=Kt[31];	Kp[37]+=Kt[32];	Kp[38]+=Kt[33];	Kp[39]+=Kt[34];	Kp[40]+=Kt[35];	Kp[47]+=Kt[38];	Kp[48]+=Kt[39];	Kp[49]+=Kt[40];	Kp[50]+=Kt[41];	Kp[51]+=Kt[42];	
				Kp[52]+=Kt[43];	Kp[53]+=Kt[44];	Kp[54]+=Kt[45];	Kp[55]+=Kt[46];	Kp[56]+=Kt[47];	Kp[57]+=Kt[48];	Kp[58]+=Kt[49];	Kp[59]+=Kt[50];	Kp[60]+=Kt[51];	Kp[61]+=Kt[52];	Kp[62]+=Kt[53];	Kp[69]+=Kt[57];	Kp[70]+=Kt[58];	Kp[71]+=Kt[59];	Kp[72]+=Kt[60];	Kp[73]+=Kt[61];	Kp[74]+=Kt[62];	Kp[75]+=Kt[63];	Kp[76]+=Kt[64];	Kp[77]+=Kt[65];	Kp[78]+=Kt[66];	Kp[79]+=Kt[67];	Kp[80]+=Kt[68];	Kp[81]+=Kt[69];	Kp[82]+=Kt[70];	Kp[83]+=Kt[71];	Kp[90]+=Kt[76];	Kp[91]+=Kt[77];	Kp[92]+=Kt[78];	Kp[93]+=Kt[79];	Kp[94]+=Kt[80];	Kp[95]+=Kt[81];	Kp[96]+=Kt[82];	Kp[97]+=Kt[83];	Kp[98]+=Kt[84];	Kp[99]+=Kt[85];	Kp[100]+=Kt[86];	Kp[101]+=Kt[87];	Kp[102]+=Kt[88];	Kp[103]+=Kt[89];	
				Kp[110]+=Kt[95];	Kp[111]+=Kt[96];	Kp[112]+=Kt[97];	Kp[113]+=Kt[98];	Kp[114]+=Kt[99];	Kp[115]+=Kt[100];	Kp[116]+=Kt[101];	Kp[117]+=Kt[102];	Kp[118]+=Kt[103];	Kp[119]+=Kt[104];	Kp[120]+=Kt[105];	Kp[121]+=Kt[106];	Kp[122]+=Kt[107];	Kp[129]+=Kt[114];	Kp[130]+=Kt[115];	Kp[131]+=Kt[116];	Kp[132]+=Kt[117];	Kp[133]+=Kt[118];	Kp[134]+=Kt[119];	Kp[135]+=Kt[120];	Kp[136]+=Kt[121];	Kp[137]+=Kt[122];	Kp[138]+=Kt[123];	Kp[139]+=Kt[124];	Kp[140]+=Kt[125];	Kp[147]+=Kt[133];	Kp[148]+=Kt[134];	Kp[149]+=Kt[135];	Kp[150]+=Kt[136];	Kp[151]+=Kt[137];	Kp[152]+=Kt[138];	Kp[153]+=Kt[139];	Kp[154]+=Kt[140];	Kp[155]+=Kt[141];	Kp[156]+=Kt[142];	Kp[157]+=Kt[143];	Kp[164]+=Kt[152];	Kp[165]+=Kt[153];	Kp[166]+=Kt[154];	Kp[167]+=Kt[155];	
				Kp[168]+=Kt[156];	Kp[169]+=Kt[157];	Kp[170]+=Kt[158];	Kp[171]+=Kt[159];	Kp[172]+=Kt[160];	Kp[173]+=Kt[161];	Kp[180]+=Kt[171];	Kp[181]+=Kt[172];	Kp[182]+=Kt[173];	Kp[183]+=Kt[174];	Kp[184]+=Kt[175];	Kp[185]+=Kt[176];	Kp[186]+=Kt[177];	Kp[187]+=Kt[178];	Kp[188]+=Kt[179];	Kp[195]+=Kt[190];	Kp[196]+=Kt[191];	Kp[197]+=Kt[192];	Kp[198]+=Kt[193];	Kp[199]+=Kt[194];	Kp[200]+=Kt[195];	Kp[201]+=Kt[196];	Kp[202]+=Kt[197];	Kp[209]+=Kt[209];	Kp[210]+=Kt[210];	Kp[211]+=Kt[211];	Kp[212]+=Kt[212];	Kp[213]+=Kt[213];	Kp[214]+=Kt[214];	Kp[215]+=Kt[215];	Kp[222]+=Kt[228];	Kp[223]+=Kt[229];	Kp[224]+=Kt[230];	Kp[225]+=Kt[231];	Kp[226]+=Kt[232];	Kp[227]+=Kt[233];	Kp[234]+=Kt[247];	Kp[235]+=Kt[248];	Kp[236]+=Kt[249];	Kp[237]+=Kt[250];	
				Kp[238]+=Kt[251];	Kp[245]+=Kt[266];	Kp[246]+=Kt[267];	Kp[247]+=Kt[268];	Kp[248]+=Kt[269];	Kp[255]+=Kt[285];	Kp[256]+=Kt[286];	Kp[257]+=Kt[287];	Kp[264]+=Kt[304];	Kp[265]+=Kt[305];	Kp[272]+=Kt[323];	

				break;
			case 2:
			
				Kp[129]+=Kt[0];	Kp[130]+=Kt[1];	Kp[131]+=Kt[2];	Kp[132]+=Kt[3];	Kp[133]+=Kt[4];	Kp[134]+=Kt[5];	Kp[135]+=Kt[6];	Kp[136]+=Kt[7];	Kp[137]+=Kt[8];	Kp[138]+=Kt[9];	Kp[139]+=Kt[10];	Kp[140]+=Kt[11];	Kp[141]+=Kt[12];	Kp[142]+=Kt[13];	Kp[143]+=Kt[14];	Kp[144]+=Kt[15];	Kp[145]+=Kt[16];	Kp[146]+=Kt[17];	Kp[147]+=Kt[19];	Kp[148]+=Kt[20];	Kp[149]+=Kt[21];	Kp[150]+=Kt[22];	Kp[151]+=Kt[23];	Kp[152]+=Kt[24];	Kp[153]+=Kt[25];	Kp[154]+=Kt[26];	Kp[155]+=Kt[27];	Kp[156]+=Kt[28];	Kp[157]+=Kt[29];	Kp[158]+=Kt[30];	Kp[159]+=Kt[31];	Kp[160]+=Kt[32];	Kp[161]+=Kt[33];	Kp[162]+=Kt[34];	Kp[163]+=Kt[35];	Kp[164]+=Kt[38];	Kp[165]+=Kt[39];	Kp[166]+=Kt[40];	Kp[167]+=Kt[41];	Kp[168]+=Kt[42];	
				Kp[169]+=Kt[43];	Kp[170]+=Kt[44];	Kp[171]+=Kt[45];	Kp[172]+=Kt[46];	Kp[173]+=Kt[47];	Kp[174]+=Kt[48];	Kp[175]+=Kt[49];	Kp[176]+=Kt[50];	Kp[177]+=Kt[51];	Kp[178]+=Kt[52];	Kp[179]+=Kt[53];	Kp[180]+=Kt[57];	Kp[181]+=Kt[58];	Kp[182]+=Kt[59];	Kp[183]+=Kt[60];	Kp[184]+=Kt[61];	Kp[185]+=Kt[62];	Kp[186]+=Kt[63];	Kp[187]+=Kt[64];	Kp[188]+=Kt[65];	Kp[189]+=Kt[66];	Kp[190]+=Kt[67];	Kp[191]+=Kt[68];	Kp[192]+=Kt[69];	Kp[193]+=Kt[70];	Kp[194]+=Kt[71];	Kp[195]+=Kt[76];	Kp[196]+=Kt[77];	Kp[197]+=Kt[78];	Kp[198]+=Kt[79];	Kp[199]+=Kt[80];	Kp[200]+=Kt[81];	Kp[201]+=Kt[82];	Kp[202]+=Kt[83];	Kp[203]+=Kt[84];	Kp[204]+=Kt[85];	Kp[205]+=Kt[86];	Kp[206]+=Kt[87];	Kp[207]+=Kt[88];	Kp[208]+=Kt[89];	
				Kp[209]+=Kt[95];	Kp[210]+=Kt[96];	Kp[211]+=Kt[97];	Kp[212]+=Kt[98];	Kp[213]+=Kt[99];	Kp[214]+=Kt[100];	Kp[215]+=Kt[101];	Kp[216]+=Kt[102];	Kp[217]+=Kt[103];	Kp[218]+=Kt[104];	Kp[219]+=Kt[105];	Kp[220]+=Kt[106];	Kp[221]+=Kt[107];	Kp[222]+=Kt[114];	Kp[223]+=Kt[115];	Kp[224]+=Kt[116];	Kp[225]+=Kt[117];	Kp[226]+=Kt[118];	Kp[227]+=Kt[119];	Kp[228]+=Kt[120];	Kp[229]+=Kt[121];	Kp[230]+=Kt[122];	Kp[231]+=Kt[123];	Kp[232]+=Kt[124];	Kp[233]+=Kt[125];	Kp[234]+=Kt[133];	Kp[235]+=Kt[134];	Kp[236]+=Kt[135];	Kp[237]+=Kt[136];	Kp[238]+=Kt[137];	Kp[239]+=Kt[138];	Kp[240]+=Kt[139];	Kp[241]+=Kt[140];	Kp[242]+=Kt[141];	Kp[243]+=Kt[142];	Kp[244]+=Kt[143];	Kp[245]+=Kt[152];	Kp[246]+=Kt[153];	Kp[247]+=Kt[154];	Kp[248]+=Kt[155];	
				Kp[249]+=Kt[156];	Kp[250]+=Kt[157];	Kp[251]+=Kt[158];	Kp[252]+=Kt[159];	Kp[253]+=Kt[160];	Kp[254]+=Kt[161];	Kp[255]+=Kt[171];	Kp[256]+=Kt[172];	Kp[257]+=Kt[173];	Kp[258]+=Kt[174];	Kp[259]+=Kt[175];	Kp[260]+=Kt[176];	Kp[261]+=Kt[177];	Kp[262]+=Kt[178];	Kp[263]+=Kt[179];	Kp[264]+=Kt[190];	Kp[265]+=Kt[191];	Kp[266]+=Kt[192];	Kp[267]+=Kt[193];	Kp[268]+=Kt[194];	Kp[269]+=Kt[195];	Kp[270]+=Kt[196];	Kp[271]+=Kt[197];	Kp[272]+=Kt[209];	Kp[273]+=Kt[210];	Kp[274]+=Kt[211];	Kp[275]+=Kt[212];	Kp[276]+=Kt[213];	Kp[277]+=Kt[214];	Kp[278]+=Kt[215];	Kp[279]+=Kt[228];	Kp[280]+=Kt[229];	Kp[281]+=Kt[230];	Kp[282]+=Kt[231];	Kp[283]+=Kt[232];	Kp[284]+=Kt[233];	Kp[285]+=Kt[247];	Kp[286]+=Kt[248];	Kp[287]+=Kt[249];	Kp[288]+=Kt[250];	
				Kp[289]+=Kt[251];	Kp[290]+=Kt[266];	Kp[291]+=Kt[267];	Kp[292]+=Kt[268];	Kp[293]+=Kt[269];	Kp[294]+=Kt[285];	Kp[295]+=Kt[286];	Kp[296]+=Kt[287];	Kp[297]+=Kt[304];	Kp[298]+=Kt[305];	Kp[299]+=Kt[323];	

				break;
			case 3:
			
				Kp[0]+=Kt[0];	Kp[1]+=Kt[1];	Kp[2]+=Kt[2];	Kp[3]+=Kt[3];	Kp[4]+=Kt[4];	Kp[5]+=Kt[5];	Kp[12]+=Kt[6];	Kp[13]+=Kt[7];	Kp[14]+=Kt[8];	Kp[15]+=Kt[9];	Kp[16]+=Kt[10];	Kp[17]+=Kt[11];	Kp[18]+=Kt[12];	Kp[19]+=Kt[13];	Kp[20]+=Kt[14];	Kp[21]+=Kt[15];	Kp[22]+=Kt[16];	Kp[23]+=Kt[17];	Kp[24]+=Kt[19];	Kp[25]+=Kt[20];	Kp[26]+=Kt[21];	Kp[27]+=Kt[22];	Kp[28]+=Kt[23];	Kp[35]+=Kt[24];	Kp[36]+=Kt[25];	Kp[37]+=Kt[26];	Kp[38]+=Kt[27];	Kp[39]+=Kt[28];	Kp[40]+=Kt[29];	Kp[41]+=Kt[30];	Kp[42]+=Kt[31];	Kp[43]+=Kt[32];	Kp[44]+=Kt[33];	Kp[45]+=Kt[34];	Kp[46]+=Kt[35];	Kp[47]+=Kt[38];	Kp[48]+=Kt[39];	Kp[49]+=Kt[40];	Kp[50]+=Kt[41];	Kp[57]+=Kt[42];	
				Kp[58]+=Kt[43];	Kp[59]+=Kt[44];	Kp[60]+=Kt[45];	Kp[61]+=Kt[46];	Kp[62]+=Kt[47];	Kp[63]+=Kt[48];	Kp[64]+=Kt[49];	Kp[65]+=Kt[50];	Kp[66]+=Kt[51];	Kp[67]+=Kt[52];	Kp[68]+=Kt[53];	Kp[69]+=Kt[57];	Kp[70]+=Kt[58];	Kp[71]+=Kt[59];	Kp[78]+=Kt[60];	Kp[79]+=Kt[61];	Kp[80]+=Kt[62];	Kp[81]+=Kt[63];	Kp[82]+=Kt[64];	Kp[83]+=Kt[65];	Kp[84]+=Kt[66];	Kp[85]+=Kt[67];	Kp[86]+=Kt[68];	Kp[87]+=Kt[69];	Kp[88]+=Kt[70];	Kp[89]+=Kt[71];	Kp[90]+=Kt[76];	Kp[91]+=Kt[77];	Kp[98]+=Kt[78];	Kp[99]+=Kt[79];	Kp[100]+=Kt[80];	Kp[101]+=Kt[81];	Kp[102]+=Kt[82];	Kp[103]+=Kt[83];	Kp[104]+=Kt[84];	Kp[105]+=Kt[85];	Kp[106]+=Kt[86];	Kp[107]+=Kt[87];	Kp[108]+=Kt[88];	Kp[109]+=Kt[89];	
				Kp[110]+=Kt[95];	Kp[117]+=Kt[96];	Kp[118]+=Kt[97];	Kp[119]+=Kt[98];	Kp[120]+=Kt[99];	Kp[121]+=Kt[100];	Kp[122]+=Kt[101];	Kp[123]+=Kt[102];	Kp[124]+=Kt[103];	Kp[125]+=Kt[104];	Kp[126]+=Kt[105];	Kp[127]+=Kt[106];	Kp[128]+=Kt[107];	Kp[222]+=Kt[114];	Kp[223]+=Kt[115];	Kp[224]+=Kt[116];	Kp[225]+=Kt[117];	Kp[226]+=Kt[118];	Kp[227]+=Kt[119];	Kp[228]+=Kt[120];	Kp[229]+=Kt[121];	Kp[230]+=Kt[122];	Kp[231]+=Kt[123];	Kp[232]+=Kt[124];	Kp[233]+=Kt[125];	Kp[234]+=Kt[133];	Kp[235]+=Kt[134];	Kp[236]+=Kt[135];	Kp[237]+=Kt[136];	Kp[238]+=Kt[137];	Kp[239]+=Kt[138];	Kp[240]+=Kt[139];	Kp[241]+=Kt[140];	Kp[242]+=Kt[141];	Kp[243]+=Kt[142];	Kp[244]+=Kt[143];	Kp[245]+=Kt[152];	Kp[246]+=Kt[153];	Kp[247]+=Kt[154];	Kp[248]+=Kt[155];	
				Kp[249]+=Kt[156];	Kp[250]+=Kt[157];	Kp[251]+=Kt[158];	Kp[252]+=Kt[159];	Kp[253]+=Kt[160];	Kp[254]+=Kt[161];	Kp[255]+=Kt[171];	Kp[256]+=Kt[172];	Kp[257]+=Kt[173];	Kp[258]+=Kt[174];	Kp[259]+=Kt[175];	Kp[260]+=Kt[176];	Kp[261]+=Kt[177];	Kp[262]+=Kt[178];	Kp[263]+=Kt[179];	Kp[264]+=Kt[190];	Kp[265]+=Kt[191];	Kp[266]+=Kt[192];	Kp[267]+=Kt[193];	Kp[268]+=Kt[194];	Kp[269]+=Kt[195];	Kp[270]+=Kt[196];	Kp[271]+=Kt[197];	Kp[272]+=Kt[209];	Kp[273]+=Kt[210];	Kp[274]+=Kt[211];	Kp[275]+=Kt[212];	Kp[276]+=Kt[213];	Kp[277]+=Kt[214];	Kp[278]+=Kt[215];	Kp[279]+=Kt[228];	Kp[280]+=Kt[229];	Kp[281]+=Kt[230];	Kp[282]+=Kt[231];	Kp[283]+=Kt[232];	Kp[284]+=Kt[233];	Kp[285]+=Kt[247];	Kp[286]+=Kt[248];	Kp[287]+=Kt[249];	Kp[288]+=Kt[250];	
				Kp[289]+=Kt[251];	Kp[290]+=Kt[266];	Kp[291]+=Kt[267];	Kp[292]+=Kt[268];	Kp[293]+=Kt[269];	Kp[294]+=Kt[285];	Kp[295]+=Kt[286];	Kp[296]+=Kt[287];	Kp[297]+=Kt[304];	Kp[298]+=Kt[305];	Kp[299]+=Kt[323];	

				break;
			default:
				FEI_COUT << " ERROR: COMPUTING ELEMENT MATRIX OF CQUAD ELEMENT FAILED" << FEI_ENDL;
				return 1;
		}
	}
     
	return 0;
}

// Build geometric stiffness matrix of the element
int CQUAD::BuildGeometricElementMatrix(Solver *solver, bool makeNegative){
    // Allocate storage
    int size = 24;
	double storage[427];
	for(int i=0;i<427;i++)storage[i] = 0;
	
	double* elemStrains = &storage[0];		// Element Strains (size 6x5 = 30)
	double* elemForces  = &storage[30];		// Element Forces  (size 6x5 = 30)
	
	double* dNdx        = &storage[60];     // Derivative of shape function of Triangle in x dir (size 3x1 = 3)
	double* dNdy        = &storage[63];     // Derivative of shape function of Triangle in y dir (size 3x1 = 3)
	
	double* LambdaT     = &storage[66];		// Contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 3x3 = 9)
	double* Nodes2d     = &storage[75];     // Location of the the three nodes in element coordinates (size 3x2 = 6)
	double* d1          = &storage[81];		// Relative nodal position 1 vector and normal vector of a Triangle (size 3)
	double* d2          = &storage[84];		// Relative nodal position 2 vector and normal vector of a Triangle (size 3)
	double* n           = &storage[87];		// Normal vector of a Triangle (size 3)
	
	double* Be1         = &storage[90];		// Contains geometric data of Triangular element (size 13)
	double* KgT         = &storage[103];     // Element geometric stiffness matrix of one Triangle (size 18x18 = 324)
    
	// Initialize vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle
    
    // Compute element strains and element forces
    ElementStrains(elemStrains, solver);
    ElementForces(elemForces, elemStrains);
    
    // Loop through four triangles
	for(std::size_t tria=0; tria<4; tria++){
		// Initialize storage to zero (except for first 60 terms which contain element forces and strains)
		for(std::size_t i=60;i<427;i++)storage[i] = 0;
		
		// Get nodes corresponding to this triangle
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute the geometric characteristics of this element
		
		// Compute normal direction of element, which is e_z
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,&LambdaT[6]); // ez=  d1 x d2
		double lN = sqrt(MatrixOperations::Dot(&LambdaT[6],&LambdaT[6],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[6],&LambdaT[6], 3); // ez = 1/norm(ez);
		
		// Compute e_x = emat-(emat'*ez)*ez;
		MatrixOperations::Scale(MatrixOperations::Dot(&LambdaT[6],MaterialDirection_,3),
								0, &LambdaT[6],&LambdaT[0], 3); // ex = dot(emat,ez)*ez;
		MatrixOperations::Add(1,-1,MaterialDirection_, &LambdaT[0],&LambdaT[0],3); //ex = emat-ex;
		lN = sqrt(MatrixOperations::Dot(&LambdaT[0],&LambdaT[0],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[0],&LambdaT[0], 3); // ex = 1/norm(ex);
		
		// Compute ey = cross(ez,ex);
		MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
		
		// Compute nodes2D
		Nodes2d[0] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[1] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[2] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[3] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[4] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[5] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[3],3);
		
		// Store coordinate differences in Be1
		Be1[0] = Nodes2d[0]-Nodes2d[2]; Be1[1] = Nodes2d[1]-Nodes2d[3]; // Be1[0] = x12;  Be1[1] = y12 
		Be1[2] = -Be1[0]; 				Be1[3] = -Be1[1]; 				// Be1[2] = x21;  Be1[3] = y21 
		Be1[4] = Nodes2d[2]-Nodes2d[4]; Be1[5] = Nodes2d[3]-Nodes2d[5]; // Be1[4] = x23;  Be1[5] = y23 
		Be1[6] = -Be1[4]; 				Be1[7] = -Be1[5]; 				// Be1[6] = x32;  Be1[7] = y32 
		Be1[8] = Nodes2d[4]-Nodes2d[0]; Be1[9] = Nodes2d[5]-Nodes2d[1]; // Be1[8] = x31;  Be1[9] = y31 
		Be1[10]= -Be1[8]; 				Be1[11]= -Be1[9]; 				// Be1[10]= x13;  Be1[11]= y13
		Be1[12]= (Be1[3]*Be1[10]-Be1[2]*Be1[11])/2; 					// Be1[12]= area
		
		// Compute dNdx and dNdy
		dNdx[0] = 0.5*Be1[5]/Be1[12];	//0.5*y23/area 
		dNdx[1] = 0.5*Be1[9]/Be1[12];	//0.5*y31/area 
		dNdx[2] = 0.5*Be1[1]/Be1[12];	//0.5*y12/area 
		
		dNdy[0] = 0.5*Be1[6]/Be1[12];	//0.5*x32/area 
		dNdy[1] = 0.5*Be1[10]/Be1[12];	//0.5*x13/area 
		dNdy[2] = 0.5*Be1[2]/Be1[12];	//0.5*x21/area 
		
		// Get normal forces from forcevector
		double &Nx=elemForces[6*tria+0], &Ny=elemForces[6*tria+1], &Nxy=elemForces[6*tria+2];	// Nx, Ny, Nxy in elemForces
		
		// Compute geometric element of this Triangle in local coordinates, 
		// using Kg=-Bg**T*[Nx,Nxy;Nxy,Ny]*Bg * area with Bg(2x18) contains dNdx and dNdy
		// NOTE: The 0.5 term is added, to ensure the contribution of each triangle is halved.
		KgT[38] = -0.5*(Nx*dNdx[0]*dNdx[0] + 2*Nxy*dNdx[0]*dNdy[0] 							+ Ny*dNdy[0]*dNdy[0])*Be1[12];
		KgT[44] = -0.5*(Nx*dNdx[0]*dNdx[1] +   Nxy*dNdx[0]*dNdy[1] + Nxy*dNdx[1]*dNdy[0] 	+ Ny*dNdy[0]*dNdy[1])*Be1[12];
		KgT[50] = -0.5*(Nx*dNdx[0]*dNdx[2] + 	 Nxy*dNdx[0]*dNdy[2] + Nxy*dNdx[2]*dNdy[0]	+ Ny*dNdy[0]*dNdy[2])*Be1[12];
		
		KgT[152]= -0.5*(Nx*dNdx[1]*dNdx[1] + 2*Nxy*dNdx[1]*dNdy[1] 							+ Ny*dNdy[1]*dNdy[1])*Be1[12];
		KgT[158]= -0.5*(Nx*dNdx[1]*dNdx[2] + 	 Nxy*dNdx[1]*dNdy[2] + Nxy*dNdx[2]*dNdy[1]	+ Ny*dNdy[1]*dNdy[2])*Be1[12];
		
		KgT[266]= -0.5*(Nx*dNdx[2]*dNdx[2] + 2*Nxy*dNdx[2]*dNdy[2] 							+ Ny*dNdy[2]*dNdy[2])*Be1[12];
		
		KgT[146] = KgT[44]; KgT[254]=KgT[50]; KgT[260]=KgT[158];
		
		// Consider shell orientation: Kg:= LambdaT*Kg*LambdaT**T
		MatrixOperations::dsyrk_f(KgT, LambdaT, KgT, 18, 3);
		
		// Make KgT its negative if so required (used for differential stiffness of prestress)
		if(makeNegative) for(std::size_t i=0;i<18*18;i++) KgT[i]*=-1;	
				
		// Finally add this stiffness contribution to the elementMatrix in PACKED storage,
		// using the relevant dof of this triangle.
		// Code has been generated by MATLAB.
		double * Kgp = &GeometricElementMatrix_[0];
		switch(tria){
			case 0:
				
				Kgp[0]+=KgT[0];	Kgp[1]+=KgT[1];	Kgp[2]+=KgT[2];	Kgp[3]+=KgT[3];	Kgp[4]+=KgT[4];	Kgp[5]+=KgT[5];	Kgp[6]+=KgT[6];	Kgp[7]+=KgT[7];	Kgp[8]+=KgT[8];	Kgp[9]+=KgT[9];	Kgp[10]+=KgT[10];	Kgp[11]+=KgT[11];	Kgp[18]+=KgT[12];	Kgp[19]+=KgT[13];	Kgp[20]+=KgT[14];	Kgp[21]+=KgT[15];	Kgp[22]+=KgT[16];	Kgp[23]+=KgT[17];	Kgp[24]+=KgT[19];	Kgp[25]+=KgT[20];	Kgp[26]+=KgT[21];	Kgp[27]+=KgT[22];	Kgp[28]+=KgT[23];	Kgp[29]+=KgT[24];	Kgp[30]+=KgT[25];	Kgp[31]+=KgT[26];	Kgp[32]+=KgT[27];	Kgp[33]+=KgT[28];	Kgp[34]+=KgT[29];	Kgp[41]+=KgT[30];	Kgp[42]+=KgT[31];	Kgp[43]+=KgT[32];	Kgp[44]+=KgT[33];	Kgp[45]+=KgT[34];	Kgp[46]+=KgT[35];	Kgp[47]+=KgT[38];	Kgp[48]+=KgT[39];	Kgp[49]+=KgT[40];	Kgp[50]+=KgT[41];	Kgp[51]+=KgT[42];	
				Kgp[52]+=KgT[43];	Kgp[53]+=KgT[44];	Kgp[54]+=KgT[45];	Kgp[55]+=KgT[46];	Kgp[56]+=KgT[47];	Kgp[63]+=KgT[48];	Kgp[64]+=KgT[49];	Kgp[65]+=KgT[50];	Kgp[66]+=KgT[51];	Kgp[67]+=KgT[52];	Kgp[68]+=KgT[53];	Kgp[69]+=KgT[57];	Kgp[70]+=KgT[58];	Kgp[71]+=KgT[59];	Kgp[72]+=KgT[60];	Kgp[73]+=KgT[61];	Kgp[74]+=KgT[62];	Kgp[75]+=KgT[63];	Kgp[76]+=KgT[64];	Kgp[77]+=KgT[65];	Kgp[84]+=KgT[66];	Kgp[85]+=KgT[67];	Kgp[86]+=KgT[68];	Kgp[87]+=KgT[69];	Kgp[88]+=KgT[70];	Kgp[89]+=KgT[71];	Kgp[90]+=KgT[76];	Kgp[91]+=KgT[77];	Kgp[92]+=KgT[78];	Kgp[93]+=KgT[79];	Kgp[94]+=KgT[80];	Kgp[95]+=KgT[81];	Kgp[96]+=KgT[82];	Kgp[97]+=KgT[83];	Kgp[104]+=KgT[84];	Kgp[105]+=KgT[85];	Kgp[106]+=KgT[86];	Kgp[107]+=KgT[87];	Kgp[108]+=KgT[88];	Kgp[109]+=KgT[89];	
				Kgp[110]+=KgT[95];	Kgp[111]+=KgT[96];	Kgp[112]+=KgT[97];	Kgp[113]+=KgT[98];	Kgp[114]+=KgT[99];	Kgp[115]+=KgT[100];	Kgp[116]+=KgT[101];	Kgp[123]+=KgT[102];	Kgp[124]+=KgT[103];	Kgp[125]+=KgT[104];	Kgp[126]+=KgT[105];	Kgp[127]+=KgT[106];	Kgp[128]+=KgT[107];	Kgp[129]+=KgT[114];	Kgp[130]+=KgT[115];	Kgp[131]+=KgT[116];	Kgp[132]+=KgT[117];	Kgp[133]+=KgT[118];	Kgp[134]+=KgT[119];	Kgp[141]+=KgT[120];	Kgp[142]+=KgT[121];	Kgp[143]+=KgT[122];	Kgp[144]+=KgT[123];	Kgp[145]+=KgT[124];	Kgp[146]+=KgT[125];	Kgp[147]+=KgT[133];	Kgp[148]+=KgT[134];	Kgp[149]+=KgT[135];	Kgp[150]+=KgT[136];	Kgp[151]+=KgT[137];	Kgp[158]+=KgT[138];	Kgp[159]+=KgT[139];	Kgp[160]+=KgT[140];	Kgp[161]+=KgT[141];	Kgp[162]+=KgT[142];	Kgp[163]+=KgT[143];	Kgp[164]+=KgT[152];	Kgp[165]+=KgT[153];	Kgp[166]+=KgT[154];	Kgp[167]+=KgT[155];	
				Kgp[174]+=KgT[156];	Kgp[175]+=KgT[157];	Kgp[176]+=KgT[158];	Kgp[177]+=KgT[159];	Kgp[178]+=KgT[160];	Kgp[179]+=KgT[161];	Kgp[180]+=KgT[171];	Kgp[181]+=KgT[172];	Kgp[182]+=KgT[173];	Kgp[189]+=KgT[174];	Kgp[190]+=KgT[175];	Kgp[191]+=KgT[176];	Kgp[192]+=KgT[177];	Kgp[193]+=KgT[178];	Kgp[194]+=KgT[179];	Kgp[195]+=KgT[190];	Kgp[196]+=KgT[191];	Kgp[203]+=KgT[192];	Kgp[204]+=KgT[193];	Kgp[205]+=KgT[194];	Kgp[206]+=KgT[195];	Kgp[207]+=KgT[196];	Kgp[208]+=KgT[197];	Kgp[209]+=KgT[209];	Kgp[216]+=KgT[210];	Kgp[217]+=KgT[211];	Kgp[218]+=KgT[212];	Kgp[219]+=KgT[213];	Kgp[220]+=KgT[214];	Kgp[221]+=KgT[215];	Kgp[279]+=KgT[228];	Kgp[280]+=KgT[229];	Kgp[281]+=KgT[230];	Kgp[282]+=KgT[231];	Kgp[283]+=KgT[232];	Kgp[284]+=KgT[233];	Kgp[285]+=KgT[247];	Kgp[286]+=KgT[248];	Kgp[287]+=KgT[249];	Kgp[288]+=KgT[250];	
				Kgp[289]+=KgT[251];	Kgp[290]+=KgT[266];	Kgp[291]+=KgT[267];	Kgp[292]+=KgT[268];	Kgp[293]+=KgT[269];	Kgp[294]+=KgT[285];	Kgp[295]+=KgT[286];	Kgp[296]+=KgT[287];	Kgp[297]+=KgT[304];	Kgp[298]+=KgT[305];	Kgp[299]+=KgT[323];	

				break;
			case 1:
			
				Kgp[0]+=KgT[0];	Kgp[1]+=KgT[1];	Kgp[2]+=KgT[2];	Kgp[3]+=KgT[3];	Kgp[4]+=KgT[4];	Kgp[5]+=KgT[5];	Kgp[6]+=KgT[6];	Kgp[7]+=KgT[7];	Kgp[8]+=KgT[8];	Kgp[9]+=KgT[9];	Kgp[10]+=KgT[10];	Kgp[11]+=KgT[11];	Kgp[12]+=KgT[12];	Kgp[13]+=KgT[13];	Kgp[14]+=KgT[14];	Kgp[15]+=KgT[15];	Kgp[16]+=KgT[16];	Kgp[17]+=KgT[17];	Kgp[24]+=KgT[19];	Kgp[25]+=KgT[20];	Kgp[26]+=KgT[21];	Kgp[27]+=KgT[22];	Kgp[28]+=KgT[23];	Kgp[29]+=KgT[24];	Kgp[30]+=KgT[25];	Kgp[31]+=KgT[26];	Kgp[32]+=KgT[27];	Kgp[33]+=KgT[28];	Kgp[34]+=KgT[29];	Kgp[35]+=KgT[30];	Kgp[36]+=KgT[31];	Kgp[37]+=KgT[32];	Kgp[38]+=KgT[33];	Kgp[39]+=KgT[34];	Kgp[40]+=KgT[35];	Kgp[47]+=KgT[38];	Kgp[48]+=KgT[39];	Kgp[49]+=KgT[40];	Kgp[50]+=KgT[41];	Kgp[51]+=KgT[42];	
				Kgp[52]+=KgT[43];	Kgp[53]+=KgT[44];	Kgp[54]+=KgT[45];	Kgp[55]+=KgT[46];	Kgp[56]+=KgT[47];	Kgp[57]+=KgT[48];	Kgp[58]+=KgT[49];	Kgp[59]+=KgT[50];	Kgp[60]+=KgT[51];	Kgp[61]+=KgT[52];	Kgp[62]+=KgT[53];	Kgp[69]+=KgT[57];	Kgp[70]+=KgT[58];	Kgp[71]+=KgT[59];	Kgp[72]+=KgT[60];	Kgp[73]+=KgT[61];	Kgp[74]+=KgT[62];	Kgp[75]+=KgT[63];	Kgp[76]+=KgT[64];	Kgp[77]+=KgT[65];	Kgp[78]+=KgT[66];	Kgp[79]+=KgT[67];	Kgp[80]+=KgT[68];	Kgp[81]+=KgT[69];	Kgp[82]+=KgT[70];	Kgp[83]+=KgT[71];	Kgp[90]+=KgT[76];	Kgp[91]+=KgT[77];	Kgp[92]+=KgT[78];	Kgp[93]+=KgT[79];	Kgp[94]+=KgT[80];	Kgp[95]+=KgT[81];	Kgp[96]+=KgT[82];	Kgp[97]+=KgT[83];	Kgp[98]+=KgT[84];	Kgp[99]+=KgT[85];	Kgp[100]+=KgT[86];	Kgp[101]+=KgT[87];	Kgp[102]+=KgT[88];	Kgp[103]+=KgT[89];	
				Kgp[110]+=KgT[95];	Kgp[111]+=KgT[96];	Kgp[112]+=KgT[97];	Kgp[113]+=KgT[98];	Kgp[114]+=KgT[99];	Kgp[115]+=KgT[100];	Kgp[116]+=KgT[101];	Kgp[117]+=KgT[102];	Kgp[118]+=KgT[103];	Kgp[119]+=KgT[104];	Kgp[120]+=KgT[105];	Kgp[121]+=KgT[106];	Kgp[122]+=KgT[107];	Kgp[129]+=KgT[114];	Kgp[130]+=KgT[115];	Kgp[131]+=KgT[116];	Kgp[132]+=KgT[117];	Kgp[133]+=KgT[118];	Kgp[134]+=KgT[119];	Kgp[135]+=KgT[120];	Kgp[136]+=KgT[121];	Kgp[137]+=KgT[122];	Kgp[138]+=KgT[123];	Kgp[139]+=KgT[124];	Kgp[140]+=KgT[125];	Kgp[147]+=KgT[133];	Kgp[148]+=KgT[134];	Kgp[149]+=KgT[135];	Kgp[150]+=KgT[136];	Kgp[151]+=KgT[137];	Kgp[152]+=KgT[138];	Kgp[153]+=KgT[139];	Kgp[154]+=KgT[140];	Kgp[155]+=KgT[141];	Kgp[156]+=KgT[142];	Kgp[157]+=KgT[143];	Kgp[164]+=KgT[152];	Kgp[165]+=KgT[153];	Kgp[166]+=KgT[154];	Kgp[167]+=KgT[155];	
				Kgp[168]+=KgT[156];	Kgp[169]+=KgT[157];	Kgp[170]+=KgT[158];	Kgp[171]+=KgT[159];	Kgp[172]+=KgT[160];	Kgp[173]+=KgT[161];	Kgp[180]+=KgT[171];	Kgp[181]+=KgT[172];	Kgp[182]+=KgT[173];	Kgp[183]+=KgT[174];	Kgp[184]+=KgT[175];	Kgp[185]+=KgT[176];	Kgp[186]+=KgT[177];	Kgp[187]+=KgT[178];	Kgp[188]+=KgT[179];	Kgp[195]+=KgT[190];	Kgp[196]+=KgT[191];	Kgp[197]+=KgT[192];	Kgp[198]+=KgT[193];	Kgp[199]+=KgT[194];	Kgp[200]+=KgT[195];	Kgp[201]+=KgT[196];	Kgp[202]+=KgT[197];	Kgp[209]+=KgT[209];	Kgp[210]+=KgT[210];	Kgp[211]+=KgT[211];	Kgp[212]+=KgT[212];	Kgp[213]+=KgT[213];	Kgp[214]+=KgT[214];	Kgp[215]+=KgT[215];	Kgp[222]+=KgT[228];	Kgp[223]+=KgT[229];	Kgp[224]+=KgT[230];	Kgp[225]+=KgT[231];	Kgp[226]+=KgT[232];	Kgp[227]+=KgT[233];	Kgp[234]+=KgT[247];	Kgp[235]+=KgT[248];	Kgp[236]+=KgT[249];	Kgp[237]+=KgT[250];	
				Kgp[238]+=KgT[251];	Kgp[245]+=KgT[266];	Kgp[246]+=KgT[267];	Kgp[247]+=KgT[268];	Kgp[248]+=KgT[269];	Kgp[255]+=KgT[285];	Kgp[256]+=KgT[286];	Kgp[257]+=KgT[287];	Kgp[264]+=KgT[304];	Kgp[265]+=KgT[305];	Kgp[272]+=KgT[323];	

				break;
			case 2:
			
				Kgp[129]+=KgT[0];	Kgp[130]+=KgT[1];	Kgp[131]+=KgT[2];	Kgp[132]+=KgT[3];	Kgp[133]+=KgT[4];	Kgp[134]+=KgT[5];	Kgp[135]+=KgT[6];	Kgp[136]+=KgT[7];	Kgp[137]+=KgT[8];	Kgp[138]+=KgT[9];	Kgp[139]+=KgT[10];	Kgp[140]+=KgT[11];	Kgp[141]+=KgT[12];	Kgp[142]+=KgT[13];	Kgp[143]+=KgT[14];	Kgp[144]+=KgT[15];	Kgp[145]+=KgT[16];	Kgp[146]+=KgT[17];	Kgp[147]+=KgT[19];	Kgp[148]+=KgT[20];	Kgp[149]+=KgT[21];	Kgp[150]+=KgT[22];	Kgp[151]+=KgT[23];	Kgp[152]+=KgT[24];	Kgp[153]+=KgT[25];	Kgp[154]+=KgT[26];	Kgp[155]+=KgT[27];	Kgp[156]+=KgT[28];	Kgp[157]+=KgT[29];	Kgp[158]+=KgT[30];	Kgp[159]+=KgT[31];	Kgp[160]+=KgT[32];	Kgp[161]+=KgT[33];	Kgp[162]+=KgT[34];	Kgp[163]+=KgT[35];	Kgp[164]+=KgT[38];	Kgp[165]+=KgT[39];	Kgp[166]+=KgT[40];	Kgp[167]+=KgT[41];	Kgp[168]+=KgT[42];	
				Kgp[169]+=KgT[43];	Kgp[170]+=KgT[44];	Kgp[171]+=KgT[45];	Kgp[172]+=KgT[46];	Kgp[173]+=KgT[47];	Kgp[174]+=KgT[48];	Kgp[175]+=KgT[49];	Kgp[176]+=KgT[50];	Kgp[177]+=KgT[51];	Kgp[178]+=KgT[52];	Kgp[179]+=KgT[53];	Kgp[180]+=KgT[57];	Kgp[181]+=KgT[58];	Kgp[182]+=KgT[59];	Kgp[183]+=KgT[60];	Kgp[184]+=KgT[61];	Kgp[185]+=KgT[62];	Kgp[186]+=KgT[63];	Kgp[187]+=KgT[64];	Kgp[188]+=KgT[65];	Kgp[189]+=KgT[66];	Kgp[190]+=KgT[67];	Kgp[191]+=KgT[68];	Kgp[192]+=KgT[69];	Kgp[193]+=KgT[70];	Kgp[194]+=KgT[71];	Kgp[195]+=KgT[76];	Kgp[196]+=KgT[77];	Kgp[197]+=KgT[78];	Kgp[198]+=KgT[79];	Kgp[199]+=KgT[80];	Kgp[200]+=KgT[81];	Kgp[201]+=KgT[82];	Kgp[202]+=KgT[83];	Kgp[203]+=KgT[84];	Kgp[204]+=KgT[85];	Kgp[205]+=KgT[86];	Kgp[206]+=KgT[87];	Kgp[207]+=KgT[88];	Kgp[208]+=KgT[89];	
				Kgp[209]+=KgT[95];	Kgp[210]+=KgT[96];	Kgp[211]+=KgT[97];	Kgp[212]+=KgT[98];	Kgp[213]+=KgT[99];	Kgp[214]+=KgT[100];	Kgp[215]+=KgT[101];	Kgp[216]+=KgT[102];	Kgp[217]+=KgT[103];	Kgp[218]+=KgT[104];	Kgp[219]+=KgT[105];	Kgp[220]+=KgT[106];	Kgp[221]+=KgT[107];	Kgp[222]+=KgT[114];	Kgp[223]+=KgT[115];	Kgp[224]+=KgT[116];	Kgp[225]+=KgT[117];	Kgp[226]+=KgT[118];	Kgp[227]+=KgT[119];	Kgp[228]+=KgT[120];	Kgp[229]+=KgT[121];	Kgp[230]+=KgT[122];	Kgp[231]+=KgT[123];	Kgp[232]+=KgT[124];	Kgp[233]+=KgT[125];	Kgp[234]+=KgT[133];	Kgp[235]+=KgT[134];	Kgp[236]+=KgT[135];	Kgp[237]+=KgT[136];	Kgp[238]+=KgT[137];	Kgp[239]+=KgT[138];	Kgp[240]+=KgT[139];	Kgp[241]+=KgT[140];	Kgp[242]+=KgT[141];	Kgp[243]+=KgT[142];	Kgp[244]+=KgT[143];	Kgp[245]+=KgT[152];	Kgp[246]+=KgT[153];	Kgp[247]+=KgT[154];	Kgp[248]+=KgT[155];	
				Kgp[249]+=KgT[156];	Kgp[250]+=KgT[157];	Kgp[251]+=KgT[158];	Kgp[252]+=KgT[159];	Kgp[253]+=KgT[160];	Kgp[254]+=KgT[161];	Kgp[255]+=KgT[171];	Kgp[256]+=KgT[172];	Kgp[257]+=KgT[173];	Kgp[258]+=KgT[174];	Kgp[259]+=KgT[175];	Kgp[260]+=KgT[176];	Kgp[261]+=KgT[177];	Kgp[262]+=KgT[178];	Kgp[263]+=KgT[179];	Kgp[264]+=KgT[190];	Kgp[265]+=KgT[191];	Kgp[266]+=KgT[192];	Kgp[267]+=KgT[193];	Kgp[268]+=KgT[194];	Kgp[269]+=KgT[195];	Kgp[270]+=KgT[196];	Kgp[271]+=KgT[197];	Kgp[272]+=KgT[209];	Kgp[273]+=KgT[210];	Kgp[274]+=KgT[211];	Kgp[275]+=KgT[212];	Kgp[276]+=KgT[213];	Kgp[277]+=KgT[214];	Kgp[278]+=KgT[215];	Kgp[279]+=KgT[228];	Kgp[280]+=KgT[229];	Kgp[281]+=KgT[230];	Kgp[282]+=KgT[231];	Kgp[283]+=KgT[232];	Kgp[284]+=KgT[233];	Kgp[285]+=KgT[247];	Kgp[286]+=KgT[248];	Kgp[287]+=KgT[249];	Kgp[288]+=KgT[250];	
				Kgp[289]+=KgT[251];	Kgp[290]+=KgT[266];	Kgp[291]+=KgT[267];	Kgp[292]+=KgT[268];	Kgp[293]+=KgT[269];	Kgp[294]+=KgT[285];	Kgp[295]+=KgT[286];	Kgp[296]+=KgT[287];	Kgp[297]+=KgT[304];	Kgp[298]+=KgT[305];	Kgp[299]+=KgT[323];	

				break;
			case 3:
			
				Kgp[0]+=KgT[0];	Kgp[1]+=KgT[1];	Kgp[2]+=KgT[2];	Kgp[3]+=KgT[3];	Kgp[4]+=KgT[4];	Kgp[5]+=KgT[5];	Kgp[12]+=KgT[6];	Kgp[13]+=KgT[7];	Kgp[14]+=KgT[8];	Kgp[15]+=KgT[9];	Kgp[16]+=KgT[10];	Kgp[17]+=KgT[11];	Kgp[18]+=KgT[12];	Kgp[19]+=KgT[13];	Kgp[20]+=KgT[14];	Kgp[21]+=KgT[15];	Kgp[22]+=KgT[16];	Kgp[23]+=KgT[17];	Kgp[24]+=KgT[19];	Kgp[25]+=KgT[20];	Kgp[26]+=KgT[21];	Kgp[27]+=KgT[22];	Kgp[28]+=KgT[23];	Kgp[35]+=KgT[24];	Kgp[36]+=KgT[25];	Kgp[37]+=KgT[26];	Kgp[38]+=KgT[27];	Kgp[39]+=KgT[28];	Kgp[40]+=KgT[29];	Kgp[41]+=KgT[30];	Kgp[42]+=KgT[31];	Kgp[43]+=KgT[32];	Kgp[44]+=KgT[33];	Kgp[45]+=KgT[34];	Kgp[46]+=KgT[35];	Kgp[47]+=KgT[38];	Kgp[48]+=KgT[39];	Kgp[49]+=KgT[40];	Kgp[50]+=KgT[41];	Kgp[57]+=KgT[42];	
				Kgp[58]+=KgT[43];	Kgp[59]+=KgT[44];	Kgp[60]+=KgT[45];	Kgp[61]+=KgT[46];	Kgp[62]+=KgT[47];	Kgp[63]+=KgT[48];	Kgp[64]+=KgT[49];	Kgp[65]+=KgT[50];	Kgp[66]+=KgT[51];	Kgp[67]+=KgT[52];	Kgp[68]+=KgT[53];	Kgp[69]+=KgT[57];	Kgp[70]+=KgT[58];	Kgp[71]+=KgT[59];	Kgp[78]+=KgT[60];	Kgp[79]+=KgT[61];	Kgp[80]+=KgT[62];	Kgp[81]+=KgT[63];	Kgp[82]+=KgT[64];	Kgp[83]+=KgT[65];	Kgp[84]+=KgT[66];	Kgp[85]+=KgT[67];	Kgp[86]+=KgT[68];	Kgp[87]+=KgT[69];	Kgp[88]+=KgT[70];	Kgp[89]+=KgT[71];	Kgp[90]+=KgT[76];	Kgp[91]+=KgT[77];	Kgp[98]+=KgT[78];	Kgp[99]+=KgT[79];	Kgp[100]+=KgT[80];	Kgp[101]+=KgT[81];	Kgp[102]+=KgT[82];	Kgp[103]+=KgT[83];	Kgp[104]+=KgT[84];	Kgp[105]+=KgT[85];	Kgp[106]+=KgT[86];	Kgp[107]+=KgT[87];	Kgp[108]+=KgT[88];	Kgp[109]+=KgT[89];	
				Kgp[110]+=KgT[95];	Kgp[117]+=KgT[96];	Kgp[118]+=KgT[97];	Kgp[119]+=KgT[98];	Kgp[120]+=KgT[99];	Kgp[121]+=KgT[100];	Kgp[122]+=KgT[101];	Kgp[123]+=KgT[102];	Kgp[124]+=KgT[103];	Kgp[125]+=KgT[104];	Kgp[126]+=KgT[105];	Kgp[127]+=KgT[106];	Kgp[128]+=KgT[107];	Kgp[222]+=KgT[114];	Kgp[223]+=KgT[115];	Kgp[224]+=KgT[116];	Kgp[225]+=KgT[117];	Kgp[226]+=KgT[118];	Kgp[227]+=KgT[119];	Kgp[228]+=KgT[120];	Kgp[229]+=KgT[121];	Kgp[230]+=KgT[122];	Kgp[231]+=KgT[123];	Kgp[232]+=KgT[124];	Kgp[233]+=KgT[125];	Kgp[234]+=KgT[133];	Kgp[235]+=KgT[134];	Kgp[236]+=KgT[135];	Kgp[237]+=KgT[136];	Kgp[238]+=KgT[137];	Kgp[239]+=KgT[138];	Kgp[240]+=KgT[139];	Kgp[241]+=KgT[140];	Kgp[242]+=KgT[141];	Kgp[243]+=KgT[142];	Kgp[244]+=KgT[143];	Kgp[245]+=KgT[152];	Kgp[246]+=KgT[153];	Kgp[247]+=KgT[154];	Kgp[248]+=KgT[155];	
				Kgp[249]+=KgT[156];	Kgp[250]+=KgT[157];	Kgp[251]+=KgT[158];	Kgp[252]+=KgT[159];	Kgp[253]+=KgT[160];	Kgp[254]+=KgT[161];	Kgp[255]+=KgT[171];	Kgp[256]+=KgT[172];	Kgp[257]+=KgT[173];	Kgp[258]+=KgT[174];	Kgp[259]+=KgT[175];	Kgp[260]+=KgT[176];	Kgp[261]+=KgT[177];	Kgp[262]+=KgT[178];	Kgp[263]+=KgT[179];	Kgp[264]+=KgT[190];	Kgp[265]+=KgT[191];	Kgp[266]+=KgT[192];	Kgp[267]+=KgT[193];	Kgp[268]+=KgT[194];	Kgp[269]+=KgT[195];	Kgp[270]+=KgT[196];	Kgp[271]+=KgT[197];	Kgp[272]+=KgT[209];	Kgp[273]+=KgT[210];	Kgp[274]+=KgT[211];	Kgp[275]+=KgT[212];	Kgp[276]+=KgT[213];	Kgp[277]+=KgT[214];	Kgp[278]+=KgT[215];	Kgp[279]+=KgT[228];	Kgp[280]+=KgT[229];	Kgp[281]+=KgT[230];	Kgp[282]+=KgT[231];	Kgp[283]+=KgT[232];	Kgp[284]+=KgT[233];	Kgp[285]+=KgT[247];	Kgp[286]+=KgT[248];	Kgp[287]+=KgT[249];	Kgp[288]+=KgT[250];	
				Kgp[289]+=KgT[251];	Kgp[290]+=KgT[266];	Kgp[291]+=KgT[267];	Kgp[292]+=KgT[268];	Kgp[293]+=KgT[269];	Kgp[294]+=KgT[285];	Kgp[295]+=KgT[286];	Kgp[296]+=KgT[287];	Kgp[297]+=KgT[304];	Kgp[298]+=KgT[305];	Kgp[299]+=KgT[323];	

				break;
			default:
				FEI_COUT << " ERROR: COMPUTING ELEMENT MATRIX OF CQUAD ELEMENT FAILED" << FEI_ENDL;
				return 1;
		}
	}
	
	return 0;
}

// Get element matrix
void CQUAD::ElementMatrix(double*& stiffMat, int& size){
	size = 24;
	MatrixOperations::PackedToConventional(ElementMatrix_,stiffMat,size);
}

// Get geometric element matrix
void CQUAD::GeometricElementMatrix(double*& geomMat, int& size){
	size = 24;
	MatrixOperations::PackedToConventional(GeometricElementMatrix_,geomMat,size);
}

// Compute Temperature Load Vector
int CQUAD::TemperatureLoadVector(double* loadVector){
	// First extract the element properties and interpolation matrix
	Property** prop = Element::Properties();
	int nProp = Element::NumProperties();
	
	// Get Element delta Temperatures
	double* temp = Element::DeltaTemperatures();
	int nTemp = Element::NumDeltaTemperatures();
	
	// Interpolate the element properties and delta Temp to the center 
	// of the triangles, by creating the interpolation matrix InterpP or
	// interpT: (actually we store P**T)
	double InterpP[4*nProp];
	if(InterpolateProperties(nProp,InterpP)==1)return 1; 
	
	double* InterpT;
	if(nTemp>0){
		InterpT = new double[4*nTemp];
		if(InterpolateProperties(nTemp,InterpT)==1)return 1;
	}
	
	// Initialize vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle

	// Allocate storage
    int size = 24;
	double storage[332];  
	double* BLm = &storage[0],  *B4m = &storage[27], *B5m = &storage[54], *B6m = &storage[81];
	double* BLb = &storage[108],*B4b = &storage[135],*B5b = &storage[162],*B6b = &storage[189];
	double* Be1 = &storage[216],*Be2 = &storage[243],*Be3 = &storage[270],*Be4 = &storage[297];
	
	double* loadVec   = &storage[27];   // Temperature loadvector for a triangle in global coordinates (size 3*6=18)
	double* loadVecE  = &storage[45];   // loadVec in element coordinate system (size 3*6=18)
	double* loadVecEm = &storage[63];   // Membranal components of loadVecE (size 3*3=9)
	double* loadVecEb = &storage[72];   // Bending components of loadVecE (size 3*3=9)
	
	// Initialize some matrices and vectors used for each of the four tia elements
	double* LambdaT  = &storage[81];   // contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 9)
	double* Nodes2d  = &storage[90];   // Location of the the three nodes in element coordinates (size 6)
	double* d1= &storage[96], *d2=&storage[99], *n=&storage[102]; // relative position vector of nodes, and normal vector
	
	// Temperature vectors and area
	double *am 		= &storage[324];			// am - Temperature A vector (3x1)
	double *bm 		= &storage[327];			// bm - Temperature B vector (3x1)
	double &dtemp	=  storage[330]; 			// dTemp
	double &area    =  storage[331];  			// element area
	
	
	// Initialize LoadVector to zero
	for(std::size_t i=0;i<size;i++)loadVector[i] = 0;
	
	// Loop through four triangles
	for(std::size_t tria=0; tria<4; tria++){
		
		// Initialize storage to zero
		for(std::size_t i=0;i<332;i++)storage[i] = 0;
		
		// Get nodes corresponding to this triangle
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute the geometric characteristics of this element
		
		// Compute normal direction of element, which is e_z
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,&LambdaT[6]); // ez=  d1 x d2
		double lN = sqrt(MatrixOperations::Dot(&LambdaT[6],&LambdaT[6],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[6],&LambdaT[6], 3); // ez = 1/norm(ez);
		
		// Compute e_x = emat-(emat'*ez)*ez;
		MatrixOperations::Scale(MatrixOperations::Dot(&LambdaT[6],MaterialDirection_,3),
								0, &LambdaT[6],&LambdaT[0], 3); // ex = dot(emat,ez)*ez;
		MatrixOperations::Add(1,-1,MaterialDirection_, &LambdaT[0],&LambdaT[0],3); //ex = emat-ex;
		lN = sqrt(MatrixOperations::Dot(&LambdaT[0],&LambdaT[0],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[0],&LambdaT[0], 3); // ex = 1/norm(ex);
		
		// Compute ey = cross(ez,ex);
		MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
	
		// Compute nodes2D
		Nodes2d[0] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[1] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[2] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[3] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[4] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[5] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[3],3);
		
		
		// Store coordinate differences in Be1
		Be1[0] = Nodes2d[0]-Nodes2d[2]; Be1[1] = Nodes2d[1]-Nodes2d[3]; // Be1[0] = x12;  Be1[1] = y12 
		Be1[2] = -Be1[0]; 				Be1[3] = -Be1[1]; 				// Be1[2] = x21;  Be1[3] = y21 
		Be1[4] = Nodes2d[2]-Nodes2d[4]; Be1[5] = Nodes2d[3]-Nodes2d[5]; // Be1[4] = x23;  Be1[5] = y23 
		Be1[6] = -Be1[4]; 				Be1[7] = -Be1[5]; 				// Be1[6] = x32;  Be1[7] = y32 
		Be1[8] = Nodes2d[4]-Nodes2d[0]; Be1[9] = Nodes2d[5]-Nodes2d[1]; // Be1[8] = x31;  Be1[9] = y31 
		Be1[10]= -Be1[8]; 				Be1[11]= -Be1[9]; 				// Be1[10]= x13;  Be1[11]= y13
		Be1[12]= (Be1[3]*Be1[10]-Be1[2]*Be1[11])/2; 					// Be1[12]= area
		Be1[13]= Be1[2]*Be1[2]+Be1[3]*Be1[3];							// Be1[13]= L21*L21
		Be1[14]= Be1[6]*Be1[6]+Be1[7]*Be1[7];							// Be1[14]= L32*L32
		Be1[15]= Be1[10]*Be1[10]+Be1[11]*Be1[11];						// Be1[15]= L13*L13
			
		// Calculate membranal B matrices
		CHK_ERR(CTRIA::FelippaTriMembrane(storage,false));
		double* Bmt = &BLm[0];
		
		// Calculate bending B matrices
		CTRIA::FelippaTriBending(storage);
		double* Bbt = &BLb[0];	
		MatrixOperations::Add(1,1./3,Bbt, B4b,Bbt,27); // Bbt+= B4b/3
		MatrixOperations::Add(1,1./3,Bbt, B5b,Bbt,27); // Bbt+= B5b/3
		MatrixOperations::Add(1,1./3,Bbt, B6b,Bbt,27); // Bbt+= B6b/3
		
		// Comute am and bm (Temperature matrices)
		if(MaterialTemperatureMatrices(am,bm,dtemp,prop,nProp,InterpP,temp,nTemp,InterpT,tria)!=0)return 1;
				
		// Compute temperature load vector for membrane and bending part (to be integrated over area and multiplied by dtemp)
		MatrixOperations::dgemv9_3_1(loadVecEm,Bmt,am);	// In-plane load
		MatrixOperations::dgemv9_3_1(loadVecEb,Bbt,bm);   // Out of plane load
		
		// Assemble Membranal and Bending components
		loadVecE[0] = loadVecEm[0];  loadVecE[1] = loadVecEm[1];  loadVecE[5] = loadVecEm[2];
		loadVecE[6] = loadVecEm[3];  loadVecE[7] = loadVecEm[4];  loadVecE[11]= loadVecEm[5];
		loadVecE[12]= loadVecEm[6];  loadVecE[13]= loadVecEm[7];  loadVecE[17]= loadVecEm[8];
		
		loadVecE[2] = loadVecEb[0];  loadVecE[3] = loadVecEb[1];  loadVecE[4] = loadVecEb[2];
		loadVecE[8] = loadVecEb[3];  loadVecE[9] = loadVecEb[4];  loadVecE[10]= loadVecEb[5];
		loadVecE[14]= loadVecEb[6];  loadVecE[15]= loadVecEb[7];  loadVecE[16]= loadVecEb[8];
		
		// Compute loadVec in global coordinates: loadVec = LambdaT*loadVecE
		for(size_t index=0; index<18; index+=3){
			MatrixOperations::dgemv3_3_1(&loadVec[index],LambdaT,&loadVecE[index]);
		}
		
		// Multiply by dtemp
		MatrixOperations::Scale(dtemp, 0, loadVec,loadVec, 18); // Ftemp *= dtemp
		
		// Put this loadVec in loadVector of complete quad
		// Code has been generated by MATLAB.
		switch(tria){
			case 0:
				
				loadVector[0]+=loadVec[0];	loadVector[1]+=loadVec[1];	loadVector[2]+=loadVec[2];	loadVector[3]+=loadVec[3];	loadVector[4]+=loadVec[4];	loadVector[5]+=loadVec[5];	loadVector[6]+=loadVec[6];	loadVector[7]+=loadVec[7];	loadVector[8]+=loadVec[8];	
				loadVector[9]+=loadVec[9];	loadVector[10]+=loadVec[10];	loadVector[11]+=loadVec[11];	loadVector[18]+=loadVec[12];	loadVector[19]+=loadVec[13];	loadVector[20]+=loadVec[14];	loadVector[21]+=loadVec[15];	loadVector[22]+=loadVec[16];	loadVector[23]+=loadVec[17];	

				break;
			case 1:
			
				loadVector[0]+=loadVec[0];	loadVector[1]+=loadVec[1];	loadVector[2]+=loadVec[2];	loadVector[3]+=loadVec[3];	loadVector[4]+=loadVec[4];	loadVector[5]+=loadVec[5];	loadVector[6]+=loadVec[6];	loadVector[7]+=loadVec[7];	loadVector[8]+=loadVec[8];	
				loadVector[9]+=loadVec[9];	loadVector[10]+=loadVec[10];	loadVector[11]+=loadVec[11];	loadVector[12]+=loadVec[12];	loadVector[13]+=loadVec[13];	loadVector[14]+=loadVec[14];	loadVector[15]+=loadVec[15];	loadVector[16]+=loadVec[16];	loadVector[17]+=loadVec[17];	

				break;
			case 2:
			
				loadVector[6]+=loadVec[0];	loadVector[7]+=loadVec[1];	loadVector[8]+=loadVec[2];	loadVector[9]+=loadVec[3];	loadVector[10]+=loadVec[4];	loadVector[11]+=loadVec[5];	loadVector[12]+=loadVec[6];	loadVector[13]+=loadVec[7];	loadVector[14]+=loadVec[8];	
				loadVector[15]+=loadVec[9];	loadVector[16]+=loadVec[10];	loadVector[17]+=loadVec[11];	loadVector[18]+=loadVec[12];	loadVector[19]+=loadVec[13];	loadVector[20]+=loadVec[14];	loadVector[21]+=loadVec[15];	loadVector[22]+=loadVec[16];	loadVector[23]+=loadVec[17];	

				break;
			case 3:
			
				loadVector[0]+=loadVec[0];	loadVector[1]+=loadVec[1];	loadVector[2]+=loadVec[2];	loadVector[3]+=loadVec[3];	loadVector[4]+=loadVec[4];	loadVector[5]+=loadVec[5];	loadVector[12]+=loadVec[6];	loadVector[13]+=loadVec[7];	loadVector[14]+=loadVec[8];	
				loadVector[15]+=loadVec[9];	loadVector[16]+=loadVec[10];	loadVector[17]+=loadVec[11];	loadVector[18]+=loadVec[12];	loadVector[19]+=loadVec[13];	loadVector[20]+=loadVec[14];	loadVector[21]+=loadVec[15];	loadVector[22]+=loadVec[16];	loadVector[23]+=loadVec[17];	

				break;
			default:
				FEI_COUT << " ERROR: COMPUTING TEMPERATURE FORCE VECTOR OF CQUAD ELEMENT FAILED" << FEI_ENDL;
				return 1;
		}
	}
	
	// Multiply by element area and divide by four, since we have four triangles
	ElementArea(area);
	MatrixOperations::Scale(0.25*area, 0, loadVector,loadVector, size); // Ftemp *= 0.25*area
	
	return 0;
}

// Compute normal direction of this element
int CQUAD::NormalDirection(double *n){
	
	// Compute (average) normal direction of element
	double d1[3] = {0,0,0}, d2[3]={0,0,0};
	MatrixOperations::Add(-1,1,Nodes_[0]->GlobalCoordinates(),
	                           Nodes_[2]->GlobalCoordinates(),d1,3); // d1 = node3-node1
	MatrixOperations::Add(-1,1,Nodes_[1]->GlobalCoordinates(),
	                           Nodes_[3]->GlobalCoordinates(),d2,3); // d2 = node4-node2
	MatrixOperations::Cross(d1,d2,n); // n=  d1 x d2
	double lN = sqrt(MatrixOperations::Dot(n,n,3));
	MatrixOperations::Scale(1/lN, 0, n,n, 3); // n = 1/norm(n);
	
	return 0;
}

// Compute area of this element
int CQUAD::ElementArea(double& area){
	// Initialize vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle
	                      
    // Compute element area by adding the area of the four triangles and divide by 2
    area = 0;
    double d1[3] = {0,0,0}, d2[3]={0,0,0}, n[3]={0,0,0};
	for(std::size_t tria=0; tria<4; tria++){
		
		// Get nodes corresponding to this triangle
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute area using cross product of two vertices
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,n); // n=  d1 x d2
		area += sqrt(MatrixOperations::Dot(n,n,3))/4;		// Divide by four as each triangle area is divided by 2
		
	}
	
	return 0;
}

// Interpolate element properties to gauss points
int CQUAD::InterpolateProperties(const int nProp, double* interpMat){
	
	// Write the interpolation matrix depending on the number of properties
	switch(nProp){
		case 0:
			FEI_COUT << " ERROR: EACH ELEMENT MUST HAVE MORE THAN" 
			<< " ZERO PROPERTIES" << FEI_ENDL;
			return 1;
		case 1:
			interpMat[0] = 1;  		// N1(-1/3,-1/3)
			interpMat[1] = 1; 		// N1( 1/3,-1/3)
			interpMat[2] = 1;  		// N1( 1/3, 1/3)
			interpMat[3] = 1; 		// N1(-1/3, 1/3)
			break;
		case 4:
			interpMat[0] = 4./9;  		// N1(-1/3,-1/3)
			interpMat[1] = 2./9; 		// N2(-1/3,-1/3)
			interpMat[2] = 1./9;  		// N3(-1/3,-1/3)
			interpMat[3] = 2./9; 		// N4(-1/3,-1/3)
			
			interpMat[4] = 2./9;  		// N1( 1/3,-1/3)
			interpMat[5] = 4./9; 		// N2( 1/3,-1/3)
			interpMat[6] = 2./9;  		// N3( 1/3,-1/3)
			interpMat[7] = 1./9; 		// N4( 1/3,-1/3)
			
			interpMat[8] = 1./9;  		// N1( 1/3, 1/3)
			interpMat[9] = 2./9; 		// N2( 1/3, 1/3)
			interpMat[10]= 4./9;  		// N3( 1/3, 1/3)
			interpMat[11]= 2./9; 		// N4( 1/3, 1/3)
			
			interpMat[12]= 2./9;  		// N1(-1/3, 1/3)
			interpMat[13]= 1./9; 		// N2(-1/3, 1/3)
			interpMat[14]= 2./9;  		// N3(-1/3, 1/3)
			interpMat[15]= 4./9; 		// N4(-1/3, 1/3)
			break;
		default:
			FEI_COUT << " ERROR: THE CQUAD ELEMENT ONLY SUPPORTS ONE OR" 
			<< " FOUR PROPERTIES PER ELEMENT" << FEI_ENDL;
			return 1;
	}
	
	return 0;
	
}

// Read CQUAD from file
int CQUAD::ReadFromFile(std::ifstream& fin, Domain& domain){
	char line[DomainConstants::MAXLINESIZE];
	
	// Read first line and parse
	fin.getline(line,DomainConstants::MAXLINESIZE);
	std::stringstream parse(line);
	
	// The line consists of CTRIA, elemID, node1,node2,node3,node4 and material direction
	char elementCard[DomainConstants::MAXNAMESIZE];
	int elemID,node1,node2,node3,node4;
	
	if(parse>>elementCard>> elemID>>node1>>node2>>node3>>node4
			>>MaterialDirection_[0]>>MaterialDirection_[1]>>MaterialDirection_[2]);
	else{
		FEI_COUT << " ERROR: COULD NOT READ ELEMENT DATA " << FEI_ENDL;
		return 1;		
	}
	
	// Check elementCard
	if(strcmp(elementCard,"CQUAD")!=0){
		FEI_COUT << " ERROR: ELEMENT CARD: " << elementCard <<
		" IS EXPECTED TO BE: CQUAD" << FEI_ENDL;
		return 1;		
	}
	
	// Add the elementID to the Local Global conversion
	int elemIDLocal = domain.ElementID_LG.size();
	domain.ElementID_LG.push_back(elemID);
	domain.ElementID_GL.insert(
		std::pair<const int,const int>(
			elemID, elemIDLocal));
			
	// Set local element id
	Element::SetLocalElementID(elemIDLocal);
	
	// Check if node1, node2, node3 and node4 are defined, and set them in Nodes_
	if(domain.NodeID_GL.find(node1) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node1 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	} else if(domain.NodeID_GL.find(node2) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node2 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	}else if(domain.NodeID_GL.find(node3) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node3 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	}else if(domain.NodeID_GL.find(node4) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node4 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	}
	Nodes_[0] = domain.NodeList[domain.NodeID_GL[node1]];
	Nodes_[1] = domain.NodeList[domain.NodeID_GL[node2]];
	Nodes_[2] = domain.NodeList[domain.NodeID_GL[node3]];
	Nodes_[3] = domain.NodeList[domain.NodeID_GL[node4]];
	
	// Set # of dof on these nodes to 6:
	domain.NodalNumDOF[domain.NodeID_GL[node1]] = 6;
	domain.NodalNumDOF[domain.NodeID_GL[node2]] = 6;
	domain.NodalNumDOF[domain.NodeID_GL[node3]] = 6;
	domain.NodalNumDOF[domain.NodeID_GL[node4]] = 6;
	
	// Normalize MaterialDirection to unit length
	double lD = sqrt(MatrixOperations::Dot(MaterialDirection_,MaterialDirection_,3));
	MatrixOperations::Scale(1/lD, 0, MaterialDirection_,MaterialDirection_, 3);
	
	// Read property lines of this element
	if(ReadPropertyLines(fin,domain,Property::Shell)!=0) return 1;
	
	// Verify if the MaterialDirection of all properies is not normal to
	// this element
	if(VerifyMaterialDirection(elemID)!=0)return 1;
	
	return 0;
}

// Verify if the MaterialDirection of all properies is not normal to
// this element
int CQUAD::VerifyMaterialDirection(int GlobalElemID){	
	
	// Initialized vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle
	
	// Loop through four triangles
	double d1[3] = {0,0,0}, d2[3]={0,0,0}, n[3]={0,0,0}; // Initialization of d1,d2,n
	for(std::size_t tria=0; tria<4; tria++){
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute normal direction of element
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,n); // n= d1 x d2
		double lN = sqrt(MatrixOperations::Dot(n,n,3));
		MatrixOperations::Scale(1/lN, 0, n,n, 3); // n = 1/norm(n);
		
		// Compute angle theta between normal direction and material direction (which may not be too small)
		double theta;
		
		theta = acos(std::fabs(MatrixOperations::Dot(n,MaterialDirection_,3)));
		if(theta<DomainConstants::ERRORANGLE){
			FEI_COUT<< "ERROR: THE MATERIAL DIRECTION OF ELEMENT " <<
			 GlobalElemID << " IS TOO CLOSE TO THE NORMAL DIRECTION TO THIS ELEMENT."
					<< FEI_ENDL;
			return 1;
		}else if(theta<DomainConstants::WARNINGANGLE){
			FEI_COUT<< "WARNING: THE MATERIAL DIRECTION OF ELEMENT " <<
			 GlobalElemID << " IS CLOSE TO THE NORMAL DIRECTION TO THIS ELEMENT.."
					<< FEI_ENDL;
		}		
	}
	
	
	return 0;
}

// Write to outputfile
int CQUAD::WriteOutput(Domain& domain, Solver *solver,
		                     std::ofstream* tempFiles, bool doTitle){
	
	OutputRequest& outReq = *domain.OutputReq;
	
	// Check whether element output needs to be printed
	if(!outReq.PrintElementOutput_) return 0;	
	
	// Check whether element strains forces, or strain energy need to be printed
	if(outReq.PrintElementStrains_ || outReq.PrintElementForces_ || outReq.PrintElementStrainEnergy_){
		// Introduce elemStrains for strains at center of four triangles and the average
		double elemStrains[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		
		// Compute element strains at center of element
		if(ElementStrains(elemStrains,solver)!=0) return 1;
		
		if(outReq.PrintElementStrains_){
			if(doTitle){
				tempFiles[0] << FEI_ENDL << "CQUAD element strains in element coordinates:" << FEI_ENDL << std::left
				 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
				 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "TID"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   EX"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   EY"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   GXY"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   KXX"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   KYY"
				 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   KXY" << FEI_ENDL;
			}
			// Loop through four triangles and print strains to file
			int width = DomainConstants::OUTPUT_INT_FIELD_WIDTH;
			
			// First triangle
			tempFiles[0] << std::left << std::setw(width) << domain.ElementID_LG[Element::LocalElementID()]
				 << std::left << std::setw(width) << 1;
				 UtilityFunctions::PrintMatrix(&elemStrains[0],1,6,tempFiles[0]);
			// Second - fourth triangles
			for(std::size_t tria=1; tria<4; tria++){
				tempFiles[0] << std::left << std::setw(width) << " " << std::left << std::setw(width) << tria+1;
				UtilityFunctions::PrintMatrix(&elemStrains[tria*6],1,6,tempFiles[0]);
			}
			// Print average strains
			tempFiles[0] << std::left << std::setw(width) << " " << std::left << std::setw(width) << "AVG";
			UtilityFunctions::PrintMatrix(&elemStrains[24],1,6,tempFiles[0]);
			
		}
		
		// Check whether element forces or strain energy need to be printed
		if(outReq.PrintElementForces_ || outReq.PrintElementStrainEnergy_){
			
			// Compute element forces at center of element
			double elemForces[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			if(ElementForces(elemForces,elemStrains)!=0) return 1;
			
			// Print element forces if requested
			if(outReq.PrintElementForces_){
				if(doTitle){
					tempFiles[1] << FEI_ENDL << "CQUAD element forces in element coordinates:" << FEI_ENDL << std::left
					 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
					 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "TID"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   NX"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   NY"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   NXY"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MXX"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MYY"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MXY" << FEI_ENDL;
				}
				
				
				// Loop through four triangles and print element forces to file
				int width = DomainConstants::OUTPUT_INT_FIELD_WIDTH;
				
				// First triangle
				tempFiles[1] << std::left << std::setw(width) << domain.ElementID_LG[Element::LocalElementID()]
					 << std::left << std::setw(width) << 1;
					 UtilityFunctions::PrintMatrix(&elemForces[0],1,6,tempFiles[1]);
				// Second - fourth triangles
				for(std::size_t tria=1; tria<4; tria++){
					tempFiles[1] << std::left << std::setw(width) << " " << std::left << std::setw(width) << tria+1;
					UtilityFunctions::PrintMatrix(&elemForces[tria*6],1,6,tempFiles[1]);
				}
				// Print average element forces
				tempFiles[1] << std::left << std::setw(width) << " " << std::left << std::setw(width) << "AVG";
				UtilityFunctions::PrintMatrix(&elemForces[24],1,6,tempFiles[1]);
			}
			
			// Print element strain energy if requested
			if(outReq.PrintElementStrainEnergy_){
				if(doTitle){
					tempFiles[2] << FEI_ENDL << "CQUAD element strain energy:" << FEI_ENDL << std::left
					 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
					 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   Energy" << FEI_ENDL;
				}
				
				
				// Compute element strain energy at center of element
				double strainEnergy = 0;
				if(ElementStrainEnergy(strainEnergy,elemForces,elemStrains)!=0) return 1;
				
				// Print the Element strain energies to the file
				int width = DomainConstants::OUTPUT_INT_FIELD_WIDTH;
				tempFiles[2] << std::left << std::setw(width) << domain.ElementID_LG[Element::LocalElementID()];
				UtilityFunctions::PrintMatrix(&strainEnergy,1,1,tempFiles[2]);
			}
		}
		
	}
	
				
	return 0;
}

// Compute material stiffness matrices integrated over element from properties
int CQUAD::MaterialMatrices(double* Am,double* Bm,double* Dm, Property** properties, const int& nProp, double* InterpP, const int& triaID){
	// Create space for packed storage matrices
	double* Amp = new double[6];
	double* Bmp = new double[6];
	double* Dmp = new double[6];
	
	// Initialize Material matrices to 0
	for(size_t i=0; i<6; i++){ Amp[i]=0; Bmp[i]=0; Dmp[i]=0;}
	
	// Loop through properties
	PSHELL* propI;
	double temp1=0, temp2=0;  // Two temporary variables for GP1 and GP2
	for(std::size_t i=0; i<nProp; i++){
		propI = static_cast<PSHELL*>(properties[i]);
		// Fill Amp, Bmp and Dmp respectively
		for(int k=0; k<6; k++){
			
			Amp[k]+= InterpP[triaID*nProp+i]*propI->A()[k]; // Contribution of k^th property to i^th entry of Amp of triaID^th triangle 
			Bmp[k]+= InterpP[triaID*nProp+i]*propI->B()[k]; // Contribution of k^th property to i^th entry of Bmp of triaID^th triangle 
			Dmp[k]+= InterpP[triaID*nProp+i]*propI->D()[k]; // Contribution of k^th property to i^th entry of Dmp of triaID^th triangle 
		
		}		
	}
	
	// Write in conventional storage
	MatrixOperations::PackedToConventional(Amp, Am, 3);
	MatrixOperations::PackedToConventional(Bmp, Bm, 3);
	MatrixOperations::PackedToConventional(Dmp, Dm, 3);

	return 0;
	
}

// Compute material stiffness matrices integrated over element from properties
int CQUAD::MaterialTemperatureMatrices(double* am,double* bm,double &dTemp, Property** properties, const int& nProp, double* InterpP, double* temp, const int& nTemp, double* InterpT, const int& triaID){
	
	// Initialize Material matrices and dTemp o 0
	for(size_t i=0; i<3; i++){ am[i]=0; bm[i]=0;}
	dTemp = 0;
	
	// Loop through properties
	PSHELL* propI;
	for(size_t i=0; i<nProp; i++){
		propI = static_cast<PSHELL*>(properties[i]);
		
		// Add the contribution of this property to am, bm
		for(size_t k=0; k<3; k++){
			am[k] += InterpP[triaID*nProp+i]*propI->a()[k]; // Contribution of i^th property to k^th entry of am of triaID^th triangle 
			bm[k] += InterpP[triaID*nProp+i]*propI->b()[k];	// Contribution of i^th property to k^th entry of bm of triaID^th triangle 
		}			
	}
	
	// Interpolate the temperature to the triangle centers
	if(nTemp>0){
		// Loop through temperatures
		for(size_t i=0; i<nTemp; i++) dTemp += InterpT[triaID*nTemp+i]*temp[i];	// Contribution of i^th property to dTemp of triaID^th triangle 
	}

	return 0;
	
}
		
// Compute element strains at center of four triangles
int CQUAD::ElementStrains(double* elemStrains, Solver *solver){
		
	// Initialize vector containing the nodes of the four triangles
	Node* triaNodes[12]={ Nodes_[0], Nodes_[1], Nodes_[3],		// First triangle
	                      Nodes_[0], Nodes_[1], Nodes_[2],		// Second triangle
	                      Nodes_[1], Nodes_[2], Nodes_[3],		// Third triangle
	                      Nodes_[0], Nodes_[2], Nodes_[3]};		// Fourth triangle

	// Allocate storage
    int size = 24;
	double storage[576];  // 567, same as tria;
	double* BLm = &storage[0],  *B4m = &storage[27], *B5m = &storage[54], *B6m = &storage[81];
	double* BLb = &storage[108],*B4b = &storage[135],*B5b = &storage[162],*B6b = &storage[189];
	double* Be1 = &storage[216],*Be2 = &storage[243],*Be3 = &storage[270],*Be4 = &storage[297];
	double* Kmm = &storage[324],*Kbb = &storage[405],*Kbm = &storage[486];
	
	double* dispVec   = &storage[27];   // Displacement Vector for a triangle (size 3*6=18)
	double* dispVecE  = &storage[45];   // Displacement vector in element coordinate system (size 3*6=18)
	double* dispVecEm = &storage[63];   // Membranal components of dispVecE (size 3*3=9)
	double* dispVecEb = &storage[72];   // Bending components of dispVecE (size 3*3=9)
	
	// Initialize some matrices and vectors used for each of the four tia elements
	double* LambdaT  = &storage[81];   // contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 9)
	double* Nodes2d  = &storage[90];   // Location of the the three nodes in element coordinates (size 6)
	double* d1= &storage[96], *d2=&storage[99], *n=&storage[102]; // relative position vector of nodes, and normal vector
	
	
	// Loop through four triangles
	for(std::size_t tria=0; tria<4; tria++){
		// Initialize storage to zero
		for(std::size_t i=0;i<567;i++)storage[i] = 0;
		
		// Get nodes corresponding to this triangle
		Node** nodesT = &triaNodes[3*tria];
		
		// Compute the geometric characteristics of this element
		
		// Compute normal direction of element, which is e_z
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[1]->GlobalCoordinates(),d1,3); // d1 = node2-node1
		MatrixOperations::Add(-1,1,nodesT[0]->GlobalCoordinates(),
								   nodesT[2]->GlobalCoordinates(),d2,3); // d2 = node3-node1
		MatrixOperations::Cross(d1,d2,&LambdaT[6]); // ez=  d1 x d2
		double lN = sqrt(MatrixOperations::Dot(&LambdaT[6],&LambdaT[6],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[6],&LambdaT[6], 3); // ez = 1/norm(ez);
		
		// Compute e_x = emat-(emat'*ez)*ez;
		MatrixOperations::Scale(MatrixOperations::Dot(&LambdaT[6],MaterialDirection_,3),
								0, &LambdaT[6],&LambdaT[0], 3); // ex = dot(emat,ez)*ez;
		MatrixOperations::Add(1,-1,MaterialDirection_, &LambdaT[0],&LambdaT[0],3); //ex = emat-ex;
		lN = sqrt(MatrixOperations::Dot(&LambdaT[0],&LambdaT[0],3));
		MatrixOperations::Scale(1/lN, 0, &LambdaT[0],&LambdaT[0], 3); // ex = 1/norm(ex);
		
		// Compute ey = cross(ez,ex);
		MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
		
		// Get Nodal displacement vector for this triangle
		nodesT[0]->Displacements(solver, &dispVec[0]);
		nodesT[1]->Displacements(solver, &dispVec[DomainConstants::DISP_FIELD_SIZE]);
		nodesT[2]->Displacements(solver, &dispVec[2*DomainConstants::DISP_FIELD_SIZE]);
		
		// Compute nodal displacement vector in element coordinates uE = Lambda*u, so uE**T = u**T*LambdaT
		for(size_t index=0; index<18; index+=3){
			MatrixOperations::dgemv1_3_3(&dispVecE[index],&dispVec[index],LambdaT);
		}
		
		// Get Membranal and Bending components
		dispVecEm[0] = dispVecE[0];  dispVecEm[1] = dispVecE[1];  dispVecEm[2] = dispVecE[5];
		dispVecEm[3] = dispVecE[6];  dispVecEm[4] = dispVecE[7];  dispVecEm[5] = dispVecE[11];
		dispVecEm[6] = dispVecE[12]; dispVecEm[7] = dispVecE[13]; dispVecEm[8] = dispVecE[17];
		
		dispVecEb[0] = dispVecE[2];  dispVecEb[1] = dispVecE[3];  dispVecEb[2] = dispVecE[4];
		dispVecEb[3] = dispVecE[8];  dispVecEb[4] = dispVecE[9];  dispVecEb[5] = dispVecE[10];
		dispVecEb[6] = dispVecE[14]; dispVecEb[7] = dispVecE[15]; dispVecEb[8] = dispVecE[16];
		
	
		// Compute nodes2D
		Nodes2d[0] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[1] = MatrixOperations::Dot(nodesT[0]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[2] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[3] = MatrixOperations::Dot(nodesT[1]->GlobalCoordinates(),&LambdaT[3],3);
		Nodes2d[4] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[0],3);
		Nodes2d[5] = MatrixOperations::Dot(nodesT[2]->GlobalCoordinates(),&LambdaT[3],3);
		
		
		// Store coordinate differences in Be1
		Be1[0] = Nodes2d[0]-Nodes2d[2]; Be1[1] = Nodes2d[1]-Nodes2d[3]; // Be1[0] = x12;  Be1[1] = y12 
		Be1[2] = -Be1[0]; 				Be1[3] = -Be1[1]; 				// Be1[2] = x21;  Be1[3] = y21 
		Be1[4] = Nodes2d[2]-Nodes2d[4]; Be1[5] = Nodes2d[3]-Nodes2d[5]; // Be1[4] = x23;  Be1[5] = y23 
		Be1[6] = -Be1[4]; 				Be1[7] = -Be1[5]; 				// Be1[6] = x32;  Be1[7] = y32 
		Be1[8] = Nodes2d[4]-Nodes2d[0]; Be1[9] = Nodes2d[5]-Nodes2d[1]; // Be1[8] = x31;  Be1[9] = y31 
		Be1[10]= -Be1[8]; 				Be1[11]= -Be1[9]; 				// Be1[10]= x13;  Be1[11]= y13
		Be1[12]= (Be1[3]*Be1[10]-Be1[2]*Be1[11])/2; 					// Be1[12]= area
		Be1[13]= Be1[2]*Be1[2]+Be1[3]*Be1[3];							// Be1[13]= L21*L21
		Be1[14]= Be1[6]*Be1[6]+Be1[7]*Be1[7];							// Be1[14]= L32*L32
		Be1[15]= Be1[10]*Be1[10]+Be1[11]*Be1[11];						// Be1[15]= L13*L13
			
		// Calculate membranal B matrices
		CHK_ERR(CTRIA::FelippaTriMembrane(storage,false));
		double* Bmt = &BLm[0];
		
		// Calculate bending B matrices
		CTRIA::FelippaTriBending(storage);
		double* Bbt = &BLb[0];	
		MatrixOperations::Add(1,1./3,Bbt, B4b,Bbt,27); // Bbt+= B4b/3
		MatrixOperations::Add(1,1./3,Bbt, B5b,Bbt,27); // Bbt+= B5b/3
		MatrixOperations::Add(1,1./3,Bbt, B6b,Bbt,27); // Bbt+= B6b/3
		
		// Compute strains for this triangle
		MatrixOperations::dgemv1_9_3(&elemStrains[tria*6],dispVecEm,Bmt);		// In-plane strains
		MatrixOperations::dgemv1_9_3(&elemStrains[tria*6+3],dispVecEb,Bbt);  // Curvatures
		
		// Add contribution of this triangle to average strains
		MatrixOperations::Add(1,0.25,&elemStrains[24], &elemStrains[tria*6],&elemStrains[24],6); // strain_avg+= 1/4 *strainT
		
	}
	
	return 0;
}

// Compute element forces at center of four triangles
int CQUAD::ElementForces(double* elemForces, const double* elemStrains){
	// First extract the element properties and interpolation matrix
	Property** prop = Element::Properties();
	int nProp = Element::NumProperties();
	
	// Get Element delta Temperatures
	double* temp = Element::DeltaTemperatures();
	int nTemp = Element::NumDeltaTemperatures();
	
	// Interpolate the element properties and delta temperatures to the 
	// center of the triangles, by creating the interpolation matrices
	// InterpP and InterpT: (actually we store IntperP**T and InterpT**T)
	double InterpP[4*nProp];
	if(InterpolateProperties(nProp,InterpP)==1)return 1; 
	
	double* InterpT;
	if(nTemp>0){
		InterpT = new double[4*nTemp];
		if(InterpolateProperties(nTemp,InterpT)==1)return 1;
	}

	// Get space for material Am,Bm,Dm matrices in conventional storage
	double *Am = new double[9], *Bm = new double[9], *Dm = new double[9];
	double *am = new double[3], *bm = new double[3], dTemp=0;
	
	// Loop through triangles
	for(std::size_t tria=0;tria<4;tria++){
		// Compute interpolated Am, Bm, Dm, am, bm, dTemp
		if(MaterialMatrices(Am,Bm,Dm,prop,nProp,InterpP,tria)!=0) return 1; 
		if(MaterialTemperatureMatrices(am,bm,dTemp,prop,nProp,InterpP,temp,nTemp,InterpT,tria)!=0)return 1;
		
		// Compute Normal loads
		MatrixOperations::dgemv(1, &elemForces[tria*6], 1, Am, &elemStrains[tria*6], 3, 3);  // N+= A*eps
		MatrixOperations::dgemv(1, &elemForces[tria*6], 1, Bm, &elemStrains[tria*6+3], 3, 3);  // N+= B*kapp
		MatrixOperations::Add(1,-dTemp,&elemForces[tria*6],am,&elemForces[tria*6],3); 	   // N-= a*dTemp
		
		// Compute Moment loads
		MatrixOperations::dgemv(1, &elemForces[tria*6+3], 1, Bm, &elemStrains[tria*6], 3, 3);  // M+= B*eps
		MatrixOperations::dgemv(1, &elemForces[tria*6+3], 1, Dm, &elemStrains[tria*6+3], 3, 3);  // M+= D*kapp
		MatrixOperations::Add(1,-dTemp,&elemForces[tria*6+3],bm,&elemForces[tria*6+3],3); 	   // M-= b*dTemp
		
		// Add contribution of this triangle to average elementForces
		MatrixOperations::Add(1,0.25,&elemForces[24], &elemForces[tria*6],&elemForces[24],6); // forc_avg+= 1/4 *forceT
	}
	
	return 0;
}

// Compute element strain energy at center of triangle
int CQUAD::ElementStrainEnergy(double& elemEnergy, double* elemForces, const double* elemStrains){
	// Create storage to store cvec=[a;b], dTemp and area
	double storage[71];
		
	double* am   = &storage[0];   		// a vector (size 3)
	double* bm   = &storage[3];   		// b vector (size 3)
	double* cvec = &storage[0];   		// cvec vector = [a;b]
	double &dTemp=  storage[6]; 		// dTemp
	double &area =  storage[7];			// area
	
	double* Am   = &storage[8];   		// A matrix (size 9)
	double* Bm   = &storage[17];   		// B matrix (size 9)
	double* Dm   = &storage[26];     	// D matrix (size 9)
	double* Cmat = &storage[35];		// Cmat matrix = [A B; B D] (size 36)
	
	// Extract the element properties and interpolation matrix
	Property** prop = Element::Properties();
	int nProp = Element::NumProperties();
	
	// Get Element delta Temperatures
	double* temp = Element::DeltaTemperatures();
	int nTemp = Element::NumDeltaTemperatures();
	
	// Interpolate the element properties and delta temperatures to the 
	// center of the triangles, by creating the interpolation matrices
	// InterpP and InterpT: (actually we store IntperP**T and InterpT**T)
	double InterpP[4*nProp];
	if(InterpolateProperties(nProp,InterpP)==1)return 1; 
	
	double* InterpT;
	if(nTemp>0){
		InterpT = new double[4*nTemp];
		if(InterpolateProperties(nTemp,InterpT)==1)return 1;
	}
	
	
	// We compute strain energy based on average strain, average element 
	// forces and total elemene area
	elemEnergy = 0.;			// Initialize to zero
	
	// Loop through triangles
	for(std::size_t tria=0;tria<4;tria++){
		// Initialize storage to zero
		for(std::size_t i=0;i<71;i++) storage[i]=0;
		
		// The first contribution to the energy is the dot product of elemForces and elemStrains
		elemEnergy+= MatrixOperations::Dot(&elemForces[tria*6],&elemStrains[tria*6],6);
		
		// Compute am,bm,dTemp
		if(MaterialTemperatureMatrices(am,bm,dTemp,prop,nProp,InterpP,temp,nTemp,InterpT,tria)!=0)return 1;
		
		// Continue if dTemp!=0
		if(dTemp!=0){
			// Second contribution: -strain**T*cvec*dTemp
			elemEnergy-= MatrixOperations::Dot(&elemStrains[tria*6],cvec,6)*dTemp;
			
			// Compute material matrices
			if(MaterialMatrices(Am,Bm,Dm,prop,nProp,InterpP,tria)!=0) return 1; 
			
			// Fill Cmat (only lower triangular part)
			Cmat[0] =Am[0];  Cmat[1] =Am[1];  Cmat[2] =Am[2];  Cmat[3] =Bm[0];  Cmat[4] =Bm[1];  Cmat[5] =Bm[2];
							 Cmat[7] =Am[4];  Cmat[8] =Am[5];  Cmat[9] =Bm[3];  Cmat[10]=Bm[4];  Cmat[11]=Bm[5];
											  Cmat[14]=Am[8];  Cmat[15]=Bm[6];  Cmat[16]=Bm[7];  Cmat[17]=Bm[8];
			
															   Cmat[21]=Dm[0];  Cmat[22]=Dm[1];  Cmat[23]=Dm[2];
																				Cmat[28]=Dm[4];  Cmat[29]=Dm[5];
																								 Cmat[35]=Dm[8];
			// Compute Cholesky decomposition L, Cmat:=L
			int info;
			MatrixOperations::dpotrf(Cmat,6,&info);
			if(info!=0){ FEI_COUT << "ERROR: LAPACK CHOLESKY DECOMPOSITION FAILED WITH ERROR: " 
			                     << info << ".";
			             return 1;}
			
			// Compute inverse of L, Cmat:= inv(L)
			MatrixOperations::dtrtri6(Cmat,Cmat);
			
			// Compute inv(L)*cvec, store in first column of Cmat, Cmat:=inv(L)*cvec
			MatrixOperations::dtrmv6(Cmat, cvec, Cmat);
			
			// Third contribution: (inv(L)*cvec)**T*(inv(L)*cvec)*dTemp^2 = cvec**T*inv(Cmat)*cvec*dTemp^2
			elemEnergy+= MatrixOperations::Dot(Cmat,Cmat,6)*dTemp*dTemp;
		}
	}
	
	
	if(ElementArea(area)!=0)return 1; 
	
	elemEnergy*=area*0.125;	// Integrate over element area and divide by 4 (since phi = 0.5*int() and the sum of the four tria's is twice the quad)
	
	return 0;
}


