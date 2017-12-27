#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "fei_base.hpp"
#include <cmath>

// Constructors:
CBEAM::CBEAM(): Element(Beam){}
CBEAM::~CBEAM(){}

// Build stiffness matrix of the element
int CBEAM::BuildElementMatrix(){
	// Allocate storage
    int size = 12;
	double storage[256];
	for(int i=0;i<256;i++)storage[i] = 0;
	double* Ce1 = &storage[0],  *Ce2 = &storage[16], *Ce3 = &storage[32], *Ce4 = &storage[48];
	double* Ce5 = &storage[52],*Ge1 = &storage[56],*Ge2 = &storage[64],*Ge3 = &storage[72];
	double* Se = &storage[74];
	double* Be1 = &storage[78],*Ke = &storage[102];
	
	double* LambdaT  = &storage[246];   // contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 9)
	double  &L       =  storage[255]; 		// Length of element
	
	// Put vector from node 1 to node 2 into LambdaT, and scale with length
	MatrixOperations::Add(-1,1,Nodes_[0]->GlobalCoordinates(),
	                           Nodes_[1]->GlobalCoordinates(),LambdaT,3);
	L = sqrt(MatrixOperations::Dot(LambdaT,LambdaT,3));   
	MatrixOperations::Scale(1/L, 0, LambdaT,LambdaT, 3);
	
	// Set z direction to the user defined one minus the component in the
	// e_x direction (in order to have them perpendicular)
	MatrixOperations::Add(1,1,&LambdaT[6], ZDirection_,&LambdaT[6],3);
	MatrixOperations::Add(-1*MatrixOperations::Dot(LambdaT,&LambdaT[6],3),
							1,LambdaT,&LambdaT[6],&LambdaT[3],3); // ez = ez-dot(ex,ez)*ex
	MatrixOperations::Scale(1/sqrt(MatrixOperations::Dot(&LambdaT[3],&LambdaT[3],3)),
	                        0, &LambdaT[3],&LambdaT[6], 3); // ez = ez/norm(ez);
	                        
	// Compute e_y = cross(e_z,e_x);
	MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
	
	// Fill material matrices by integrating material stiffnesses over
	// element using the gauss points:
	// Ce1=C00, Ce2=C01, Ce3=C11, Ce4=c0, Ce5=c1 Ge1=G0, Ge2=G1, Ge3=g0, Se=S
	if(MaterialMatrices(Ce1,Ce2,Ce3,Ge1,Ge2,Se,Ce4,Ce5,Ge3)!=0)return 1;
	
	// Do cholesky factorization on Ce3==C11
	MatrixOperations::dpotrf4(Ce3,Ce3); // Ce3:= L, with C11 = L*L**T
	
	// Compute inv(L) ==inv(Ce3)
	MatrixOperations::dtrtri4(Ce3,Ce3);  // Ce3:= inv(L)
	
	// Compute C01p=inv(L)*C01 == Ce3*Ce2 and store in Ce2
	MatrixOperations::dtrmm4(Ce3, Ce2, Ce2); // Ce2:=inv(L)*C01 (==CO1p)
	
	// Compute C00-C01p**T*C01p == Ce1-Ce2**T*Ce2 and store in Ce1
	MatrixOperations::dsyrd4(Ce1, Ce2, Ce1); // Ce1:=C00-C01*inv(C11)*C01
	
	// Compute B3-2/L*G1, by scaling Ge2 and adding B3 manually
	MatrixOperations::Scale(-2/L, 0, Ge2,Ge2, 8);
	Ge2[5]+=1./3; Ge2[6]-=1./3;			// Ge2:=B3-2/L*G1;
	
	// Compute Gp=(B3-2/L*G1)*inv(L)**T == Ge2*Ce3**T, and store in Ge2
	MatrixOperations::dtrmm2_4t(Ge2, Ce3, Ge2); // Ge2 := Gp;
	
	// Compute S+L^2/4*Gp*Gp**T == Se+L^2/4*Ge2*Ge2**T and store in Se
	MatrixOperations::dsyr2_4(Se, L*L/4, Ge2, Se); // Se:= inv(Cs)
	
	// Compute Ls s.t. Ls*Ls**T = inv(Cs) and store in Se
	MatrixOperations::dpotrf2(Se, Se); // Se:= Ls, s.t. Ls*Ls**T = inv(Cs)
	MatrixOperations::dtrtri2(Se, Se); // Se:= inv(Ls);
	
	// Compute Gp*C01p == (B3-2/L*G1)*inv(C11)*C01 == Ge2*Ce2,
	// and store in Ce3 (only first half of Ce3 is occupied)
	MatrixOperations::dgemm2_4(Ge2, Ce2, Ce3); //Ce3:=Gp*C01p
	
	// Compute Gstar = 0.5*Gp*C01p+G0/L == 0.5*Ce3+Ge1/L, store in Ge1
	// C:= alpha*A + beta*B, with A,B,C(pxq) with pxq=M
	MatrixOperations::Add(0.5, 1/L, Ce3,Ge1, Ge1, 8); // Ge1:=Gstar
	
	// Compute Bstar = Gstar*B0+B2 == Be1, using hard programming
	Be1[0] =-Ge1[0];	Be1[1] =-Ge1[1];		Be1[2] =-1/L;
						Be1[5] =-1/L;			Be1[6] =-Ge1[2];		Be1[7] =-Ge1[3], 
	Be1[8] =-Ge1[4];	Be1[9] =-Ge1[5]+0.5;	Be1[10]=-Ge1[6]-0.5;	Be1[11]=-Ge1[7],
	Be1[12]= Ge1[0];	Be1[13]= Ge1[1];		Be1[14]= 1/L; 
						Be1[17]= 1/L;			Be1[18]= Ge1[2];		Be1[19]= Ge1[3];
	Be1[20]= Ge1[4];	Be1[21]= Ge1[5]+0.5;	Be1[22]= Ge1[6]-0.5;	Be1[23]= Ge1[7];
	 
	 // Compute Be1 = inv(Ls)*Bstar == Se*Be1;
     MatrixOperations::dtrmm2_4L(Se, Be1, Be1); //Be1:=inv(Ls)*Bstar
     MatrixOperations::dtrmm2_4L(Se, &Be1[8], &Be1[8]); //Be1:=inv(Ls)*Bstar
     MatrixOperations::dtrmm2_4L(Se, &Be1[16], &Be1[16]); //Be1:=inv(Ls)*Bstar
     
     // Compute part of Ke = B0**T*(C00-C01*inv(C11)*C01)*B0 using
     // hard coding
     Ke[0] =  Ce1[0]; Ke[3]= Ce1[1]; Ke[4] = Ce1[2]; Ke[5] = Ce1[3];
     Ke[6] = -Ce1[0]; Ke[9]=-Ce1[1]; Ke[10]=-Ce1[2]; Ke[11]=-Ce1[3];
     
     Ke[36] =  Ce1[4]; Ke[39]= Ce1[5]; Ke[40] = Ce1[6]; Ke[41] = Ce1[7];
     Ke[42] = -Ce1[4]; Ke[45]=-Ce1[5]; Ke[46] =-Ce1[6]; Ke[47] =-Ce1[7];
     
     Ke[48]= Ce1[8]; Ke[51]= Ce1[9]; Ke[52] = Ce1[10]; Ke[53] = Ce1[11];
     Ke[54]=-Ce1[8]; Ke[57]=-Ce1[9]; Ke[58] =-Ce1[10]; Ke[59] =-Ce1[11];
     
     Ke[60]= Ce1[12]; Ke[63]= Ce1[13]; Ke[64]= Ce1[14]; Ke[65]= Ce1[15];
     Ke[66]=-Ce1[12]; Ke[69]=-Ce1[13]; Ke[70]=-Ce1[14]; Ke[71]=-Ce1[15];
     
     Ke[72]= -Ce1[0]; Ke[75]=-Ce1[1]; Ke[76] =-Ce1[2]; Ke[77] =-Ce1[3];
     Ke[78]=  Ce1[0]; Ke[81]= Ce1[1]; Ke[82] = Ce1[2]; Ke[83] = Ce1[3];
     
     Ke[108]= -Ce1[4]; Ke[111]=-Ce1[5]; Ke[112]=-Ce1[6]; Ke[113]=-Ce1[7];
     Ke[114]=  Ce1[4]; Ke[117]= Ce1[5]; Ke[118]= Ce1[6]; Ke[119]= Ce1[7];
     
     Ke[120]=-Ce1[8]; Ke[123]=-Ce1[9]; Ke[124]=-Ce1[10]; Ke[125]=-Ce1[11];
     Ke[126]= Ce1[8]; Ke[129]= Ce1[9]; Ke[130]= Ce1[10]; Ke[131]= Ce1[11];
     
     Ke[132]=-Ce1[12]; Ke[135]=-Ce1[13]; Ke[136]=-Ce1[14]; Ke[137]=-Ce1[15];
     Ke[138]= Ce1[12]; Ke[141]= Ce1[13]; Ke[142]= Ce1[14]; Ke[143]= Ce1[15];
     
     // Update Ke:= 1/L*Ke + L*(inv(Ls)*Bstar)**T*(inv(Ls)*Bstar)
     // == 1/L*Ke + L* Be1**T*Be1
     MatrixOperations::dsyrk_nt(1/L, Ke, L, Be1,2, 12);
     
     // Consider bar orientation: Ke:= LambdaT*Ke*LambdaT**T
     MatrixOperations::dsyrk_f(Ke, LambdaT, Ke, size, 3);
    
	// Finally, store Ke in packed storage
	MatrixOperations::ConventionalToPacked(Ke,ElementMatrix_,size);
	
	//FEI_COUT << " Stiffness matrix: " << FEI_ENDL; UtilityFunctions::PrintMatrix(Ke,size,size);
    
	return 0;
}

// Build geometric stiffness matrix of the element
int CBEAM::BuildGeometricElementMatrix(Solver *solver, bool makeNegative){
    int size = 12;
    
    // Set Kg to be identity matrix
    double *Kg =  new double[size*size];
    for(std::size_t i=0;i<size*size;i++) Kg[i]=0;
    for(std::size_t i=0;i<size;i++) Kg[i*size+i]=1;
    // Make Kg its negative if so required (used for differential stiffness of prestress)
    if(makeNegative) for(std::size_t i=0;i<size*size;i++) Kg[i]*=-1;
	MatrixOperations::ConventionalToPacked(Kg,GeometricElementMatrix_,size);
	
	return 0;
}

// Get element matrix
void CBEAM::ElementMatrix(double*& stiffMat, int& size){
	size = 12;
	MatrixOperations::PackedToConventional(ElementMatrix_,stiffMat,size);
}

// Get geometric element matrix
void CBEAM::GeometricElementMatrix(double*& geomMat, int& size){
	size = 12;
	MatrixOperations::PackedToConventional(GeometricElementMatrix_,geomMat,size);
}

// Compute Temperature Load Vector
int CBEAM::TemperatureLoadVector(double* loadVector){
	int size = 12;
	
	// First check if temperature is applied
	if(Element::NumDeltaTemperatures()==0){
		for(std::size_t i=0;i<size;i++) loadVector[i]=0;
	}
	
	// Allocate storage
    double storage[172];
	for(int i=0;i<172;i++)storage[i] = 0;
	double* Ce1 = &storage[0],  *Ce2 = &storage[16], *Ce3 = &storage[32], *Ce4 = &storage[48];
	double* Ce5 = &storage[52],*Ge1 = &storage[56],*Ge2 = &storage[64],*Ge3 = &storage[72];
	double* Se = &storage[74];
	double* Be1 = &storage[78],*Be2 = &storage[126];
	
	double* LambdaT      = &storage[162];   // contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 9)
	double  &L       =  storage[171]; 		// Length of element
	
	// Compute the geometric characteristics of this beam element
		
	// Put vector from node 1 to node 2 into LambdaT, and scale with length
	MatrixOperations::Add(-1,1,Nodes_[0]->GlobalCoordinates(),
	                           Nodes_[1]->GlobalCoordinates(),LambdaT,3);
	L = sqrt(MatrixOperations::Dot(LambdaT,LambdaT,3));   
	MatrixOperations::Scale(1/L, 0, LambdaT,LambdaT, 3);
	
	// Set z direction to the user defined one minus the component in the
	// e_x direction (in order to have them perpendicular)
	MatrixOperations::Add(1,1,&LambdaT[6], ZDirection_,&LambdaT[6],3);
	MatrixOperations::Add(-1*MatrixOperations::Dot(LambdaT,&LambdaT[6],3),
							1,LambdaT,&LambdaT[6],&LambdaT[3],3); // ez = ez-dot(ex,ez)*ex
	MatrixOperations::Scale(1/sqrt(MatrixOperations::Dot(&LambdaT[3],&LambdaT[3],3)),
	                        0, &LambdaT[3],&LambdaT[6], 3); // ez = ez/norm(ez);
	                        
	// Compute e_y = cross(e_z,e_x);
	MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
	
	// Fill material matrices by integrating material stiffnesses over
	// element using the gauss points:
	// Ce1=C00, Ce2=C01, Ce3=C11, Ce4=c0, Ce5=c1 Ge1=G0, Ge2=G1, Ge3=g0, Se=S
	if(MaterialMatrices(Ce1,Ce2,Ce3,Ge1,Ge2,Se,Ce4,Ce5,Ge3)!=0)return 1;
	
	// Do cholesky factorization on Ce3==C11
	MatrixOperations::dpotrf4(Ce3,Ce3); // Ce3:= L, with C11 = L*L**T
	
	// Compute inv(L) ==inv(Ce3)
	MatrixOperations::dtrtri4(Ce3,Ce3);  // Ce3:= inv(L)
	
	// Compute C01p=inv(L)*C01 == Ce3*Ce2 and store in Ce2
	MatrixOperations::dtrmm4(Ce3, Ce2, Ce2); // Ce2:=inv(L)*C01 (==CO1p)
	
	// Compute c1p=inv(L)*c1 ==Ce3*Ce5 and store in Ce5
	MatrixOperations::dtrmv4(Ce3, Ce5, Ce5); // Ce5:=inv(L)*c1 (==c1p)
	
	// Compute B3-2/L*G1, by scaling Ge2 and adding B3 manually
	MatrixOperations::Scale(-2/L, 0, Ge2,Ge2, 8);
	Ge2[5]+=1./3; Ge2[6]-=1./3;			// Ge2:=B3-2/L*G1;
	
	// Compute Gp=(B3-2/L*G1)*inv(L)**T == Ge2*Ce3**T, and store in Ge2
	MatrixOperations::dtrmm2_4t(Ge2, Ce3, Ge2); // Ge2 := Gp;
	
	// Compute S+L^2/4*Gp*Gp**T == Se+L^2/4*Ge2*Ge2**T and store in Se
	MatrixOperations::dsyr2_4(Se, L*L/4, Ge2, Se); // Se:= inv(Cs)
	
	// Compute Ls s.t. Ls*Ls**T = inv(Cs) and store in Se
	MatrixOperations::dpotrf2(Se, Se); // Se:= Ls, s.t. Ls*Ls**T = inv(Cs)
	MatrixOperations::dtrtri2(Se, Se); // Se:= inv(Ls);
	
	// Compute Gp*C01p == (B3-2/L*G1)*inv(C11)*C01 == Ge2*Ce2,
	// and store in Ce3 (only first half of Ce3 is occupied)
	MatrixOperations::dgemm2_4(Ge2, Ce2, Ce3); //Ce3:=Gp*C01p
	
	// Compute C01p**T*c1p == C01*inv(C11)*c1 == Ce2**T*Ce5
	// and store in Ce2 (only first column of Ce2 is occupied)
	MatrixOperations::dgemv4_4_1t(Ce2, Ce2, Ce5); //Ce2(1:4):=Ce2**T*Ce5
	
	// Compute Gp*c1p == (B3-2/L*G1)*inv(C11)*c1 == Ge2*Ce5,
	// and store in Ge2 (only first column of Ge2 is occupied)
	MatrixOperations::dgemv2_4_1(Ge2, Ge2, Ce5); //Ge2:=Gp*c1p
	
	// Compute Gstar = 0.5*Gp*C01p+G0/L == 0.5*Ce3+Ge1/L, store in Ge1
	MatrixOperations::Add(0.5, 1/L, Ce3,Ge1, Ge1, 8); // Ge1:=Gstar
	
	// Compute gstar = 0.5*Gp*c1p+g0/L == 0.5*Ge2+Ge3/L, store in Ge3
	MatrixOperations::Add(0.5, 1/L, Ge2,Ge3, Ge3, 2); // Ge3:=gstar
	
	// Compute Ge2 = inv(Ls)*Gstar = Se*Ge1
    MatrixOperations::dtrmm2_4L(Se, Ge1, Ge2); //Ge2 := inv(Ls)*Gstar
    
    // Compute Ge3 = inv(Ls)*gstar = Se*Ge3
    MatrixOperations::dtrmv2(Se, Ge3, Ge3); //Ge3 := inv(Ls)*gstar
    
    // Compute Ce1 = 1/L*Ce1 + L * Ge2**T*Ge2
    MatrixOperations::dsyrk_nt(1/L, Ce1, L, Ge2,2, 4); //Ce1:= 1/L*Ce1+L*Ge2**T*Ge2;
    
    // Compute Ce5 = Ge2**T*Ge3 == Gstar**T*Cs*gstar
    MatrixOperations::dgemv2_4_1t(Ce5, Ge2, Ge3);     //Ce5:=Gstar**T*Cs*gstar 
    
    // Compute Ge2 = inv(Ls)**T*inv(Ls)*Gstar == Cs*Gstar == Se**T*Ge2
    MatrixOperations::dtrmm2_4tL(Se,Ge2,Ge2);
    
    // Compute Ge3 = inv(Ls)**T*inv(Ls)*gstar == Cs*gstar == Se**T*Ge3
    MatrixOperations::dtrmv2t(Se,Ge3,Ge3);
    
    // Compute Ce4 = -c0 + C01*inv(C11)*c1 = -Ce4+Ce2
    MatrixOperations::Add(-1,1,Ce4,Ce2,Ce4,4);   // Ce4:= -c0 + C01*inv(C11)*c1
    
    // Compute Ge3 = -L*Cs*gstar = -L*Ge3
    MatrixOperations::Scale(-L,0,Ge3,Ge3,2);      	// Ge3:=x2
    
    // Compute Ce5 = -L^2*Gstar**T*Cs*gstar+(-c0+C01*inv(C11)*c1) = -L^2*Ce5+Ce4 == x3
    MatrixOperations::Add(-L*L,1,Ce5,Ce4,Ce5,4);	// Ce5:=x3
    
    // Compute loadVector in element coordinates and store in 
    // Be1 = -(B0**T*x3+L*B2**T*x2). Use hardcoding
    Be1[0] = Ce5[0]; 				Be1[1] = Ge3[0]; 
    Be1[2] = Ge3[1]; 				Be1[3] = Ce5[1];
    Be1[4] = Ce5[2]-Ge3[1]*0.5*L; 	Be1[5] = Ce5[3]+Ge3[0]*0.5*L;
    Be1[6] =-Ce5[0]; 				Be1[7] =-Ge3[0]; 
    Be1[8] =-Ge3[1]; 				Be1[9] =-Ce5[1];
    Be1[10]=-Ce5[2]-Ge3[1]*0.5*L; 	Be1[11]=-Ce5[3]+Ge3[0]*0.5*L;
    
    // Compute loadVector in global coordinates ad store in: loadVector = LambdaT*loadVecE
	for(size_t index=0; index<size; index+=3){
		MatrixOperations::dgemv3_3_1(&loadVector[index],LambdaT,&Be1[index]);
	}	
	
	return 0;
}

// Interpolate element properties to gauss points
int CBEAM::InterpolateProperties(const int nProp, double* interpMat){
	double xi = 1/sqrt(3);
	// Write the interpolation matrix depending on the number of properties
	switch(nProp){
		case 0:
			FEI_COUT << " ERROR: EACH ELEMENT MUST HAVE MORE THAN" 
			<< " ZERO PROPERTIES" << FEI_ENDL;
			return 1;
		case 1:
			interpMat[0] = 1;  		// N1(xi1)
			interpMat[1] = 1; 		// N1(xi2)
			break;
		case 2:
			interpMat[0] = 0.5*(1+xi);  		// N1(xi1)
			interpMat[1] = 0.5*(1-xi);  		// N2(xi1)
			
			interpMat[2] = 0.5*(1-xi);  		// N1(xi2)
			interpMat[3] = 0.5*(1+xi);  		// N2(xi2)
			break;
		case 3:
			interpMat[0] = 0.5*(xi*xi+xi);  		// N1(xi1)
			interpMat[1] = 1-xi*xi;					// N2(xi1)
			interpMat[2] = 0.5*(xi*xi-xi);  		// N3(xi1)
			
			interpMat[3] = 0.5*(xi*xi-xi);  		// N1(xi2)
			interpMat[4] = 1-xi*xi;					// N2(xi2)
			interpMat[5] = 0.5*(xi*xi+xi);  		// N3(xi2)
					
			break;
		
		case 4:
			interpMat[0] = (9*xi*xi*xi+9*xi*xi-xi-1)/16;  		// N1(xi1)
			interpMat[1] = (-27*xi*xi*xi-9*xi*xi+27*xi+9)/16;	// N2(xi1)
			interpMat[2] = (27*xi*xi*xi-9*xi*xi-27*xi+9)/16;  	// N3(xi1)
			interpMat[3] = (-9*xi*xi*xi+9*xi*xi+xi-1)/16;  		// N4(xi1)
			
			interpMat[4] = (-9*xi*xi*xi+9*xi*xi+xi-1)/16;  		// N1(xi2)
			interpMat[5] = (27*xi*xi*xi-9*xi*xi-27*xi+9)/16;	// N2(xi2)
			interpMat[6] = (-27*xi*xi*xi-9*xi*xi+27*xi+9)/16;  	// N3(xi2)
			interpMat[7] = (9*xi*xi*xi+9*xi*xi-xi-1)/16;  		// N4(xi2)
			break;
		default:
			FEI_COUT << " ERROR: THE CBEAM ELEMENT SUPPORTS UP TO ONLY" 
			<< " FOUR PROPERTIES PER ELEMENT" << FEI_ENDL;
			return 1;
	}
	return 0;
	
}

// Read CBEAM from file
int CBEAM::ReadFromFile(std::ifstream& fin, Domain& domain){

	char line[DomainConstants::MAXLINESIZE];
	
	// Read first line and parse
	fin.getline(line,DomainConstants::MAXLINESIZE);
	std::stringstream parse(line);
	
	// The line consists of CBEAM, elemID, node1,node2, zX,zY,zZ
	char elementCard[DomainConstants::MAXNAMESIZE];
	int elemID,node1,node2;
	double zX,zY,zZ;
	
	if(parse>>elementCard>> elemID>>node1>>node2>>zX>>zY>>zZ);
	else{
		FEI_COUT << " ERROR: COULD NOT READ ELEMENT DATA " << FEI_ENDL;
		return 1;		
	}
	
	// Check elementCard
	if(strcmp(elementCard,"CBEAM")!=0){
		FEI_COUT << " ERROR: ELEMENT CARD: " << elementCard <<
		" IS EXPECTED TO BE: CBEAM" << FEI_ENDL;
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
	
	// Check if node1 and node2 are defined, and set them in Nodes_
	if(domain.NodeID_GL.find(node1) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node1 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	} else if(domain.NodeID_GL.find(node2) == domain.NodeID_GL.end()) {
	   FEI_COUT << " ERROR: NODE " << node2 << " IN ELEMENT " << elemID 
	            << " UNDEFINED" << FEI_ENDL;
		return 1;
	}
	Nodes_[0] = domain.NodeList[domain.NodeID_GL[node1]];
	Nodes_[1] = domain.NodeList[domain.NodeID_GL[node2]];
	
	// Set # of dof on these nodes to 6:
	domain.NodalNumDOF[domain.NodeID_GL[node1]] = 6;
	domain.NodalNumDOF[domain.NodeID_GL[node2]] = 6;
	
	
	// Set ZDirection and normalize to unit lenght
	ZDirection_[0] = zX;
	ZDirection_[1] = zY;
	ZDirection_[2] = zZ;
	double lZ = sqrt(MatrixOperations::Dot(ZDirection_,ZDirection_,3));
	MatrixOperations::Scale(1/lZ, 0, ZDirection_,ZDirection_, 3);
	
	// Verify if this ZDirection_ is not parallel to the Xdir
	double xDir[] = {0,0,0};
	MatrixOperations::Add(-1,1,Nodes_[0]->GlobalCoordinates(),
	                           Nodes_[1]->GlobalCoordinates(),xDir,3);
	MatrixOperations::Scale(1/sqrt(MatrixOperations::Dot(xDir,xDir,3)), 
	                  0, xDir,xDir, 3);
	                  
	double theta = acos(std::fabs(MatrixOperations::Dot(xDir,ZDirection_,3)));
	if(theta<DomainConstants::ERRORANGLE){
		FEI_COUT<< "ERROR: THE NORMAL DIRECTION OF ELEMENT " <<
		 elemID << " IS TOO CLOSE TO ITS LONGITUDINAL DIRECTION:"
			    << FEI_ENDL;
		return 1;
	}else if(theta<DomainConstants::WARNINGANGLE){
		FEI_COUT<< "WARNING: THE NORMAL DIRECTION OF ELEMENT " <<
		 elemID << " IS CLOSE TO ITS LONGITUDINAL DIRECTION "
			    << FEI_ENDL;
	}
	
	
	// Read property lines of this element
	if(ReadPropertyLines(fin,domain,Property::Beam)!=0) return 1;

	return 0;
}

// Write to outputfile
int CBEAM::WriteOutput(Domain& domain, Solver *solver,
		                     std::ofstream* tempFiles, bool doTitle){
								 
	OutputRequest& outReq = *domain.OutputReq;
	
	// Check whether element output needs to be printed
	if(!outReq.PrintElementOutput_) return 0;	
	
	// If no return, then continue:
	if(outReq.PrintElementForces_){		
		if(doTitle){
			tempFiles[0] << FEI_ENDL << "CBEAM element forces in element coordinates:" << FEI_ENDL << std::left
			 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
			 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "GP"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   NX"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   VY"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   VZ"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MX"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MY"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   MZ" << FEI_ENDL;
		}
		
		// Introduce ElementForce vector for both Gauss Points and average (i.e. 18x1)
		double elemForces[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		
		// Compute element forces
		if(ElementForces(elemForces,solver)!=0) return 1;
		
		
		// Print the Element forces to the file
		int width = DomainConstants::OUTPUT_INT_FIELD_WIDTH;
		tempFiles[0] << std::left << std::setw(width) << domain.ElementID_LG[Element::LocalElementID()]
			 << std::left << std::setw(width) << 1;
		UtilityFunctions::PrintMatrix(&elemForces[0],1,6,tempFiles[0]);
		tempFiles[0] << std::left << std::setw(width) << " "
			 << std::left << std::setw(width) << 2;
		UtilityFunctions::PrintMatrix(&elemForces[6],1,6,tempFiles[0]);
		tempFiles[0] << std::left << std::setw(width) << " "
			 << std::left << std::setw(width) << "AVG";
		UtilityFunctions::PrintMatrix(&elemForces[12],1,6,tempFiles[0]);
	}
	if(outReq.PrintElementStrainEnergy_){
		if(doTitle){
			tempFiles[1] << FEI_ENDL << "CBEAM element strain energy:" << FEI_ENDL << std::left
			 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   Energy" << FEI_ENDL;
		}
		
		// Compute element strain energy of element
		double strainEnergy = 0;
		if(ElementStrainEnergy(strainEnergy,solver)!=0) return 1;
		
		// Print the Element strain energies to the file
		int width = DomainConstants::OUTPUT_INT_FIELD_WIDTH;
		tempFiles[1] << std::left << std::setw(width) << domain.ElementID_LG[Element::LocalElementID()];
		UtilityFunctions::PrintMatrix(&strainEnergy,1,1,tempFiles[1]);
		
	}	
						
	return 0;
}

// Compute material stiffness matrices integrated over element from properties
int CBEAM::MaterialMatrices(double* C00,double* C01,double* C11,double* G0,double* G1,double* S,double* c0, double* c1, double* g0){
	// Get Element Properties
	Property** prop = Element::Properties();
	int nProp = Element::NumProperties();
	
	// Get Element delta Temperatures
	double* dtemp = Element::DeltaTemperatures();
	int nTemp = Element::NumDeltaTemperatures();
	
	// Interpolate the element properties to the gauss points, by
	// creating the interpolation matrix P: (actually we store P**T)
	double InterpP[2*nProp];
	if(InterpolateProperties(nProp,InterpP)==1)return 1;
	
	// Loop through properties
	PBEAM* propI;
	double temp1=0, temp2=0;  // Two temporary variables for GP1 and GP2
	for(int i=0; i<nProp; i++){
		propI = static_cast<PBEAM*>(prop[i]);
		// Fill Ce1, Ce2, and Ce3 with C00, C01 and C11 respectively
		for(int k=0; k<16; k++){
			// The contribution of this property to GP1 and GP2:
			temp1 = InterpP[i]*propI->BendingStiffness()[k];			// GP1
			temp2 = InterpP[nProp+i]*propI->BendingStiffness()[k];	// GP2
			
			// Add the contribution of this property to Ce1, Ce2, and 
			// Ce3, which are int(p dxi), int(p*xi dxi) and 
			//int(p*xi^2 dxi). Use numerical integration with GP1 and 
			// GP2 with w1=w2=1, and xi1=-1/sqrt(3)=-xi2
			C00[k] += 0.5*(temp1+temp2);				
			C01[k] += 0.5*(-temp1+temp2)/sqrt(3);		
			C11[k] += 0.5*(temp1+temp2)/3;				
		}
		// Fill Ge1, Ge2 with G0 and G1 respectively
		for(int k=0; k<8; k++){
			// The contribution of this property to GP1 and GP2:
			temp1 = InterpP[i]*propI->CouplingStiffness()[k];		// GP1
			temp2 = InterpP[nProp+i]*propI->CouplingStiffness()[k];	// GP2
			
			// Add the contribution of this property to Ce1
			G0[k] += 0.5*(temp1+temp2);				
			G1[k] += 0.5*(-temp1+temp2)/sqrt(3);		
		}
		// Fill Se with S
		for(int k=0; k<4; k++){
			// The contribution of this property to GP1 and GP2:
			temp1 = InterpP[i]*propI->ShearStiffness()[k];			// GP1
			temp2 = InterpP[nProp+i]*propI->ShearStiffness()[k];	// GP2
			
			// Add the contribution of this property to Ce1
			S[k] += 0.5*(temp1+temp2);				
		}				
	}
	
	if(nTemp>0){
		double cGP1[4] = {0,0,0,0};
		double cGP2[4] = {0,0,0,0};
		double gGP1[2] = {0,0};
		double gGP2[2] = {0,0};
		double dTGP1 = 0, dTGP2=0;
		
		
		// Loop through properties
		for(int i=0; i<nProp; i++){
			propI = static_cast<PBEAM*>(prop[i]);
			// Fill Ce4, and Ce5 with c0 and c1 respectively
			for(int k=0; k<4; k++){
				// The contribution of this property to GP1 and GP2:
				temp1 = InterpP[i]*propI->BendingExpansion()[k];			// GP1
				temp2 = InterpP[nProp+i]*propI->BendingExpansion()[k];	// GP2
				
				cGP1[k] += temp1;
				cGP2[k] += temp2;			
			}
			// Fill Ge3 with g0
			for(int k=0; k<2; k++){
				// The contribution of this property to GP1 and GP2:
				temp1 = InterpP[i]*propI->ShearExpansion()[k];		 	// GP1
				temp2 = InterpP[nProp+i]*propI->ShearExpansion()[k];	// GP2
				
				gGP1[k] += temp1;
				gGP2[k] += temp2;							
			}				
		}
	
		// Interpolate the element dTemp
		double InterpT[2*nTemp];
		if(InterpolateProperties(nTemp,InterpT)==1)return 1;
		
		// Loop through temperatures
		for(int i=0; i<nTemp; i++){		
			dTGP1 += InterpT[i]*dtemp[i];
			dTGP2 += InterpT[nTemp+i]*dtemp[i];				
		}
		
		MatrixOperations::Add(0.5*dTGP1,0.5*dTGP2,cGP1,cGP2,c0,4);
		MatrixOperations::Add(-0.5*dTGP1/sqrt(3),0.5*dTGP2/sqrt(3),cGP1,cGP2,c1,4);
		MatrixOperations::Add(0.5*dTGP1,0.5*dTGP2,gGP1,gGP2,g0,2);
		
	}
	
	return 0;
	
}

// Compute element forces at gauss points
int CBEAM::ElementForces(double * elemForces, Solver *solver){
	// Allocate storage
    int size = 12;
	double storage[184];
	for(int i=0;i<184;i++)storage[i] = 0;
	double* Ce1 = &storage[0],  *Ce2 = &storage[16], *Ce3 = &storage[32], *Ce4 = &storage[48];
	double* Ce5 = &storage[52],*Ge1 = &storage[56],*Ge2 = &storage[64],*Ge3 = &storage[72];
	double* Se = &storage[74];
	double* Be1 = &storage[78],*Be2 = &storage[126];
	
	double* dispVecE = &storage[150];  // Nodal displacement vector in element coordinates (size 12)
	double* dispVec  = &storage[162]; 	// Nodal displacement vector in global coordinates (size 12)
	double* LambdaT  = &storage[174];   // contains ex,ey,ez, i.e. LambdaT == Lambda**T (size 9)
	double  &L       =  storage[183]; 		// Length of element
	
	// Compute the geometric characteristics of this beam element
		
	// Put vector from node 1 to node 2 into LambdaT, and scale with length
	MatrixOperations::Add(-1,1,Nodes_[0]->GlobalCoordinates(),
	                           Nodes_[1]->GlobalCoordinates(),LambdaT,3);
	L = sqrt(MatrixOperations::Dot(LambdaT,LambdaT,3));   
	MatrixOperations::Scale(1/L, 0, LambdaT,LambdaT, 3);
	
	// Set z direction to the user defined one minus the component in the
	// e_x direction (in order to have them perpendicular)
	MatrixOperations::Add(1,1,&LambdaT[6], ZDirection_,&LambdaT[6],3);
	MatrixOperations::Add(-1*MatrixOperations::Dot(LambdaT,&LambdaT[6],3),
							1,LambdaT,&LambdaT[6],&LambdaT[3],3); // ez = ez-dot(ex,ez)*ex
	MatrixOperations::Scale(1/sqrt(MatrixOperations::Dot(&LambdaT[3],&LambdaT[3],3)),
	                        0, &LambdaT[3],&LambdaT[6], 3); // ez = ez/norm(ez);
	                        
	// Compute e_y = cross(e_z,e_x);
	MatrixOperations::Cross(&LambdaT[6],&LambdaT[0], &LambdaT[3]);
	
	// Get Nodal displacement vector
	Nodes_[0]->Displacements(solver, &dispVec[0]);
	Nodes_[1]->Displacements(solver, &dispVec[DomainConstants::DISP_FIELD_SIZE]);
	
	// Compute nodal displacement vector in element coordinates uE = Lambda*u, so uE**T = u**T*LambdaT
	for(size_t index=0; index<size; index+=3){
		MatrixOperations::dgemv1_3_3(&dispVecE[index],&dispVec[index],LambdaT);
	}
	
	// Fill material matrices by integrating material stiffnesses over
	// element using the gauss points:
	// Ce1=C00, Ce2=C01, Ce3=C11, Ce4=c0, Ce5=c1 Ge1=G0, Ge2=G1, Ge3=g0, Se=S
	if(MaterialMatrices(Ce1,Ce2,Ce3,Ge1,Ge2,Se,Ce4,Ce5,Ge3)!=0)return 1;
	
	// Do cholesky factorization on Ce3==C11
	MatrixOperations::dpotrf4(Ce3,Ce3); // Ce3:= L, with C11 = L*L**T
	
	// Compute inv(L) ==inv(Ce3)
	MatrixOperations::dtrtri4(Ce3,Ce3);  // Ce3:= inv(L)
	
	// Compute C01p=inv(L)*C01 == Ce3*Ce2 and store in Ce2
	MatrixOperations::dtrmm4(Ce3, Ce2, Ce2); // Ce2:=inv(L)*C01 (==CO1p)
	
	// Compute c1p=inv(L)*c1 ==Ce3*Ce5 and store in Ce5
	MatrixOperations::dtrmv4(Ce3, Ce5, Ce5); // Ce5:=inv(L)*c1 (==c1p)
	
	// Compute C00-C01p**T*C01p == Ce1-Ce2**T*Ce2 and store in Ce1
	MatrixOperations::dsyrd4(Ce1, Ce2, Ce1); // Ce1:=C00-C01*inv(C11)*C01
	
	// Compute B3-2/L*G1, by scaling Ge2 and adding B3 manually
	MatrixOperations::Scale(-2/L, 0, Ge2,Ge2, 8);
	Ge2[5]+=1./3; Ge2[6]-=1./3;			// Ge2:=B3-2/L*G1;
	
	// Compute Gp=(B3-2/L*G1)*inv(L)**T == Ge2*Ce3**T, and store in Ge2
	MatrixOperations::dtrmm2_4t(Ge2, Ce3, Ge2); // Ge2 := Gp;
	
	// Compute S+L^2/4*Gp*Gp**T == Se+L^2/4*Ge2*Ge2**T and store in Se
	MatrixOperations::dsyr2_4(Se, L*L/4, Ge2, Se); // Se:= inv(Cs)
	
	// Compute Ls s.t. Ls*Ls**T = inv(Cs) and store in Se
	MatrixOperations::dpotrf2(Se, Se); // Se:= Ls, s.t. Ls*Ls**T = inv(Cs)
	MatrixOperations::dtrtri2(Se, Se); // Se:= inv(Ls);
	
	// Compute Gp*C01p == (B3-2/L*G1)*inv(C11)*C01 == Ge2*Ce2,
	// and store in Ce3 (only first half of Ce3 is occupied)
	MatrixOperations::dgemm2_4(Ge2, Ce2, Ce3); //Ce3:=Gp*C01p
	
	// Compute C01p**T*c1p == C01*inv(C11)*c1 == Ce2**T*Ce5
	// and store in Ce2 (only first column of Ce2 is occupied)
	MatrixOperations::dgemv4_4_1t(Ce2, Ce2, Ce5); //Ce2(1:4):=Ce2**T*Ce5
	
	// Compute Gp*c1p == (B3-2/L*G1)*inv(C11)*c1 == Ge2*Ce5,
	// and store in Ge2 (only first column of Ge2 is occupied)
	MatrixOperations::dgemv2_4_1(Ge2, Ge2, Ce5); //Ge2:=Gp*c1p
	
	// Compute Gstar = 0.5*Gp*C01p+G0/L == 0.5*Ce3+Ge1/L, store in Ge1
	// C:= alpha*A + beta*B, with A,B,C(pxq) with pxq=M
	MatrixOperations::Add(0.5, 1/L, Ce3,Ge1, Ge1, 8); // Ge1:=Gstar
	
	// Compute gstar = 0.5*Gp*c1p+g0/L == 0.5*Ge2+Ge3/L, store in Ge3
	MatrixOperations::Add(0.5, 1/L, Ge2,Ge3, Ge3, 2); // Ge3:=gstar
	
	// Compute Ge2 = inv(Ls)*Gstar = Se*Ge1
    MatrixOperations::dtrmm2_4L(Se, Ge1, Ge2); //Ge2 := inv(Ls)*Gstar
    
    // Compute Ge3 = inv(Ls)*gstar = Se*Ge3
    MatrixOperations::dtrmv2(Se, Ge3, Ge3); //Ge3 := inv(Ls)*gstar
    
    // Compute Ce1 = 1/L*Ce1 + L * Ge2**T*Ge2
    MatrixOperations::dsyrk_nt(1/L, Ce1, L, Ge2,2, 4); //Ce1:= 1/L*Ce1+L*Ge2**T*Ge2;
    
    // Compute Ce5 = Ge2**T*Ge3 == Gstar**T*Cs*gstar
    MatrixOperations::dgemv2_4_1t(Ce5, Ge2, Ge3);     //Ce5:=Gstar**T*Cs*gstar 
    
    // Compute Ge2 = inv(Ls)**T*inv(Ls)*Gstar == Cs*Gstar == Se**T*Ge2
    MatrixOperations::dtrmm2_4tL(Se,Ge2,Ge2);
    
    // Compute Ge3 = inv(Ls)**T*inv(Ls)*gstar == Cs*gstar == Se**T*Ge3
    MatrixOperations::dtrmv2t(Se,Ge3,Ge3);
    
    // Compute Ce4 = -c0 + C01*inv(C11)*c1 = -Ce4+Ce2
    MatrixOperations::Add(-1,1,Ce4,Ce2,Ce4,4);   // Ce4:= -c0 + C01*inv(C11)*c1
    
    // Compute Ge3 = -L*Cs*gstar = -L*Ge3
    MatrixOperations::Scale(-L,0,Ge3,Ge3,2);      	// Ge3:=x2
    
    // Compute Ce5 = -L^2*Gstar**T*Cs*gstar+(-c0+C01*inv(C11)*c1) = -L^2*Ce5+Ce4 == x3
    MatrixOperations::Add(-L*L,1,Ce5,Ce4,Ce5,4);	// Ce5:=x3
    
    // Compute Se = inv(Ls)**T*inv(Ls) == Cs
    MatrixOperations::dtrmtrm2(Se, Se); // Se :=Cs
    
    // Compute Be1 = L*(Cs*Gstar)**T*B2 by hard coding    
    Be1[4] =-Ge2[0];		Be1[5] =-Ge2[2];		Be1[6] =-Ge2[4];		Be1[7] =-Ge2[6];
	Be1[8] =-Ge2[1];		Be1[9] =-Ge2[3];		Be1[10]=-Ge2[5];		Be1[11]=-Ge2[7];
	Be1[16]= Ge2[1]*L/2;	Be1[17]= Ge2[3]*L/2;	Be1[18]= Ge2[5]*L/2;	Be1[19]= Ge2[7]*L/2;
	Be1[20]=-Ge2[0]*L/2;	Be1[21]=-Ge2[2]*L/2;	Be1[22]=-Ge2[4]*L/2;	Be1[23]=-Ge2[6]*L/2;
	Be1[28]= Ge2[0];		Be1[29]= Ge2[2];		Be1[30]= Ge2[4];		Be1[31]= Ge2[6];
	Be1[32]= Ge2[1];		Be1[33]= Ge2[3];		Be1[34]= Ge2[5];		Be1[35]= Ge2[7];
	Be1[40]= Ge2[1]*L/2;	Be1[41]= Ge2[3]*L/2;	Be1[42]= Ge2[5]*L/2;	Be1[43]= Ge2[7]*L/2;
	Be1[44]=-Ge2[0]*L/2;	Be1[45]=-Ge2[2]*L/2;	Be1[46]=-Ge2[4]*L/2;	Be1[47]=-Ge2[6]*L/2;
   				  
	// Update Be1 = Be1 + Ce1*B0 by hard coding
	Be1[0] -= Ce1[0];  Be1[1] -= Ce1[1];  Be1[2] -= Ce1[2];  Be1[3] -= Ce1[3];
	Be1[12]-= Ce1[4];  Be1[13]-= Ce1[5];  Be1[14]-= Ce1[6];  Be1[15]-= Ce1[7];
    Be1[16]-= Ce1[8];  Be1[17]-= Ce1[9];  Be1[18]-= Ce1[10]; Be1[19]-= Ce1[11];
    Be1[20]-= Ce1[12]; Be1[21]-= Ce1[13]; Be1[22]-= Ce1[14]; Be1[23]-= Ce1[15];
    Be1[24]+= Ce1[0];  Be1[25]+= Ce1[1];  Be1[26]+= Ce1[2];  Be1[27]+= Ce1[3];
	Be1[36]+= Ce1[4];  Be1[37]+= Ce1[5];  Be1[38]+= Ce1[6];  Be1[39]+= Ce1[7];
    Be1[40]+= Ce1[8];  Be1[41]+= Ce1[9];  Be1[42]+= Ce1[10]; Be1[43]+= Ce1[11];
    Be1[44]+= Ce1[12]; Be1[45]+= Ce1[13]; Be1[46]+= Ce1[14]; Be1[47]+= Ce1[15];
    
    // Compute Be2 = Cs*Gstar*B0+Cs*B2 using hard coding in two steps
    Be2[0] =-Ge2[0];	Be2[1] =-Ge2[1];	Be2[6] =-Ge2[2]; 	Be2[7] =-Ge2[3];
	Be2[8] =-Ge2[4];	Be2[9] =-Ge2[5];	Be2[10]=-Ge2[6];	Be2[11]=-Ge2[7];
	Be2[12]= Ge2[0];	Be2[13]= Ge2[1];	Be2[18]= Ge2[2];	Be2[19]= Ge2[3];
	Be2[20]= Ge2[4];	Be2[21]= Ge2[5];	Be2[22]= Ge2[6];	Be2[23]= Ge2[7];
   	  
	Be2[2]-=Se[0]/L; Be2[3]-=Se[1]/L; Be2[4]-=Se[2]/L; Be2[5]-=Se[3]/L;
	Be2[8]+=Se[2]/2; Be2[9]+=Se[3]/2; Be2[10]-=Se[0]/2; Be2[11]-=Se[1]/2;
	Be2[14]+=Se[0]/L; Be2[15]+=Se[1]/L; Be2[16]+=Se[2]/L; Be2[17]+=Se[3]/L;
	Be2[20]+=Se[2]/2; Be2[21]+=Se[3]/2; Be2[22]-=Se[0]/2; Be2[23]-=Se[1]/2;
	
	// Compute N0, V, and N1 and store in Ce1,Ce2,and Ce3 respectively
	
	// Ce1 = Be1*nodDisp == N0  (only first 4 entries of Ce1)
    MatrixOperations::dgemv(0, Ce1, 1,Be1, dispVecE, 4, 12);  //Ce1:=N0 (without temperature effects)
    MatrixOperations::Add(1,1,Ce1,Ce5,Ce1,4);				 //Ce1:=N0 (with temperature efects: Ce1+=x3))
    
    // Ce3 = Be2*nodDisp == V(only first 2 entries of Ce3);
    MatrixOperations::dgemv(0, Ce3, 1,Be2, dispVecE, 2, 12);  	 // Ce3:=V (without temperature effects)
    MatrixOperations::Add(1,1,Ce3,Ge3,Ce3,2);					 // Ce3:=V (with temperature efects: Ce3+=x2)
    
    // Ce2 = 3*L/2*B3**T*V == N1 (only first 4 entries of Ce2, hard coding)
    Ce2[0] = 0; Ce2[1] = 0; Ce2[2] = Ce3[1]*L/2; Ce2[3] = -Ce3[0]*L/2; // Ce2:= N1
    
    // Finally fill in the beam element force vector at the Gauss Points
    // Use that N1[0] = N1[1] = 0 always, and N = N0+xi*N1; xi1= -sqrt(3)/3, xi2 = sqrt(3)/3
    double xip = sqrt(3)/3;
    elemForces[0] = Ce1[0]; elemForces[1] = Ce3[0]; elemForces[2] = Ce3[1];
    elemForces[3] = Ce1[1]; elemForces[4] = Ce1[2]-xip*Ce2[2]; elemForces[5] = Ce1[3]-xip*Ce2[3];
    elemForces[6] = Ce1[0]; elemForces[7] = Ce3[0]; elemForces[8] = Ce3[1];
    elemForces[9] = Ce1[1]; elemForces[10] = Ce1[2]+xip*Ce2[2]; elemForces[11] = Ce1[3]+xip*Ce2[3];
    
    // Add average values
    elemForces[12] = 0.5*(elemForces[0]+elemForces[6]);  	elemForces[13] = 0.5*(elemForces[1]+elemForces[7]);
    elemForces[14] = 0.5*(elemForces[2]+elemForces[8]);  	elemForces[15] = 0.5*(elemForces[3]+elemForces[9]);
    elemForces[16] = 0.5*(elemForces[4]+elemForces[10]);  	elemForces[17] = 0.5*(elemForces[5]+elemForces[11]);
		
	return 0;
}
	
// Compute element strain energy in beam element
int CBEAM::ElementStrainEnergy(double& elemEnergy, Solver *solver){
	// Initialize parameters
	int size = 12;
	double* elemMat = new double[size*size];
	double* tempLoadVector = new double[size];
	double* dispVec = new double[size];
	
	// Extract stiffness matrix
	ElementMatrix(elemMat, size);
	
	// Compute Temperature load vector
	TemperatureLoadVector(tempLoadVector);
	
	// Extract nodal displacement vector
	Nodes_[0]->Displacements(solver, &dispVec[0]);
	Nodes_[1]->Displacements(solver, &dispVec[DomainConstants::DISP_FIELD_SIZE]);
	
	// Compute first part of strain energy, i.e. 0.5*u**T*K*u
	elemEnergy = 0.5*MatrixOperations::dgemc(elemMat,dispVec,size);
	
	// Add second part: -u**T*Ftemp
	elemEnergy -= MatrixOperations::Dot(dispVec,tempLoadVector,size);
	
	return 0;
	
}
