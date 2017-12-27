#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 
#include "fei_base.hpp"


// Constructors
SPC::SPC(){
}

SPC::SPC(Node* spcNode,const char* preDOF, double preVal) :
	Node_(spcNode), PrescribedValue_(preVal){
	// Initialize int PrescribedDOF_ to zero:
	for(int i=0;i<6;i++){PrescribedDOF_[i]=0;}
		
	// Interpret the char preDOF, which states the prescribed dof with
	// values from 1-6
	for (const char* it = preDOF; *it; ++it){
		if(((int)(*it)-1-48)>-1 && ((int)(*it)-1-48)<6)
			PrescribedDOF_[(int)(*it)-1-48]=1;  // -1 because 1-6 in file
												// -48 because '0' = 48
		}
	
}

// Initialize SPC
int SPC::Initialize(fei::MatrixGraph* matrixGraph){
	/*
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int LocalNodeID = Node_->LocalNodeID();
	int offsetOfSlave = 0;
	double weights[] = {0,0,0,0,0,0};
	int offsetIntoSlaveField;
	
	for(int offsetIntoSlaveField=0; offsetIntoSlaveField<6; 
		offsetIntoSlaveField++){
		if(PrescribedDOF_[offsetIntoSlaveField]==1){
			weights[offsetIntoSlaveField]=1;
			if(matrixGraph->initSlaveConstraint(1,
							  &nodeTypeID,
							  &LocalNodeID,
							  &dispFieldID,
							  offsetOfSlave,
							  offsetIntoSlaveField,
							  weights,
							  PrescribedValue_) !=0 ) {
			  FEI_COUT << " ERROR: INITIALIZING SPC FAILED" << FEI_ENDL;
			  return 1;}
			weights[offsetIntoSlaveField]=0;
		}
	}
	*/
	
	/* Uncomment this for first try
	int  numIDs = 1;
	int* nodeTypes= new int[2]; nodeTypes[0]= 0; nodeTypes[1] = 0;
	int* nodeIDs  = new int[2]; nodeIDs[0] = 1; nodeIDs[1] = 2;
	int* fieldIDs = new int[2]; fieldIDs[0] = 0; fieldIDs[1] = 0;
	int offsetOfSlave = 0;
	int offsetIntoSlaveField = 1;
	
	double* weights = new double[12];
	for(int i=0; i<12; i++) weights[i] = 0;
	weights[1] = 0;
	weights[7] = 1;
	
	double rhsValue = -0.0;
	
	int ierr = matrixGraph->initSlaveConstraint(numIDs,
							  nodeTypes,
							  nodeIDs,
							  fieldIDs,
							  offsetOfSlave,
							  offsetIntoSlaveField,
							  weights,
							  rhsValue);
	
	
	
	delete [] fieldIDs;
	delete [] nodeTypes;
	delete [] nodeIDs;
	delete [] weights;
	
	FEI_COUT << " I initialized an SPC with ierr = : " << ierr << FEI_ENDL;
	
	*/
	
	/*
	// Try lagrange multipliers (this works)
	int constraintID = 0;
	int constraintTypeID = 2;
	int numIDs = 1;
	int idType = DomainConstants::NODE_TYPE_ID;
	int nodeID = 0;
	int fieldID = DomainConstants::DISP_FIELD_ID;
	int ierr = matrixGraph->initLagrangeConstraint(constraintID,
											    constraintTypeID,
											    numIDs,
											    &idType,
											    &nodeID,
											    &fieldID);
	FEI_COUT << " I initialized an SPC with ierr = : " << ierr << FEI_ENDL;										   
	*/
	
	return 0;
}

// Apply SPC
int SPC::Apply(fei::LinearSystem* linSys){
	
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int LocalNodeID = Node_->LocalNodeID();
	int offsetIntoField;
	
	// Loop trhough the 6 dof of 
	for(int offsetIntoField=0; offsetIntoField<6; offsetIntoField++){
		// Check if this dof is prescribed
		if(PrescribedDOF_[offsetIntoField]==1){
			if(linSys->loadEssentialBCs(1,
							  &LocalNodeID,
							  nodeTypeID,
							  dispFieldID,
							  &offsetIntoField,
							  &PrescribedValue_) !=0 ) {
			  FEI_COUT << " ERROR: APPLYING SPC FAILED" << FEI_ENDL;
			  return 1;}
		}
	}
	
	return 0;
}
