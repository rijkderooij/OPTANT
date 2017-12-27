#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 
#include "fei_base.hpp"



LOAD::LOAD(){}

LOAD::LOAD(Node* node, double fx, double fy, double fz, double mx,
	 double my, double mz): Node_(node) {
		 
	// Initialize LoadVector_
	LoadVector_[0]=fx; LoadVector_[1]=fy; LoadVector_[2]=fz;
	LoadVector_[3]=mx; LoadVector_[4]=my; LoadVector_[5]=mz;
}
	 
LOAD::LOAD(Node* node, double* loadvec): Node_(node) {
	// Assign loadvec to LoadVector_
	for(int i=0;i<6;i++)LoadVector_[i]=loadvec[i];
}

// Apply LOAD
int LOAD::Apply(fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec){
		
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int dispFieldSize = DomainConstants::DISP_FIELD_SIZE;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int LocalNodeID = Node_->LocalNodeID();
	std::vector<int> DOFIndices(dispFieldSize);
	
	// Compute the DOFIndices from the LocalNodeID
	CHK_ERR(matrixGraph->getRowSpace()->getGlobalIndices(1,
													&LocalNodeID,
													nodeTypeID,
													dispFieldID,
													&DOFIndices[0]));	
	// Sum load in forceVector
	CHK_ERR( rhsVec->sumIn(dispFieldSize, &DOFIndices[0], LoadVector_, 0));
	
	
	return 0;
}

