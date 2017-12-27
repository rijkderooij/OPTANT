#include "../include/PCH_OPTANT.h"
#include "fei_iostream.hpp" 


TEMP::TEMP(){
	FEI_COUT << " I created a TEMP" << FEI_ENDL;
}

TEMP::TEMP(Element* element, int NumTemperatures, double* deltaTemp)
			: Element_(element), 
			  NumDeltaTemperatures_(NumTemperatures) {
	
	// Assign deltaTemp to DeltaTemperatures_
	DeltaTemperatures_ = new double[NumDeltaTemperatures_];
	
	for(std::size_t i=0;i<NumDeltaTemperatures_;i++)
			DeltaTemperatures_[i]=deltaTemp[i];
}

// Apply TEMP
int TEMP::Apply(Domain& domain, fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec){
	// Initialize variables
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int dispFieldSize = DomainConstants::DISP_FIELD_SIZE;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int LocalNodeID = -1; 
	std::vector<int> DOFIndices(dispFieldSize);
	
	// Get numNodes and initialize loadVector
	int numNodes = domain.NodesPerElementBlock[domain.ElementBlockID[Element_->LocalElementID()]];
	Node** nodes =  Element_->Nodes();
	double* loadVector = new double[numNodes*dispFieldSize];
	
	// Set this temp load in the element definition
	Element_->SetNumDeltaTemperatures(NumDeltaTemperatures_);
	for(std::size_t ind=0;ind<NumDeltaTemperatures_;ind++)
			Element_->AddDeltaTemperature(ind,DeltaTemperatures_[ind]);
	
	// Compute temperature load vector for this element
	Element_->TemperatureLoadVector(loadVector);
	
	// Loop through the nodes and add nodal force
	for (std::size_t nid=0; nid<numNodes; nid++){
		LocalNodeID = nodes[nid]->LocalNodeID();
		
		// Compute the DOFIndices from the LocalNodeID
		CHK_ERR(matrixGraph->getRowSpace()->getGlobalIndices(1,
														&LocalNodeID,
														nodeTypeID,
														dispFieldID,
														&DOFIndices[0]));	
		// Sum load corresponding to this node to forceVector
		CHK_ERR( rhsVec->sumIn(dispFieldSize, &DOFIndices[0], &loadVector[nid*dispFieldSize], 0));
	}
	
	return 0;
}

