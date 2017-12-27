#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "fei_base.hpp"


PLOAD::PLOAD(){}

PLOAD::PLOAD(int numPressElements, Element** pressElem, double pressure)
		:  NumPressElements_(numPressElements), Pressure_(pressure){
	
	// Initialize array
	PressElements_ 			= new Element*[NumPressElements_];
	
	// Copy input to PLOAD variable
	for(std::size_t i=0;i<NumPressElements_;i++){
		PressElements_[i] = pressElem[i];	
	}
	
}
	 
PLOAD::~PLOAD(){
	delete[] PressElements_;
}


// Apply PLOAD
int PLOAD::Apply(Domain& domain, fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec){
		
	// Initialize
	double normalVec[3] = {0,0,0};			// Unit normal vector
	double area = 0.;						// Element area
	int numNodes = 0;
	Node** nodes;
	double loadVector[6] ={0,0,0,0,0,0};	// Load Vector
	
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int dispFieldSize = DomainConstants::DISP_FIELD_SIZE;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	int LocalNodeID = -1; 
	std::vector<int> DOFIndices(dispFieldSize);
	
	// Loop through the elements
	for (std::size_t id=0; id<NumPressElements_; id++){
		
		// Obtain the element type and check whether it is a shell	
		if(PressElements_[id]->ElementType() != Element::Shell){
			FEI_COUT << "WARNING: PLOAD CANNOT BE APPLIED TO ELEMENT" 
			         << domain.ElementID_LG[PressElements_[id]->LocalElementID()] 
			         << "AS THIS IS NOT SHELL ELEMENT" << FEI_ENDL;
			return 0;
		}
		
		// Compute normal direction and element area
		PressElements_[id]->NormalDirection(normalVec);
		PressElements_[id]->ElementArea(area);
				
		// Compute number of nodes and the nodes
		numNodes = domain.NodesPerElementBlock[domain.ElementBlockID[PressElements_[id]->LocalElementID()]];
		nodes = PressElements_[id]->Nodes();
		
		// Fill the load vector for each of the nodes
		MatrixOperations::Scale(Pressure_*area/numNodes, 0, normalVec,loadVector, 3); // loadVector = p*A/nNodes* normalVec
		
		// Loop through the nodes and add nodal force
		for (std::size_t nid=0; nid<numNodes; nid++){
			LocalNodeID = nodes[nid]->LocalNodeID();
			
			// Compute the DOFIndices from the LocalNodeID
			CHK_ERR(matrixGraph->getRowSpace()->getGlobalIndices(1,
															&LocalNodeID,
															nodeTypeID,
															dispFieldID,
															&DOFIndices[0]));	
			// Sum load in forceVector
			CHK_ERR( rhsVec->sumIn(dispFieldSize, &DOFIndices[0], loadVector, 0));
			
		}
	}
	
	return 0;
}

