#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include <unordered_map>


MPC::MPC(){
	std::cout << " I created a MPC" << std::endl;
}

MPC::MPC(int numNodes, Node** mpcNodes, int* prescribedDOF, double* weights, double rhsConstant)
		: NumNodes_(numNodes), RHSConstant_(rhsConstant){
	
	// Initialize arrays
	MPCNodes_ 			= new Node*[NumNodes_];
	PrescribedDOF_ 		= new int[NumNodes_];
	Weights_ 			= new double[NumNodes_];
	
	// Copy input to MPC variables
	for(std::size_t i=0;i<NumNodes_;i++){
		MPCNodes_[i] 			= mpcNodes[i];
		PrescribedDOF_[i] 		= prescribedDOF[i];
		Weights_[i] 			= weights[i];		
	}
	
	// Compute derived parameters:
	nodeIDs_ = new int[NumNodes_]; 			// Stores all nodeIDs
	uniqueNodeIDsCounter_ = 0;			 	// Counter for the number of unique nodeIDs
	

	// Fill nodeIDs using mpcNodes, and create map to link nodeIDs to uniqueNodeIDsCounter
	for (std::size_t i=0; i<NumNodes_;i++){
		nodeIDs_[i] = MPCNodes_[i]->LocalNodeID();
		
		// Check whether this nodeID has been used before already
		if(!nodeMap_.count(nodeIDs_[i])){			// This nodeID has not been used before
			uniqueNodeIDs_.push_back(nodeIDs_[i]); 				
			nodeMap_.insert(std::pair<int,int>(nodeIDs_[i],uniqueNodeIDsCounter_));
			uniqueNodeIDsCounter_++;
		}
	}
	
}

MPC::~MPC(){
	delete[] MPCNodes_;
}

// Initialize MPC
int MPC::Initialize(fei::MatrixGraph* matrixGraph){
	// NOTE: WE SOLVE: u_slave = sum(weights_m*u_m)-rhsConstant
	
		
	// Compute additional required parameters based on total number of unique nodeIDs
	std::vector<int> nodeTypes(uniqueNodeIDsCounter_,DomainConstants::NODE_TYPE_ID);
	std::vector<int> fieldIDs(uniqueNodeIDsCounter_,DomainConstants::DISP_FIELD_ID);
	std::vector<double> weights(uniqueNodeIDsCounter_*DomainConstants::DISP_FIELD_SIZE,0);  // Initialize all weigts to zero
	
	// Fill the weights vector for all NumNodes_
	for(size_t i=0; i<NumNodes_; i++){
		
		// Add weight at following location in weights:
		// Local nodeID related to this terms:			nodeIDs_[i]
		// Map this nodeIDs to the uniqueNodeIDs entry: nodeMap_[nodeIDs_[i]]
		// Find related dof for this term: 				nodeIDs_[i]]*DomainConstants::DISP_FIELD_SIZE+PrescribedDOF_[i]-1
		// The -1 in prev. relation is to translate the 1-6 dof to 0-5 entry id.
		weights[nodeMap_[nodeIDs_[i]]*DomainConstants::DISP_FIELD_SIZE+PrescribedDOF_[i]-1]=Weights_[i]; 
		
	}
	
	// Initialize the slave constraint
	int ierr = matrixGraph->initSlaveConstraint(uniqueNodeIDsCounter_, 		// Number of unique nodeIDs
								  &*nodeTypes.begin(),						// Node types
								  &*uniqueNodeIDs_.begin(),					// Unique nodeIDs
								  &*fieldIDs.begin(),						// Field IDs
								  0,										// Offset of slave node (first nodeID)
								  PrescribedDOF_[0]-1,						// DOF of slave node (first nodeID),
								  &*weights.begin(),						// Weights of complete equation
								  RHSConstant_);							// Constant at RHS of equation
	
	
	/*
	FEI_COUT << "uniqueNodeIDsCounter:   " << uniqueNodeIDsCounter_ << FEI_ENDL;
	FEI_COUT << "nodeIDs:       "; UtilityFunctions::PrintArray(nodeIDs_,NumNodes_);
	FEI_COUT << "uniqueNodeIDs: "; UtilityFunctions::PrintArray(&*uniqueNodeIDs_.begin(),uniqueNodeIDsCounter_);
	
	FEI_COUT << "Elements in nodeMap " << FEI_ENDL;
	for (std::unordered_map<int, int>::iterator it = nodeMap_.begin();
		it != nodeMap_.end(); ++it)
		FEI_COUT << "  [" << (*it).first << ", " << (*it).second << "]" << FEI_ENDL;
	
	
	FEI_COUT << "Weights_:      "; UtilityFunctions::PrintArray(Weights_,NumNodes_);
	FEI_COUT << "weights:       "; UtilityFunctions::PrintArray(&*weights.begin(),uniqueNodeIDsCounter_*DomainConstants::DISP_FIELD_SIZE);
	*/
	return ierr;
	
}
