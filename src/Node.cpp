#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp" 
#include "arssym.h"

// Constructors
Node::Node(){
}
Node::Node(double x, double y, double z){
	GlobalCoordinates_[0] = x;
	GlobalCoordinates_[1] = y;
	GlobalCoordinates_[2] = z;
}
	
Node::Node(double* globCoord){
	GlobalCoordinates_[0] = globCoord[0];
	GlobalCoordinates_[1] = globCoord[1];
	GlobalCoordinates_[2] = globCoord[2];
}

// Write output
int Node::WriteOutput(Domain& domain, Solver* solver, std::ofstream* tempFiles, bool doTitle){
	
	// Check whether nodal output needs to be printed
	if(!domain.OutputReq->PrintNodalOutput_)return 0;
	
	// If no return, then print:
	if(doTitle)
	tempFiles[0] <<FEI_ENDL<< "Nodal displacements:" << FEI_ENDL << std::left
		 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UX"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UY"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UZ"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RX"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RY"
		 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RZ" << FEI_ENDL;
	
	// Get nodal displacements
	std::vector<double> DOFValues(DomainConstants::DISP_FIELD_SIZE);
	Displacements(solver,&*DOFValues.begin());
	
												
	// Print the DOFValues to the file
	tempFiles[0] << std::left << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH) 
		 << domain.NodeID_LG[LocalNodeID_];
	UtilityFunctions::PrintMatrix(&DOFValues[0],1,DomainConstants::DISP_FIELD_SIZE,tempFiles[0]);
	
	return 0;
}

// Write Linear Buckling Output
int Node::WriteLinearBucklingOutput(Domain& domain, LinearBucklingSolver* solver, std::ofstream* tempFiles, bool doTitle){
	// Check whether nodal output needs to be printed
	if(!domain.OutputReq->PrintNodalBucklingOutput_)return 0;
	
	int nModes = solver->NumConvergedModes();
	
	// If no return, then print:
	if(doTitle){
		double buckLoad = 0;
		
		// First temp file contains list with buckling load factors
		tempFiles[0] <<FEI_ENDL<< "Buckling Load Factors:" << FEI_ENDL << std::left
			 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
			 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   BLF" << FEI_ENDL;
		
		// Loop through modes and write buckling load
		for(std::size_t fid=0; fid<nModes; fid++){
			// Compute buckLoad (note that solver returns 1/buckLoad in reversed order)
			buckLoad = 1./(solver->ArpackProblem()->Eigenvalue(nModes-fid-1));
			tempFiles[0] << std::left << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH) << fid+1;
			UtilityFunctions::PrintMatrix(&buckLoad,1,1,tempFiles[0]);
		}
			 
		// The other temp files contain buckling modes
		for(std::size_t fid=0; fid<nModes; fid++){
			buckLoad = 1./(solver->ArpackProblem()->Eigenvalue(nModes-fid-1));
			
			tempFiles[fid+1] <<FEI_ENDL<< "Buckling Mode " << fid+1 << 
					" with BLF = " << buckLoad << FEI_ENDL << std::left
							 << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH)  << "#"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UX"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UY"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   UZ"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RX"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RY"
							 << std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)  << "   RZ" << FEI_ENDL;
		}
	}
	
	// Loop through modes and the nodal displacement at this node
	std::vector<double> DOFValues(DomainConstants::DISP_FIELD_SIZE);
	
	for(std::size_t fid=0; fid<nModes; fid++){
		// Get buckling mode displacements
		BucklingModeDisplacements(solver,nModes-fid-1,&*DOFValues.begin());
		
		// Print the DOFValues to the file
		tempFiles[fid+1] << std::left << std::setw(DomainConstants::OUTPUT_INT_FIELD_WIDTH) 
			 << domain.NodeID_LG[LocalNodeID_];
		UtilityFunctions::PrintMatrix(&DOFValues[0],1,DomainConstants::DISP_FIELD_SIZE,tempFiles[fid+1]);
	}
	
	return 0;
}

// Extract nodal displacements
int Node::Displacements(Solver* solver, double* DOFValues){
	
	// First get the indices of the dof at this node
	std::vector<int> DOFIndices(DomainConstants::DISP_FIELD_SIZE);
	
	// Compute the DOFIndices from the LocalNodeID
	CHK_ERR(solver->NodeSpace()->getGlobalIndices(1,
												&LocalNodeID_,
												DomainConstants::NODE_TYPE_ID,
												DomainConstants::DISP_FIELD_ID,
												&DOFIndices[0]));
												
	// Copy the values into DOFValues
	CHK_ERR(solver->DisplacementVector()->copyOut(DomainConstants::DISP_FIELD_SIZE,
												&DOFIndices[0],
												&DOFValues[0]));
												
	return 0;
}

// Extract buckling mode displacements
int Node::BucklingModeDisplacements(LinearBucklingSolver* solver, int mode, double* DOFValues){
	
	// First get the indices of the dof at this node
	std::vector<int> DOFIndices(DomainConstants::DISP_FIELD_SIZE);
	
	// Compute the DOFIndices from the LocalNodeID
	CHK_ERR(solver->NodeSpace()->getGlobalIndices(1,
												&LocalNodeID_,
												DomainConstants::NODE_TYPE_ID,
												DomainConstants::DISP_FIELD_ID,
												&DOFIndices[0]));
	
	for(std::size_t dof=0;dof<DomainConstants::DISP_FIELD_SIZE;dof++){
		DOFValues[dof] = solver->ArpackProblem()->Eigenvector(mode,DOFIndices[dof]);
	}
												
	return 0;
}



// Setters
void Node::SetGlobalCoordinates(double x, double y, double z){
	GlobalCoordinates_[0] = x;
	GlobalCoordinates_[1] = y;
	GlobalCoordinates_[2] = z;
}
