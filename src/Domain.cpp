#include "fei_iostream.hpp" 
#include "fei_fstream.hpp"
#include "fei_base.hpp"
#include "fei_ErrMacros.hpp"


// Include own classes
#include "../include/PCH_OPTANT.h"

// Constructors
Domain::Domain(): NumTotalElements(0), NumLocalElements(0), 
				  NumTotalNodes(0) {}

Domain::Domain(const char* filename) : FileName(filename),
		NumTotalElements(0), NumLocalElements(0), NumTotalNodes(0), ElementMatricesBuilt(false) {}



// ReadInputFile
int Domain::ReadInputFile(){
	
	// Read the input file
	if(UtilityFunctions::read_input_file(FileName,*this)){
		FEI_COUT << " Error in reading input file " << FEI_ENDL;
		return 1;
	}
	return 0;
}

// Initialize Element Connectivities
int Domain::InitializeElementConnectivities(fei::MatrixGraph* matrixGraph){
	FEI_COUT << "Initializing Element Connectivities... " << FEI_ENDL;
	
	// First check whether elements exist (if not we are on another process)
	if(ElementList.size()==0){
		FEI_COUT << " No Elements" << FEI_ENDL;
		return 0;
	}
	
	// Extract dispFieldID and nodeTypeID
	int dispFieldID = DomainConstants::DISP_FIELD_ID;
	int nodeTypeID = DomainConstants::NODE_TYPE_ID;
	
	
	// Inititalize connectivity blocks, and set the index as the block ID
	int patternID, numElemPerBlock, numNodesPerElem;
	for(std::size_t blockID=0; blockID!=NumElementsBlock.size(); ++blockID){
		numNodesPerElem = NodesPerElementBlock[blockID];
		numElemPerBlock = NumElementsBlock[blockID];
		
		// Define pattern
		patternID = matrixGraph->definePattern(numNodesPerElem,
			       nodeTypeID, dispFieldID);
			       
	    // Initialize connectivity block
	    CHK_ERR( matrixGraph->initConnectivityBlock(blockID, numElemPerBlock, patternID));	     
	}
	
	// Compute the maximum number of nodes per element
	int maxNodesPerElem = *(std::max_element(NodesPerElementBlock.begin(), NodesPerElementBlock.end()));
	int* nodeIDs = new int[maxNodesPerElem];
	
	// Iterate through all elements and initialize connectivity
	int elemID = -1, blockID = -1;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		
		elemID = (*itrE)->LocalElementID();
		blockID = ElementBlockID[elemID];
		
		// Loop through nodes and get nodeIDs
		for(size_t i=0; i!=NodesPerElementBlock[blockID]; ++i)
			nodeIDs[i] = (*itrE)->Nodes()[i]->LocalNodeID();
		
		// Initialize connectivity
		CHK_ERR( matrixGraph->initConnectivity(blockID, elemID, nodeIDs) );
	}
	
	
	delete [] nodeIDs;
	FEI_COUT << " Done!" << FEI_ENDL;
	return 0;
}

// Initialize Constraints
int Domain::InitializeConstraints(fei::MatrixGraph* matrixGraph, int globalLoadCaseID){
	FEI_COUT << "Initializing Constraints... " << FEI_ENDL;
	
	LoadCase* loadcase = LoadCaseList[globalLoadCaseID];
	
	/*
	// Check if loadCases do exist on this process
	if(LoadCaseList.size() == 0){
		FEI_COUT << " No LoadCases on this process" << FEI_ENDL;
		return 0;
	}
	// Check whether this loadCaseID is in the LoadCaseList
	if(loadCaseID<0 || loadCaseID >= LoadCaseList.size()){
		FEI_COUT << "ERROR: LOADCASE ID OUT OF RANGE" << FEI_ENDL;
		return 1;}
		
			 
	// Get iterator for this loadcase
	std::vector<LoadCase*>::iterator itrLC = LoadCaseList.begin();
	std::advance(itrLC,loadCaseID);
	*/
	
	// Loop through the SPC in the SPCset of this load case and initialize
	std::multimap<const int, SPC*>::iterator itrS;
	for(itrS =SPCList.equal_range(loadcase->SPCSetID()).first; 
		itrS!=SPCList.equal_range(loadcase->SPCSetID()).second; ++itrS){
		// Initialize:
		if((*itrS).second->Initialize(matrixGraph)==1) return 1;
	}
	
	// Loop throught the MPC in the MPCset of this load case and initialize
	std::multimap<const int, MPC*>::iterator itrM;
	for(itrM =MPCList.equal_range(loadcase->MPCSetID()).first; 
		itrM!=MPCList.equal_range(loadcase->MPCSetID()).second; ++itrM){
		// Apply:
		if((*itrM).second->Initialize(matrixGraph)==1) return 1;
	}
	
	FEI_COUT << " Done!" << FEI_ENDL;
	return 0;
}

// Apply Constraints
int Domain::ApplyConstraints(fei::LinearSystem* linSys,  int globalLoadCaseID){
	FEI_COUT << "Applying Constraints... " << FEI_ENDL;
	
	LoadCase* loadcase = LoadCaseList[globalLoadCaseID];
	
	// Loop throught the SPC in the SPCset of this load case and apply
	std::multimap<const int, SPC*>::iterator itrS;
	for(itrS =SPCList.equal_range(loadcase->SPCSetID()).first; 
		itrS!=SPCList.equal_range(loadcase->SPCSetID()).second; ++itrS){
		// Apply:
		if((*itrS).second->Apply(linSys)==1) return 1;
	}
	
	// NOTE: Here we should go through all nodes with only 3 active dof,
	// and we should constrain their rotations to zero!!
	
	FEI_COUT << " Done!" << FEI_ENDL;
	 return 0;
}

// Apply Loads
int Domain::ApplyLoads(fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec, int globalLoadCaseID){
	FEI_COUT << "Applying Loads... " << FEI_ENDL;
	
	LoadCase* loadcase = LoadCaseList[globalLoadCaseID];
	
	// Loop through the LOAD in the LOADset of this load case and apply
	std::multimap<const int, LOAD*>::iterator itrL;
	for(itrL =LOADList.equal_range(loadcase->LOADSetID()).first; 
		itrL!=LOADList.equal_range(loadcase->LOADSetID()).second; ++itrL){
		// Initialize:
		if((*itrL).second->Apply(matrixGraph,rhsVec)==1) return 1;
	}
	
	// Loop through the PLOAD in the PLOADset of this load case and apply
	std::multimap<const int, PLOAD*>::iterator itrPL;
	for(itrPL =PLOADList.equal_range(loadcase->PLOADSetID()).first; 
		itrPL!=PLOADList.equal_range(loadcase->PLOADSetID()).second; ++itrPL){
		// Initialize:
		if((*itrPL).second->Apply(*this,matrixGraph,rhsVec)==1) return 1;
	}
	
	// Loop through the TEMP in the TEMPset of this load case and apply
	std::multimap<const int, TEMP*>::iterator itrT;
	for(itrT =TEMPList.equal_range(loadcase->TEMPSetID()).first; 
		itrT!=TEMPList.equal_range(loadcase->TEMPSetID()).second; ++itrT){
		// Initialize:
		if((*itrT).second->Apply(*this,matrixGraph,rhsVec)==1) return 1;
	}
	
	FEI_COUT << " Done!" << FEI_ENDL;
	 return 0;
 }
 
// Reset all element delta temperatures
int Domain::ResetElementTemperatures(){
	
	// Loop through elements and reset temperatures
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
			
		(*itrE)->ResetDeltaTemperatures();
	}
	
	return 0;
}

// Build element matrices
int Domain::BuildElementMatrices(){
	FEI_COUT << " Building Element Matrices... " << FEI_ENDL;
	
	// First compute the total number of entries in all element matrices
	// in packed storage
	ElemMatPackStorageSize = 0;
	int nDOF=0;  // number of dof per element in a block of elements
	for(std::size_t blockID=0; blockID!=NumElementsBlock.size(); ++blockID){
		nDOF = NodesPerElementBlock[blockID]*DomainConstants::DISP_FIELD_SIZE;
		
		// Number of entries per element is nDOF*(nDOF+1)/2:
		ElemMatPackStorageSize+= NumElementsBlock[blockID]*nDOF*(nDOF+1)/2;
	}
	
	FEI_COUT <<"Total Num Entries in packed storage: " << ElemMatPackStorageSize << FEI_ENDL;
	
	// Resize vector containing all element matrices in packed storage:
	AllElementMatrices.assign(ElemMatPackStorageSize,0);
	
	// Loop trhough all elements and set pointer to its element
	double *matPointer = &AllElementMatrices[0];
	int blockID, elemID=-1;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		// Set pointer:
		(*itrE)->SetElementMatrixLocation(matPointer);
		// Increment pointer
		elemID++;
		blockID = ElementBlockID[elemID];
		nDOF = NodesPerElementBlock[blockID]*DomainConstants::DISP_FIELD_SIZE;
		matPointer += (nDOF*(nDOF+1)/2);	
	}
	
	// Loop through elements and build
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE)
		if((*itrE)->BuildElementMatrix()==1) return 1;
	
	// Set ElementMatricesBuilt to true;
	ElementMatricesBuilt = true;
	
	FEI_COUT << " Done... " << FEI_ENDL;
	return 0;
}

// Build geometric element matrices
int Domain::BuildGeometricElementMatrices(Solver *solver, bool makeNegative){
	FEI_COUT << " Building Geometric Element Matrices... " << FEI_ENDL;
	
	// First compute the total number of entries in all element matrices
	// in packed storage
	ElemMatPackStorageSize = 0;
	int nDOF=0;  // number of dof per element in a block of elements
	for(std::size_t blockID=0; blockID!=NumElementsBlock.size(); ++blockID){
		nDOF = NodesPerElementBlock[blockID]*DomainConstants::DISP_FIELD_SIZE;
		
		// Number of entries per element is nDOF*(nDOF+1)/2:
		ElemMatPackStorageSize+= NumElementsBlock[blockID]*nDOF*(nDOF+1)/2;
	}
	
	FEI_COUT <<"Total Num Entries in packed storage: " << ElemMatPackStorageSize << FEI_ENDL;
	
	// Resize vector containing all geometric element matrices in packed storage:
	AllGeometricElementMatrices.assign(ElemMatPackStorageSize,0);
	
	// Loop through all elements and set pointer to its geometric stiffness element
	double *matPointer = &AllGeometricElementMatrices[0];
	int blockID, elemID=-1;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		// Set pointer:
		(*itrE)->SetGeometricElementMatrixLocation(matPointer);
		// Increment pointer
		elemID++;
		blockID = ElementBlockID[elemID];
		nDOF = NodesPerElementBlock[blockID]*DomainConstants::DISP_FIELD_SIZE;
		matPointer += (nDOF*(nDOF+1)/2);	
	}
	
	// Loop through elements and build
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE)
		if((*itrE)->BuildGeometricElementMatrix(solver,makeNegative)==1) return 1;
	
	FEI_COUT << " Done... " << FEI_ENDL;
	
	return 0;
}

// Assemble Stiffness Matrix
int Domain::AssembleStiffnessMatrix(fei::MatrixGraph* matrixGraph,
									fei::Matrix* mat, bool useGeomStiffness){
	FEI_COUT << "Assembling Global Stiffness Matrix ..." << FEI_ENDL;
	
	// Compute the maximum number of dof per element
	int maxDOFPerElem =(*(std::max_element(NodesPerElementBlock.begin(), 
						NodesPerElementBlock.end())))
						*DomainConstants::DISP_FIELD_SIZE;
	
	// Initialize element matrix
	double* elemMat = new double[maxDOFPerElem*maxDOFPerElem];
	double** elemMat2D = new double*[maxDOFPerElem];
	int matSize = 0;		//matrix Size, i.e. elemMat(matSize*matSize)
	std::vector<int> indices(maxDOFPerElem);
	
	bool block_matrix = mat->usingBlockEntryStorage();
	if (block_matrix) {
		mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::POINT_ENTRY_GRAPH);
    }
	
	// Loop through all elements and extract the element
	int elemID = -1,blockID;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		// Lood elementMatrix and size:
		if(useGeomStiffness)
			(*itrE)->GeometricElementMatrix(elemMat,matSize);
		else
			(*itrE)->ElementMatrix(elemMat,matSize);
		
		// Fill elemMat2D
		for(size_t i=0;i<matSize;++i){
			elemMat2D[i] = &(elemMat[i*matSize]);
		}
			
		// Compute elementID and blockID
		elemID++;
		blockID = ElementBlockID[elemID];	
		
		// Get connectivity indices of this element
		CHK_ERR(mat->getMatrixGraph()->getConnectivityIndices(blockID, elemID, matSize,&indices[0],matSize));	
		
		// Sum element Matrix:
		CHK_ERR( mat->sumIn(matSize, &indices[0], matSize, &indices[0],
                        elemMat2D, FEI_DENSE_COL) );
		
	}
	
	if (block_matrix) {
    mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::BLOCK_ENTRY_GRAPH);
	}
	
	
	delete[] elemMat;
	delete[] elemMat2D;
	
	
	FEI_COUT << " Done!" << FEI_ENDL;
	
	return 0;
}
/*
// Assemble Geometric Stiffness Matrix
int Domain::AssembleGeometricStiffnessMatrix(fei::MatrixGraph* matrixGraph,
									fei::Matrix* mat){
	FEI_COUT << "Assembling Global Geometric Stiffness Matrix ..." << FEI_ENDL;
	
	// Compute the maximum number of dof per element
	int maxDOFPerElem =(*(std::max_element(NodesPerElementBlock.begin(), 
						NodesPerElementBlock.end())))
						*DomainConstants::DISP_FIELD_SIZE;
	
	// Initialize element matrix
	double* elemGeomMat = new double[maxDOFPerElem*maxDOFPerElem];
	double** elemGeomMat2D = new double*[maxDOFPerElem];
	int matSize = 0;		//matrix Size, i.e. elemMat(matSize*matSize)
	std::vector<int> indices(maxDOFPerElem);
	
	bool block_matrix = mat->usingBlockEntryStorage();
	if (block_matrix) {
		mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::POINT_ENTRY_GRAPH);
    }
	
	// Loop through all elements and extract the element
	int elemID = -1,blockID;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		// Load GeometricElementMatrix and size:
		(*itrE)->GeometricElementMatrix(elemGeomMat,matSize);
		
		// Fill elemMat2D
		for(size_t i=0;i<matSize;++i){
			elemGeomMat2D[i] = &(elemGeomMat[i*matSize]);
		}
			
		// Compute elementID and blockID
		elemID++;
		blockID = ElementBlockID[elemID];	
		
		// Get connectivity indices of this element
		CHK_ERR(mat->getMatrixGraph()->getConnectivityIndices(blockID, elemID, matSize,&indices[0],matSize));	
		
		// Sum element Matrix:
		CHK_ERR( mat->sumIn(matSize, &indices[0], matSize, &indices[0],
                        elemGeomMat2D, FEI_DENSE_COL) );
		
	}
	
	if (block_matrix) {
    mat->getMatrixGraph()->setIndicesMode(fei::MatrixGraph::BLOCK_ENTRY_GRAPH);
	}
	
	delete[] elemGeomMat;
	delete[] elemGeomMat2D;
	
	FEI_COUT << " Done!" << FEI_ENDL;
	
	return 0;
}
*/

// Write Output file
int Domain::WriteOutput(Solver* solver,
							const char* outputFile, int globalLoadCaseID){
	
	//LoadCase* loadcase = LoadCaseList[globalLoadCaseID];
	
	// Create output stream
	std::ofstream fout;
	
	// Check whether the output file has been opened already for a previous
	// load case. If not, write new file; if yes, append to file;
	if(!OutputReq->OutputFileAlreadyOpened_){
		fout.open(outputFile, std::ofstream::out);
		OutputReq->OutputFileAlreadyOpened_ = true;
	}else fout.open(outputFile, std::ofstream::app);
	
	// Check how many temporary files need to be generated
	int numTempFiles = OutputReq->numNodalOutputs_;
	if(OutputReq->numElementOutputs_>OutputReq->numNodalOutputs_) numTempFiles = OutputReq->numElementOutputs_;
	FEI_COUT << " The number of temporary files is: " << numTempFiles << FEI_ENDL;
	
	// Create vector of temporary files and open them
	std::ofstream* tempFiles = new std::ofstream[numTempFiles];
	const char* tempFilePath = "temp/tempOutput";
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFiles[id].open(name.c_str(),std::ofstream::out);
		tempFiles[id] <<" "; 
	}
	/*
	// Write general load case info
	fout << "------ LOADCASE " << std::setw(3)
	     << globalLoadCaseID << " ------" << FEI_ENDL
	     << "SPC set:   " << loadcase->SPCSetID()  << FEI_ENDL
	     << "MPC set:   " << loadcase->MPCSetID()  << FEI_ENDL
	     << "LOAD set:  " << loadcase->LOADSetID() << FEI_ENDL
	     << "PLOAD set: " << loadcase->PLOADSetID() << FEI_ENDL
	     << "TEMP set:  " << loadcase->TEMPSetID() << FEI_ENDL << FEI_ENDL;
	*/
	// Write modelInfo
	OutputReq->PrintModelInfo(*this,fout);
	
	
	// Iterate through nodes and write output
	bool doTitle = true;
	for (std::vector<Node*>::iterator itrN = NodeList.begin();
		itrN != NodeList.end(); ++itrN){
		
		CHK_ERR((*itrN)->WriteOutput(*this,solver,tempFiles,doTitle));
		doTitle = false;
	}
	
	// Add all temp files to the outputfile
	CHK_ERR(HandleTempFiles(fout,tempFiles,numTempFiles,tempFilePath));
	
	
	// Iterate through elements and write output
	doTitle = false;
	int blockID = -1, elemID = -1;
	for (std::vector<Element*>::iterator itrE = ElementList.begin();
		itrE != ElementList.end(); ++itrE){
		
		elemID++;
		if(ElementBlockID[elemID] != blockID){
			// First handle the current temporary files by adding them to outputfile
			//CHK_ERR(HandleTempFiles(fout,tempFiles,numTempFiles,tempFilePath));
			
			// Update blockID and doTitle to continue to next block of elements
			blockID = ElementBlockID[elemID];
			doTitle = true;
		}else doTitle = false;
		
		CHK_ERR((*itrE)->WriteOutput(*this,solver,tempFiles,doTitle));
	}
	
	// Add all temp files to the outputfile
	CHK_ERR(HandleTempFiles(fout,tempFiles,numTempFiles,tempFilePath));
	
	// Close and delete temporary files
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFiles[id].close();
		
		if( remove( name.c_str() ) != 0 )
			FEI_COUT << "ERROR: DELETING TEMPORARY FILE FAILED" << FEI_ENDL;
	}
	
	fout << FEI_ENDL << FEI_ENDL;
	fout.close();
	
	
	/*
	// Test delta Temperatures for first element
	std::vector<Element*>::iterator itrE = ElementList.begin();
	(*itrE)->SetNumDeltaTemperatures(2);
	(*itrE)->AddDeltaTemperature(0,5);
	(*itrE)->AddDeltaTemperature(1,10);
	
	// Print contents of DeltaTemperatures
	FEI_COUT << "DeltaTemperature 1: "; UtilityFunctions::PrintArray((*itrE)->DeltaTemperatures(),(*itrE)->NumDeltaTemperatures());
	
	// Reset and print
	(*itrE)->ResetDeltaTemperatures();
	FEI_COUT << "DeltaTemperature 2: "; UtilityFunctions::PrintArray((*itrE)->DeltaTemperatures(),(*itrE)->NumDeltaTemperatures());

	// Give new delta temp vector
	(*itrE)->SetNumDeltaTemperatures(3);
	(*itrE)->AddDeltaTemperature(0,1);
	(*itrE)->AddDeltaTemperature(1,1.3);
	(*itrE)->AddDeltaTemperature(2,25.5);
	FEI_COUT << "DeltaTemperature 3: "; UtilityFunctions::PrintArray((*itrE)->DeltaTemperatures(),(*itrE)->NumDeltaTemperatures());

	
	// Extract stiffness matrix for second element
	
	// Compute the maximum number of dof per element
	int maxDOFPerElem =(*(std::max_element(NodesPerElementBlock.begin(), 
						NodesPerElementBlock.end())))
						*DomainConstants::DISP_FIELD_SIZE;
	FEI_COUT << " The maximum number of dof per element is: " << 
	          maxDOFPerElem << FEI_ENDL;
	
	std::vector<Element*>::iterator itrE = ElementList.begin(); itrE++;
	int size;
	double* stiffmat = new double[maxDOFPerElem*maxDOFPerElem];
	
	(*itrE)->ElementMatrix(stiffmat,size);
	FEI_COUT<< "stiffmat:" << FEI_ENDL; UtilityFunctions::PrintMatrix(stiffmat,size,size);
	
	delete[] stiffmat;
	
	//FEI_COUT<< "Value of stiffmat[6] is: " << stiffmat[6]<<  FEI_ENDL;
	
	// Write output for this element
	//(*itrE)->WriteOutput(*this,outReq,fout);
	*/
	
	
	
	
	
	return 0;
}

// Write Output for Linear Buckling
int Domain::WriteLinearBucklingOutput(LinearBucklingSolver* solver,
							const char* outputFile){
	
	// Create output stream
	std::ofstream fout;
	
	// Check whether the output file has been opened already for a previous
	// load case. If not, write new file; if yes, append to file;
	if(!OutputReq->OutputFileAlreadyOpened_){
		fout.open(outputFile, std::ofstream::out);
		OutputReq->OutputFileAlreadyOpened_ = true;
	}else fout.open(outputFile, std::ofstream::app);
	
	
	
	// Check how many temporary files need to be generated (nConvergedModes+1)
	int numTempFiles = solver->NumConvergedModes()+1;
	FEI_COUT << " The number of temporary files for buckling is: " << numTempFiles << FEI_ENDL;
	
	
	
	// Create vector of temporary files and open them
	std::ofstream* tempFiles = new std::ofstream[numTempFiles];
	const char* tempFilePath = "temp/tempOutput";
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFiles[id].open(name.c_str(),std::ofstream::out);
		tempFiles[id] <<" "; 
	}
	
	// Iterate through nodes and write output
	bool doTitle = true;
	for (std::vector<Node*>::iterator itrN = NodeList.begin();
		itrN != NodeList.end(); ++itrN){
		
		CHK_ERR((*itrN)->WriteLinearBucklingOutput(*this,solver,tempFiles,doTitle));
		doTitle = false;
	}
	
	// Add all temp files to the outputfile
	CHK_ERR(HandleTempFiles(fout,tempFiles,numTempFiles,tempFilePath));
	
	
	// Close and delete temporary files
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFiles[id].close();
		
		if( remove( name.c_str() ) != 0 )
			FEI_COUT << "ERROR: DELETING TEMPORARY FILE FAILED" << FEI_ENDL;
	}
	
	fout << FEI_ENDL << FEI_ENDL;
	fout.close();
	
	
	return 0;
}

// Append tempfiles to output file
int Domain::HandleTempFiles(std::ofstream& fout, std::ofstream* tempFiles, const int& numTempFiles, const char* tempFilePath){
	// Close all tempFiles
	for(std::size_t id=0; id<numTempFiles; id++){
		tempFiles[id].close();
	}
	
	// Introduce the tempfiles, but now as ifstream
	std::ifstream tempFileIn;
	
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFileIn.open(name.c_str());
		
		// Add this file to fout
		fout << tempFileIn.rdbuf();
		
		// Close the file
		tempFileIn.close();
	}
	
	// Open the temporary files again
	for(std::size_t id=0; id<numTempFiles; id++){
		std::string name = tempFilePath+std::to_string(id);
		tempFiles[id].open(name.c_str(),std::ofstream::out);
		tempFiles[id] << " "; 
	}
	
	return 0;
}
		
