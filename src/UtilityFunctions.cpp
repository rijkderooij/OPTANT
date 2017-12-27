#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp"
#include "fei_fstream.hpp"
#include "fei_base.hpp"
#include <cmath>


namespace UtilityFunctions {
	using std::cout;
	using std::endl;

// Get id of argument
int whichArg(int argc, const char*const* argv, const char* findarg)
{
  for(int i=0; i<argc; i++) {
    if (argv[i] != NULL) {
      if (strcmp(findarg, argv[i]) == 0) return(i);
    }
  }
  return(-1);
}

// Get filename given in execute command
int getFileName(int argc, char** argv, const char* argTag, char* destination){
	std::string filename;
	const char* result;

	int dashIarg = whichArg(argc, argv, argTag);

	if (dashIarg < 0) {
		fei::console_out() << "fei_test_utils::construct_filename: argument '" 
			<< argTag <<"' not found." << FEI_ENDL;
		return(-1);
	}

	if (argc > dashIarg+1) {
	if (argv[dashIarg+1] != 0) {
	  filename = argv[dashIarg+1];
	}
	}
	
	strcpy(destination,filename.c_str());
	return 0;
}


// Set of functions to read input file
int read_input_file(const char* filename, Domain & domain){
	int ierr = 0;
	
	// Create ifstream
	std::ifstream fin(filename, std::ifstream::in);
	
	// Read materials
	ierr += ReadMaterials(fin, domain);
	
	// Read properties
	ierr += ReadProperties(fin,domain);
	
	// Read nodes
	ierr += ReadNodes(fin,domain);
	
	// Read elements
	ierr += ReadElements(fin,domain);
	
	// Read SPCs
	ierr += ReadSPC(fin,domain);
	
	// Read MPCs
	ierr += ReadMPC(fin,domain);
	
	// Read LOADs
	ierr += ReadLOAD(fin,domain);
	
	// Read PLOADs
	ierr += ReadPLOAD(fin,domain);
	
	// Read TEMPs
	ierr += ReadTEMP(fin,domain);
	
	// Read LoadCases
	ierr += ReadLoadCases(fin,domain);
	
	// Read solution type
	ierr += ReadSolutionType(fin, domain);
	
	return ierr;
}

// Read SolutionType
int ReadSolutionType(std::ifstream& fin, Domain &domain){
	
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char solTypeStart[] = "$Solution0";
	char solTypeEnd[]   = "$Solution1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find solTypeStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,solTypeStart)==0){
			break;
		}
	}
	
	// Check if solTypeStart was found:
	if(strcmp(line,solTypeStart)!=0){
			FEI_COUT << "WARNING: NO SOLUTION TYPE IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// Declare variables
	char tempCard[DomainConstants::MAXNAMESIZE]; // Should be SOL 
	int solutionType; 
	
	// First line gives solution type
	fin.getline(line,DomainConstants::MAXLINESIZE);
	std::stringstream parse(line);
	if(parse >> tempCard >> solutionType);
		else{FEI_COUT << " ERROR: IN SOLUTION TYPE DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;}
	
	// Check tempCard
	if(strcmp(tempCard,"SOL")!=0){
		FEI_COUT << " ERROR: SOLUTION TYPE CARD: " << tempCard <<
		" IS EXPECTED TO BE: SOL" << FEI_ENDL;
		return 1;		
	}
	
	// Set solutionType to domain
	domain.SolutionType = (Domain::SolType)solutionType;
	
	switch(solutionType){
		case Domain::LinearStatic:
			// Nothing has to be done
			FEI_COUT << "This is a Linear Static Solution type: " << domain.SolutionType << FEI_ENDL;
			break;
		case Domain::LinearBuckling:
			{
				// The next line contains the number of required modes
				FEI_COUT << "This is a Linear Buckling Solution type: " << domain.SolutionType << FEI_ENDL;
				int nModes;
				
				// ------------- Read next line: nMODES --------------//
				fin.getline(line,DomainConstants::MAXLINESIZE);
				parse.str(line);
				
				if(parse >> tempCard >> nModes);
					else{FEI_COUT << " ERROR: IN SOLUTION TYPE DESCRIPTION FORMAT" << FEI_ENDL;
						 return 1;}
				
				// Check tempCard
				if(strcmp(tempCard,"NMODES")!=0){
					FEI_COUT << " ERROR: CARD: " << tempCard <<
					" IS EXPECTED TO BE: NMODES" << FEI_ENDL;
					return 1;		
				}
				
				// Check whether nModes is 1 or larger
				if(nModes<1){
					FEI_COUT << " ERROR: THE NUMBER OF DESIRED BUCKLING MODES SHOULD BE GREATER THAN 0" << FEI_ENDL;
					return 1;
				}
				
				// Set nmodes to SolverParams
				fei::Param iparam("nModes", nModes);
				domain.SolverParams.add(iparam,false);
				
				// Check nModes
				int nModesTest = 0;
				domain.SolverParams.getIntParamValue("nModes", nModesTest); 
				FEI_COUT << "The number of modes is: " << nModesTest << FEI_ENDL;
				
				
				// ------------ Read next line: preStress ------------//
				fin.getline(line,DomainConstants::MAXLINESIZE);
				parse.str(line);
				
				// Check for end statement
				if(strcmp(line,solTypeEnd)==0)return 0;
				
				int preStress = 0;
				if(parse >> tempCard >> preStress);
					else{FEI_COUT << " ERROR: IN SOLUTION TYPE DESCRIPTION FORMAT" << FEI_ENDL;
						 return 1;}
				
				// Check tempCard
				if(strcmp(tempCard,"PRESTRESS")!=0){
					FEI_COUT << " ERROR: CARD: " << tempCard <<
					" IS EXPECTED TO BE: PRESTRESS" << FEI_ENDL;
					return 1;		
				}
				
				// Check whether preStress loadCase is defined
				if(domain.LoadCaseList.find(preStress) == domain.LoadCaseList.end()) {
				   FEI_COUT << " ERROR: THE PRESTRESS LOADCASE " << preStress 
							<< " IS UNDEFINED" << FEI_ENDL;
					return 1;
				}
				
				// Set prestress to SolverParams
				fei::Param pparam("preStress", preStress);
				domain.SolverParams.add(pparam,false);
				
				// Check nModes
				int preStressTest = 0;
				domain.SolverParams.getIntParamValue("preStress", preStressTest); 
				FEI_COUT << "The prestress loadcase is: " << preStressTest << FEI_ENDL;
				
				break;
			}
		default:
			FEI_COUT << " ERROR: THE PROVIDED SOLUTION TYPE IS NOT SUPPORTED" << FEI_ENDL;
			return 1;
	}
	
	// The final line should be solTypeEnd:
	fin.getline(line,DomainConstants::MAXLINESIZE);
	if(strcmp(line,solTypeEnd)!=0){
		FEI_COUT << "ERROR: THE END STATEMENT " << solTypeEnd << " IS EXPECTED" << FEI_ENDL;
		return 1;
	}
	return 0;
}

// Read Materials
int ReadMaterials(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char matStart[] = "$Materials0";
	char matEnd[]   = "$Materials1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find matStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,matStart)==0){
			break;
		}
	}
	// Check if matStart was found:
	if(strcmp(line,matStart)!=0){
			FEI_COUT << "WARNING: NO MATERIALS IN INPUT FILE" << FEI_ENDL;
			return 0;}
	
	
	// The first line always contains: matID matType numLines matName
	int matID,matType,numLines;
	char matName[DomainConstants::MAXNAMESIZE];
	
	// Loop through the materials, read the first line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for matEnd string
		if(strcmp(line,matEnd)==0){
			break;
		}
		
		// If break was not applied, another material is defined:
		std::stringstream parse(line);
		parse >> matID >> matType >> numLines >> matName;
		
		// Compute lenght of propArray based on matType
		int propArrayLength;
		switch (matType){
			case 0: propArrayLength = 6; break;
			default: propArrayLength = 0; break;
		}
		double propArray[propArrayLength];
		
		// Next loop through the following numLines lines and fill propArray
		int count = 0;
		for (int i=0; i<numLines; i++){
			// Get next line and parse it
			fin.getline(line,DomainConstants::MAXLINESIZE);	
			std::stringstream parse(line);
			
			// Substitute in propArray;
			while(parse>>propArray[count]) count++;
		}
		
		// Check if count equals propArrayLength
		if(count!=propArrayLength){
			FEI_COUT << "ERROR IN MATERIAL INPUT PARAMETERS!!!" << FEI_ENDL;
			return 1;}
		
		// Create new Material
		domain.MaterialList.insert(
			std::pair<const int, Material*>(
				matID, new Material(matName, matType,propArray)));
				
	}
	
	return 0;
}

// Read Properties
int ReadProperties(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char propStart[] = "$Properties0";
	char propEnd[]   = "$Properties1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find propStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,propStart)==0){
			break;
		}
	}
	
	// Check if propStart was found:
	if(strcmp(line,propStart)!=0){
			FEI_COUT << "WARNING: NO PROPERTIES IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// The first line always contains: propID propOption numLines propCard
	int propID, propOption,numLines;
	char propCard[DomainConstants::MAXNAMESIZE];
	
	// Loop through the properties, read the first line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for propEnd string
		if(strcmp(line,propEnd)==0){
			break;
		}
		
		// If break was not applied, another property is defined:
		std::stringstream parse(line);
		parse >> propID >> propOption >> numLines >> propCard;
		
		
		// Switch between propertyCards to create property
		if(strcmp(propCard,"PBEAM")==0){
			// Create new PBEAM property
			domain.PropertyList.push_back(new PBEAM());
		}else if(strcmp(propCard,"PSHELL")==0){
			// Create new PSHELL property
			domain.PropertyList.push_back(new PSHELL());
		}else{
			FEI_COUT << "ERROR: PROP CARD " << propCard 
			         << " NOT SUPPORTED!!" <<FEI_ENDL;
			return 1;
		}
		
		// Add the propID to the Local Global conversion
		domain.PropertyID_LG.push_back(propID);
		domain.PropertyID_GL.insert(
		std::pair<const int,const int>(
			propID, domain.PropertyID_LG.size()-1));
			
		// Read the following numLines lines for this property
		if(domain.PropertyList.back()
				 ->ReadFromFile(fin,propOption, numLines,domain)){
			FEI_COUT << "ERROR: READING PROP " << propID 
			         << " FAILED!!" <<FEI_ENDL;
			return 1;
		}
			
			
				
	}
	
	
	return 0;
}

// Read Nodes
int ReadNodes(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char nodeStart[] = "$Nodes0";
	char nodeEnd[]   = "$Nodes1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find nodeStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,nodeStart)==0){
			break;
		}
	}
	
	// Check if nodeStart was found:
	if(strcmp(line,nodeStart)!=0){
			FEI_COUT << "WARNING: NO NODES IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// First line is total number of nodes
	int nTotalNodes=0;
	
	// Read line and parse
	fin.getline(line,DomainConstants::MAXLINESIZE);
	std::stringstream parse(line);
	if(strcmp(line,nodeEnd)!=0 && parse>>nTotalNodes);
	else{ FEI_COUT << "ERROR: # OF NODES COULD NOT BE READ" << FEI_ENDL;
			return 1;}
	domain.NumTotalNodes=nTotalNodes;
	
	// Declare variables
	double x,y,z;
	int GlobalNodeID; //global nodeid;
	
	// First resize NodeList and initialize with 0's
	domain.NodeList.resize(nTotalNodes,0);
	
	// Then create all Nodes
	for(int i=0; i<nTotalNodes; i++){
		domain.NodeList[i] = new Node();
	}
	
	// Loop through the nTotalNodes lines and fill the nodes:
	for(int i=0; i<nTotalNodes; i++){
		// Read line and parse
		fin.getline(line,DomainConstants::MAXLINESIZE);
		std::stringstream parse(line);
		
		// Fill gid, x,y,z
		if(parse >> GlobalNodeID >> x >> y >> z){
			
			domain.NodeList[i]->SetGlobalCoordinates(x,y,z);
			domain.NodeList[i]->SetLocalNodeID(i);
			
			
			// Add the nodeID to the Local Global conversion
			domain.NodeID_LG.push_back(GlobalNodeID);
			domain.NodeID_GL.insert(
				std::pair<const int,const int>(
					GlobalNodeID, i));
					
			// Set default number of DOF on this node to 3:
			domain.NodalNumDOF.push_back(3);	
		}else{
			FEI_COUT << "ERROR: # OF NODES INCONSISTENT" << 
			         " or WRONG NODE DESCRIPTION" << FEI_ENDL;
			return 1;}
		
	}
	// Next line should be nodeEnd[]:
	fin.getline(line,DomainConstants::MAXLINESIZE);
	if(strcmp(line,nodeEnd)!=0){
		FEI_COUT << "ERROR: # OF NODES INCONSISTENT" << FEI_ENDL;
		return 1;
	}
	
	
	return 0;
}

// Read Elements
int ReadElements(std::ifstream& fin, Domain &domain){

	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char elementStart[] = "$Elements0";
	char elementEnd[]   = "$Elements1";
	char elementType[]  = "$ElementType";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find elementStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,elementStart)==0){
			break;
		}
	}
	
	// Check if elementStart was found:
	if(strcmp(line,elementStart)!=0){
			FEI_COUT << "WARNING: NO ELEMENTS IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	// Next the first line after elementShart should contain elementType[]		
	fin.getline(line,DomainConstants::MAXLINESIZE);
	if(strcmp(line,elementType)!=0){
			FEI_COUT << "ERROR: EXPECTED: " << elementType <<
			" IN ELEMENT DECLARATION "<< FEI_ENDL;
			return 1;}
	
	// Loop through element types, but first declare some variables
	// Size of elementList so far:
	int numTotElemPrev = domain.ElementList.size(); 
			
	// The first line after elementType contains the element type and 
	// and number of elements in this type;
	char elementCard[DomainConstants::MAXNAMESIZE];
	int numElem=0, NumNodePerElem=0;
	while(1){	
		
		// Read line and parse
		fin.getline(line,DomainConstants::MAXLINESIZE);
		std::stringstream parse(line);
		if(parse >> elementCard >> numElem);
		else {
			FEI_COUT << "ERROR: COULD NOT READ ELEMENT CARD " <<
				" AND NUMBER OF ELEMENTS "<< FEI_ENDL;
			return 1;}
	
		
		// First resize ElementList and initialize with 0
		domain.ElementList.resize(numTotElemPrev+numElem,0);
		
		// The read elements from file and add to ElementList
		// Create empty elements of element type
		if(strcmp(elementCard,"CBEAM")==0){
			for(int i=0; i<numElem; i++){
				domain.ElementList[numTotElemPrev+i] =new CBEAM();
			}
			NumNodePerElem=2;
		}else if(strcmp(elementCard,"CTRIA")==0){
			for(int i=0; i<numElem; i++){
				domain.ElementList[numTotElemPrev+i] =new CTRIA();
			}
			NumNodePerElem=3;
		}else if(strcmp(elementCard,"CQUAD")==0){
			for(int i=0; i<numElem; i++){
				domain.ElementList[numTotElemPrev+i] =new CQUAD();
			}
			NumNodePerElem=4;
		}else{
			FEI_COUT << "ERROR: ELEMENT CARD " << elementCard 
			         << " NOT SUPPORTED!!" <<FEI_ENDL;
			return 1;
		}
		
		// Read these elements from file
		for(int i=0; i<numElem; i++){
			if(domain.ElementList[numTotElemPrev+i]->
											ReadFromFile(fin, domain))
				return 1;
		}
		
		// Resize ElementBlockID and initialize with the current lenght 
		// of NumElementsBlock
		domain.ElementBlockID.resize(numTotElemPrev+numElem,domain.NumElementsBlock.size());
		
		// This new element card will give a new pattern and block in the
		// calculations: add the info:
		domain.NodesPerElementBlock.push_back(NumNodePerElem);
		domain.NumElementsBlock.push_back(numElem);
		
		// Increment NumTotalElements
		domain.NumTotalElements = domain.ElementList.size();
		numTotElemPrev = domain.ElementList.size();
		
		// The next line should contain elementType[] or elementEnd[].
		fin.getline(line,DomainConstants::MAXLINESIZE);
		if(strcmp(line,elementType)==0){
			FEI_COUT << "LET'S GO TO THE NEXT ELEMENT TYPE" << FEI_ENDL;
		}else if(strcmp(line,elementEnd)==0){
			FEI_COUT << "WE GOT ALL ELEMENTS CORRECTLY" << FEI_ENDL;
			break;
		}else{
			FEI_COUT << "ERROR IN READING ELEMENTS, EXPECTED: " << 
			 elementType << " OR " << elementEnd << FEI_ENDL;
			return 1;
		}
	}
		
	return 0;
}

// Read SPCs
int ReadSPC(std::ifstream& fin, Domain &domain){
	
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char spcStart[] = "$SPC0";
	char spcEnd[]   = "$SPC1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find spcStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,spcStart)==0){
			break;
		}
	}
	
	// Check if spcStart was found:
	if(strcmp(line,spcStart)!=0){
			FEI_COUT << "WARNING: NO SPCs IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// Declare variables
	char PrescribedDOF[DomainConstants::MAXNAMESIZE];
	char spcCard[DomainConstants::MAXNAMESIZE]; // Should be SPC
	double PrescribedValue;
	int GlobalSPCsetID; 
	int GlobalNodeID; // nodeid;
	Node* node;       // constrained node
	
	// Loop through the SPCs, and read the line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for spcEnd string
		if(strcmp(line,spcEnd)==0){
			break;
		}
		
		// If break was not applied, another SPC is defined:
		std::stringstream parse(line);
		if(parse >> spcCard >> GlobalSPCsetID >>
			 GlobalNodeID >> PrescribedDOF >> PrescribedValue);
		else{FEI_COUT << " ERROR: IN SPC DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;}
		
		// Check spcCard
		if(strcmp(spcCard,"SPC")!=0){
			FEI_COUT << " ERROR: SPC CARD: " << spcCard <<
			" IS EXPECTED TO BE: SPC" << FEI_ENDL;
			return 1;		
		}
		
		// Compute localNodeID
		if(domain.NodeID_GL.find(GlobalNodeID) == domain.NodeID_GL.end()) {
		   FEI_COUT << " ERROR: NODE " << GlobalNodeID << " IN SPC " 
					<< GlobalSPCsetID << " UNDEFINED" << FEI_ENDL;
			return 1;
		}
		node = domain.NodeList[domain.NodeID_GL[GlobalNodeID]];
		
		// Create new SPC and add to list:
		domain.SPCList.insert(
			std::pair<const int, SPC*>(
				GlobalSPCsetID, 
				new SPC(node, PrescribedDOF,PrescribedValue)));
	}
	
	return 0;
}

// Read MPCs
int ReadMPC(std::ifstream& fin, Domain &domain){
	
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char mpcStart[] = "$MPC0";
	char mpcEnd[]   = "$MPC1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find mpcStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,mpcStart)==0){
			break;
		}
	}
	
	// Check if mpcStart was found:
	if(strcmp(line,mpcStart)!=0){
			FEI_COUT << "WARNING: NO MPCs IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	// Variables
	double temp;
	char mpcCard[DomainConstants::MAXNAMESIZE]; // Should be MPC
	double rhsConstant;				// Constant at rhs of MPC
	int GlobalMPCsetID; 
	int GlobalNodeID; 				// nodeid;
	int numNodes = 0;				// Number of nodes per MPC
			
	// Loop through the MPCs, and read the line
	int doContinue = 1;				// 1 - other mpcs, 0 - no more mpcs, i.e. stop
	while(doContinue){
		// Store current position in file, to return to after number of nodes have been counted
		int posReturn = fin.tellg();
		
		// Read first line of MPC definition
		if(fin.getline(line,DomainConstants::MAXLINESIZE));else break;
		
		// Find the number of nodes by looping through the next few lines
		numNodes = 1; 		// Count the slave nodes already
						
		while(fin.getline(line,DomainConstants::MAXLINESIZE)){
			// Check whether this line already defines a new MPC, or contains the mpcEnd string
			if(strcmp(line,mpcEnd)==0){ 		// Check for mpcEnd string
				doContinue = 0;
				break;
			}else if(strncmp(line,"MPC",3)==0){ 	// Check whether this line already defines the next MPC
				break;
			}
			
			// If break is not applied, this line constains 1 or two master nodes, represented by three values:
			std::stringstream parse(line);
			int countLine = 0;
			while(parse>>temp>>temp>>temp){
				if(countLine<2){
					numNodes++;
					countLine++;
				}else{
					FEI_COUT << " ERROR: MPC DESCRIPTION FORMAT ALLOWS " 
					         << "ONLY TWO MASTER TERMS PER LINE" << FEI_ENDL;
					return 1;
				}
			}	
		}
		
		// Check numNodes>1
		if(!(numNodes>1)){
			FEI_COUT << " ERROR: MPC DEFINITIONS REQUIRE AT LEAST ONE MASTER TERM" << FEI_ENDL;
			return 1;
		}
		
		// Return to position in file before counting started.
		fin.seekg(posReturn);
		
		
		// -----------------------------------------------------------//
		
		// Initialize variables		
		Node** mpcNodes = new Node*[numNodes];
		double* weights = new double[numNodes];			// Weights of all terms in MPC
		int* prescribedDOF = new int[numNodes];  		// DOF of all terms in MPC
		
		// Read first line:
		fin.getline(line,DomainConstants::MAXLINESIZE);
		std::stringstream parse(line);
		
		if(parse >> mpcCard >> GlobalMPCsetID >>
			 GlobalNodeID >> prescribedDOF[0] >> rhsConstant){
				// Initialize weight of slave node to 1
				weights[0] = 1;
				// Compute localNodeID of this slave node
				if(domain.NodeID_GL.find(GlobalNodeID) == domain.NodeID_GL.end()) {
				   FEI_COUT << " ERROR: NODE " << GlobalNodeID << " IN MPC " 
							<< GlobalMPCsetID << " UNDEFINED" << FEI_ENDL;
					return 1;
				}
				mpcNodes[0] = domain.NodeList[domain.NodeID_GL[GlobalNodeID]];
			}
		else{FEI_COUT << " ERROR: IN MPC DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;}
		
		// Check mpcCard
		if(strcmp(mpcCard,"MPC")!=0){
			FEI_COUT << " ERROR: MPC CARD: " << mpcCard <<
			" IS EXPECTED TO BE: MPC" << FEI_ENDL;
			return 1;		
		}
		
		
		// The number of lines that follow after the first line is ceil((numNodes-1)/2), which is numNodes/2 for integers
		int count = 1; 						// Counter for nodes, starts at 1 since slave node has been determined already.
		for(std::size_t i=0; i<(numNodes/2); i++){
			// Read the line
			fin.getline(line,DomainConstants::MAXLINESIZE);
			std::stringstream parse(line);
			
			// Read master nodes on this line
			int countLine = 0;
			while(countLine<2 && count<numNodes){
				countLine++;
				
				if(parse >> GlobalNodeID >> prescribedDOF[count] >> weights[count]){
					// Compute localNodeID of this master node
					if(domain.NodeID_GL.find(GlobalNodeID) == domain.NodeID_GL.end()) {
					   FEI_COUT << " ERROR: NODE " << GlobalNodeID << " IN MPC " 
								<< GlobalMPCsetID << " UNDEFINED" << FEI_ENDL;
						return 1;
					}
					mpcNodes[count] = domain.NodeList[domain.NodeID_GL[GlobalNodeID]];
					
					// Increment the node count
					count++;
				}
				else{FEI_COUT << " ERROR: IN MPC DESCRIPTION FORMAT" << FEI_ENDL;
					 return 1;}
				 
			}
		}
		
		// Create new MPC and add to list:
		domain.MPCList.insert(
			std::pair<const int, MPC*>(
				GlobalMPCsetID, 
				new MPC(numNodes,mpcNodes, prescribedDOF,weights,rhsConstant)));
		
		// Delete the pointer-pointer array
		delete[] mpcNodes;
		
	}
	
	return 0;
}

// Read LOADs
int ReadLOAD(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char loadStart[] = "$LOAD0";
	char loadEnd[]   = "$LOAD1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find loadStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,loadStart)==0){
			break;
		}
	}
	
	// Check if loadStart was found:
	if(strcmp(line,loadStart)!=0){
			FEI_COUT << "WARNING: NO LOADs IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// Declare variables
	char loadCard[DomainConstants::MAXNAMESIZE]; // Should be LOAD
	double LoadVec[6];
	int GlobalLOADsetID; 
	int GlobalNodeID; // nodeid;
	Node* node;
	
	// Loop through the LOADs, and read the line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for loadEnd string
		if(strcmp(line,loadEnd)==0){
			break;
		}
		
		// If break was not applied, another LOAD is defined:
		std::stringstream parse(line);
		if(parse >> loadCard >> GlobalLOADsetID >> GlobalNodeID >> 
		     LoadVec[0] >> LoadVec[1] >> LoadVec[2] >> LoadVec[3] >> 
		     LoadVec[4] >> LoadVec[5]);
		else{FEI_COUT << " ERROR: IN LOAD DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;} 
		
		// Check loadCard
		if(strcmp(loadCard,"LOAD")!=0){
			FEI_COUT << " ERROR: LOAD CARD: " << loadCard <<
			" IS EXPECTED TO BE: LOAD" << FEI_ENDL;
			return 1;		
		}
		
		// Compute localNodeID
		if(domain.NodeID_GL.find(GlobalNodeID) == domain.NodeID_GL.end()) {
		   FEI_COUT << " ERROR: NODE " << GlobalNodeID << " IN LOAD " 
					<< GlobalLOADsetID << " UNDEFINED" << FEI_ENDL;
			return 1;
		}
		node = domain.NodeList[domain.NodeID_GL[GlobalNodeID]];
		
		// Create new LOAD and add to list:
		domain.LOADList.insert(
			std::pair<const int, LOAD*>(
				GlobalLOADsetID, 
				new LOAD(node, LoadVec)));
	}
	
	return 0;
}

// Read PLOADs
int ReadPLOAD(std::ifstream& fin, Domain &domain){
	
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char ploadStart[] = "$PLOAD0";
	char ploadEnd[]   = "$PLOAD1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find ploadStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,ploadStart)==0){
			break;
		}
	}
	
	// Check if ploadStart was found:
	if(strcmp(line,ploadStart)!=0){
			FEI_COUT << "WARNING: NO PLOADs IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	// Variables
	double temp;
	char tempChar[DomainConstants::MAXNAMESIZE];
	char ploadCard[DomainConstants::MAXNAMESIZE]; // Should be PLOAD
	double pressure;				// Pressure of PLOAD
	int GlobalPLOADsetID; 
	int GlobalElementID;
	int numElem = 0;				// Number of elements per PLOAD
	int count=0;					// Counter for reading elements
			
	// Loop through the PLOADs, and read the line
	int doContinue = 1;				// 1 - other ploads, 0 - no more ploads, i.e. stop
	while(doContinue){
		// Store current position in file, to return to after number of elements has been counted
		int posReturn = fin.tellg();
		
		// Read first line of PLOAD definition
		if(fin.getline(line,DomainConstants::MAXLINESIZE));else break;
		
		// Count number of elements on this first line
		numElem = 0; 		// Initialized
		std::stringstream parseF(line);
		if(parseF>>tempChar>>temp>>temp); 	// First three entries of first line
		else{FEI_COUT << " ERROR: IN PLOAD DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;
		}
			 
		int countLine = 0;
		while(parseF>>temp){
			if(countLine<8){
				numElem++;
				countLine++;
			}else{
				FEI_COUT << " ERROR: PLOAD DESCRIPTION FORMAT ALLOWS " 
						 << "ONLY EIGHT ELEMENTS PER LINE" << FEI_ENDL;
				return 1;
			}
		}	
		
		// Find the number of elements in next lines of input file
						
		while(fin.getline(line,DomainConstants::MAXLINESIZE)){
			// Check whether this line already defines a new MPC, or contains the mpcEnd string
			if(strcmp(line,ploadEnd)==0){ 		// Check for ploadEnd string
				doContinue = 0;
				break;
			}else if(strncmp(line,"PLOAD",5)==0){ 	// Check whether this line already defines the next PLOAD
				break;
			}
			
			// If break is not applied, this line constains 1 or more elements (up to 8), represented by its global id:
			std::stringstream parse(line);
			int countLine = 0;
			while(parse>>temp){
				if(countLine<8){
					numElem++;
					countLine++;
				}else{
					FEI_COUT << " ERROR: PLOAD DESCRIPTION FORMAT ALLOWS " 
					         << "ONLY EIGHT ELEMENTS PER LINE" << FEI_ENDL;
					return 1;
				}
			}	
		}
		
		// Check numNodes>1
		if(!(numElem>0)){
			FEI_COUT << " ERROR: PLOAD DEFINITIONS REQUIRE AT LEAST ONE ELEMENT" << FEI_ENDL;
			return 1;
		}
		
		// Return to position in file before counting started.
		fin.seekg(posReturn);
		
		
		// -----------------------------------------------------------//
		
		// Initialize variables		
		Element** pressElem = new Element*[numElem]; 			// Array of element pointers
		
		// Read first line:
		fin.getline(line,DomainConstants::MAXLINESIZE);
		std::stringstream parse(line);
		
		if(parse >> ploadCard >> GlobalPLOADsetID >> pressure );
		else{FEI_COUT << " ERROR: IN PLOAD DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;}
		
		// Check ploadCard
		if(strcmp(ploadCard,"PLOAD")!=0){
			FEI_COUT << " ERROR: PLOAD CARD: " << ploadCard <<
			" IS EXPECTED TO BE: PLOAD" << FEI_ENDL;
			return 1;		
		}
		
		// Read the element ids in this line and following lines.
		// The number of lines (including first) is ceil(numElem/8)
		count =0;
		while(count<numElem){
						
			// Check whether a new line has to be read
			if(count>0 && (count%8)==0){
				fin.getline(line,DomainConstants::MAXLINESIZE);
				parse.str(line);
			}
			
			// Read global element id
			if(parse >> GlobalElementID){
				// Compute localElementID of this master node
				if(domain.ElementID_GL.find(GlobalElementID) == domain.ElementID_GL.end()) {
				   FEI_COUT << " ERROR: ELEMENT " << GlobalElementID << " IN PLOAD " 
							<< GlobalPLOADsetID << " UNDEFINED" << FEI_ENDL;
					return 1;
				}
				pressElem[count] = domain.ElementList[domain.ElementID_GL[GlobalElementID]];
				
				// Increment the element count
				count++;
			}
			else{FEI_COUT << " ERROR: IN PLOAD DESCRIPTION FORMAT" << FEI_ENDL;
				 return 1;}
				 
			
		}
		
		// Create new PLOAD and add to list:
		domain.PLOADList.insert(
			std::pair<const int, PLOAD*>(
				GlobalPLOADsetID, 
				new PLOAD(numElem,pressElem, pressure)));
		
		// Delete the pointer-pointer array
		delete[] pressElem;
	}
	
	return 0;
}

// Read TEMPs
int ReadTEMP(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char tempStart[] = "$TEMP0";
	char tempEnd[]   = "$TEMP1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find tempStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,tempStart)==0){
			break;
		}
	}
	
	// Check if tempStart was found:
	if(strcmp(line,tempStart)!=0){
			FEI_COUT << "WARNING: NO TEMPs IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// Declare variables
	const int maxTempPerElem = 6;				 // Maximum number of temperatures that can be defined per element
	int numTempPerElem = 0;						 // Number of delta temperatures defined for the element
	char tempCard[DomainConstants::MAXNAMESIZE]; // Should be TEMP
	double DeltaTemp[maxTempPerElem];
	double temp = 0;							 // Temporary variable to store deltaTemp
	int GlobalTEMPsetID; 
	int GlobalElementID; 						 // elementid;
	Element* element;
	
	// Loop through the TEMPs, and read the line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for tempEnd string
		if(strcmp(line,tempEnd)==0){
			break;
		}
		
		// If break was not applied, another TEMP is defined:
		std::stringstream parse(line);
		
		// This line first contains tempCard, GlobalTEMPsetID and GlobalElementID
		if(parse >> tempCard >> GlobalTEMPsetID >> GlobalElementID);
		else{FEI_COUT << " ERROR: IN TEMP DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;} 
		
		// Check tempCard
		if(strcmp(tempCard,"TEMP")!=0){
			FEI_COUT << " ERROR: TEMP CARD: " << tempCard <<
			" IS EXPECTED TO BE: TEMP" << FEI_ENDL;
			return 1;		
		}
		
		// Compute localElementID
		if(domain.ElementID_GL.find(GlobalElementID) == domain.ElementID_GL.end()) {
		   FEI_COUT << " ERROR: ELEMENT " << GlobalElementID << " IN TEMP " 
					<< GlobalTEMPsetID << " UNDEFINED" << FEI_ENDL;
			return 1;
		}
		element = domain.ElementList[domain.ElementID_GL[GlobalElementID]];
		
		// Read the number of delta temperature and temperatures
		numTempPerElem=0;
		
		while(parse>>temp){
			if(numTempPerElem>=maxTempPerElem){
				FEI_COUT << " ERROR: THE MAXIMUM NUMBER OF TEMPERATURES" 
				         << " PER ELEMENT ("<< maxTempPerElem << ") IS" 
				         << " EXCEEDED IN TEMP " << GlobalTEMPsetID << FEI_ENDL;
				return 1;
			}
			
			DeltaTemp[numTempPerElem]=temp;
			numTempPerElem++;
		}
		
		// Create new TEMP and add to list:
		domain.TEMPList.insert(
			std::pair<const int, TEMP*>(
				GlobalTEMPsetID, 
				new TEMP(element, numTempPerElem, DeltaTemp)));
	
	}
	
	return 0;
}

// Read LoadCases
int ReadLoadCases(std::ifstream& fin, Domain &domain){
	
	// Go to beginning of file
	fin.clear();
	fin.seekg(0, std::ios::beg);
	
	char lcStart[] = "$LoadCases0";
	char lcEnd[]   = "$LoadCases1";
	
	char line[DomainConstants::MAXLINESIZE];

	// Go through file to find lcStart
	while (fin.getline(line,DomainConstants::MAXLINESIZE)){
		if(strcmp(line,lcStart)==0){
			break;
		}
	}
	
	// Check if lcStart was found:
	if(strcmp(line,lcStart)!=0){
			FEI_COUT << "WARNING: NO LOADCASES IN INPUT FILE" << FEI_ENDL;
			return 0;}
			
	
	// Declare variables
	int gid; // global id 
	int spc, mpc, load, pload, temp; // SPCset, MPCset, LOADset, PLOADset, TEMPset
	
	// Loop through the LoadCases, and read the line
	while(fin.getline(line,DomainConstants::MAXLINESIZE)){
		// Check for lcEnd string
		if(strcmp(line,lcEnd)==0){
			break;
		}
		
		// If break was not applied, another LoadCase is defined:
		std::stringstream parse(line);
		if(parse >> gid >> spc >> mpc >> load >> pload >> temp);
		else{FEI_COUT << " ERROR: IN LOADCASE DESCRIPTION FORMAT" << FEI_ENDL;
			 return 1;}
			 
			 
		// Check if the input spc,mps,load,temp have been defined
		if(spc>0 && domain.SPCList.find(spc) == domain.SPCList.end()) {
		   FEI_COUT << " ERROR: SPC " << spc << " IN LOADCASE " 
					<< gid << " UNDEFINED" << FEI_ENDL;
			return 1;}	 
		if(mpc>0 && domain.MPCList.find(mpc) == domain.MPCList.end()) {
		   FEI_COUT << " ERROR: MPC " << mpc << " IN LOADCASE " 
					<< gid << " UNDEFINED" << FEI_ENDL;
			return 1;}	
		if(load>0 && domain.LOADList.find(load) == domain.LOADList.end()) {
		   FEI_COUT << " ERROR: LOAD " << load << " IN LOADCASE " 
					<< gid << " UNDEFINED" << FEI_ENDL;
			return 1;}	 
		if(pload>0 && domain.PLOADList.find(pload) == domain.PLOADList.end()) {
		   FEI_COUT << " ERROR: PLOAD " << pload << " IN LOADCASE " 
					<< gid << " UNDEFINED" << FEI_ENDL;
			return 1;}	 
		if(temp>0 && domain.TEMPList.find(temp) == domain.TEMPList.end()) {
		   FEI_COUT << " ERROR: TEMP " << temp << " IN LOADCASE " 
					<< gid << " UNDEFINED" << FEI_ENDL;
			return 1;}	 

		// Create new LoadCase and add to list:
		domain.LoadCaseList.insert(
			std::pair<const int, LoadCase*>(
				gid, new LoadCase(spc,mpc,load,pload,temp)));
		
	}
	
	return 0;
}


// Read file with solver parameter
int read_param_file(const char* filename, MPI_Comm comm,
							std::vector<std::string>& file_contents){
  int localProc =0;
#ifndef FEI_SER
  int numProcs = 1;
  MPI_Comm_rank(comm, &localProc);
  MPI_Comm_size(comm, &numProcs);
#endif

  if (localProc == 0) {
    read_file_lines_into_strings(filename, file_contents);

#ifndef FEI_SER
    if (numProcs > 1) {
      int num = file_contents.size();

      MPI_Bcast(&num, 1, MPI_INT, 0, comm);

      for(int j=0; j<num; ++j) {
        const char* cptr = file_contents[j].c_str();
        int length = strlen(cptr) + 1;
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        MPI_Bcast((void*)cptr, length, MPI_CHAR, 0, comm);
      }
    }
#endif
  }
#ifndef FEI_SER
  else {//localProc != 0
    int num = -1;
    MPI_Bcast(&num, 1, MPI_INT, 0, comm);

    file_contents.resize(0);

    for(int j=0; j<num; ++j) {
      int length = 0;
      MPI_Bcast(&length, 1, MPI_INT, 0, comm);
      char* newstring = new char[length];
      MPI_Bcast(newstring, length, MPI_CHAR, 0, comm);
      std::string strg(newstring);
      file_contents.push_back(strg);
      delete [] newstring;
    }
  }
#endif
}

// Read file lines, containing a key and value, into a string vector;
void read_file_lines_into_strings(const char* filename,
				   std::vector<std::string>& file_contents){
  FEI_IFSTREAM infile(filename);
  if (!infile) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei_test_utils::read_file_lines_into_strings ERROR, couldn't open file '"
         << filename << "'";
    throw std::runtime_error(osstr.str());
  }

  file_contents.clear();

  std::string line;

  getline(infile, line);
  while(!infile.eof()) {
    file_contents.push_back(line);

    line.resize(0);
    getline(infile, line);
  }
}


// Print array
//template <class T>
void PrintArray(double* array, int size){
	for(int i=0; i<size; i++){
		FEI_COUT << array[i] << "  ";
	}
	FEI_COUT << FEI_ENDL;
}
void PrintArray(int* array, int size){
	for(int i=0; i<size; i++){
		FEI_COUT << array[i] << "  ";
	}
	FEI_COUT << FEI_ENDL;
}
void PrintArray(double* array, int size,std::ofstream& fout){
	for(int i=0; i<size; i++){
		fout << array[i] << "  ";
	}
	fout << FEI_ENDL;
}
void PrintArray(int* array, int size,std::ofstream& fout){
	for(int i=0; i<size; i++){
		fout << array[i] << "  ";
	}
	fout << FEI_ENDL;
}

// Print matrix of size nrow*ncol
void PrintMatrix(double* a, int nrow,int ncol){
	if(ncol==-1)ncol=nrow;
	
	// Loop through rows
	for (int i=0;i<nrow;i++){
		// Loop through columns
		for (int j=0; j<ncol;j++){
			FEI_COUT << a[j*nrow+i] << "   ";
		}
		FEI_COUT<<FEI_ENDL;
	}
}
void PrintMatrix(double* a, int nrow,int ncol,std::ofstream& fout){
	if(ncol==-1)ncol=nrow;
	
	fout << std::scientific;
	fout << std::setprecision(DomainConstants::OUTPUT_PRECISION);
	fout << std::right; // right float
	
	// Loop through rows
	for (int i=0;i<nrow;i++){
		// Loop through columns
		for (int j=0; j<ncol;j++){
			fout<< std::setw(DomainConstants::OUTPUT_DOUBLE_FIELD_WIDTH)
			    << a[j*nrow+i];
		}
		fout<<FEI_ENDL;
	}
}




}
