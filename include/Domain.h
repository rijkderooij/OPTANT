#ifndef DOMAIN_H
#define DOMAIN_H


#include "fei_iostream.hpp" 
#include "fei_base.hpp"

// Forward declarations
class Material;
class Section;
class Property;
class Node;
class Element;
class LoadCase;
class LOAD;
class PLOAD;
class SPC;
class MPC;
class TEMP;
class OutputRequest;
class Solver;
class LinearBucklingSolver;


class Domain{
	
	public:
		enum SolType{ 
				LinearStatic=10,
				LinearBuckling=11
		};
	
		// Constructor
		Domain();
		Domain(const char* filename);
		
		// Functions
		int ReadInputFile();
		int BuildElementMatrices();
		int BuildGeometricElementMatrices(Solver* solver, bool makeNegative = false);
		int AssembleStiffnessMatrix(fei::MatrixGraph* matrixGraph,
									fei::Matrix* mat, bool useGeomStiffness = false);
		int InitializeElementConnectivities(fei::MatrixGraph* matrixGraph);	
		int InitializeConstraints(fei::MatrixGraph* matrixGraph, int globalLoadCaseID);				
		int ApplyConstraints(fei::LinearSystem* linSys, int globalLoadCaseID);
		int ApplyLoads(fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec, int globalLoadCaseID);
		int ResetElementTemperatures();
		int WriteOutput(Solver* solver, const char* outputFile, int globalLoadCaseID);
		int WriteLinearBucklingOutput(LinearBucklingSolver* solver, const char* outputFile);
		
		// Lists
		std::map<int, Material*> MaterialList;
		std::vector<Section*> SectionList;
		std::vector<Property*> PropertyList;
		std::vector<Node*> NodeList;
		std::vector<Element*> ElementList;	
		std::map<const int, LoadCase*> LoadCaseList;
		std::multimap<const int, SPC*>SPCList;
		std::multimap<const int, MPC*> MPCList;
		std::multimap<const int, LOAD*> LOADList;
		std::multimap<const int, PLOAD*> PLOADList;
		std::multimap<const int, TEMP*> TEMPList;
		
		
		// Conversions from local to global ID and vice versa
		std::vector<int> NodeID_LG, ElementID_LG, SectionID_LG, PropertyID_LG;
		std::map<int,int> NodeID_GL, ElementID_GL, SectionID_GL, PropertyID_GL;
		
		// Number of DOF per node
		std::vector<int> NodalNumDOF;
		
		// Vectors required to define patterns and blocks
		std::vector<int> NodesPerElementBlock; // Length=Num Total Blocks
		std::vector<int> NumElementsBlock;   // Length=Num Total Blocks
		std::vector<int> ElementBlockID;     // Length=NumTotalElements
		int ElemMatPackStorageSize; // Number of entries in all element matrices in packed storage
		std::vector<double> AllElementMatrices; // Large vector containing all element matrices in packed storage
		std::vector<double> AllGeometricElementMatrices; // Large vector containing all geometric element matrices in packed storage
		
		// Solver parameters
		fei::ParameterSet SolverParams;
		
		// OutputRequest
		OutputRequest* OutputReq;
		
		// General variables
		Domain::SolType SolutionType;
		int NumTotalElements, NumLocalElements;
		int NumTotalNodes;
		bool ElementMatricesBuilt;
		
		
	private:
		// Filename
		const char* FileName;
		
		
		// Functions
		
		// Function to add append tempfiles to output file
		int HandleTempFiles(std::ofstream& fout, std::ofstream* tempFiles, const int& numTempFiles, const char* tempFilePath);
		
	
	
};
#endif
