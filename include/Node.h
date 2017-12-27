#ifndef NODE_H
#define NODE_H

#include "fei_iostream.hpp" 

// Forward Declarations
class Domain;
class OutputRequest;
class Solver;
class LinearBucklingSolver;

class Node{
	
	public:
		// Constructor
		Node();
		Node(double x, double y, double z);
		Node(double* globCoord);
		
		// Functions
		int WriteOutput(Domain& domain, Solver* solver, std::ofstream* tempFiles, bool doTitle=false);
		int WriteLinearBucklingOutput(Domain& domain, LinearBucklingSolver* solver, std::ofstream* tempFiles, bool doTitle=false);
		int Displacements(Solver* solver, double* DOFValues);		// Extract nodal displacements
		int BucklingModeDisplacements(LinearBucklingSolver* solver, int mode, double* DOFValues);	// Extract buckling mode displacements
		
		// Getters
		int LocalNodeID(){return LocalNodeID_;};
		double* GlobalCoordinates() {return GlobalCoordinates_;};
		
		// Setters
		void SetLocalNodeID(int id){LocalNodeID_=id;};
		void SetGlobalCoordinates(double x, double y, double z);
	
	private:
		double GlobalCoordinates_[3];
		int LocalNodeID_;
	
};
#endif
