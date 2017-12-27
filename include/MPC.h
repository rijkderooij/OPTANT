#ifndef MPC_H
#define MPC_H

#include "fei_base.hpp"
#include <unordered_map>


class Node;

class MPC{
	
	public:
		MPC();
		MPC(int numNodes, Node** mpcNodes, int* prescribedDOF, double* weights, double rhsConstant);
		~MPC();
		
		int Initialize(fei::MatrixGraph* matrixGraph);
		
		// Getters
		int NumNodes(){return NumNodes_;};
		Node**  MPCNodes(){return MPCNodes_;};
		int* PrescribedDOF(){return PrescribedDOF_;};
		double* Weights(){return Weights_;};
		double RHSConstant(){return RHSConstant_;};
		
	private:
		// Variables
		int NumNodes_; 							// Number of nodes (including slave node)
		Node** MPCNodes_;
			int* PrescribedDOF_;
		double* Weights_;
		double RHSConstant_;
		
		// Derived variables
		int* nodeIDs_; 							// Stores all nodeIDs
		std::vector<int> uniqueNodeIDs_;        // Store the set of unique nodeIDs
		std::unordered_map<int,int> nodeMap_;	// Relates the nodeIDs to the index at which they are located in uniqueNodeIDs;
		int uniqueNodeIDsCounter_;				// Counter for the number of unique nodeIDs
		
	
	
};
#endif
