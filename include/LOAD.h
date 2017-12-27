#ifndef LOAD_H
#define LOAD_H

#include "fei_base.hpp"

// Forward declarations:
class Domain;
class Node;

class LOAD{
	
	public:
		LOAD();
		LOAD(Node* node, double fx, double fy, double fz, double mx,
			 double my, double mz);
		LOAD(Node* node, double* loadvec);
		
		int Initialize(Domain& domain){return 0;};
		int Apply(fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec);
		
		// Getters
		Node* GetNode(){return Node_;};
		double* LoadVector(){return LoadVector_;};
		
	private:
		// Variables
		Node* Node_;
		double LoadVector_[6];
		
	
};
#endif
