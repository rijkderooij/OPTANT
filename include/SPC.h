#ifndef SPC_H
#define SPC_H

#include "fei_base.hpp"


class Domain;
class Node;

class SPC{
	
	public:
		SPC();
		SPC(Node* spcNode,const char* preDOF, double preVal);
		
		int Initialize(fei::MatrixGraph* matrixGraph);
		int Apply(fei::LinearSystem* linSys);
		
		// Getters
		Node*  GetNode(){return Node_;};
		int* PrescribedDOF(){return PrescribedDOF_;};
		double  PrescribedValue(){return PrescribedValue_;};
		
	private:
		// Variables
		Node* Node_;
		int PrescribedDOF_[6];
		double PrescribedValue_;
	
	
};
#endif
