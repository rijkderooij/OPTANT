#ifndef PLOAD_H
#define PLOAD_H

#include "fei_base.hpp"

// Forward declarations:
class Domain;
class Element;

class PLOAD{
	
	public:
		PLOAD();
		PLOAD(int numPressElements, Element** pressElem, double pressure);
		~PLOAD();
		
		int Initialize(Domain& domain){return 0;};
		int Apply(Domain& domain, fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec);
		
		// Getters
		int NumPressElements(){return NumPressElements_;};
		Element**  PressElements(){return PressElements_;};
		double Pressure(){return Pressure_;};
		
	private:
		// Variables
		int NumPressElements_;				// Number of elements for this PLOAD
		Element** PressElements_;
		double Pressure_;
		
	
};
#endif
