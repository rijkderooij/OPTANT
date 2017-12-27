#ifndef TEMP_H
#define TEMP_H

// Forward declarations:
class Domain;
class Element;

class TEMP{
	
	public:
		TEMP();
		TEMP(Element* element, int NumTemperatures, double* deltaTemp);
		
		int Initialize(Domain& domain){return 0;};
		int Apply(Domain& domain, fei::MatrixGraph* matrixGraph,fei::Vector* rhsVec);
		
		// Getters
		Element* GetElement(){return Element_;};
		int NumDeltaTemperatures(){return NumDeltaTemperatures_;};
		double* DeltaTemperatures(){return DeltaTemperatures_;};
		
	private:
		// Variables
		Element* Element_;
		int NumDeltaTemperatures_;
		double* DeltaTemperatures_;
};
#endif
