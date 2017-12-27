#ifndef CBEAM_H
#define CBEAM_H

#include "Element.h"

#include "fei_iostream.hpp" 

//Forward declarations
class Domain;
class Node;
class OutputRequest;
class Solver;


class CBEAM : public Element{
	
	public:		
		// Constructor
		CBEAM();
		~CBEAM();
		
		// StiffnessMatrix
		int BuildElementMatrix();
		int BuildGeometricElementMatrix(Solver *solver, bool makeNegative = false);
		void ElementMatrix(double*& stiffMat, int& size);
		void GeometricElementMatrix(double*& geomMat, int& size);
		int TemperatureLoadVector(double* loadVector);
		
		// Functions for element loads (i.e. PLOAD)
		int NormalDirection(double *n){return 0;};
		int ElementArea(double& area){return 0;};
		
		// Functions
		int ReadFromFile(std::ifstream& fin, Domain& domain);
		int WriteOutput(Domain& domain, Solver *solver,
		                     std::ofstream* tempFiles, bool doTitle=false);
		
		
		// Getters
		Node** Nodes(){return Nodes_;};
		double* ZDirection(){return ZDirection_;};
		
		
	private:
		//Variables
		Node* Nodes_[2];
		double ZDirection_[3];
		
		// Functions
		
		// Compute material stiffness matrices integrated over element from properties
		int MaterialMatrices(double* C00,double* C01,double* C11,double* G0,double* G1,double* S,double* c0, double* c1, double* g0);
		
		// Interpolate element properties to gauss points
		int InterpolateProperties(const int nProp, double* interpMat);
		
		// Compute element forces at gauss points
		int ElementForces(double* elemForces, Solver *solver);
		
		// Compute element strain energy in beam element
		int ElementStrainEnergy(double& elemEnergy, Solver *solver);
		
		
	
};
#endif
