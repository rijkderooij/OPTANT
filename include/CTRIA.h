#ifndef CTRIA_H
#define CTRIA_H

#include "Element.h"

#include "fei_iostream.hpp" 

//Forward declarations
class Domain;
class Node;
class OutputRequest;
class Solver;


class CTRIA : public Element{
	
	public:		
		// Constructor
		CTRIA();
		~CTRIA();
		
		// StiffnessMatrix
		int BuildElementMatrix();
		int BuildGeometricElementMatrix(Solver *solver, bool makeNegative = false);
		void ElementMatrix(double*& stiffMat, int& size);
		void GeometricElementMatrix(double*& geomMat, int& size);
		int TemperatureLoadVector(double* loadVector);
		
		// Calculate membranal and bending B matrices according to Felippa
		// These functions are static, since they are also called by the CQUAD element
		static int FelippaTriMembrane(double* storage, bool forStiffnessMatrix);
		static void FelippaTriBending(double* storage);
		
		// Functions for element loads (i.e. PLOAD)
		int NormalDirection(double *n);
		int ElementArea(double& area);
		
		// Functions
		int ReadFromFile(std::ifstream& fin, Domain& domain);
		int WriteOutput(Domain& domain, Solver* solver,
		                     std::ofstream* tempFiles, bool doTitle=false);
		
		
		// Getters
		double* MaterialDirection() {return MaterialDirection_;};
		Node** Nodes(){return Nodes_;};
		
		
	private:
		//Variables
		double MaterialDirection_[3]; // 1-direction of element property
		Node* Nodes_[3];
		
		// Functions
		
		// Verify if the MaterialDirection of all properies is not normal to
		// this element
		int VerifyMaterialDirection(int GlobalElemID);
		
		// Compute material stiffness matrices integrated over element from properties
		int MaterialMatrices(double* Am,double* Bm,double* Dm);
		int MaterialTemperatureMatrices(double* am,double* bm, double& dTemp);
		
		// Interpolate element properties to gauss points
		int InterpolateProperties(const int nProp, double* interpMat);
		
		// Compute element strains at center of element
		int ElementStrains(double* elemStrains, Solver* solver);
		
		// Compute element forces at center of element
		int ElementForces(double* elemForces, const double* elemStrains);
		
		// Compute element strain energy at center of triangle
		int ElementStrainEnergy(double& elemEnergy, double* elemForces, const double* elemStrains);
	
};
#endif
