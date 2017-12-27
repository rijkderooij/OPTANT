#ifndef ELEMENT_H
#define ELEMENT_H

#include "fei_iostream.hpp" 

#include "Property.h"

//Forward declarations
class Domain;
class Node;
class OutputRequest;
class Solver;

class Element{
	
	public:
		enum Type{ 
				Rod=0,
				Beam=1,
				Shell=2,
				Solid=3
		};
		
		Element(Element::Type et);
		~Element();
		
		// StiffnessMatrix
		virtual int BuildElementMatrix() = 0;
		virtual int BuildGeometricElementMatrix(Solver *solver, bool makeNegative = false) =0;
		virtual void ElementMatrix(double*& stiffMat, int& size)=0;
		virtual void GeometricElementMatrix(double*& geomMat, int& size)=0;
		virtual int TemperatureLoadVector(double* loadVector)=0;
		
		// Interpolate element properties to e.g. gauss points
		virtual int InterpolateProperties(const int nProp, double* interpMat) = 0;
		
		// Functions for element loads (i.e. PLOAD)
		virtual int NormalDirection(double *n) = 0;
		virtual int ElementArea(double& area) = 0;
		
		
		
		// General Functions
		int LocalElementID(){return LocalElementID_;};
		int NumProperties() {return NumProperties_;};
		Property** Properties() {return Properties_;};
		int NumDeltaTemperatures() {return NumDeltaTemperatures_;};
		double* DeltaTemperatures() {return DeltaTemperatures_;};
		virtual Node** Nodes() = 0;
		
		void SetElementMatrixLocation(double* loc){ElementMatrix_=loc;};
		void SetGeometricElementMatrixLocation(double* loc){GeometricElementMatrix_=loc;};
		void SetNumProperties(int np);
		void SetNumDeltaTemperatures(int nt);
		void ResetDeltaTemperatures();
		void SetLocalElementID(int id){LocalElementID_=id;};
		int AddProperty(int index, int GlobalPropID, Domain& domain);
		int AddDeltaTemperature(int index, double dTemp);
		
		virtual int ReadFromFile(std::ifstream& fin, Domain& domain) = 0;
		int ReadPropertyLines(std::ifstream& fin, Domain& domain, Property::Type propType);
		
		virtual int WriteOutput(Domain& domain, Solver *solver,
		                     std::ofstream* tempFiles, bool doTitle=false)=0;
		
		Element::Type ElementType(){return ElementType_;};
		
		
	protected:
	
		double* ElementMatrix_;  			// Pointer to the ElementMatrix in packed storage on the heap
		double* GeometricElementMatrix_;  	// Pointer to the GeometricElementMatrix in packed storage on the heap
		
		int LocalElementID_;
		int NumProperties_; // Total number of properties in the element
		Property** Properties_; // Array of property pointers
		int NumDeltaTemperatures_; // Total number of applied delta temperatures in the element (0 by default, changes per load case) 
		double* DeltaTemperatures_; // Applied delta temperatures 
		
		Element::Type ElementType_;
		

};
#endif
