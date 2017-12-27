#ifndef PSHELL_H
#define PSHELL_H

#include "Property.h"

#include "fei_iostream.hpp" 

//Forward declarations
class Domain;

class PSHELL : public Property{
	
	public:		
		// Constructor
		PSHELL();
		
		// Functions
		int ReadFromFile(std::ifstream& fin, int PropOption,
						   int NumLines, Domain& domain);
		
		// Getters		
		double* A() {return A_;};
		double* B() {return B_;};
		double* D() {return D_;};
		
		double* a() {return a_;};
		double* b() {return b_;};
		
		
		
		
	private:
		//Variables
		
		double A_[6]; // 3x3 A matrix in packed storage
		double B_[6]; // 3x3 B matrix in packed storage
		double D_[6]; // 3x3 D matrix in packed storage
		
		double a_[3]; // 3x1 a vector s.t. [Nx;Ny;Nxy]-=a*dTemp
		double b_[3]; // 3x1 b vector s.t. [Mx;My;Mxy]-=b*dTemp
		
		
	
};
#endif
