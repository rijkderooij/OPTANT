#ifndef ARPACKOPERATIONS_H
#define ARPACKOPERATIONS_H

#include "fei_base.hpp"
#include "cholmod.h"
class Epetra_MultiVector;
class Epetra_CrsMatrix;

typedef fei::SharedPtr<fei::Factory> feiFactoryPtr;
typedef fei::SharedPtr<fei::Matrix> feiMatrixPtr;
typedef fei::SharedPtr<fei::Vector> feiVectorPtr;

class ArpackOperations
{
    public:
        ArpackOperations(feiFactoryPtr fact, feiMatrixPtr stiffMat, feiMatrixPtr geomStiffMat, int ndof,bool verbose);
		~ArpackOperations();
		
		
		// Operator B : w= B*v
		void OpB(double* v, double*w);
		
		// Operator A : w= A*v
		void OpA(double* v);
		
		// Operator BiA: w = inv(B)*A*v
		void OpBiA(double* v, double* w);
		
		// Set EpetraMultivector into fei::Vector. Based on get_Epetra_MultiVector in fei_Trilinos_Helpers.hpp
		void set_Epetra_MultiVector(fei::Vector* feivec, Epetra_MultiVector* epvec);
		
		// Getters
		bool StiffMatIsPosDef(){return StiffMatIsPosDef_;};
		
    
    private:
		feiMatrixPtr StiffnessMatrix_;
		feiMatrixPtr GeometricStiffnessMatrix_;
		
		feiVectorPtr x_;
		feiVectorPtr y_;
		
		feiFactoryPtr Factory_;
		
		// Epetra parameters
		Epetra_CrsMatrix* EpetraStiffnessMatrix_;
		Epetra_CrsMatrix* EpetraGeometricStiffnessMatrix_;
		Epetra_MultiVector* EpetraSolution_;
		Epetra_MultiVector* EpetraForce_;
		
		// Cholmod parameters
		cholmod_common  Chol_;	 			 // Cholmod
		cholmod_factor *CholFactorization_;	 // Factorization of stiffnessMatrix
		cholmod_dense  *CholSolution_;		 // Solution vector of Cholmod system
		cholmod_dense  *CholForce_;			 // Force vector for Cholmod system
		bool StiffMatIsPosDef_;				 // Boolean to check whether stiffmat is positive definite
		
		
		int ndof_;
		bool mpcPresent_;
	
	



};

#endif // ARPACKOPERATIONS_H

