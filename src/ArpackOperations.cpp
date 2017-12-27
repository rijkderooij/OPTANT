#include "../include/PCH_OPTANT.h"

#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include <fei_Trilinos_Helpers.hpp>


// Constructor
ArpackOperations::ArpackOperations(feiFactoryPtr fact, feiMatrixPtr stiffMat, feiMatrixPtr geomStiffMat,int ndof, bool verbose){
	ndof_ = ndof;
	StiffMatIsPosDef_ = false;  // Initialize the positive definiteness to be false
	mpcPresent_ = true;			// Initialize the presence of mpc's to be true
	
	StiffnessMatrix_ = stiffMat;
	GeometricStiffnessMatrix_ = geomStiffMat;
	
	Factory_ = fact;
	
	// Assign x_ and y_
	x_ = Factory_->createVector(StiffnessMatrix_->getMatrixGraph(),true);
	y_ = Factory_->createVector(StiffnessMatrix_->getMatrixGraph());
	
	// Start Cholmod
	cholmod_start(&Chol_);
	
	if(verbose){
		
		FEI_COUT << "Factorizing Global Effective Stiffness Matrix..." << FEI_ENDL;
		
		// Get StiffnessMatrix_ as an epetra_Crsmatrix												
		EpetraStiffnessMatrix_ = Trilinos_Helpers::get_Epetra_CrsMatrix(StiffnessMatrix_.get());
		EpetraGeometricStiffnessMatrix_ = Trilinos_Helpers::get_Epetra_CrsMatrix(GeometricStiffnessMatrix_.get());
		EpetraSolution_ = Trilinos_Helpers::get_Epetra_MultiVector(x_.get(), true);
		EpetraForce_    = Trilinos_Helpers::get_Epetra_MultiVector(y_.get(), false);
		
		// Write EpetraStiffnessMatrix_ in CCS format
		int nrows = EpetraStiffnessMatrix_->NumGlobalRows();
		if(nrows==ndof_) mpcPresent_ = false;
		int nnz = EpetraStiffnessMatrix_->NumGlobalNonzeros();  // Total number of nonzeros
		double* AThis = new double[nnz];
		int* iThis = new int[nnz];
		int* pThis = new int[nrows+1];
		
		// Loop through the rows of EpetraStiffnessMatrix_ and fill these vectors
		int count = 0;
		int numEntries = 0;
		int length = nnz;
		
		for (int row = 0; row<nrows; row++){
			EpetraStiffnessMatrix_->ExtractGlobalRowCopy(row,length,numEntries,&AThis[count],&iThis[count]);
			pThis[row] = count;
			count+=numEntries;
			length-=numEntries;
		}
		pThis[nrows] = count;	
	
		// Start cholmod
		cholmod_sparse *A ;
		cholmod_dense *x;
		
		// Allocate space for stiffness matrix
		A = cholmod_allocate_sparse(nrows,nrows,nnz,true,true,-1,CHOLMOD_REAL,&Chol_);
		
		// Fill A by setting the pointers	
		A->i = iThis;
		A->p = pThis;
		A->x = AThis;
		
		// Factorize stiffness matrix and store in CholFactorization_
		CholFactorization_ = cholmod_analyze (A, &Chol_) ; 					// analyze
		cholmod_factorize (A, CholFactorization_, &Chol_) ; 					// factorize
		
		// Allocate space for CholSolution_ and CholForce_ 
		CholForce_    = cholmod_allocate_dense (nrows, 1, nrows,A->xtype, &Chol_); 
		CholSolution_ = cholmod_allocate_dense (nrows, 1, nrows,A->xtype, &Chol_); 
		
		/*
		// Print StiffnessMatrix_, StiffnessMatrixNoMPC_, and GeometricStiffnessMatrixNoMPC_
		double* coefs = new double[ndof_];
		int* ind = new int[ndof_];
		
		for(std::size_t id=0;id<ndof_;id++){coefs[id]=0; ind[id]=-1;}
		FEI_COUT << FEI_ENDL<< FEI_ENDL;
		for(std::size_t row=0; row<nrows;row++){
			EpetraStiffnessMatrix_->ExtractGlobalRowCopy(row,ndof_,numEntries,coefs,ind);
			
			for(std::size_t id=0;id<numEntries;id++){
				if(coefs[id]< -0.01 || coefs[id]> 0.01)
				FEI_COUT << "K(" <<row<<"," <<ind[id] << ") = " << coefs[id] << "    ";
				
			}
		}
		FEI_COUT << FEI_ENDL<< FEI_ENDL;
		for(std::size_t row=0; row<nrows;row++){
			
			EpetraGeometricStiffnessMatrix_->ExtractGlobalRowCopy(row,ndof_,numEntries,coefs,ind);
			
			for(std::size_t id=0;id<numEntries;id++){
				if(coefs[id]< -0.01 || coefs[id]> 0.01)
				FEI_COUT << "Kg(" <<row<<"," <<ind[id] << ") = " << coefs[id] << "    ";
				
			}
		}
		FEI_COUT << FEI_ENDL<< FEI_ENDL;
		//*/
		
		FEI_COUT << " Done... " << FEI_ENDL;
		
	}

}

// Destructor
ArpackOperations::~ArpackOperations(){
	cholmod_finish(&Chol_);
}
	
// Operator B : w= StiffnessMatrix_*v
void ArpackOperations::OpB(double* v, double*w){
	// Put v into x_
	for(int id=0;id<ndof_;id++){
		x_->copyIn(1, &id, &v[id], 0);
	}
	
	// Put x_ as the EpetraSolution_ vector
	EpetraSolution_ = Trilinos_Helpers::get_Epetra_MultiVector(x_.get(), true);
		
	// Compute EpetraForce_ = EpetraStiffnessMatrix_*EpetraSolution_										
	if(EpetraStiffnessMatrix_->Apply(*EpetraSolution_,*EpetraForce_)!=0) 
		FEI_COUT << "ERROR: FUNCTION EpetraStiffnessMatrix_ FAILED " << FEI_ENDL;
		
	// Set EpetraForce_ into y_
	set_Epetra_MultiVector(y_.get(),EpetraForce_);
	
	// Put y_ into w
	for(int id=0;id<ndof_;id++){
		y_->copyOut(1, &id, &w[id], 0);
	}

}


void ArpackOperations::OpA(double* v){
	// Put v into x_
	for(int id=0;id<ndof_;id++){
		x_->copyIn(1, &id, &v[id], 0);
	}
	
	// Put x_ as the EpetraSolution_ vector
	EpetraSolution_ = Trilinos_Helpers::get_Epetra_MultiVector(x_.get(), true);
	
	// Compute EpetraForce_ = EpetraGeometricStiffnessMatrix_*EpetraSolution_	
	if(EpetraGeometricStiffnessMatrix_->Apply(*EpetraSolution_,*EpetraForce_)!=0) 
		FEI_COUT << "ERROR: FUNCTION EpetraGeometricStiffnessMatrix_ FAILED " << FEI_ENDL;
	
	// Set EpetraForce_ into y_
	set_Epetra_MultiVector(y_.get(),EpetraForce_);																	
}


// Operator BiA: w = inv(B)*A*v
void ArpackOperations::OpBiA(double* v, double* w){
	// Compute A*v and store in y_, with A=GeometricStiffnessMatrix_
	OpA(v);
	
	// Set y_ as force vector of Cholmod system
	EpetraForce_->ExtractCopy((double*)CholForce_->x,EpetraForce_->GlobalLength());
	
	// Solve StiffnessMatrix_*CholSolution_ = CholForce_
	cholmod_free_dense(&CholSolution_,&Chol_);				// First clear the solution
	CholSolution_ = cholmod_solve(CHOLMOD_A, CholFactorization_, CholForce_, &Chol_); 
	
	if(!mpcPresent_){		
		// Put CholSolution_ into w
		for(int id=0;id<ndof_;id++){
			w[id] = ((double*)CholSolution_->x)[id];
		}
		
	}
	else{
		// Set CholSolution_ to EpetraForce_
		for(int id=0;id<CholSolution_->nrow;id++){
			EpetraForce_->ReplaceGlobalValue(id,0,((double*)CholSolution_->x)[id]);
		}
			
		// Set EpetraForce_ to y_
		set_Epetra_MultiVector(y_.get(),EpetraForce_);
		
		// Put y_ into w
		for(int id=0;id<ndof_;id++){
			y_->copyOut(1, &id, &w[id], 0);
		}
	}
}

// Set EpetraMultivector into fei::Vector. Based on get_Epetra_MultiVector in fei_Trilinos_Helpers.hpp
void ArpackOperations::set_Epetra_MultiVector(fei::Vector* feivec, Epetra_MultiVector* epvec)
{
  fei::Vector* vecptr = feivec;
  fei::VectorReducer* feireducer = dynamic_cast<fei::VectorReducer*>(feivec);
  if (feireducer != NULL) vecptr = feireducer->getTargetVector().get();

  fei::Vector_Impl<Epetra_MultiVector>* fei_epetra_vec =
    dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>*>(vecptr);

  if (fei_epetra_vec == NULL) {
    throw std::runtime_error("failed to obtain Epetra_MultiVector from fei::Vector.");
  }

  if (fei_epetra_vec != NULL) {
    fei_epetra_vec->setUnderlyingVector(epvec);
  }
}
