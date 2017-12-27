#include "../include/PCH_OPTANT.h"

#include "fei_iostream.hpp"
#include <Epetra_BLAS.h>
#include <Epetra_LAPACK.h>

namespace MatrixOperations {
// Variables
Epetra_BLAS blas;
Epetra_LAPACK lapack;




// BLAS and LaPACK operations customized for own use

// Cholesky decomposition of A(4x4) or(2x2) matrix A=L*L**T, A=spd
// L may be stored on position of A
void dpotrf4(const double *A, double* L){
	L[0] = sqrt(A[0]);
	L[1] = A[1]/L[0];
	L[2] = A[2]/L[0];
	L[3] = A[3]/L[0];
	L[5] = sqrt(A[5]-L[1]*L[1]);
	L[6] = (A[6]-L[2]*L[1])/L[5];
	L[7] = (A[7]-L[3]*L[1])/L[5];
	L[10] = sqrt(A[10]-L[2]*L[2]-L[6]*L[6]);
	L[11] = (A[11]-L[3]*L[2]-L[7]*L[6])/L[10];
	L[15] = sqrt(A[15]-L[3]*L[3]-L[7]*L[7]-L[11]*L[11]);
}
void dpotrf2(const double *A, double* L){
	L[0] = sqrt(A[0]);
	L[1] = A[1]/L[0];
	L[3] = sqrt(A[3]-L[1]*L[1]);
}
void dpotrf(double *A, const int N, int* info){ // L is stored on location of A
	lapack.POTRF('L',N,A,N,info);
}

// B = L**(-1), with L(4x4) lower triangular matrix in conventional storage
// B may be stored on position of L
void dtrtri6(const double *L, double* B){
	B[0] = 1/L[0];
	B[1] = (-L[1]*B[0])/L[7];
	B[2] = (-L[2]*B[0] - L[8]*B[1])/L[14];
	B[3] = (-L[3]*B[0] - L[9]*B[1] - L[15]*B[2])/L[21];
	B[4] = (-L[4]*B[0] - L[10]*B[1] - L[16]*B[2] - L[22]*B[3])/L[28];
	B[5] = (-L[5]*B[0] - L[11]*B[1] - L[17]*B[2] - L[23]*B[3] - L[29]*B[4])/L[35];
	
	B[7] = 1/L[7];
	B[8] = (-L[8] *B[7])/L[14];
	B[9] = (-L[9] *B[7] - L[15]*B[8])/L[21];
	B[10]= (-L[10]*B[7] - L[16]*B[8] - L[22]*B[9])/L[28];
	B[11]= (-L[11]*B[7] - L[17]*B[8] - L[23]*B[9] - L[29]*B[10])/L[35];
	
	B[14] = 1/L[14];
	B[15] = (-L[15]*B[14])/L[21];
	B[16] = (-L[16]*B[14] - L[22]*B[15])/L[28];
	B[17] = (-L[17]*B[14] - L[23]*B[15] - L[29]*B[16])/L[35];
	
	B[21] = 1/L[21];
	B[22] = (-L[22]*B[21])/L[28];
	B[23] = (-L[23]*B[21] - L[29]*B[22])/L[35];
	
	B[28] = 1/L[28];
	B[29] = (-L[29]*B[28])/L[35];
	
	B[35] = 1/L[35];
	
}
void dtrtri4(const double *L, double* B){
	B[0] = 1/L[0];
	B[5] = 1/L[5];
	B[10] = 1/L[10];
	B[15] = 1/L[15];
	
	B[1] = -L[1]*B[0]*B[5];
	B[6] = -L[6]*B[5]*B[10];
	B[11] = -L[11]*B[10]*B[15];
	
	B[3] = -L[3]*B[0]*B[15] - L[7]*B[1]*B[15] - L[2]*B[11]*B[0] + B[1]*B[6]*B[11]/B[5]/B[10]; 
	
	B[2] = -L[2]*B[0]*B[10] + B[6]*B[1]/B[5];
	B[7] = -L[7]*B[5]*B[15] + B[11]*B[6]/B[10];
	
}
void dtrtri2(const double *L, double* B){
	B[0] = 1/L[0];
	B[3] = 1/L[3];
	
	B[1] = -L[1]*B[0]*B[3];
}

// B = A^(-1) with A(2x2)
// B may be stored on position of A
void dpotri2(const double *A, double* B){
	B[2] = A[0];
	B[0] = A[3];
	B[3] = B[2];
	B[2] = B[0]*B[3]-A[1]*A[1];
	B[0] /= B[2];
	B[3] /= B[2];
	B[1] /= -1*B[2];
	B[2] = B[1];
}

// B = L*A, with L lower triangular in conventional storage, A(4x4)
// B may be stored on position of A
void dtrmm4(const double *L, const double* A, double* B){
	B[15] = L[3]*A[12]+L[7]*A[13]+L[11]*A[14]+L[15]*A[15];
	B[14] = L[2]*A[12]+L[6]*A[13]+L[10]*A[14];
	B[13] = L[1]*A[12]+L[5]*A[13];
	B[12] = L[0]*A[12];
	B[11] = L[3]*A[8]+L[7]*A[9]+L[11]*A[10]+L[15]*A[11];
	B[10] = L[2]*A[8]+L[6]*A[9]+L[10]*A[10];
	B[9] =  L[1]*A[8]+L[5]*A[9];
	B[8] =  L[0]*A[8];
	B[7] =  L[3]*A[4]+L[7]*A[5]+L[11]*A[6]+L[15]*A[7];
	B[6] =  L[2]*A[4]+L[6]*A[5]+L[10]*A[6];
	B[5] =  L[1]*A[4]+L[5]*A[5];
	B[4] =  L[0]*A[4];
	B[3] =  L[3]*A[0]+L[7]*A[1]+L[11]*A[2]+L[15]*A[3];
	B[2] =  L[2]*A[0]+L[6]*A[1]+L[10]*A[2];
	B[1] =  L[1]*A[0]+L[5]*A[1];
	B[0] =  L[0]*A[0];
}

// B = A*L or B=L*A, with L lower triangular in conventional storage, A(2x4)
// B may be stored on position of A
void dtrmm2_4(const double* A, const double *L, double* B){  // A*L(4x4)
	B[0] = L[0]*A[0]+L[1]*A[2]+L[2]*A[4]+L[3]*A[6];
	B[2] =           L[5]*A[2]+L[6]*A[4]+L[7]*A[6];
	B[4] =                    L[10]*A[4]+L[11]*A[6];
	B[6] =                               L[15]*A[6];
	B[1] = L[0]*A[1]+L[1]*A[3]+L[2]*A[5]+L[3]*A[7];
	B[3] =           L[5]*A[3]+L[6]*A[5]+L[7]*A[7];
	B[5] =                    L[10]*A[5]+L[11]*A[7];
	B[7] =                               L[15]*A[7];
}
void dtrmm2_4L(const double *L, const double* A, double* B){ // L(2x2)*A
	B[1] = A[0]*L[1] + A[1]*L[3];
	B[0] = A[0]*L[0];
	B[3] = A[2]*L[1] + A[3]*L[3];
	B[2] = A[2]*L[0];
	B[5] = A[4]*L[1] + A[5]*L[3];
	B[4] = A[4]*L[0];
	B[7] = A[6]*L[1] + A[7]*L[3];
	B[6] = A[6]*L[0];
}

// C = A*B, A(2x4), B(4x4)
// C MAYNOT be stored on position of A or B
void dgemm2_4(const double* A, const double *B, double* C){
	C[0] = B[0]*A[0]+B[1]*A[2]+B[2]*A[4]+B[3]*A[6];
	C[2] = B[4]*A[0]+B[5]*A[2]+B[6]*A[4]+B[7]*A[6];
	C[4] = B[8]*A[0]+B[9]*A[2]+B[10]*A[4]+B[11]*A[6];
	C[6] = B[12]*A[0]+B[13]*A[2]+B[14]*A[4]+B[15]*A[6];
	C[1] = B[0]*A[1]+B[1]*A[3]+B[2]*A[5]+B[3]*A[7];
	C[3] = B[4]*A[1]+B[5]*A[3]+B[6]*A[5]+B[7]*A[7];
	C[5] = B[8]*A[1]+B[9]*A[3]+B[10]*A[5]+B[11]*A[7];
	C[7] = B[12]*A[1]+B[13]*A[3]+B[14]*A[5]+B[15]*A[7];

}

// B = L**T*A, with L lower triangular in conventional storage,
// B may be stored on position of A
void dtrmm4tL(const double *L, const double* A, double* B){ // L(4x4), A(4x4)
	B[0] = L[0]*A[0]+L[1]*A[1]+L[2]*A[2]+L[3]*A[3];
	B[1] =           L[5]*A[1]+L[6]*A[2]+L[7]*A[3];
	B[2] =                    L[10]*A[2]+L[11]*A[3];
	B[3] =                               L[15]*A[3];
	B[4] = L[0]*A[4]+L[1]*A[5]+L[2]*A[6]+L[3]*A[7];
	B[5] =           L[5]*A[5]+L[6]*A[6]+L[7]*A[7];
	B[6] =                    L[10]*A[6]+L[11]*A[7];
	B[7] =                               L[15]*A[7];
	B[8] = L[0]*A[8]+L[1]*A[9]+L[2]*A[10]+L[3]*A[11];
	B[9] =           L[5]*A[9]+L[6]*A[10]+L[7]*A[11];
	B[10] =                   L[10]*A[10]+L[11]*A[11];
	B[11] =                               L[15]*A[11];
	B[12] = L[0]*A[12]+L[1]*A[13]+L[2]*A[14]+L[3]*A[15];
	B[13] =            L[5]*A[13]+L[6]*A[14]+L[7]*A[15];
	B[14] =                     L[10]*A[14]+L[11]*A[15];
	B[15] =                                 L[15]*A[15];
}
void dtrmm2_4tL(const double *L, const double* A, double* B){ // L(2x2), A(2x4)
	B[0] = L[0]*A[0]+L[1]*A[1];
	B[1] =           L[3]*A[1];
	B[2] = L[0]*A[2]+L[1]*A[3];
	B[3] =           L[3]*A[3];
	B[4] = L[0]*A[4]+L[1]*A[5];
	B[5] =           L[3]*A[5];
	B[6] = L[0]*A[6]+L[1]*A[7];
	B[7] =           L[3]*A[7];

}


// B = A*L**T, with L lower triangular in conventional storage, A(2x4)
// B may be stored on position of A
void dtrmm2_4t(const double *A, const double* L, double* B){
	B[6] = L[3]*A[0]+L[7]*A[2]+L[11]*A[4]+L[15]*A[6];
	B[4] = L[2]*A[0]+L[6]*A[2]+L[10]*A[4];
	B[2] = L[1]*A[0]+L[5]*A[2];
	B[0] = L[0]*A[0];
	
	B[7] = L[3]*A[1]+L[7]*A[3]+L[11]*A[5]+L[15]*A[7];
	B[5] = L[2]*A[1]+L[6]*A[3]+L[10]*A[5];
	B[3] = L[1]*A[1]+L[5]*A[3];
	B[1] = L[0]*A[1];
}

// B = L**T*L, with L lower triangular in conventional storage
// B may be stored on position of L
void dtrmtrm2(const double *L, double* B){ // L(2x2)
	B[0] = L[0]*L[0]+L[1]*L[1];
	B[1] = L[1]*L[3];
	B[2] = B[1];
	B[3] = L[3]*L[3];
}


// C = A-B**T*B, with A(4x4) = A**T, i.e. full rank downdate
// C may be stored on position of A
void dsyrd4(const double* A, const double* B, double* C){
	C[0] = A[0]-(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]+B[3]*B[3]);
	C[1] = A[1]-(B[0]*B[4]+B[1]*B[5]+B[2]*B[6]+B[3]*B[7]);
	C[2] = A[2]-(B[0]*B[8]+B[1]*B[9]+B[2]*B[10]+B[3]*B[11]);
	C[3] = A[3]-(B[0]*B[12]+B[1]*B[13]+B[2]*B[14]+B[3]*B[15]);
	
	C[4] = C[1];
	C[5] = A[5]-(B[4]*B[4]+B[5]*B[5]+B[6]*B[6]+B[7]*B[7]);
	C[6] = A[6]-(B[4]*B[8]+B[5]*B[9]+B[6]*B[10]+B[7]*B[11]);
	C[7] = A[7]-(B[4]*B[12]+B[5]*B[13]+B[6]*B[14]+B[7]*B[15]);
	
	C[8] = C[2];
	C[9] = C[6];
	C[10] = A[10]-(B[8]*B[8]+B[9]*B[9]+B[10]*B[10]+B[11]*B[11]);
	C[11] = A[11]-(B[8]*B[12]+B[9]*B[13]+B[10]*B[14]+B[11]*B[15]);
	
	C[12] = C[3];
	C[13] = C[7];
	C[14] = C[11];
	C[15] = A[15]-(B[12]*B[12]+B[13]*B[13]+B[14]*B[14]+B[15]*B[15]);	
	
}

// C = A+beta*B*B**T, with A(2x2) = A**T, B(2x4)
// C may be storedd on position of A
void dsyr2_4(const double* A, const double beta, const double* B, double* C){
	C[0] = A[0] + beta*(B[0]*B[0]+B[2]*B[2]+B[4]*B[4]+B[6]*B[6]);
	C[1] = A[1] + beta*(B[0]*B[1]+B[2]*B[3]+B[4]*B[5]+B[6]*B[7]);
	C[2] = C[1];
	C[3] = A[3] + beta*(B[1]*B[1]+B[3]*B[3]+B[5]*B[5]+B[7]*B[7]);
}

// A := alpha*A+beta*B**T*B, with A(nxn) = A**T, B(kxn)
void dsyrk_nt(const double alpha, double* A, const double beta, const double* B,const int k, const int n){
	blas.SYRK('L','T',n,k,beta,B,k,alpha,A,n);
	// Only lower part is filled, i.e. fill the rest
	for(int i=0;i<n;i++){
		for(int j=0; j<i; j++){
			A[i*n+j] = A[j*n+i];
		}	
	}
}

// C = Btilde*A*Btilde**T, with Btilde being a diagonal block matrix
// with diagonal blocks B(mxm). A(nxn), m<n, n%m=0, A==A**T, C(nxn) (full update matrix)
// C may be stored on position of A
void dsyrk_f(const double* A, const double* B, double* C, const int n, const int m){
	// The computation is a two-step process, and we need to generate an
	// array to store the intermediate result.
	double* D = new double[n*n];
	for(int i=0;i<(n*n);i++)D[i]=0;
	int ki0=0, li0=0,index=0;
	
	// Compute D = Btilde*A
	for(int i=0; i<n; i++){
		ki0 = i/m*m;
		li0 = i%m;
		index = i;
		for(int j=0; j<n; j++){
			for(int ij=0; ij<m; ij++){
				D[index]+= A[ki0+ij]*B[li0+ij*m];
			}
			ki0+=n;
			index+=n;
		}
	}
	
	for(int i=0;i<(n*n);i++)C[i]=0;
	// Compute C = D*Btilde**T
	for(int j=0; j<n; j++){
		ki0 = j/m*m*n+j;
		for(int i=j; i<n; i++){
			index = i+n*j;
			li0 = j%m;
			for(int ij=0; ij<m; ij++){
				C[index]+= D[ki0+ij*n]*B[li0+ij*m];
			}
			ki0++;
		}
	}
	
	// Compute symmetric part
	for(int j=0; j<n; j++){
		for(int i=0; i<j; i++){
			C[i+n*j] = C[n*i+j];
		}
	}	
	
	delete [] D;
}

// C := alpha*A*B; A(mxp), B(pxn), C(mxn);
// C MAYNOT be stored on position of A or B
void dgemm(const double ALPHA, const double *A, const double *B, double *C, const int M, const int P, const int N){
	blas.GEMM('N','N',M,N,P,ALPHA,A,M,B,P,0,C,M);
}


// K:= K + alpha* B*A*C**T, with A(mxm)==A**T; B(nxm),C(nxm) in conv. storage
// K MAYNOT be stored on either position A,B,C
void dsymm_bact(const double alpha, const double* B, const double *A, const double *C, double *K, const int M, const int N){
	int indexA = 0;
	for(size_t m1=0; m1<M; m1++){
		for(size_t m2=0; m2<M; m2++){
			dger_n(A[indexA]*alpha,&B[m1*N],&C[m2*N],K,N,N);
			indexA++;
		}		
	}
}


// Matrix vector multiplication

// y = beta*y + alpha*A*x;  A(mxn)
void dgemv(const double beta, double *y, const double alpha, const double* A, const double* x, const int m, const int n){
	blas.GEMV('N',m,n,alpha,A,m,x,beta,y);
}

// y = x**T*A;  x**T(1x3), A(3x3)
// y may be stored on first column of A
void dgemv1_3_3(double *y, const double* x, const double* A){ // x(1x3), A(3x3)
	y[0] = x[0]*A[0]+x[1]*A[1]+x[2]*A[2];
	y[1] = x[0]*A[3]+x[1]*A[4]+x[2]*A[5];
	y[2] = x[0]*A[6]+x[1]*A[7]+x[2]*A[8];	
}
void dgemv1_9_3(double *y, const double* x, const double* A){ // x(1x9), A(9x3)
	y[0] = x[0]*A[0] +x[1]*A[1] +x[2]*A[2] +x[3]*A[3] +x[4]*A[4] +x[5]*A[5]+ x[6]*A[6] +x[7]*A[7] +x[8]*A[8];
	y[1] = x[0]*A[9] +x[1]*A[10]+x[2]*A[11]+x[3]*A[12]+x[4]*A[13]+x[5]*A[14]+x[6]*A[15]+x[7]*A[16]+x[8]*A[17];
	y[2] = x[0]*A[18]+x[1]*A[19]+x[2]*A[20]+x[3]*A[21]+x[4]*A[22]+x[5]*A[23]+x[6]*A[24]+x[7]*A[25]+x[8]*A[26];
}

// y = A*x;
// y may be stored on first column of A
void dgemv3_3_1(double *y, const double* A, const double* x){ // A(3x3), x(3x1)
	y[0] = A[0]*x[0] + A[3]*x[1] + A[6]*x[2];
	y[1] = A[1]*x[0] + A[4]*x[1] + A[7]*x[2];
	y[2] = A[2]*x[0] + A[5]*x[1] + A[8]*x[2];
}
void dgemv2_4_1(double *y, const double* A, const double* x){ // A(2x4), x(4x1)
	y[0] = A[0]*x[0] + A[2]*x[1] + A[4] *x[2]+ A[6]*x[3];
	y[1] = A[1]*x[0] + A[3]*x[1] + A[5] *x[2]+ A[7]*x[3];
}
void dgemv4_4_1(double *y, const double* A, const double* x){ // A(4x4), x(4x1)
	y[0] = A[0]*x[0] + A[4]*x[1] + A[8] *x[2]+ A[12]*x[3];
	y[1] = A[1]*x[0] + A[5]*x[1] + A[9] *x[2]+ A[13]*x[3];
	y[2] = A[2]*x[0] + A[6]*x[1] + A[10]*x[2]+ A[14]*x[3];
	y[3] = A[3]*x[0] + A[7]*x[1] + A[11]*x[2]+ A[15]*x[3];
}
void dgemv9_3_1(double *y, const double* A, const double* x){ // A(9x3), x(3x1)
	y[0] = A[0]*x[0] + A[9] *x[1] + A[18]*x[2];
	y[1] = A[1]*x[0] + A[10]*x[1] + A[19]*x[2];
	y[2] = A[2]*x[0] + A[11]*x[1] + A[20]*x[2];
	y[3] = A[3]*x[0] + A[12]*x[1] + A[21]*x[2];
	y[4] = A[4]*x[0] + A[13]*x[1] + A[22]*x[2];
	y[5] = A[5]*x[0] + A[14]*x[1] + A[23]*x[2];
	y[6] = A[6]*x[0] + A[15]*x[1] + A[24]*x[2];
	y[7] = A[7]*x[0] + A[16]*x[1] + A[25]*x[2];
	y[8] = A[8]*x[0] + A[17]*x[1] + A[26]*x[2];
}

// y = A**T*x;
// y may be stored on first column of A
void dgemv4_4_1t(double *y, const double* A, const double* x){ // A(4x4), x(4x1)
	y[0] = A[0] *x[0] + A[1] *x[1] + A[2] *x[2]+ A[3] *x[3];
	y[1] = A[4] *x[0] + A[5] *x[1] + A[6] *x[2]+ A[7] *x[3];
	y[2] = A[8] *x[0] + A[9] *x[1] + A[10]*x[2]+ A[11]*x[3];
	y[3] = A[12]*x[0] + A[13]*x[1] + A[14]*x[2]+ A[15]*x[3];
}
void dgemv2_4_1t(double *y, const double* A, const double* x){ // A(2x4), x(2x1)
	y[0] = A[0]*x[0] + A[1]*x[1] ;
	y[1] = A[2]*x[0] + A[3]*x[1] ;
	y[2] = A[4]*x[0] + A[5]*x[1] ;
	y[3] = A[6]*x[0] + A[7]*x[1] ;
}


// y = L*x with L(2x2), L(4x4) or L(6x6) lower triangular in conventional storage, x(6x1)
// y may be stored on position of x or on first column of L
void dtrmv2(const double *L, const double* x, double* y){
	y[1] = L[1]*x[0] + L[3]*x[1];
	y[0] = L[0]*x[0];
}
void dtrmv4(const double *L, const double* x, double* y){
	y[3] = L[3]*x[0] + L[7]*x[1] + L[11]*x[2] + L[15]*x[3];
	y[2] = L[2]*x[0] + L[6]*x[1] + L[10]*x[2];
	y[1] = L[1]*x[0] + L[5]*x[1];
	y[0] = L[0]*x[0];
}
void dtrmv6(const double *L, const double* x, double* y){
	y[5] = L[5]*x[0] + L[11]*x[1] + L[17]*x[2] + L[23]*x[3] + L[29]*x[4] + L[35]*x[5];
	y[4] = L[4]*x[0] + L[10]*x[1] + L[16]*x[2] + L[22]*x[3] + L[28]*x[4];
	y[3] = L[3]*x[0] + L[ 9]*x[1] + L[15]*x[2] + L[21]*x[3];
	y[2] = L[2]*x[0] + L[ 8]*x[1] + L[14]*x[2];
	y[1] = L[1]*x[0] + L[ 7]*x[1];
	y[0] = L[0]*x[0];
}

// y = L**T*x, with L lower triangular in conventional storage,
// y may be stored on position of x
void dtrmv2t(const double *L, const double* x, double* y){ // L(2x2), x(2x1)
	y[0] = L[0]*x[0]+L[1]*x[1];
	y[1] =           L[3]*x[1];
}
void dtrmv4t(const double *L, const double* x, double* y){ // L(4x4), x(4x1)
	y[0] = L[0]*x[0]+L[1]*x[1]+L[2]*x[2]+L[3]*x[3];
	y[1] =           L[5]*x[1]+L[6]*x[2]+L[7]*x[3];
	y[2] =                    L[10]*x[2]+L[11]*x[3];
	y[3] =                               L[15]*x[3];
}


// Vector matrix vector multiplication
// a = x**T*A*x, with A(nxn) = A**T, x(nx1) i.e. double matrix contraction
double dgemc(const double* A, const double* x, const int n){
	double a=0;
	for(size_t i=0; i<n; i++){
		a+=A[i*(n+1)]*x[i]*x[i];
		
		for(size_t j=i+1; j<n; j++){
			a+=2*A[i*n+j]*x[i]*x[j];
		}
	}
	return a;
}

// Vector vector multiplication

// B = alpha*x*y' + B; x(mx1), y(nx1), A(mxn)
void dger_n(const double alpha, const double *x, const double *y, double *B, const int M, const int N){
	for(size_t n=0; n<N; n++){
		Add(alpha*y[n],1,x,&B[M*n],&B[M*n],M);
	}
}



// Matrix-matrix multiplication



// C := alpha*A*B; A(mxm), A=A**T, B(mxn), C(mxn); (SIDE='L')
// C := alpha*B*A; A(mxm), A=A**T, B(nxm), C(nxm); (SIDE='R')
void SymMultiply(const double ALPHA, const double *A, const double *B, double *C, const int M, const int N,const char SIDE){
	if(SIDE=='L')
		blas.SYMM(SIDE,'L',M,N,ALPHA,A,M,B,M,0,C,M);
	else if(SIDE=='R')
		blas.SYMM(SIDE,'L',M,N,ALPHA,A,N,B,M,0,C,M);
}


// Matrix addition

// C:= alpha*A + beta*B, with A,B,C(pxq) with pxq=M
void Add(const double alpha, const double beta, const double *A,const double *B, double *C, const int M){
	for (int i=0;i<M;i++) C[i] = alpha*A[i]+beta*B[i];
}


// Matrix scaling

// C:= alpha*A + beta; with A,C(pxq), with pxq=M
void Scale(const double alpha, const double beta, const double *A,double *C, const int M){
	for (int i=0;i<M;i++) C[i] = alpha*A[i]+beta;
}


// Dot product dot(A,B);
double Dot(const double *A,const double *B, const int M){
	double dot = 0;
	for (int i=0;i<M;i++) dot+= A[i]*B[i];
	return dot;
}

// Cross product C=cross(A,B), with A,B,C (3x1);
void Cross(const double *A,const double *B, double* C){
	C[0] = A[1]*B[2]-A[2]*B[1];
	C[1] =-A[0]*B[2]+A[2]*B[0];
	C[2] = A[0]*B[1]-A[1]*B[0];
}



// Conventional to packed storage and other way around

// C(nxn) - conventional storage, P(n*(n+1)/2) - packed storage
void ConventionalToPacked(const double* C, double* P, const int n){
	int indexP=0;
	for(int j=0; j<n; j++){
		for(int i=j; i<n; i++){
		P[indexP] = C[i+n*j];
		indexP++;
		}
	}
}
void PackedToConventional(const double* P, double* C, const int n){
	int indexP=0;
	for(int j=0; j<n; j++){
		for(int i=j; i<n; i++){
		C[i+n*j]=P[indexP];
		C[j+n*i]=P[indexP];
		indexP++;
		}
	}
}

}
