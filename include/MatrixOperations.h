#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

namespace MatrixOperations {



// BLAS and LaPACK perations customized for own use

// Cholesky decomposition of 4x4 (or2x2) matrix A=L*L**T
// L may be stored on position of A
void dpotrf4(const double *A, double* L);
void dpotrf2(const double *A, double* L);
void dpotrf(double *A, const int N, int* info); // L is stored on location of A

// B = L**(-1), with L(4x4) lower triangular matrix in conventional storage
// B may be stored on position of L
void dtrtri6(const double *L, double* B);
void dtrtri4(const double *L, double* B);
void dtrtri2(const double *L, double* B);

// B = A^(-1) with A(2x2) , A=A**T
// B may be stored on position of A
void dpotri2(const double *A, double* B);

// B = L*A, with L lower triangular in conventional storage, A(4x4)
// B may be stored on position of A
void dtrmm4(const double *L, const double* A, double* B);

// B = A*L or B=L*A, with L lower triangular in conventional storage, A(2x4)
// B may be stored on position of A
void dtrmm2_4(const double* A, const double *L, double* B);  // A*L(4x4)
void dtrmm2_4L(const double *L, const double* A, double* B); // L(2x2)*A

// C = A*B, A(2x4), B(4x4)
// C MAYNOT be stored on position of A or B
void dgemm2_4(const double* A, const double *B, double* C);

// B = L**T*A, with L lower triangular in conventional storage, A(4x4)
// B may be stored on position of A
void dtrmm4tL(const double *L, const double* A, double* B); // L(4x4), A(4x4)
void dtrmm2_4tL(const double *L, const double* A, double* B); // L(2x2), A(2x4)

// B = A*L**T, with L lower triangular in conventional storage, A(2x4)
// B may be stored on position of A
void dtrmm2_4t(const double *A, const double* L, double* B);

// B = L**T*L, with L lower triangular in conventional storage
// B may be stored on position of L
void dtrmtrm2(const double *L, double* B); // L(2x2)

// C = A-B**T*B, with A(4x4) = A**T, B(4x4) i.e. full rank downdate
// C may be stored on position of A
void dsyrd4(const double* A, const double* B, double* C);

// C = A+beta*B*B**T, with A(2x2) = A**T, B(2x4)
// C may be stored on position of A
void dsyr2_4(const double* A, const double beta, const double* B, double* C);

// A = alpha*A+beta*B**T*B, with A(nxn) = A**T, B(kxn)
void dsyrk_nt(const double alpha, double* A, const double beta, const double* B,const int k, const int n);

// C = Btilde*A*Btilde**T, with Btilde being a diagonal block matrix
// with diagonal blocks B(mxm). A(nxn), m<n, n%m=0, A==A**T, C(nxn) (full update matrix)
// C may be stored on position of A
void dsyrk_f(const double* A, const double* B, double* C, const int n, const int m);

// C := alpha*A*B; A(mxp), B(pxn), C(mxn);
// C MAYNOT be stored on position of A or B
void dgemm(const double ALPHA, const double *A, const double *B, double *C, const int M, const int P, const int N);

// K:= K + B*A*C**T, with A(mxm)==A**T; B(nxm),C(nxm) in conv. storage
// K MAYNOT be stored on either position A,B,C
void dsymm_bact(const double alpha, const double* B, const double *A, const double *C, double *K, const int m, const int n);


// Matrix vector multiplication

// y = beta*y + alpha*A*x;  A(mxn)
void dgemv(const double beta, double *y, const double alpha, const double* A, const double* x, const int m, const int n);

// y = x**T*A;  x(1x3), A(3x3)
// y may be stored on first column of A
void dgemv1_3_3(double *y, const double* x, const double* A); // x(1x3), A(3x3)
void dgemv1_9_3(double *y, const double* x, const double* A); // x(1x9), A(9x3)

// y = A*x;  
// y may be stored on first column of A
void dgemv3_3_1(double *y, const double* A, const double* x); // A(3x3), x(3x1)
void dgemv2_4_1(double *y, const double* A, const double* x); // A(2x4), x(4x1)
void dgemv4_4_1(double *y, const double* A, const double* x); // A(4x4), x(4x1)
void dgemv9_3_1(double *y, const double* A, const double* x); // A(9x3), x(3x1)

// y = A**T*x;
// y may be stored on first column of A
void dgemv4_4_1t(double *y, const double* A, const double* x); // A(4x4), x(4x1)
void dgemv2_4_1t(double *y, const double* A, const double* x); // A(2x4), x(2x1)

// y = L*x with L(4x4) or L(6x6) lower triangular in conventional storage, x(6x1)
// y may be stored on position of x or on first column of L
void dtrmv2(const double *L, const double* x, double* y);
void dtrmv4(const double *L, const double* x, double* y);
void dtrmv6(const double *L, const double* x, double* y);

// y = L**T*x, with L lower triangular in conventional storage,
// y may be stored on position of x
void dtrmv2t(const double *L, const double* x, double* y); // L(2x2), x(2x1)
void dtrmv4t(const double *L, const double* x, double* y); // L(4x4), x(4x1)


// Vector matrix vector multiplication
// a = x**T*A*x, with A(nxn) = A**T, x(nx1) i.e. double matrix contraction
double dgemc(const double* A, const double* x, const int n);



// Vector vector multiplication

// A = alpha*x*y' + A; x(mx1), y(nx1)
void dger_n(const double alpha, const double *x, const double *y, double *A, const int m, const int n);


// C := alpha*A*B; A(mxm), A=A**T, B(mxn), C(mxn); (SIDE='L')
// C := alpha*B*A; A(mxm), A=A**T, B(nxm), C(nxm); (SIDE='R')
void SymMultiply(const double ALPHA, const double *A, const double *B, double *C, const int M, const int N,const char SIDE='L');

// Matrix addition

// C:= alpha*A + beta*B, with A,B,C(pxq) with pxq=M
// C may be stored on both position A,B
void Add(const double alpha, const double beta, const double *A,const double *B, double *C, const int M);


// Matrix scaling

// C:= alpha*A + beta; with A,C(pxq), with pxq=M
void Scale(const double alpha, const double beta, const double *A,double *C, const int M);

// Dot product dot(A,B);
double Dot(const double *A,const double *B, const int M);

// Cross product C=cross(A,B), with A,B,C (3x1);
void Cross(const double *A,const double *B, double* C);


// Conventional to packed storage and other way around

// C(nxn) - conventional storage, P(n*(n+1)/2) - packed storage
void ConventionalToPacked(const double* C, double* P, const int n);
void PackedToConventional(const double* P, double* C, const int n);

}

#endif 
