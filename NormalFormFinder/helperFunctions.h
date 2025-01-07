#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include "../containers/Polynomial.h"
#include "../containers/PolynomialMatrix.h"
#include "../containers/PairMap.h"

Polynomial<capd::Complex> getTaylorSeries(const CMap &function, int degree);

Polynomial<capd::Complex> projP(const Polynomial<capd::Complex> &poly, int upToDegree = -1);

Polynomial<capd::Complex> projR(const Polynomial<capd::Complex> &poly, int upToDegree = -1);

CVector gamma(int p, int q, capd::Complex lambda1, capd::Complex lambda2);

PairMap<Polynomial<capd::Complex>> pqCoefficients(const Polynomial<capd::Complex> &poly, int upToDegree);

Polynomial<capd::Complex> operatorL(const Polynomial<capd::Complex> Psi, const Polynomial<capd::Complex> &N, const CMatrix &lambda);

PolynomialMatrix<capd::Complex, 4> D(const Polynomial<capd::Complex> &F);

#endif