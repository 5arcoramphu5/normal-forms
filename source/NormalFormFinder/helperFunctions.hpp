#pragma once

#include "../typedefs.h"
#include "../containers/Polynomial.hpp"
#include "../containers/PolynomialMatrix.hpp"
#include "../containers/PairMap.hpp"

Polynomial<capd::Complex> projP(const Polynomial<capd::Complex> &poly, int upToDegree = -1);

Polynomial<capd::Complex> projR(const Polynomial<capd::Complex> &poly, int upToDegree = -1);

CVector gamma(int p, int q, capd::Complex lambda1, capd::Complex lambda2);

Polynomial<capd::Complex> operatorL(const Polynomial<capd::Complex> Psi, const Polynomial<capd::Complex> &N, const CMatrix &lambda);

PolynomialMatrix<capd::Complex, 4> D(const Polynomial<capd::Complex> &F);