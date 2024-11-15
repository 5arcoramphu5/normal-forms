#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include "../Polynomials/PolynomialOf4Variables.h"
#include "capd/diffAlgebra/lib.h"

void getLinearPartWithReminder(const capd::DJet &taylor, capd::DMatrix *linearPart, PolynomialOf4Variables4 *reminder);

capd::DJet getTaylorSeries(const capd::DMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

PolynomialOf4Variables4 proj_P(const PolynomialOf4Variables4 &poly);

PolynomialOf4Variables4 proj_R(const PolynomialOf4Variables4 &poly);

#endif