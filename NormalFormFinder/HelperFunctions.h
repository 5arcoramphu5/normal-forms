#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include "../Polynomials/PolynomialOf4Variables.h"

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, PolynomialOf4Variables4 *reminder);

CJet getTaylorSeries(const CMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

PolynomialOf4Variables4 proj_P(const PolynomialOf4Variables4 &poly);

PolynomialOf4Variables4 proj_R(const PolynomialOf4Variables4 &poly);

#endif