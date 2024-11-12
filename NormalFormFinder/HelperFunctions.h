#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"
#include "capd/capdlib.h"
#include "../Polynomials/PolynomialOf4Variables.h"

void getLinearPartWithReminder(const capd::DJet &taylor, capd::DMatrix *linearPart, PolynomialOf4Variables4 *reminder);

capd::DJet getTaylorSeries(const capd::DMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

#endif