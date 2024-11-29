#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, CJet *reminder);

CJet getTaylorSeries(const CMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

CJet projP(const CJet &poly);

CJet projR(const CJet &poly);

#endif