#ifndef _HELPER_FUNCTIONS_
#define _HELPER_FUNCTIONS_

#include "../typedefs.h"

void getLinearPartWithReminder(const CJet &taylor, CMatrix *linearPart, CJet *reminder);

CJet getTaylorSeries(const CMap &function, int degree);

CVector getEigenvalues(const capd::DMatrix &matrix);

CJet projP(const CJet &poly);

CJet projR(const CJet &poly);

std::string toString(CJet polynomial, std::string var1 = "x1", std::string var2 = "x2", std::string var3 = "x3", std::string var4 = "x4");

#endif