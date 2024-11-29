#ifndef _DEBUG_UTILS_
#define _DEBUG_UTILS_

#include "capd/capdlib.h"
#include "../typedefs.h"
#include "../NormalFormFinder/PseudoNormalForm.h"

std::string toString(CJet polynomial, std::string var1 = "x1", std::string var2 = "x2", std::string var3 = "x3", std::string var4 = "x4");

void checkPseudoNormalCondition(const PseudoNormalForm &normalForm);

#endif