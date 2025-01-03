#ifndef _DEBUG_UTILS_
#define _DEBUG_UTILS_

#include "capd/capdlib.h"
#include "../typedefs.h"
#include "../NormalFormFinder/PseudoNormalForm.h"

const std::string defaultVars[] = {"x1", "x2", "x3", "x4"};
std::string toString(CJet polynomial, const std::string vars[] = defaultVars);

#endif