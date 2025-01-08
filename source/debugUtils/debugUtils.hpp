#pragma once

#include "capd/capdlib.h"
#include "../typedefs.h"
#include "../NormalFormFinder/PseudoNormalForm.h"
#include "../containers/Polynomial.hpp"

const std::string defaultVars[] = {"x1", "x2", "x3", "x4"};

template<ArithmeticType Coeff>
std::string toString(Polynomial<Coeff> polynomial, const std::string vars[] = defaultVars);

template<ArithmeticType Coeff>
std::string toCoefficientString(Polynomial<Coeff> polynomial);

#include "debugUtils.tpp"