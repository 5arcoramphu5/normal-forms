#ifndef _TYPEDEFS_
#define _TYPEDEFS_

#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/fields/lib.h"

#include "capd/map/Map.hpp"

typedef capd::vectalg::Vector<capd::Complex, 0> CVector;
typedef capd::vectalg::Matrix<capd::Complex, 0, 0> CMatrix;
typedef capd::map::Map<CMatrix> CMap;
typedef capd::diffAlgebra::Jet<CMatrix, 0> CJet;

#endif