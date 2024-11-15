#ifndef _TYPEDEFS_
#define _TYPEDEFS_

#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/fields/lib.h"

#include "capd/map/Map.hpp"

#define DIMENSION 4
typedef capd::vectalg::Vector<capd::Complex, DIMENSION> CVector;
typedef capd::vectalg::Matrix<capd::Complex, DIMENSION, DIMENSION> CMatrix;
typedef capd::map::Map<CMatrix> CMap; // TODO: use it

#endif