#ifndef _TYPEDEFS_
#define _TYPEDEFS_

#include "capd/capdlib.h"
#include "capd/capdAlglib.h"
#include "capd/map/Map.hpp"

#define DIMENSION 4
typedef capd::vectalg::Vector<capd::Complex, DIMENSION> CVector;
typedef capd::vectalg::Matrix<capd::Complex, DIMENSION, DIMENSION> CMatrix;
typedef capd::map::Map<CMatrix> CMap; // TODO: use it

#endif