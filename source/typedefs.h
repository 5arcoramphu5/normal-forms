#pragma once

#include "capd/vectalg/lib.h"
#include "capd/map/lib.h"
#include "capd/fields/lib.h"

#include "capd/map/Map.hpp"

template<typename Coeff>
using Vector = capd::vectalg::Vector<Coeff, 0>;

using CVector = Vector<capd::Complex>;

template<typename Coeff>
using Matrix = capd::vectalg::Matrix<Coeff, 0, 0>;

using CMatrix = Matrix<capd::Complex>;

template<typename Coeff>
using Map = capd::map::Map<Matrix<Coeff>>;

using CMap = Map<capd::Complex>;

typedef capd::diffAlgebra::Jet<CMatrix, 0> CJet;
typedef capd::vectalg::ColumnVector<capd::Complex, 0> CColumnVector;