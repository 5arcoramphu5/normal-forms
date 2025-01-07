#ifndef _POLYNOMIAL_MATRIX_H_
#define _POLYNOMIAL_MATRIX_H_

#include "Polynomial.h"

template<int N>
class PolynomialMatrix : public std::array<std::array<Polynomial, N>, N>
{
    int matrixDegree;

    public:
    PolynomialMatrix(int degree) : matrixDegree(degree)
    { 
        for(auto &array : *this)
            for(auto &jet : array)
                jet = Polynomial(1, N, degree);
    }

    int degree() const { return matrixDegree; }
};

Polynomial operator*(const PolynomialMatrix<4> &jetMatrix, const Polynomial &jet);

#endif
