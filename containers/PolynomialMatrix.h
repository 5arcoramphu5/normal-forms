#pragma once

#include "Polynomial.h"

template<ArithmeticType Coeff, int N>
class PolynomialMatrix : public std::array<std::array<Polynomial<Coeff>, N>, N>
{
    int matrixDegree;

    public:
    PolynomialMatrix(int degree) : matrixDegree(degree)
    { 
        for(auto &array : *this)
            for(auto &jet : array)
                jet = Polynomial<Coeff>(1, N, degree);
    }

    int degree() const { return matrixDegree; }
};

template<ArithmeticType Coeff, int N>
Polynomial<Coeff> operator*(const PolynomialMatrix<Coeff, N> &jetMatrix, const Polynomial<Coeff> &jet);