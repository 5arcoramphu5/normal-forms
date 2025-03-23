#pragma once

#include "../templateUtils.hpp"
#include "../typedefs.h"
#include "../containers/Polynomial.hpp"

#include "capd/autodiff/NodeType.h"

using MapFunction = void (*)(capd::autodiff::Node, capd::autodiff::Node[], int, capd::autodiff::Node[], int, capd::autodiff::Node[], int);

// class used to diagonalize 4x4
template<ArithmeticType Coeff>
class Diagonalization
{
    public:
        const Matrix<Coeff> lambda;

        Diagonalization(MapFunction f, int noParams, const Vector<Coeff> &p, const Matrix<Coeff> &J, const Matrix<Coeff> &invJ, const Matrix<Coeff> &lambda, int maxDerivative);

        Polynomial<Coeff> getDiagonalizedTaylorSeries(int degree) const;

        Vector<Coeff> toDiag(const Vector<Coeff> &vector);
        Vector<Coeff> toOriginal(const Vector<Coeff> &vector);

        void setParameter(int index, Coeff value);

    private:

        const MapFunction f;
        const Vector<Coeff> p;
        const Matrix<Coeff> J;
        const Matrix<Coeff> invJ;
        const int maxDerivative;
        int noParams;

        Map<Coeff> diagonalizedF;

        // function f(p + invJ*x) that can be used in capd::map::Map<...> 
        static void functionWithSubstitution(MapFunction f, int noParams, capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParamsInner);
};

#include "Diagonalization.tpp"
