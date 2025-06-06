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

        Diagonalization();
        Diagonalization(const Diagonalization &other);
        Diagonalization& operator=(const Diagonalization& other);  

        Diagonalization(MapFunction f, int noParams, const Vector<Coeff> &p, const Matrix<Coeff> &J, const Matrix<Coeff> &invJ, const Matrix<Coeff> &lambda, int maxDerivative);

        Polynomial<Coeff> getDiagonalizedTaylorSeries(int degree) const;

        Vector<Coeff> toDiag(const Vector<Coeff> &vector) const;
        Vector<Coeff> toOriginal(const Vector<Coeff> &vector) const;

        void setParameter(int index, Coeff value);

        inline const Matrix<Coeff>& getLambda() const { return lambda; }
        inline const Vector<Coeff>& getP() const { return p; }
        inline const Matrix<Coeff>& getJ() const { return J; }
        inline const Matrix<Coeff>& getinvJ() const { return invJ; }

        Polynomial<Coeff> polynomialComposition(const Polynomial<Coeff> &poly) const;
        Polynomial<Coeff> polynomialCompositionWithReminder(const Polynomial<Coeff> &poly) const;

    private:

        MapFunction f;
        Vector<Coeff> p;
        Matrix<Coeff> J;
        Matrix<Coeff> invJ;
        int maxDerivative;
        int noParams;
        Matrix<Coeff> lambda;

        Map<Coeff> diagonalizedF;

        void setDiagonalizedFParameters();
        // function f(p + invJ*x) that can be used in capd::map::Map<...> 
        static void functionWithSubstitution(MapFunction f, int noParams, capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParamsInner);
};

#include "Diagonalization.tpp"
