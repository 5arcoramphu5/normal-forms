#include "Diagonalization.hpp"

template <ArithmeticType Coeff>
Diagonalization<Coeff>::Diagonalization() {}

template <ArithmeticType Coeff>
Diagonalization<Coeff>::Diagonalization(MapFunction f, int noParams, const Vector<Coeff> &p, const Matrix<Coeff> &J, const Matrix<Coeff> &invJ, const Matrix<Coeff> &lambda, int maxDerivative) 
    : f(f), p(p), J(J), invJ(invJ), lambda(lambda), maxDerivative(maxDerivative), noParams(noParams),
    diagonalizedF(
        [f, noParams](capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int n){
            functionWithSubstitution(f, noParams, t, in, dimIn, out, dimOut, params, n);
        }, 
        4, 4, noParams+20, maxDerivative)
{
    setDiagonalizedFParameters();
}

template <ArithmeticType Coeff>
Diagonalization<Coeff>& Diagonalization<Coeff>::operator=(const Diagonalization<Coeff> &other)
{
    f = other.f;
    p = other.p;
    J = other.J;
    invJ = other.invJ;
    lambda = other.lambda;
    maxDerivative = other.maxDerivative;
    noParams = other.noParams;
    diagonalizedF = Map<Coeff>(
        [this](capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int n){
            functionWithSubstitution(this->f, this->noParams, t, in, dimIn, out, dimOut, params, n);
        }, 
        4, 4, noParams+20, maxDerivative);
    diagonalizedF.setDegree(maxDerivative);

    setDiagonalizedFParameters();
    for(int i = 0; i < noParams; ++i)
        diagonalizedF.setParameter(i, other.diagonalizedF.getParameter(i));

    return *this;
}

template <ArithmeticType Coeff>
Diagonalization<Coeff>::Diagonalization(const Diagonalization &other) 
    : f(other.f), p(other.p), J(other.J), invJ(other.invJ), lambda(other.lambda), maxDerivative(other.maxDerivative), noParams(other.noParams),
    diagonalizedF(
        [this](capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int n){
            functionWithSubstitution(this->f, this->noParams, t, in, dimIn, out, dimOut, params, n);
        }, 
        4, 4, noParams+20, maxDerivative)
{
    setDiagonalizedFParameters();
    for(int i = 0; i < noParams; ++i)
        diagonalizedF.setParameter(i, other.diagonalizedF.getParameter(i));
}

template <ArithmeticType Coeff>
void Diagonalization<Coeff>::setDiagonalizedFParameters()
{
    for(int i = 0; i < 4; ++i)
        diagonalizedF.setParameter(noParams+i, p[i]);
    
    int i = noParams+4;
    for(int x = 0; x < 4; ++x)
        for(int y = 0; y < 4; ++y)
        {
            diagonalizedF.setParameter(i, invJ[x][y]);
            i++;
        }
}

template <ArithmeticType Coeff>
inline Polynomial<Coeff> Diagonalization<Coeff>::getDiagonalizedTaylorSeries(int degree) const
{
    Polynomial<Coeff> taylor(4, 4, degree);
    diagonalizedF(CVector({0, 0, 0, 0}), taylor);
    taylor = J * taylor;
    return taylor;
}

template <ArithmeticType Coeff>
void Diagonalization<Coeff>::functionWithSubstitution(MapFunction f, int noParams, capd::autodiff::Node t, capd::autodiff::Node in[], int /*dimIn*/, capd::autodiff::Node out[], int /*dimOut*/, capd::autodiff::Node params[], int /*noParamsInner*/)
{
    // params 0 - noParams-1 - parameters of f
    // params noParams - noParams+3 - point p
    // params noParams+4 - noParams+19 - invJ
    capd::autodiff::Node newIn[4]; // p + invJ*in
    newIn[0] = params[noParams] + params[noParams+4] * in[0] + params[noParams+5] * in[1] + params[noParams+6] * in[2] + params[noParams+7] * in[3];
    newIn[1] = params[noParams+1] + params[noParams+8] * in[0] + params[noParams+9] * in[1] + params[noParams+10] * in[2] + params[noParams+11] * in[3];
    newIn[2] = params[noParams+2] + params[noParams+12] * in[0] + params[noParams+13] * in[1] + params[noParams+14] * in[2] + params[noParams+15] * in[3];
    newIn[3] = params[noParams+3] + params[noParams+16] * in[0] + params[noParams+17] * in[1] + params[noParams+18] * in[2] + params[noParams+19] * in[3];
    f(t, newIn, 4, out, 4, params, noParams);
}

template <ArithmeticType Coeff>
Vector<Coeff> Diagonalization<Coeff>::toDiag(const Vector<Coeff> &vector)
{
    return J * (vector - p);
}

template <ArithmeticType Coeff>
Vector<Coeff> Diagonalization<Coeff>::toOriginal(const Vector<Coeff> &vector)
{
    return p + invJ * vector;
}

template <ArithmeticType Coeff>
void Diagonalization<Coeff>::setParameter(int index, Coeff value)
{
    diagonalizedF.setParameter(index, value);
}