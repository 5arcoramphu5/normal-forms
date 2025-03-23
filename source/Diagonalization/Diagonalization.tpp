#include "Diagonalization.hpp"

template <ArithmeticType Coeff>
Diagonalization<Coeff>::Diagonalization(MapFunction f, const Vector<Coeff> &p, const Matrix<Coeff> &J, const Matrix<Coeff> &invJ, const Matrix<Coeff> &lambda, int maxDerivative) 
    : f(f), p(p), J(J), invJ(invJ), lambda(lambda), maxDerivative(maxDerivative), 
    diagonalizedF(
        [f](capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams){
            functionWithSubstitution(f, t, in, dimIn, out, dimOut, params, noParams);
        }, 
        4, 4, /* TODO */ 21, maxDerivative)
{
    diagonalizedF.setParameter(0, 0.5); // TODO: delete, mu

    for(int i = 0; i < 4; ++i)
        diagonalizedF.setParameter(i+1, p[i]);
    
    int i = 5;
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
void Diagonalization<Coeff>::functionWithSubstitution(MapFunction f, capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams)
{
    // params 0 - mu // TODO: other parameters
    // params 1-4 - point p
    // params 5-20 - invJ
    capd::autodiff::Node newIn[4]; // p + invJ*in
    newIn[0] = params[1] + params[5] * in[0] + params[6] * in[1] + params[7] * in[2] + params[8] * in[3];
    newIn[1] = params[2] + params[9] * in[0] + params[10] * in[1] + params[11] * in[2] + params[12] * in[3];
    newIn[2] = params[3] + params[13] * in[0] + params[14] * in[1] + params[15] * in[2] + params[16] * in[3];
    newIn[3] = params[4] + params[17] * in[0] + params[18] * in[1] + params[19] * in[2] + params[20] * in[3];
    f(t, newIn, 4, out, 4, params, 1);
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