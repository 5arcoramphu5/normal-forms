#include <iostream>

#include "source/typedefs.h"
#include "source/NormalFormFinder/NormalFormFinder.hpp"
using namespace std;
using namespace capd;
using capd::autodiff::Node;

#define MAX_DERIVATIVE 10
#define METHOD_DEGREE 5

#define LOGGER Logger<Diagnostic, SymbolicPolynomialPrinting, 5>

void diagonalVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
    out[0] = params[0]*in[0] + 2*in[0]*in[2]*in[3];
    out[1] = params[1]*in[1] + 5*in[1]*in[3];
    out[2] = params[2]*in[2] + params[4]*in[1]*in[1];
    out[3] = params[3]*in[3] + in[0]*in[0]*in[0] + in[1]*in[2];
}

PseudoNormalForm diagonal_matrix_test()
{
    CVector p({0, 0, 0, 0});
    CMatrix lambda({ {Complex(1, 1), 0, 0, 0}, {0, Complex(-1, -1), 0, 0}, {0, 0, Complex(1, -1), 0}, {0, 0, 0, Complex(-1, 1)} });
    CMatrix J = CMatrix::Identity(4);
    CMatrix invJ = CMatrix::Identity(4);

    Diagonalization<Complex> diagonalization(diagonalVectorField, 5, p, J, invJ, lambda, MAX_DERIVATIVE);
    diagonalization.setParameter(0, Complex(1, 1));
    diagonalization.setParameter(1, Complex(-1, -1));
    diagonalization.setParameter(2, Complex(1, -1));
    diagonalization.setParameter(3, Complex(-1, 1));
    diagonalization.setParameter(4, 1i);

    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, diagonalization);
    return finder.calculatePseudoNormalForm();
}

void pcr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{       
    Node mu = params[0]; // relative mass of Jupiter
    Node mj = 1 - mu; // relative mass of the Sun

    Node xMu = in[0] + mu;
    Node xMj = in[0] - mj;
    Node xMuSquare = xMu^2;
    Node xMjSquare = xMj^2;
    Node ySquare = in[1]^2;
    Node factor1 = mj / ((xMuSquare + ySquare)^1.5);
    Node factor2 = mu / ((xMjSquare + ySquare)^1.5);

    out[0] = in[2];
    out[1] = in[3];
    out[2] = in[0] - xMu*factor1 - xMj*factor2 + 2*in[3];
    out[3] = in[1] * (1 - factor1 - factor2) - 2*in[2];
}

// may not work because of CAPD sqrt
PseudoNormalForm PCR3BP_test() 
{
    CVector p({0, 0.866025403784438646763723170753, 0,  0}); // L4

    CMatrix lambda({ 
        {Complex(0.632075, 0.94843), 0, 0, 0}, 
        {0, Complex(-0.632075, -0.94843), 0, 0}, 
        {0, 0, Complex(0.632075, -0.94843), 0}, 
        {0, 0, 0, Complex(-0.632075, 0.94843)} });

    CMatrix J({ 
        {Complex(-0.677045, -0.120902), Complex(0, 1.56774), Complex(-0.417702, -0.958065), Complex(-0.660842, 0.440414)},
        {Complex(0.677045, 0.120902), Complex(0, 1.56774), Complex(-0.417702, -0.958065), Complex(0.660842, -0.440414)},
        {Complex(-0.677045, 0.120902), Complex(0, -1.56774), Complex(-0.417702, 0.958065), Complex(-0.660842, -0.440414)},
        {Complex(0.677045, -0.120902), Complex(0, -1.56774), Complex(-0.417702, 0.958065), Complex(0.660842, 0.440414)} 
    });

    CMatrix invJ({ 
        {Complex(-0.29122, 0.436975), Complex(0.29122, -0.436975), Complex(-0.29122, -0.436975), Complex(0.29122, 0.436975)},
        {Complex(-0.365758, -0.159465), Complex(-0.365758, -0.159465), Complex(-0.365758, 0.159465), Complex(-0.365758, 0.159465)},
        {-0.598513, -0.598513, -0.598513, -0.598513},
        {Complex(-0.0799453, -0.44769), Complex(0.0799453, 0.44769), Complex(-0.0799453, 0.44769), Complex(0.0799453, -0.44769)}
    });

    Diagonalization<Complex> diagonalization(pcr3bpVectorField, 1, p, J, invJ, lambda, MAX_DERIVATIVE);
    diagonalization.setParameter(0, 0.5); // mu parameter

    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, diagonalization);
    return finder.calculatePseudoNormalForm();
}

int main()
{
    auto result1 = diagonal_matrix_test();
    auto result2 = PCR3BP_test();

    return 0;
}