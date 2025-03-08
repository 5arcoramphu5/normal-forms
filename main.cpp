#include <iostream>

#include "source/typedefs.h"
#include "source/NormalFormFinder/NormalFormFinder.hpp"
using namespace std;
using namespace capd;
using capd::autodiff::Node;

#define MAX_DERIVATIVE 10
#define METHOD_DEGREE 4

#define LOGGER Logger<Diagnostic, SymbolicPolynomialPrinting, 13>

void diagonal_matrix_test()
{
    CMap f(
        "par: p1, p2, p3, p4, a1, a2, a3;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "p1*x1 + a1*x1*x3*x4,"
            "p2*x2 + a2*x2*x4,"
            "p3*x3 + a3*x2*x2,"
            "p4*x4 + x1*x1*x1 + x2*x3;", 
        MAX_DERIVATIVE);

    f.setParameter("p1", Complex(1, 1));
    f.setParameter("p2", Complex(-1, -1));
    f.setParameter("p3", Complex(1, -1));
    f.setParameter("p4", Complex(-1, 1));
    f.setParameter("a1", 2);
    f.setParameter("a2", 5);
    f.setParameter("a3", 1i);

    CVector p({0, 0, 0, 0});
    CMatrix lambda({ {Complex(1, 1), 0, 0, 0}, {0, Complex(-1, -1), 0, 0}, {0, 0, Complex(1, -1), 0}, {0, 0, 0, Complex(-1, 1)} });
    CMatrix J({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });
    CMatrix invJ({ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} });

    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, f, p, lambda, J, invJ);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();
}

void henon_heiles_test()
{
    CMap f(
        "par: a1, lambda;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "x2,"
            "-x1 - a1*lambda*x1*x3,"
            "x4,"
            "-x3 - lambda*(x1^2 - x3^2);", 
        MAX_DERIVATIVE);

    f.setParameter("a1", 2);
    f.setParameter("lambda", 1);

    CMatrix lambda({ 
        {Complex(0, 1), 0, 0, 0}, 
        {0, Complex(0, -1), 0, 0}, 
        {0, 0, Complex(0, 1), 0}, 
        {0, 0, 0, Complex(0, -1)} });

    CMatrix J({ 
        {Complex(0, 0.5), Complex(0.5, 0), 0, 0} ,
        {Complex(0, -0.5), Complex(0.5, 0), 0, 0}, 
        {0, 0, Complex(0, 0.5), Complex(0.5, 0)},
        {0, 0, Complex(0, -0.5), Complex(0.5, 0)} });

    CMatrix invJ({ 
        {Complex(0, -1), Complex(0, 1), 0, 0}, 
        {1, 1, 0, 0}, 
        {0, 0, Complex(0, -1), Complex(0, 1)}, 
        {0, 0, 1, 1} });

    CVector p({0, 0, 0, 0});

    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, f, p, lambda, J, invJ);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();
}

void pcr3bpVectorField(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
    Node mu = params[0]; // relative mass of Jupiter
    Node mj = 1 - mu; // relative mass of the Sun

    // in[0] = x, in[1] = vx, in[2] = y, in[3] = vy
    Node xMu = in[0] + mu;
    Node xMj = in[0] - mj;
    Node ySquare = in[2]*in[2];
    Node xMuSquareYSquare = xMu*xMu + ySquare;
    Node xMjSquareYSquare = xMj*xMu + ySquare;
    Node factor1 = mj / ( xMuSquareYSquare * sqrt(xMuSquareYSquare) );
    Node factor2 = mu / ( xMjSquareYSquare * sqrt(xMjSquareYSquare) );

    out[0] = in[1];
    out[1] = in[0] - xMu*factor1 - xMj*factor2 + 2*in[3];
    out[2] = in[3];
    out[3] = in[2]*(1 - factor1 - factor2) - 2*in[1];
}

// not working because of CAPD sqrt
void PCR3BP_test() 
{
    int dim=4, noParam=1;
    CMap f(pcr3bpVectorField, dim, dim, noParam, MAX_DERIVATIVE);
    f.setParameter(0, 0.5); // mu parameter

    CVector p({0, 0, 0.866025403784439, 0}); // L4

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

    NormalFormFinder<LOGGER> finder(METHOD_DEGREE, f, p, lambda, J, invJ);
    PseudoNormalForm normalForm = finder.calculatePseudoNormalForm();
}

int main()
{
    // diagonal_matrix_test();
    // henon_heiles_test();
    PCR3BP_test();

    return 0;
}