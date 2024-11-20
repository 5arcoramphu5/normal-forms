#include <iostream>

#include "typedefs.h"
#include "NormalFormFinder/NormalFormFinder.h" 
#include "NormalFormFinder/PseudoNormalForm.h"
#include "Polynomials/PolynomialOf4Variables.h"
using namespace std;
using namespace capd;
using namespace capd::fields;

int main()
{
    CVector x({1.5, 2, 1, 0.75});
    int maxDerivativeOrder = 5;
    CMap f(
        "par: p1, p2, p3, p4, a1, a2, a3;"
        "var: x1, x2, x3, x4;"
        "fun:"
            "p1*x1 + a1*x1*x3*x4,"
            "p2*x2 + a2*x1*x4,"
            "p3*x3 + a3*x2*x2,"
            "p4*x4 + x1*x1*x1 + x2*x3;", 
        maxDerivativeOrder);

    f.setParameter("p1", 1 + 1i);
    f.setParameter("p2", -1 - 1i);
    f.setParameter("p3", 1 - 1i);
    f.setParameter("p4", -1 + 1i);
    f.setParameter("a1", 2);
    f.setParameter("a2", 2);
    f.setParameter("a3", 1i);

    NormalFormFinder finder(3);
    finder.calculatePseudoNormalForm(f, x);

    return 0;
}