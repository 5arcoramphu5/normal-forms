#include <iostream>
#include "capd/capdlib.h"

#include "NormalFormFinder/NormalFormFinder.h" 
#include "NormalFormFinder/PseudoNormalForm.h"
#include "Polynomials/PolynomialOf4Variables.h"
using namespace std;
using namespace capd;

int main()
{
    DVector x({1.5, 2, 1, 0.75});
    int maxDerivativeOrder = 5;
    DMap f("par: s;   var: x1, x2, x3, x4;   fun: x1, s*x2 + x1*x1*x1, x3 + x4*x4*x1, x4;", maxDerivativeOrder);
    f.setParameter("s", 2);

    NormalFormFinder finder(3);
    finder.calculatePseudoNormalForm(f, x);

    return 0;
}