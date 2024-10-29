#include <iostream>
#include "capd/capdlib.h"
#include "NormalFormFinder/NormalFormFinder.h" 
#include "NormalFormFinder/PseudoNormalForm.h"

using namespace std;
using namespace capd;

int main()
{
    NormalFormFinder nff(5);
    nff.CalculatePseudoNormalForm();

    return 0;
}