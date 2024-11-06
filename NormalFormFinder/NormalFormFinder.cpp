#include "NormalFormFinder.h"
#include "PseudoNormalForm.h"

PseudoNormalForm NormalFormFinder::calculatePseudoNormalForm(const DMatrix &lambda, const DMap &f)
{
    // TODO: move system to set point to X=0

    // TODO: diagonalize matrix

    // TODO: check if saddle-focus or saddle-center in X=0

    PseudoNormalForm result;
    for(int i = 0; i < degree; ++i)
    {
        nextIteration(result);
    }
    return result;
}

void NormalFormFinder::nextIteration(PseudoNormalForm &result)
{
    
}