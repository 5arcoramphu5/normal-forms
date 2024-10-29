#include "NormalFormFinder.h"
#include "PseudoNormalForm.h"

PseudoNormalForm NormalFormFinder::CalculatePseudoNormalForm()
{
    PseudoNormalForm result;
    for(int i = 0; i < degree; ++i)
    {
        NextIteration(result);
    }
    return result;
}

void NormalFormFinder::NextIteration(PseudoNormalForm &result)
{

}