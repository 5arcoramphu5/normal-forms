#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

#include "../typedefs.h"
class PseudoNormalForm;

class NormalFormFinder
{
    private:
        int degree;

        void nextIteration(PseudoNormalForm &result);

    public:
        NormalFormFinder(int _degree) : degree(_degree) {}

        PseudoNormalForm calculatePseudoNormalForm(const CMap &f, const CVector &point);
};

#endif