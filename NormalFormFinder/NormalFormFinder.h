#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

class PseudoNormalForm;

class NormalFormFinder
{
    private:
        int degree;

        void NextIteration(PseudoNormalForm &result);

    public:
        NormalFormFinder(int _degree) : degree(_degree) {}

        PseudoNormalForm CalculatePseudoNormalForm();
};

#endif