#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

#include "../typedefs.h"
class PseudoNormalForm;

class NormalFormFinder
{
    private:
        enum PointType{Unsupported, SaddleFocus, SaddleCenter};

        const int degree; // number of iterations
        CMap f; // input function

        // variables used in computations:
        CJet taylorSeries;
        CMatrix linearPart;
        CJet reminder;

        capd::Complex lambda1, lambda2;

        PseudoNormalForm getInitialValues();
        void nextIteration(PseudoNormalForm *normalForm);
        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex* lambda1, capd::Complex* lambda2);

    public:
        NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint);

        PseudoNormalForm calculatePseudoNormalForm();
};

#endif