#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

#include "../typedefs.h"
class PseudoNormalForm;

class NormalFormFinder
{
    private:
        enum PointType{Unsupported, SaddleFocus, SaddleCenter};

        const int degree;
        CMap f;

        void nextIteration(PseudoNormalForm *result);
        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex* lambda1, capd::Complex* lambda2);

    public:
        NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint);

        PseudoNormalForm calculatePseudoNormalForm();
};

#endif