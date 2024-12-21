#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

#include "../typedefs.h"
#include "../debugUtils/logging.h"
class PseudoNormalForm;

template<LoggerType Logger = Logger<VerbosityLevel::None>>
class NormalFormFinder
{
    private:
        enum PointType{Unsupported, SaddleFocus, SaddleCenter};

        const int degree; // number of iterations
        CMap f; // input function

        // variables used in computations:
        CJet taylorSeries;
        CMatrix linearPart;
        CJet f_reminder;

        capd::Complex lambda1, lambda2;
        CJet a1_reminder;
        CJet a2_reminder;

        int iterations;

        PseudoNormalForm getInitialNormalFormValues();
        void setInitialValues();
        void nextIteration(PseudoNormalForm *normalForm);

        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex* lambda1, capd::Complex* lambda2);

        // solves equation of type: L(R(Psi)) = R(H)
        void solveFirstEquation(CJet &Psi, const CJet &H);
        // solves equation of type: N + B = P(H)
        void solveSecondEquation(CJet &N, CJet &B, const CJet &H);

        template<VerbosityLevel MessageVerbosity>
        static void log(std::string message)
        { 
            Logger::template print<MessageVerbosity>(message); 
        }

    public:
        NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint);

        PseudoNormalForm calculatePseudoNormalForm();
};

#endif