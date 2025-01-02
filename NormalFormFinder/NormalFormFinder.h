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

        // system: X' = F(X) 
        // F(X) = J f(p + J^-1 X)
        // DF = lambda = J Df J^-1
        // where:
        // * f(p) = 0, it implies F(0) = 0
        // * lambda is diagonal
        const CMap f;
        const CVector p;
        const CMatrix J;
        const CMatrix invJ;
        const CMatrix lambda;

        // variables used in computations:
        CJet F_taylorSeries;
        CJet F_reminder;
        capd::Complex lambda1, lambda2;

        CJet a1_reminder;
        CJet a2_reminder;

        int iterations;

        PseudoNormalForm getInitialNormalFormValues();
        void setInitialValues();
        void nextIteration(PseudoNormalForm &normalForm);

        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex &lambda1, capd::Complex &lambda2);

        // solves equation of type: L(R(Psi)) = R(H)
        void solveFirstEquation(CJet &Psi, const CJet &H);
        // solves equation of type: N + B = P(H)
        void solveSecondEquation(CJet &N, CJet &B, const CJet &H);

        void checkFirstEquation(const CJet &Psi, const CJet &H, const CJet &N);

        void checkSecondEquation(const CJet &N, const CJet &B, const CJet &H);

        template<VerbosityLevel MessageVerbosity, Streamable MessageType>
        static inline void log(MessageType message)
        { 
            Logger::template print<MessageVerbosity>(message); 
        }

    public:
        NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint, const CMatrix &diagonalDerivative, const CMatrix &diagonalizationMatrix, const CMatrix diagonalizationMatrixInverse);

        PseudoNormalForm calculatePseudoNormalForm();
};

#endif