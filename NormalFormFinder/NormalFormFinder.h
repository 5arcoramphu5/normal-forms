#ifndef _NORMAL_FORM_FINDER_H_
#define _NORMAL_FORM_FINDER_H_

#include "../typedefs.h"
#include "../debugUtils/logging.h"
#include "../containers/Polynomial.h"

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
        // * lambda a diagonal matrix of form:  lambda1     0           0           0
        //                                      0           -lambda1    0           0
        //                                      0           0           lambda2     0
        //                                      0           0           0           -lambda2
        const CMap f;
        const CVector p;
        const CMatrix J;
        const CMatrix invJ;
        const CMatrix lambda;

        // variables used in computations:
        Polynomial F_taylorSeries;
        Polynomial F_reminder;
        capd::Complex lambda1, lambda2;

        Polynomial a1_reminder;
        Polynomial a2_reminder;

        int iterations;

        PseudoNormalForm getInitialNormalFormValues();
        void setInitialValues();
        void nextIteration(PseudoNormalForm &normalForm);

        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex &lambda1, capd::Complex &lambda2);

        // solves equation of type: L(R(Psi)) = R(H)
        void solveFirstEquation(Polynomial &Psi, const Polynomial &H);
        // solves equation of type: N + B = P(H)
        void solveSecondEquation(Polynomial &N, Polynomial &B, const Polynomial &H);

        void checkFirstEquation(const Polynomial &Psi, const Polynomial &H, const Polynomial &N);

        void checkSecondEquation(const Polynomial &N, const Polynomial &B, const Polynomial &H);

        void checkNormalFormEquality(const PseudoNormalForm &normalForm);

        template<VerbosityLevel MessageVerbosity, Streamable MessageType>
        static inline void log(MessageType message)
        { 
            Logger::template log<MessageVerbosity>(message); 
        }

    public:
        NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint, const CMatrix &diagonalDerivative, const CMatrix &diagonalizationMatrix, const CMatrix diagonalizationMatrixInverse);

        PseudoNormalForm calculatePseudoNormalForm();
};

#endif