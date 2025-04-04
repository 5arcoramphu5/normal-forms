#pragma once

#include "../typedefs.h"
#include "../logging/logging.hpp"
#include "../containers/Polynomial.hpp"
#include "../Diagonalization/Diagonalization.hpp"
#include "PseudoNormalForm.h"

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
        // * lambda - a diagonal matrix of form:    lambda1     0           0           0
        //                                          0           -lambda1    0           0
        //                                          0           0           lambda2     0
        //                                          0           0           0           -lambda2
        const Diagonalization<capd::Complex> diagonalization;

        // variables used in computations:
        capd::Complex lambda1, lambda2;
        int iterations;

        PseudoNormalForm getInitialNormalFormValues();
        void nextIteration(PseudoNormalForm &normalForm);

        PointType getPointType(const CMatrix &diagonalMatrix, capd::Complex &lambda1, capd::Complex &lambda2);

        // solves equation of type: L(R(Psi)) = R(H)
        void solveFirstEquation(Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &a1, const Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H);
        // solves equation of type: N + B = P(H)
        void solveSecondEquation(Polynomial<capd::Complex> &N, Polynomial<capd::Complex> &B, Polynomial<capd::Complex> &a1, Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H);

        void checkFirstEquation(const Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &H, const Polynomial<capd::Complex> &N);

        void checkSecondEquation(const Polynomial<capd::Complex> &N, const Polynomial<capd::Complex> &B, const Polynomial<capd::Complex> &H);

        void checkNormalFormEquality(const PseudoNormalForm &normalForm);

        template<VerbosityLevel MessageVerbosity, Streamable... MessageTypes>
        static inline void log(MessageTypes... message)
        { 
            Logger::template log<MessageVerbosity>(message...); 
        }

    public:
        NormalFormFinder(int _degree, const Diagonalization<capd::Complex> &diagonalization);

        PseudoNormalForm calculatePseudoNormalForm();
};

#include "NormalFormFinder.tpp"