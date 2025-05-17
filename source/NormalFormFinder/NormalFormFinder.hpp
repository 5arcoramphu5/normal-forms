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

        // input system: X' = F(X) 
        // F(X) = J f(p + J^-1 X)
        // DF = lambda = J Df J^-1
        // where:
        // * f(p) = 0, it implies F(0) = 0
        // * lambda - a diagonal matrix of form:    lambda1     0           0           0
        //                                          0           -lambda1    0           0
        //                                          0           0           lambda2     0
        //                                          0           0           0           -lambda2

        // variables used in computations:
        PointType pointType;
        capd::Complex lambda1, lambda2;
        int iterations;

        void setPointTypeAndLambdas(const CMatrix &lambda);
        PseudoNormalForm getInitialNormalFormValues(const Diagonalization<capd::Complex> &diagonalization);
        void nextIteration(PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diagonalization);


        // solves equation of type: L(R(Psi)) = R(H)
        void solveFirstEquation(Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &a1, const Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H);
        // solves equation of type: N + B = P(H)
        void solveSecondEquation(Polynomial<capd::Complex> &N, Polynomial<capd::Complex> &B, Polynomial<capd::Complex> &a1, Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H);

        void checkFirstEquation(const Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &H, const Polynomial<capd::Complex> &N, const Diagonalization<capd::Complex> &diagonalization);

        void checkSecondEquation(const Polynomial<capd::Complex> &N, const Polynomial<capd::Complex> &B, const Polynomial<capd::Complex> &H);

        void checkNormalFormEquality(const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diagonalization);

        template<VerbosityLevel MessageVerbosity, Streamable... MessageTypes>
        static inline void log(MessageTypes... message)
        { 
            Logger::template log<MessageVerbosity>(message...); 
        }

    public:
        NormalFormFinder(int _degree);

        PseudoNormalForm calculatePseudoNormalForm(const Diagonalization<capd::Complex> &diagonalization);
};

#include "NormalFormFinder.tpp"