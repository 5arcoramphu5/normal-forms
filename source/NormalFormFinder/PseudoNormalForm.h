#pragma once

#include "../containers/Polynomial.hpp"
#include "../logging/logging.hpp"

template<LoggerType Logger>
class NormalFormFinder;

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree) : phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree), a1(1, 2, degree), a2(1, 2, degree) {}

        Polynomial<capd::Complex> phi;   
        Polynomial<capd::Complex> n;
        Polynomial<capd::Complex> b; 
        Polynomial<capd::Complex> a1;
        Polynomial<capd::Complex> a2;

    public:

        // transformation
        const Polynomial<capd::Complex>& getPhi() const
        { return phi; }
        
        // normal form
        const Polynomial<capd::Complex>& getN() const
        { return n; }

        // reminder term
        const Polynomial<capd::Complex>& getB() const
        { return b; }


        // solution of the system (only for B=0)
        CVector solution(double t, CVector initialPoint) const
        {
            CVector constPart({initialPoint[0]*initialPoint[1], initialPoint[2]*initialPoint[3]}); // value of xi*eta, mu*nu - constant with respect to t
            auto a1_const = a1(constPart)[0];
            auto a2_const = a2(constPart)[0];
            return CVector(
                {
                    exp(a1_const*t) * initialPoint[0],
                    exp(-a1_const*t) * initialPoint[1],
                    exp(a2_const*t) * initialPoint[2],
                    exp(-a2_const*t) * initialPoint[3]
                }
            );
        }

    template<LoggerType Logger>
    friend class NormalFormFinder;
};