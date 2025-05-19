#pragma once

#include "../containers/Polynomial.hpp"
#include "../logging/logging.hpp"

template<LoggerType Logger>
class NormalFormFinder;

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree) : phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree), a1(1, 2, degree), a2(1, 2, degree) {}
        PseudoNormalForm(const Polynomial<capd::Complex> &phi, const Polynomial<capd::Complex> &n, const Polynomial<capd::Complex> &b, const Polynomial<capd::Complex> &a1, const Polynomial<capd::Complex> &a2) : phi(phi), n(n), b(b), a1(a1), a2(a2) {}
        
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
        CVector solution(double t, CVector initialPoint) const;

        void serialize(std::ostream &stream) const;
        static PseudoNormalForm deserialize(std::istream &stream);

    template<LoggerType Logger>
    friend class NormalFormFinder;
};

#include "PseudoNormalForm.tpp"