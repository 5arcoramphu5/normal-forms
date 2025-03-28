#pragma once

#include "../containers/Polynomial.hpp"
#include "../logging/logging.hpp"

template<LoggerType Logger>
class NormalFormFinder;

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree) : phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree) {}

        Polynomial<capd::Complex> phi;   
        Polynomial<capd::Complex> n;
        Polynomial<capd::Complex> b; 

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

    template<LoggerType Logger>
    friend class NormalFormFinder;
};