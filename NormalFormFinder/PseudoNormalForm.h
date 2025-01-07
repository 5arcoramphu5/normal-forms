#ifndef _PSEUDO_NORMAL_FORM_H_
#define _PSEUDO_NORMAL_FORM_H_

#include "NormalFormFinder.h"
#include "../containers/Polynomial.h"

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree, Polynomial<capd::Complex> fTaylorSeries) : F(fTaylorSeries), phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree) {}

        Polynomial<capd::Complex> phi;   
        Polynomial<capd::Complex> n;
        Polynomial<capd::Complex> b; 

    public:
        const Polynomial<capd::Complex> F; // Taylor series expansions of the input function

        // transformation
        Polynomial<capd::Complex> getPhi() const
        { return phi; }
        
        // normal form
        Polynomial<capd::Complex> getN() const
        { return n; }

        // reminder term
        Polynomial<capd::Complex> getB() const
        { return b; }

    template<LoggerType Logger>
    friend class NormalFormFinder;
};

#endif