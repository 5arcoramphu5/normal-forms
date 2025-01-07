#ifndef _PSEUDO_NORMAL_FORM_H_
#define _PSEUDO_NORMAL_FORM_H_

#include "NormalFormFinder.h"
#include "../containers/Polynomial.h"

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree, Polynomial fTaylorSeries) : F(fTaylorSeries), phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree) {}

        Polynomial phi;   
        Polynomial n;
        Polynomial b; 

    public:
        const Polynomial F; // Taylor series expansions of the input function

        // transformation
        Polynomial getPhi() const
        { return phi; }
        
        // normal form
        Polynomial getN() const
        { return n; }

        // reminder term
        Polynomial getB() const
        { return b; }

    template<LoggerType Logger>
    friend class NormalFormFinder;
};

#endif