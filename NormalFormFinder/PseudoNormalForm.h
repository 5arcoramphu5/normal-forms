#ifndef _PSEUDO_NORMAL_FORM_H_
#define _PSEUDO_NORMAL_FORM_H_

#include "NormalFormFinder.h"
#include "../typedefs.h"

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree, CJet fTaylorSeries) : F(fTaylorSeries), phi(4, 4, degree), n(4, 4, degree), b(4, 4, degree) {}

        CJet phi;   
        CJet n;
        CJet b; 

    public:
        const CJet F; // Taylor series expansions of the input function

        // transformation
        CJet getPhi() const
        { return phi; }
        
        // normal form
        CJet getN() const
        { return n; }

        // reminder term
        CJet getB() const
        { return b; }

    template<LoggerType Logger>
    friend class NormalFormFinder;
};

#endif