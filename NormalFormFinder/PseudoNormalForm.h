#ifndef _PSEUDO_NORMAL_FORM_H_
#define _PSEUDO_NORMAL_FORM_H_

#include "NormalFormFinder.h"
#include "../typedefs.h"

class PseudoNormalForm
{
    private:
        PseudoNormalForm(int degree, CJet fTaylorSeries) : F(fTaylorSeries), Phi(4, degree), N(4, degree), B(4, degree), iteration(0) {}

        int iteration;

        CJet Phi;   
        CJet N;
        CJet B; 

    public:
        const CJet F; // Taylor series expansions of the input function

        // transformation
        CJet getPhi() const
        { return Phi; }
        
        // normal form
        CJet getN() const
        { return N; }

        // reminder term
        CJet getB() const
        { return B; }

    friend NormalFormFinder;
};

#endif