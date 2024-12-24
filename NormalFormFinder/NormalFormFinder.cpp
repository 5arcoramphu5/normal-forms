#include "NormalFormFinder.h"
#include "PseudoNormalForm.h"
#include "helperFunctions.h"
#include "../debugUtils/debugUtils.h"

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

#include <vector>

using namespace capd;
using namespace std;

template<LoggerType Logger>
NormalFormFinder<Logger>::NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint) : degree(_degree), f(_f)
{
    if(f.dimension() != 4 || f.imageDimension() != 4)
        throw runtime_error("Dimensions not supported.");
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::calculatePseudoNormalForm()
{
    taylorSeries = getTaylorSeries(f, degree);

    CMatrix _linearPart;
    CJet _reminder(4, 4, degree);
    getLinearPartWithReminder(taylorSeries, &_linearPart, &_reminder);

    linearPart = _linearPart;
    f_reminder = _reminder;

    PointType pointType = getPointType(linearPart, &lambda1, &lambda2);
    if(pointType == PointType::Unsupported)
        throw new runtime_error("Type of equilibrium point not supported.");

    PseudoNormalForm normalForm = getInitialNormalFormValues();
    setInitialValues();

    for(int i = 0; i < degree; ++i)
    {
        nextIteration(&normalForm);
        
        log<VerbosityLevel::Diagnostic>("------------------------");
        log<VerbosityLevel::Minimal>("Phi:\n" + toString(normalForm.getPhi()));
        log<VerbosityLevel::Minimal>("N:\n" + toString(normalForm.getN()));
        log<VerbosityLevel::Minimal>("B:\n" + toString(normalForm.getB()));
    }
    return normalForm;
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::getInitialNormalFormValues()
{
    PseudoNormalForm normalForm(degree, taylorSeries);

    for(int i = 0; i < 4; ++i)
    {
        int indexArr[4] = {0, 0, 0, 0};
        indexArr[i] = 1;

        normalForm.phi(i, Multiindex(4, indexArr)) = 1; // Phi(1) = Id
        normalForm.n(i, Multiindex(4, indexArr)) = linearPart(i+1, i+1); // N(1) = linear part (diagonal)
    }

    return normalForm;
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::setInitialValues()
{
    iterations = 1;

    a1_reminder = CJet(1, 2, degree+1);
    a2_reminder = CJet(1, 2, degree+1);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::nextIteration(PseudoNormalForm *normalForm)
{
    log<VerbosityLevel::Minimal>("----- iteration: " + to_string(iterations) + " -----");
        
    CJet FPhi(4, 4, f_reminder.degree() * normalForm->phi.degree());
    substitutionPowerSeries(f_reminder, normalForm->phi, FPhi, false);
    
    solveFirstEquation(normalForm->phi, FPhi);
    checkFirstEquation(normalForm->phi, FPhi, normalForm->n);

    solveSecondEquation(normalForm->n, normalForm->b, FPhi);

    iterations++;
}

template<LoggerType Logger>
typename NormalFormFinder<Logger>::PointType NormalFormFinder<Logger>::getPointType(const CMatrix &diagonalMatrix, Complex* lambda1, Complex* lambda2)
{
    Complex eigenValues[4];
    for(int i = 0; i < 4; ++i) 
        eigenValues[i] = diagonalMatrix[i][i];

    // find pairs lambda, -lambda
    for(int i = 1; i < 4; ++i)
    {
        if(eigenValues[0] == -eigenValues[i])
        {
            std::vector<int> left;
            for(int j = 1; j < 4; ++j)
                if(j != i)
                    left.push_back(j);

            if(eigenValues[left[0]] == -eigenValues[left[1]])
            {
                *lambda1 = eigenValues[0];
                *lambda2 = eigenValues[left[0]];

                if((lambda1->real() == 0 && lambda2->imag() == 0) || (lambda1->imag() == 0 && lambda2->real() == 0))
                    return PointType::SaddleCenter;

                return PointType::SaddleFocus;
            }
        }
    }

    return PointType::Unsupported;
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveFirstEquation(CJet &Psi, const CJet &H)
{
    CJet RH = projR(H, iterations+1);
    Psi = CJet(4, 4, degree);

    auto pq_coeffs = pqCoefficients(RH, iterations+1);

    for(auto [pq, h_pq] : pq_coeffs)
    {
        auto _gamma = gamma(pq.first, pq.second, lambda1, lambda2);
        CJet denominator(4, 2, iterations+1); // gamma + pa_1 + qa_2

        for(int i = 0; i < 4; ++i)
        {
            for(int deg = 0; deg <= denominator.degree(); ++deg)
            {
                Multiindex ind({deg, 0});
                do 
                {
                    denominator(i, ind) = Complex(pq.first, 0) * a1_reminder(0, ind) + Complex(pq.second, 0) * a2_reminder(0, ind);
                }while(ind.hasNext());
            }

            denominator(i, Multiindex({0, 0})) += _gamma[i];
        }

        CJet psi_pq = polyDivision(h_pq, denominator);

        for(int deg = 0; deg <= psi_pq.degree(); ++deg) // fill Psi
        {
            Multiindex ind({deg, 0});
            do 
            {
                if(ind[0] + pq.first < 0 || ind[1] + pq.second < 0)
                    continue; // TODO
                    
                Multiindex psiIndex({ind[0] + pq.first, ind[0], ind[1] + pq.second, ind[1]});
                if(2*ind[0] + pq.first + 2*ind[1] + pq.second <= Psi.degree())
                    Psi(psiIndex) += psi_pq(ind);
            }while(ind.hasNext());
        }
    }
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkFirstEquation(const CJet &Psi, const CJet &H, const CJet &N)
{
    auto LHS = upToDegree(operatorL(projR(Psi), N, linearPart), iterations+1);
    auto RHS = upToDegree(projR(H), iterations+1);
    log<VerbosityLevel::Diagnostic>("first equation (LHS - RHS):");
    log<VerbosityLevel::Diagnostic>(toString(jetSubstraction(LHS, RHS)));
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveSecondEquation(CJet &N, CJet &B, const CJet &H)
{
    // get h1, h2, h3, h4 - parts of P projection of H
    CJet h[4]; 
    for(int i = 0; i < 4; ++i)
        h[i] = CJet(1, 2, iterations+1);

    for(int i = 0; i < iterations+1; ++i)
        for(int j = 0; 2*i+2*j < iterations+1; ++j)
        {   
            Multiindex index({i, j});
            if(j > 0)
            {
                h[0](0, index) = H(0, Multiindex({i+1, i, j, j}));
                h[1](0, index) = H(1, Multiindex({i, i+1, j, j}));
            }
            if(i > 0)
            {
                h[2](0, index) = H(2, Multiindex({i, i, j+1, j}));
                h[3](0, index) = H(3, Multiindex({i, i, j, j+1}));
            }
        }

    // calculate a1, a2, b1, b2
    CJet b1(1, 2, iterations+1), b2(1, 2, iterations+1);

    for(int deg = 0; deg <= iterations+1; ++deg)
    {
        Multiindex index({deg, 0});
        do 
        {
            a1_reminder(0, index) = (h[0](0, index) - h[1](0, index)) / Complex(2, 0);
            a2_reminder(0, index) = (h[2](0, index) - h[3](0, index)) / Complex(2, 0);
            b1(0, index) = (h[0](0, index) + h[1](0, index)) / Complex(2, 0);
            b2(0, index) = (h[2](0, index) - h[3](0, index)) / Complex(2, 0);
            
        }while(index.hasNext());
    }

    // calculate N and B
    for(int deg = 0; deg <= iterations+1 && 2*deg + 1 <= N.degree(); ++deg)
    {
        Multiindex index({deg, 0});
        do 
        {
            N(0, Multiindex({index[0]+1, index[0], index[1], index[1]})) = a1_reminder(0, index);
            N(1, Multiindex({index[0], index[0]+1, index[1], index[1]})) = -a1_reminder(0, index);
            N(2, Multiindex({index[0], index[0], index[1]+1, index[1]})) = a2_reminder(0, index);
            N(3, Multiindex({index[0], index[0], index[1], index[1]+1})) = -a2_reminder(0, index);

            B(0, Multiindex({index[0]+1, index[0], index[1], index[1]})) = b1(0, index);
            B(1, Multiindex({index[0], index[0]+1, index[1], index[1]})) = b1(0, index);
            B(2, Multiindex({index[0], index[0], index[1]+1, index[1]})) = b2(0, index);
            B(3, Multiindex({index[0], index[0], index[1], index[1]+1})) = b2(0, index);
            
        }while(index.hasNext());
    }

}
