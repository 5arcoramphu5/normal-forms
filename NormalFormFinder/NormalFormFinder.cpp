#include "NormalFormFinder.h"
#include "PseudoNormalForm.h"
#include "helperFunctions.h"
#include "../debugUtils/debugUtils.h"

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

#include <vector>

using namespace capd;
using namespace std;

NormalFormFinder::NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint) : degree(_degree), f(_f)
{
    if(f.dimension() != 4 || f.imageDimension() != 4)
        throw runtime_error("Dimensions not supported.");
}

PseudoNormalForm NormalFormFinder::calculatePseudoNormalForm()
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
        cout << "----- iteration: " << i << " -----" << endl;
        
        nextIteration(&normalForm);
        
        // debug
        cout << "Phi:\n" << toString(normalForm.getPhi()) << endl;
        cout << "N:\n" << toString(normalForm.getN()) << endl;
        cout << "B:\n" << toString(normalForm.getB()) << endl;
    }
    return normalForm;
}

PseudoNormalForm NormalFormFinder::getInitialNormalFormValues()
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

void NormalFormFinder::setInitialValues()
{
    iterations = 0;

    a1_reminder = CJet(1, 2, degree);
    a2_reminder = CJet(1, 2, degree);
}

void NormalFormFinder::nextIteration(PseudoNormalForm *normalForm)
{
    iterations++;
    CJet FPhi(4, 4, f_reminder.degree() * normalForm->phi.degree());
    substitutionPowerSeries(f_reminder, normalForm->phi, FPhi, false);

    solveFirstEquation(normalForm->phi, FPhi);
    solveSecondEquation(normalForm->n, normalForm->b, FPhi);
}

NormalFormFinder::PointType NormalFormFinder::getPointType(const CMatrix &diagonalMatrix, Complex* lambda1, Complex* lambda2)
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

void NormalFormFinder::solveFirstEquation(CJet &Psi, const CJet &H)
{
    auto RH = projR(H, iterations+1);
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

void NormalFormFinder::solveSecondEquation(CJet &N, CJet &B, const CJet &H)
{
}
