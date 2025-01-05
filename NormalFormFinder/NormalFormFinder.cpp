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
NormalFormFinder<Logger>::NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint, const CMatrix &diagonalDerivative, const CMatrix &diagonalizationMatrix, const CMatrix diagonalizationMatrixInverse)
    : degree(_degree), f(_f), p(fixedPoint), lambda(diagonalDerivative), J(diagonalizationMatrix), invJ(diagonalizationMatrixInverse)
{
    if(f.dimension() != 4 || f.imageDimension() != 4)
        throw runtime_error("Dimensions not supported.");
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::calculatePseudoNormalForm()
{
    CJet f_taylorSeries = getTaylorSeries(f, degree+1);

    auto invJ_P = jetAddition(invJ, p, degree+1); // F(X) = J f(p + J^-1 X)
    F_taylorSeries = CJet(4, 4, invJ_P.degree() * f_taylorSeries.degree());
    substitutionPowerSeries(f_taylorSeries, invJ_P, F_taylorSeries, false);
    F_taylorSeries = J * F_taylorSeries;

    log<Debug>("F:");
    log<Debug>(toString(F_taylorSeries));

    F_reminder = fromToDegree(F_taylorSeries, 2, degree+1);
    log<Debug>("F reminder:");
    log<Debug>(toString(F_reminder));

    PointType pointType = getPointType(lambda, lambda1, lambda2);
    if(pointType == PointType::Unsupported)
        throw runtime_error("Type of equilibrium point not supported.");
    if(pointType == PointType::SaddleCenter) // TODO: to be deleted
        throw runtime_error("Case not implemented yet.");

    log<Debug>("lambda: ");
    log<Debug>(lambda);
    log<Debug>("lambda1, lambda2: ");
    log<Debug>(lambda1);
    log<Debug>(lambda2);

    log<Diagnostic>("Point type: " + to_string(pointType));
    log<Diagnostic>("");

    PseudoNormalForm normalForm = getInitialNormalFormValues();
    setInitialValues();

    for(iterations = 1; iterations <= degree; ++iterations)
    {
        log<Minimal>("--------- iteration: " + to_string(iterations) + " ---------");

        nextIteration(normalForm);
        
        log<Diagnostic>("--------------------------------");
        log<Minimal>("Phi:\n" + toString(normalForm.getPhi()));
        log<Minimal>("N:\n" + toString(normalForm.getN()));
        log<Minimal>("B:\n" + toString(normalForm.getB()));
        log<Diagnostic>("a1:\n" + toString(a1_reminder));
        log<Diagnostic>("a2:\n" + toString(a2_reminder));
    }

    return normalForm;
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::getInitialNormalFormValues()
{
    PseudoNormalForm normalForm(degree+1, F_taylorSeries);

    for(int i = 0; i < 4; ++i)
    {
        int indexArr[4] = {0, 0, 0, 0};
        indexArr[i] = 1;
        
        normalForm.phi(i, Multiindex(4, indexArr)) = 1; // Phi(1) = Id
        normalForm.n(i, Multiindex(4, indexArr)) = lambda[i][i]; // N(1) = linear part (diagonal)
    }

    return normalForm;
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::setInitialValues()
{
    a1_reminder = CJet(1, 2, degree+1);
    a2_reminder = CJet(1, 2, degree+1);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::nextIteration(PseudoNormalForm &normalForm)
{       
    CJet FPhi(4, 4, normalForm.phi.degree());
    substitutionPowerSeries(F_reminder, normalForm.phi, FPhi, false);
    
    solveFirstEquation(normalForm.phi, FPhi);
    checkFirstEquation(normalForm.phi, FPhi, normalForm.n);

    solveSecondEquation(normalForm.n, normalForm.b, FPhi);
    checkSecondEquation(normalForm.n, normalForm.b, FPhi);
}

template<LoggerType Logger>
typename NormalFormFinder<Logger>::PointType NormalFormFinder<Logger>::getPointType(const CMatrix &diagonalMatrix, Complex &lambda1, Complex &lambda2)
{
    // find pairs lambda, -lambda on the diagonal
    for(int i = 1; i < 4; ++i)
    {
        if(diagonalMatrix[0][0] == -diagonalMatrix[i][i])
        {
            std::vector<int> left;
            for(int j = 1; j < 4; ++j)
                if(j != i)
                    left.push_back(j);

            if(diagonalMatrix[left[0]][left[0]] == -diagonalMatrix[left[1]][left[1]])
            {
                lambda1 = diagonalMatrix[0][0];
                lambda2 = diagonalMatrix[left[0]][left[0]];

                if((lambda1.real() == 0 && lambda2.imag() == 0) || (lambda1.imag() == 0 && lambda2.real() == 0))
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

    auto pq_coeffs = pqCoefficients(RH, iterations+1);

    for(auto& [pq, h_pq] : pq_coeffs)
    {
        auto& [p, q] = pq; 

        // calculate only iterations+1 degree
        // 2*deg + p + q == iterations+1
        if( (iterations+1 - p - q) % 2 == 1 ) continue;

        auto _gamma = gamma(p, q, lambda1, lambda2);
        CJet denominator(4, 2, iterations+1); // gamma + pa_1 + qa_2

        for(int i = 0; i < 4; ++i)
        {
            for(int deg = 0; deg <= denominator.degree(); ++deg)
            {
                Multiindex ind({deg, 0});
                do 
                {
                    denominator(i, ind) = Complex(p, 0) * a1_reminder(0, ind) + Complex(q, 0) * a2_reminder(0, ind);
                }while(ind.hasNext());
            }

            denominator(i, Multiindex({0, 0})) += _gamma[i];
        }

        int deg = (iterations+1 - p - q) / 2;


        CJet psi_pq = polyDivision(h_pq, denominator, deg);

        Multiindex ind({deg, 0});
        do 
        {
            if(ind[0] + p < 0 || ind[1] + q < 0)
                continue; // TODO
                
            Multiindex psiIndex({ind[0] + p, ind[0], ind[1] + q, ind[1]});
            if(2*ind[0] + p + 2*ind[1] + q <= Psi.degree())
            {
                for(int i = 0; i < 4; ++i)
                    Psi(i, psiIndex) += psi_pq(i, ind);
            }
        }while(ind.hasNext());
    }
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkFirstEquation(const CJet &Psi, const CJet &H, const CJet &N)
{
    auto LHS = fromToDegree(operatorL(projR(reminderPart(Psi)), N, lambda), 0, iterations+1);
    auto RHS = fromToDegree(projR(H), 0, iterations+1);
    log<Diagnostic>("first equation (LHS - RHS):");
    log<Diagnostic>(toString(jetSubstraction(LHS, RHS)));
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveSecondEquation(CJet &N, CJet &B, const CJet &H)
{
    if(iterations%2 == 1)
        return;

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
            b2(0, index) = (h[2](0, index) + h[3](0, index)) / Complex(2, 0);
            
        }while(index.hasNext());
    }

    // calculate N and B (only iterations+1 degree)
    // 2*deg + 1 == iterations+1
    int deg = iterations/2;

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

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkSecondEquation(const CJet &N, const CJet &B, const CJet &H)
{
    auto LHS = fromToDegree(jetAddition(reminderPart(N), reminderPart(B)), 0, iterations+1);
    auto RHS = fromToDegree(projP(H), 0, iterations+1);
    log<Diagnostic>("second equation (LHS - RHS):");
    log<Diagnostic>(toString(jetSubstraction(LHS, RHS)));
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkNormalFormEquality(const PseudoNormalForm &normalForm)
{
    log<Diagnostic>("Normal form condition LHS:");
    auto LHS = jetAddition(multiply(D(normalForm.phi), normalForm.n), reminderPart(normalForm.b));
    // log<Diagnostic>(toString(LHS));
    // log<Diagnostic>("F: " + toString(F_taylorSeries));
    // log<Diagnostic>("Phi: " + toString(normalForm.phi));
}
