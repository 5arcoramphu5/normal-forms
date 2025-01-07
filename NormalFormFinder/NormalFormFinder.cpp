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
    Polynomial f_taylorSeries = getTaylorSeries(f, degree+1);
    
    auto invJ_P = toPolynomial(invJ, p, degree+1); // F(X) = J f(p + J^-1 X)
    F_taylorSeries = Polynomial(4, 4, invJ_P.degree() * f_taylorSeries.degree());
    polynomialComposition(f_taylorSeries, invJ_P, F_taylorSeries);
    F_taylorSeries = J * F_taylorSeries;

    log<Debug>("F:");
    log<Debug>(toString(F_taylorSeries));

    F_reminder = F_taylorSeries.fromToDegree(2, degree+1);
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
    a1_reminder = Polynomial(1, 2, degree+1);
    a2_reminder = Polynomial(1, 2, degree+1);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::nextIteration(PseudoNormalForm &normalForm)
{       
    Polynomial FPhi(4, 4, normalForm.phi.degree());
    polynomialComposition(F_reminder, normalForm.phi, FPhi);
    
    solveFirstEquation(normalForm.phi, FPhi);
    checkFirstEquation(normalForm.phi, FPhi, normalForm.n);

    solveSecondEquation(normalForm.n, normalForm.b, FPhi);
    checkSecondEquation(normalForm.n, normalForm.b, FPhi);

    // checkNormalFormEquality(normalForm);
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
void NormalFormFinder<Logger>::solveFirstEquation(Polynomial &Psi, const Polynomial &H)
{
    Polynomial RH = projR(H, iterations+1);

    auto pq_coeffs = pqCoefficients(RH, iterations+1);

    for(auto& [pq, h_pq] : pq_coeffs)
    {
        auto& [p, q] = pq; 

        // calculate only iterations+1 degree
        // 2*deg + p + q == iterations+1
        if( (iterations+1 - p - q) % 2 == 1 ) continue;

        auto _gamma = gamma(p, q, lambda1, lambda2);
        Polynomial denominator(4, 2, iterations+1); // gamma + pa_1 + qa_2

        Multiindex zero({0, 0});
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

            denominator(i, zero) += _gamma[i];
        }

        int deg = (iterations+1 - p - q) / 2;

        Polynomial psi_pq = polynomialDivision(h_pq, denominator, deg);

        Multiindex ind({deg, 0});
        do 
        {
            if(ind[0] + p < 0 || ind[1] + q < 0)
                continue;
                
            Multiindex psiIndex({ind[0] + p, ind[0], ind[1] + q, ind[1]});
            
            for(int i = 0; i < 4; ++i)
            {
                if(denominator(i ,zero) != Complex(0, 0))
                    Psi(i, psiIndex) += psi_pq(i, ind);
            }
        }while(ind.hasNext());
    }
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkFirstEquation(const Polynomial &Psi, const Polynomial &H, const Polynomial &N)
{
    auto LHS = operatorL(projR(Psi.reminderPart()), N, lambda).fromToDegree(0, iterations+1);
    auto RHS = projR(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("first equation (LHS - RHS):");
    log<Diagnostic>(toString(LHS - RHS));
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveSecondEquation(Polynomial &N, Polynomial &B, const Polynomial &H)
{
    if(iterations%2 == 1)
        return;

    // get h1, h2, h3, h4 - parts of P projection of H
    Polynomial h[4]; 
    for(int i = 0; i < 4; ++i)
        h[i] = Polynomial(1, 2, iterations+1);

    for(int i = 0; i < iterations+1; ++i)
        for(int j = 0; 2*i+2*j < iterations+1; ++j)
        {   
            Multiindex index({i, j});
            h[0](0, index) = H(0, Multiindex({i+1, i, j, j}));
            h[1](0, index) = H(1, Multiindex({i, i+1, j, j}));
            h[2](0, index) = H(2, Multiindex({i, i, j+1, j}));
            h[3](0, index) = H(3, Multiindex({i, i, j, j+1}));
        }

    // calculate a1, a2, b1, b2
    Polynomial b1(1, 2, iterations+1), b2(1, 2, iterations+1);

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
void NormalFormFinder<Logger>::checkSecondEquation(const Polynomial &N, const Polynomial &B, const Polynomial &H)
{
    auto LHS = (N.reminderPart() + B.reminderPart()).fromToDegree(0, iterations+1);
    auto RHS = projP(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("second equation (LHS - RHS):");
    log<Diagnostic>(toString(LHS - RHS));
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkNormalFormEquality(const PseudoNormalForm &normalForm)
{
    log<Diagnostic>("Normal form condition LHS:");
    auto LHS = D(normalForm.phi) * normalForm.n + normalForm.b.reminderPart();

    Polynomial RHS(4, 4, normalForm.phi.degree());
    polynomialComposition(F_reminder, normalForm.phi, RHS);

    log<Diagnostic>(toString(LHS - RHS));
}
