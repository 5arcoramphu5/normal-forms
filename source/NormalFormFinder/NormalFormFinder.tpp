#pragma once

#include "NormalFormFinder.hpp"
#include "PseudoNormalForm.h"
#include "helperFunctions.hpp"
#include "../logging/logging.hpp"
#include "../containers/Polynomial.hpp"
#include "../containers/PolynomialMatrix.hpp"

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

#include <vector>

template<LoggerType Logger>
NormalFormFinder<Logger>::NormalFormFinder(int _degree, const CMap &_f, const CVector &fixedPoint, const CMatrix &diagonalDerivative, const CMatrix &diagonalizationMatrix, const CMatrix diagonalizationMatrixInverse)
    : degree(_degree), f(_f), p(fixedPoint), lambda(diagonalDerivative), J(diagonalizationMatrix), invJ(diagonalizationMatrixInverse)
{
    if(f.dimension() != 4 || f.imageDimension() != 4)
        throw std::runtime_error("Dimensions not supported.");
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::calculatePseudoNormalForm()
{
    Polynomial<capd::Complex> f_taylorSeries = getTaylorSeries(f, degree+1);
    
    auto invJ_P = toPolynomial(invJ, p, degree+1); // F(X) = J f(p + J^-1 X)
    F_taylorSeries = Polynomial<capd::Complex>(4, 4, invJ_P.degree() * f_taylorSeries.degree());
    polynomialComposition(f_taylorSeries, invJ_P, F_taylorSeries);
    F_taylorSeries = J * F_taylorSeries;

    log<Debug>("F:\n", F_taylorSeries);

    F_reminder = F_taylorSeries.fromToDegree(2, degree+1);
    log<Debug>("F reminder:\n", F_reminder);

    PointType pointType = getPointType(lambda, lambda1, lambda2);
    if(pointType == PointType::Unsupported)
        throw std::runtime_error("Type of equilibrium point not supported.");
        
    if(pointType == PointType::SaddleCenter) // TODO: to be deleted
        throw std::runtime_error("Case not implemented yet.");

    log<Debug>("lambda:\n", lambda);
    log<Debug>("lambda1, lambda2: ", lambda1, lambda2);

    log<Diagnostic>("Point type: ", pointType, "\n");

    PseudoNormalForm normalForm = getInitialNormalFormValues();
    setInitialValues();

    for(iterations = 1; iterations <= degree; ++iterations)
    {
        log<Minimal>("--------- iteration:", iterations, "---------");

        nextIteration(normalForm);
        
        log<Diagnostic>("--------------------------------");
        log<Minimal>("Phi:\n", normalForm.getPhi());
        log<Minimal>("N:\n", normalForm.getN());
        log<Minimal>("B:\n", normalForm.getB());
        log<Debug>("a1:\n", a1_reminder);
        log<Debug>("a2:\n", a2_reminder);
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
        
        normalForm.phi(i, capd::Multiindex(4, indexArr)) = 1; // Phi(1) = Id
        normalForm.n(i, capd::Multiindex(4, indexArr)) = lambda[i][i]; // N(1) = linear part (diagonal)
    }

    return normalForm;
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::setInitialValues()
{
    a1_reminder = Polynomial<capd::Complex>(1, 2, degree+1);
    a2_reminder = Polynomial<capd::Complex>(1, 2, degree+1);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::nextIteration(PseudoNormalForm &normalForm)
{       
    Polynomial<capd::Complex> FPhi(4, 4, normalForm.phi.degree());
    polynomialComposition(F_reminder, normalForm.phi, FPhi);
    
    solveFirstEquation(normalForm.phi, FPhi);
    Logger::template enableIf<Diagnostic>( [normalForm, FPhi, this] () 
    { 
        this->checkFirstEquation(normalForm.phi, FPhi, normalForm.n); 
    } );
    
    solveSecondEquation(normalForm.n, normalForm.b, FPhi);
    Logger::template enableIf<Diagnostic>( [normalForm, FPhi, this] () 
    { 
        this->checkSecondEquation(normalForm.n, normalForm.b, FPhi); 
    } );

    Logger::template enableIf<Diagnostic>( [normalForm, this] () 
    { 
        this->checkNormalFormEquality(normalForm); 
    } );
}

template<LoggerType Logger>
typename NormalFormFinder<Logger>::PointType NormalFormFinder<Logger>::getPointType(const CMatrix &diagonalMatrix, capd::Complex &lambda1, capd::Complex &lambda2)
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
void NormalFormFinder<Logger>::solveFirstEquation(Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &H)
{
    Polynomial<capd::Complex> RH = projR(H, iterations+1);

    auto pq_coeffs = pqCoefficients(RH, iterations+1);

    for(auto& [pq, h_pq] : pq_coeffs)
    {
        auto& [p, q] = pq; 

        // calculate only iterations+1 degree
        // 2*deg + p + q == iterations+1
        if( (iterations+1 - p - q) % 2 == 1 ) continue;

        auto _gamma = gamma(p, q, lambda1, lambda2);
        Polynomial<capd::Complex> denominator(4, 2, iterations+1); // gamma + pa_1 + qa_2

        capd::Multiindex zero({0, 0});
        for(int i = 0; i < 4; ++i)
        {
            for(int deg = 0; deg <= denominator.degree(); ++deg)
            {
                capd::Multiindex ind({deg, 0});
                do 
                {
                    denominator(i, ind) = capd::Complex(p, 0) * a1_reminder(0, ind) + capd::Complex(q, 0) * a2_reminder(0, ind);
                }while(ind.hasNext());
            }

            denominator(i, zero) += _gamma[i];
        }

        int deg = (iterations+1 - p - q) / 2;

        Polynomial<capd::Complex> psi_pq = polynomialDivision(h_pq, denominator, deg);

        capd::Multiindex ind({deg, 0});
        do 
        {
            if(ind[0] + p < 0 || ind[1] + q < 0)
                continue;
                
            capd::Multiindex psiIndex({ind[0] + p, ind[0], ind[1] + q, ind[1]});
            
            for(int i = 0; i < 4; ++i)
            {
                if(denominator(i ,zero) != capd::Complex(0, 0))
                    Psi(i, psiIndex) += psi_pq(i, ind);
            }
        }while(ind.hasNext());
    }
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkFirstEquation(const Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &H, const Polynomial<capd::Complex> &N)
{
    auto LHS = operatorL(projR(Psi.reminderPart()), N, lambda).fromToDegree(0, iterations+1);
    auto RHS = projR(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("first equation (LHS - RHS):\n", LHS - RHS);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveSecondEquation(Polynomial<capd::Complex> &N, Polynomial<capd::Complex> &B, const Polynomial<capd::Complex> &H)
{
    if(iterations%2 == 1)
        return;

    // get h1, h2, h3, h4 - parts of P projection of H
    Polynomial<capd::Complex> h[4]; 
    for(int i = 0; i < 4; ++i)
        h[i] = Polynomial<capd::Complex>(1, 2, iterations+1);

    for(int i = 0; i < iterations+1; ++i)
        for(int j = 0; 2*i+2*j < iterations+1; ++j)
        {   
            capd::Multiindex index({i, j});
            h[0](0, index) = H(0, capd::Multiindex({i+1, i, j, j}));
            h[1](0, index) = H(1, capd::Multiindex({i, i+1, j, j}));
            h[2](0, index) = H(2, capd::Multiindex({i, i, j+1, j}));
            h[3](0, index) = H(3, capd::Multiindex({i, i, j, j+1}));
        }

    // calculate a1, a2, b1, b2
    Polynomial<capd::Complex> b1(1, 2, iterations+1), b2(1, 2, iterations+1);

    for(int deg = 0; deg <= iterations+1; ++deg)
    {
        capd::Multiindex index({deg, 0});
        do 
        {
            a1_reminder(0, index) = (h[0](0, index) - h[1](0, index)) / capd::Complex(2, 0);
            a2_reminder(0, index) = (h[2](0, index) - h[3](0, index)) / capd::Complex(2, 0);
            b1(0, index) = (h[0](0, index) + h[1](0, index)) / capd::Complex(2, 0);
            b2(0, index) = (h[2](0, index) + h[3](0, index)) / capd::Complex(2, 0);
            
        }while(index.hasNext());
    }

    // calculate N and B (only iterations+1 degree)
    // 2*deg + 1 == iterations+1
    int deg = iterations/2;

    capd::Multiindex index({deg, 0});
    do 
    {
        N(0, capd::Multiindex({index[0]+1, index[0], index[1], index[1]})) = a1_reminder(0, index);
        N(1, capd::Multiindex({index[0], index[0]+1, index[1], index[1]})) = -a1_reminder(0, index);
        N(2, capd::Multiindex({index[0], index[0], index[1]+1, index[1]})) = a2_reminder(0, index);
        N(3, capd::Multiindex({index[0], index[0], index[1], index[1]+1})) = -a2_reminder(0, index);

        B(0, capd::Multiindex({index[0]+1, index[0], index[1], index[1]})) = b1(0, index);
        B(1, capd::Multiindex({index[0], index[0]+1, index[1], index[1]})) = b1(0, index);
        B(2, capd::Multiindex({index[0], index[0], index[1]+1, index[1]})) = b2(0, index);
        B(3, capd::Multiindex({index[0], index[0], index[1], index[1]+1})) = b2(0, index);
        
    }while(index.hasNext());
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkSecondEquation(const Polynomial<capd::Complex> &N, const Polynomial<capd::Complex> &B, const Polynomial<capd::Complex> &H)
{
    auto LHS = (N.reminderPart() + B.reminderPart()).fromToDegree(0, iterations+1);
    auto RHS = projP(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("second equation (LHS - RHS):\n", LHS - RHS);
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkNormalFormEquality(const PseudoNormalForm &normalForm)
{
    auto LHS = D(normalForm.phi) * normalForm.n + normalForm.b.reminderPart();
    Polynomial<capd::Complex> RHS(4, 4, normalForm.phi.degree());
    polynomialComposition(F_reminder, normalForm.phi, RHS);

    log<Diagnostic>("Normal form condition (LHS - RHS):\n", LHS - RHS);
}
