#pragma once

#include "NormalFormFinder.hpp"
#include "helperFunctions.hpp"
#include "../containers/PolynomialMatrix.hpp"

#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

#include <vector>

template<LoggerType Logger>
NormalFormFinder<Logger>::NormalFormFinder(int _degree)
    : degree(_degree), iterations(0)
    {}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::calculatePseudoNormalForm(const Diagonalization<capd::Complex> &diagonalization)
{
    setPointTypeAndLambdas(diagonalization.getLambda());

    if(pointType == PointType::Unsupported)
        throw std::runtime_error("Type of equilibrium point not supported.");
        
    if(pointType == PointType::SaddleCenter)
        throw std::runtime_error("Case not implemented yet.");

    PseudoNormalForm normalForm = getInitialNormalFormValues(diagonalization);

    for(iterations = 1; iterations <= degree; ++iterations)
    {
        log<ProgressIndication>("--------- iteration:", iterations, "/", degree, "---------");

        nextIteration(normalForm, diagonalization);
        
        log<Diagnostic>("--------------------------------");
        log<Minimal>("Phi:\n", normalForm.getPhi());
        log<Minimal>("N:\n", normalForm.getN());
        log<Minimal>("B:\n", normalForm.getB());
    }

    return normalForm;
}

template<LoggerType Logger>
PseudoNormalForm NormalFormFinder<Logger>::getInitialNormalFormValues(const Diagonalization<capd::Complex> &diagonalization)
{
    PseudoNormalForm normalForm(degree+1);

    capd::Multiindex zero({0, 0});
    normalForm.a1(0, zero) = lambda1;
    normalForm.a2(0, zero) = lambda2;

    for(int i = 0; i < 4; ++i)
    {
        capd::Multiindex index({0, 0, 0, 0});
        index[i] = 1;

        normalForm.phi(i, index) = 1; // Phi(1) = Id
        normalForm.n(i, index) = diagonalization.getLambda()[i][i]; // N(1) = linear part (diagonal)
    }
    return normalForm;
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::nextIteration(PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diagonalization)
{
    auto FremPhi = diagonalization.polynomialCompositionWithReminder(normalForm.phi);
    
    solveFirstEquation(normalForm.phi, normalForm.a1, normalForm.a2, FremPhi);
    Logger::template enableIf<Diagnostic>( [normalForm, FremPhi, diagonalization, this] () 
    { 
        this->checkFirstEquation(normalForm.phi, FremPhi, normalForm.n, diagonalization); 
    } );
    
    solveSecondEquation(normalForm.n, normalForm.b, normalForm.a1, normalForm.a2, FremPhi);
    Logger::template enableIf<Diagnostic>( [normalForm, FremPhi, this] () 
    { 
        this->checkSecondEquation(normalForm.n, normalForm.b, FremPhi); 
    } );

    Logger::template enableIf<Diagnostic>( [normalForm, diagonalization, this] () 
    { 
        this->checkNormalFormEquality(normalForm, diagonalization); 
    } );
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::setPointTypeAndLambdas(const CMatrix &lambda)
{
    if(lambda[0][0] != -lambda[1][1] || lambda[2][2] != -lambda[3][3])
    {
        pointType = PointType::Unsupported;
        return;
    }
    
    lambda1 = lambda[0][0];
    lambda2 = lambda[2][2];
    
    if( (lambda1.real() == 0 && lambda1.imag() != 0 && lambda2.imag() == 0 && lambda2.real() != 0) || 
        (lambda1.imag() == 0 && lambda1.real() != 0 && lambda2.real() == 0 && lambda2.imag() != 0) )
    {
        pointType = PointType::SaddleCenter;
    }
    else if( (conj(lambda1) == lambda2 || conj(lambda1) == -lambda2) && lambda1.real() != 0 && lambda1.imag() != 0)
    {
        pointType = PointType::SaddleFocus;
    }
    else
    {
        pointType = PointType::Unsupported;
    }
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveFirstEquation(Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &a1, const Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H)
{
    Polynomial<capd::Complex> RH = projR(H, iterations+1);

    for(int deg1 = 0; deg1 <= iterations+1; ++deg1)
    {
        capd::Multiindex index({deg1, 0, 0, 0});
        do
        {
            auto coeffVector = RH(index);

            int p = index[0] - index[1]; // j - k
            int q = index[2] - index[3]; // l - m
            
            // calculate only iterations+1 degree
            // 2*deg + p + q == iterations+1
            if( (iterations+1 - p - q) % 2 == 1 ) continue;

            // prepare p, q monomial
            capd::Multiindex coefficientIndex({index[1], index[3]}); // (k, m)
            int polyDeg = index[1] + index[3];
            Polynomial<capd::Complex> pqMono(4, 2, polyDeg);
            pqMono(coefficientIndex) = coeffVector;
            
            // prepare the denominator
            auto _gamma = gamma(p, q, lambda1, lambda2);
            Polynomial<capd::Complex> denominator(4, 2, iterations+1); // gamma + pa_1 + qa_2

            capd::Multiindex zero({0, 0});
            for(int i = 0; i < 4; ++i)
            {
                for(int deg2 = 1; deg2 <= denominator.degree(); ++deg2)
                {
                    capd::Multiindex ind({deg2, 0});
                    do 
                    {
                        denominator(i, ind) = capd::Complex(p, 0) * a1(0, ind) + capd::Complex(q, 0) * a2(0, ind);
                    }while(ind.hasNext());
                }

                denominator(i, zero) += _gamma[i];
            }

            // calculate psi_pq
            int currentDegree = (iterations+1 - p - q) / 2;
            Polynomial<capd::Complex> psi_pq = polynomialDivision(pqMono, denominator, currentDegree);

            // get result from p, q form
            capd::Multiindex ind({currentDegree, 0});
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
            }
            while(ind.hasNext());
        }
        while(index.hasNext());
    }
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkFirstEquation(const Polynomial<capd::Complex> &Psi, const Polynomial<capd::Complex> &H, const Polynomial<capd::Complex> &N, const Diagonalization<capd::Complex> &diagonalization)
{
    auto LHS = operatorL(projR(Psi.reminderPart()), N, diagonalization.getLambda()).fromToDegree(0, iterations+1);
    auto RHS = projR(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("first equation (LHS - RHS):\n", LHS - RHS);
}

template<LoggerType Logger>
void NormalFormFinder<Logger>::solveSecondEquation(Polynomial<capd::Complex> &N, Polynomial<capd::Complex> &B, Polynomial<capd::Complex> &a1, Polynomial<capd::Complex> &a2, const Polynomial<capd::Complex> &H)
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

    // calculate only iterations+1 degree
    // 2*deg + 1 == iterations+1
    int deg = iterations/2;
    
    capd::Multiindex index({deg, 0});
    do 
    {
        a1(0, index) = (h[0](0, index) - h[1](0, index)) / capd::Complex(2, 0);
        a2(0, index) = (h[2](0, index) - h[3](0, index)) / capd::Complex(2, 0);
        capd::Complex b1 = (h[0](0, index) + h[1](0, index)) / capd::Complex(2, 0);
        capd::Complex b2 = (h[2](0, index) + h[3](0, index)) / capd::Complex(2, 0);

        N(0, capd::Multiindex({index[0]+1, index[0], index[1], index[1]})) = a1(0, index);
        N(1, capd::Multiindex({index[0], index[0]+1, index[1], index[1]})) = -a1(0, index);
        N(2, capd::Multiindex({index[0], index[0], index[1]+1, index[1]})) = a2(0, index);
        N(3, capd::Multiindex({index[0], index[0], index[1], index[1]+1})) = -a2(0, index);

        B(0, capd::Multiindex({index[0]+1, index[0], index[1], index[1]})) = b1;
        B(1, capd::Multiindex({index[0], index[0]+1, index[1], index[1]})) = b1;
        B(2, capd::Multiindex({index[0], index[0], index[1]+1, index[1]})) = b2;
        B(3, capd::Multiindex({index[0], index[0], index[1], index[1]+1})) = b2;
        
    }while(index.hasNext());
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkSecondEquation(const Polynomial<capd::Complex> &N, const Polynomial<capd::Complex> &B, const Polynomial<capd::Complex> &H)
{
    auto LHS = (N.reminderPart() + B).fromToDegree(0, iterations+1);
    auto RHS = projP(H).fromToDegree(0, iterations+1);
    log<Diagnostic>("second equation (LHS - RHS):\n", LHS - RHS);
}

template <LoggerType Logger>
void NormalFormFinder<Logger>::checkNormalFormEquality(const PseudoNormalForm &normalForm, const Diagonalization<capd::Complex> &diagonalization)
{
    auto LHS = D(normalForm.phi) * normalForm.n + normalForm.b;
    auto RHS = diagonalization.polynomialComposition(normalForm.phi);
    log<Diagnostic>("Normal form condition (LHS - RHS):\n", (LHS - RHS).fromToDegree(0, iterations+1));
}
