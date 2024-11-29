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

    // TODO: move system to set point to X=0
}

PseudoNormalForm NormalFormFinder::calculatePseudoNormalForm()
{
    auto taylorSeries = getTaylorSeries(f, degree);

    CMatrix linearPart;
    CJet reminder(4, 4, degree);
    getLinearPartWithReminder(taylorSeries, &linearPart, &reminder);

    // TODO: diagonalize matrix

    Complex lambda1, lambda2;
    PointType pointType = getPointType(linearPart, &lambda1, &lambda2);
    if(pointType == PointType::Unsupported)
        throw new runtime_error("Type of equilibrium point not supported.");

    PseudoNormalForm result(degree, taylorSeries);
    for(int i = 0; i < degree; ++i)
    {
        nextIteration(&result);
    }
    return result;
}

void NormalFormFinder::nextIteration(PseudoNormalForm *result)
{
    
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
