#include "../typedefs.h"
#include "NormalFormFinder.h"
#include "PseudoNormalForm.h"
#include "HelperFunctions.h"

#include "capd/capdlib.h"
#include "capd/capdAlglib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

using namespace capd;
using namespace std;

PseudoNormalForm NormalFormFinder::calculatePseudoNormalForm(const capd::DMap &f, const capd::DVector &point)
{
    if(f.dimension() != 4 || f.imageDimension() != 4)
        throw runtime_error("Dimensions not supported.");

    // TODO: move system to set point to X=0

    // TODO: diagonalize matrix

    // TODO: check if saddle-focus or saddle-center in X=0

    auto taylorSeries = getTaylorSeries(f, degree);
    cout << taylorSeries.toString() << endl;

    DMatrix lambda;
    PolynomialOf4Variables4 reminder(degree);
    getLinearPartWithReminder(taylorSeries, &lambda, &reminder);

    auto eigenValues = getEigenvalues(lambda);
    
    cout << "rozbicie na czesci:" << endl;
    cout << lambda << endl;
    cout << reminder.toString("x", "y", "z", "u") << endl;

    cout << "wartosci wlasne:" << endl;
    cout << eigenValues << endl;

    PseudoNormalForm result;
    for(int i = 0; i < degree; ++i)
    {
        nextIteration(result);
    }
    return result;
}

void NormalFormFinder::nextIteration(PseudoNormalForm &result)
{
    
}