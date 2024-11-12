#ifndef _POLYNOMIAL_OF_4_VARIABLES_H_
#define _POLYNOMIAL_OF_4_VARIABLES_H_

#include <vector>
#include <string>
#include "capd/capdlib.h"

typedef std::vector<std::vector<double>> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

// polynomial R^4 -> R
class PolynomialOf4Variables
{
    private:
        int degree;
        int futureDegree;
        Vector4D coefficients;

        void validate_indices(const capd::Multiindex &index) const;

    public:
        PolynomialOf4Variables(int _degree = 0, double fillWith = 0);

        PolynomialOf4Variables(int _degree, int _futureDegree, double fillWith = 0);

        std::string toString(std::string var1, std::string var2, std::string var3,std::string var4) const;

        int getDegree() const { return degree; };

        // make degree greater by 1
        void extend(double fillWith = 0);

        void set_coeff(const capd::Multiindex &index, double value);

        double get_coeff(const capd::Multiindex &index) const;
};

// polynomial R^4 -> R4
class PolynomialOf4Variables4
{
    private:
        PolynomialOf4Variables subfunctions[4];
    public:
        PolynomialOf4Variables4(int _degree, double fillWith = 0);

        PolynomialOf4Variables4(int _degree, int _futureDegree, double fillWith = 0);

        PolynomialOf4Variables& operator()(int index);

        std::string toString(std::string var1, std::string var2, std::string var3,std::string var4) const;
};

#endif