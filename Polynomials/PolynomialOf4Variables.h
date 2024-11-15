#ifndef _POLYNOMIAL_OF_4_VARIABLES_H_
#define _POLYNOMIAL_OF_4_VARIABLES_H_

#include <vector>
#include <string>
#include "capd/vectalg/lib.h"

typedef std::vector<std::vector<double>> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

// polynomial C^4 -> C
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

        void setCoeff(const capd::Multiindex &index, double value);

        double getCoeff(const capd::Multiindex &index) const;

        void setCoeff(int i, int j, int k, int l, double value);

        double getCoeff(int i, int j, int k, int l) const;

        PolynomialOf4Variables operator-(PolynomialOf4Variables const& obj) const;
};

// polynomial R^4 -> R4
class PolynomialOf4Variables4
{
    private:
        PolynomialOf4Variables subfunctions[4];
    public:
        PolynomialOf4Variables4(int _degree = 0, double fillWith = 0);

        PolynomialOf4Variables4(int _degree, int _futureDegree, double fillWith = 0);

        PolynomialOf4Variables& operator()(int index);

        const PolynomialOf4Variables& operator()(int index) const;

        std::string toString(std::string var1, std::string var2, std::string var3,std::string var4) const;

        int getDegree() const { return subfunctions[0].getDegree(); };
};

#endif