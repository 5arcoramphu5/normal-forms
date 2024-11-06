#ifndef _POLYNOMIAL_OF_4_VARIABLES_H_
#define _POLYNOMIAL_OF_4_VARIABLES_H_

#include <vector>
#include <string>

typedef std::vector<std::vector<double>> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

class PolynomialOf4Variables
{
    private:
        int degree;
        int futureDegree;
        Vector4D coefficients;

        void validate_indices(int i, int j, int k, int l) const;

    public:
        PolynomialOf4Variables(int _degree, double fillWith = 0);

        PolynomialOf4Variables(int _degree, int _futureDegree, double fillWith = 0);

        std::string toString(std::string var1, std::string var2, std::string var3,std::string var4) const;

        int getDegree() const { return degree; };

        // make degree greater by 1
        void extend(double fillWith = 0);

        void set_coeff(int i, int j, int k, int l, double value);

        double get_coeff(int i, int j, int k, int l) const;
};

#endif