#include "PolynomialOf4Variables.h"

#include <iostream>
#include <sstream>
#include "capd/capdlib.h"

using namespace std;
using namespace capd;

PolynomialOf4Variables::PolynomialOf4Variables(int _degree, Complex fillWith) : PolynomialOf4Variables(_degree, _degree, fillWith) {}

PolynomialOf4Variables::PolynomialOf4Variables(int _degree, int _futureDegree, Complex fillWith) : degree(_degree), futureDegree(_futureDegree)
{
    coefficients.reserve(futureDegree+1);
    for(int i = 0; i <= degree; ++i)
    {
        coefficients.emplace_back();
        coefficients.back().reserve(futureDegree+1-i);
        for(int j = 0; j <= degree-i; ++j)
        {
            coefficients.back().emplace_back();
            coefficients.back().back().reserve(futureDegree+1-i-j);
            for(int k = 0; k <= degree-i-j; ++k)
            {
                coefficients.back().back().emplace_back(degree+1-i-j-k, fillWith);
                coefficients.back().back().back().reserve(futureDegree+1-i-j-k);
            }
        }
    }
}

void PolynomialOf4Variables::extend(Complex fillWith)
{
    int degreeLeft = futureDegree-degree;
    if(degreeLeft < 0) degreeLeft = 0;

    { // i = degree+1, j=k=l=0
        coefficients.emplace_back();
        coefficients.back().reserve(degreeLeft);
        coefficients.back().emplace_back();
        coefficients.back().back().reserve(degreeLeft);
        coefficients.back().back().emplace_back(1, fillWith);
        coefficients.back().back().back().reserve(degreeLeft);
    }

    // i < degree+1
    for(int i = 0; i <= degree; ++i)
    {
        { // j = degree+1-i, k=l=0
            coefficients[i].emplace_back();
            coefficients[i].back().reserve(degreeLeft);
            coefficients[i].back().emplace_back(1, fillWith);
            coefficients[i].back().back().reserve(degreeLeft);
        }

        // j < degree+1-i
        for(int j = 0; j <= degree-i; ++j)
        {
            { // k = degree+1-i-j, l=0
                coefficients[i][j].emplace_back(1, fillWith);
                coefficients[i][j].back().reserve(degreeLeft);
            }

            // k < degree+1-i-j
            for(int k = 0; k <= degree-i-j; ++k)
                coefficients[i][j][k].push_back(fillWith);
        }
    }

    degree++;
}

string PolynomialOf4Variables::toString(string var1, string var2, string var3, string var4) const
{
    const Multiindex zero({0, 0, 0, 0});

    stringstream ss;
    for(int deg = 0; deg <= degree; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            Complex coeff = getCoeff(index);
            if(coeff != Complex(0, 0))
            {
                if(index != zero)
                {
                    if(coeff != Complex(1, 0))
                    {
                        if(coeff == Complex(-1, 0)) ss << "-";
                        else ss << coeff;
                    }
                }
                else ss << coeff;

                if(index[0] > 0) ss << var1<<(index[0] == 1 ? "" : "^"+to_string(index[0]));
                if(index[1] > 0) ss << var2<<(index[1] == 1 ? "" : "^"+to_string(index[1]));
                if(index[2] > 0) ss << var3<<(index[2] == 1 ? "" : "^"+to_string(index[2]));
                if(index[3] > 0) ss << var4<<(index[3] == 1 ? "" : "^"+to_string(index[3]));

                ss << " + ";
            }
        }while(index.hasNext());
    }

    string str = ss.str();
    if(str.length() == 0)
        str = "0";

    return str.substr(0, str.length()-3);
}

void PolynomialOf4Variables::setCoeff(const Multiindex &index, Complex value)
{
    validate_indices(index);
    coefficients[index[0]][index[1]][index[2]][index[3]] = value;
}

Complex PolynomialOf4Variables::getCoeff(const Multiindex &index) const
{
    validate_indices(index);
    return coefficients[index[0]][index[1]][index[2]][index[3]];
}

void PolynomialOf4Variables::setCoeff(int i, int j, int k, int l, Complex value)
{
    setCoeff(Multiindex({i, j, k, l}), value);
}

Complex PolynomialOf4Variables::getCoeff(int i, int j, int k, int l) const
{
    return getCoeff(Multiindex({i, j, k, l}));
}

PolynomialOf4Variables PolynomialOf4Variables::operator-(PolynomialOf4Variables const &obj) const
{
    PolynomialOf4Variables result(degree);
    for(int deg = 0; deg <= degree; ++deg)
    {
        Multiindex index({deg, 0, 0, 0});
        do
        {
            result.setCoeff(index, getCoeff(index) - obj.getCoeff(index));
        }while(index.hasNext());
    }
    return result;
}

void PolynomialOf4Variables::validate_indices(const Multiindex &index) const
{
    if(index.dimension() != 4)
        throw runtime_error("Indices invalid for polynomial of 4 variables");

    if(index[0]+index[1]+index[2]+index[3] > degree)
        throw runtime_error("Indices invalid for polynomial of degree "+to_string(degree)+".");
}

PolynomialOf4Variables4::PolynomialOf4Variables4(int _degree, Complex fillWith)
{
    for(int i = 0; i < 4; ++i)
        subfunctions[i] = PolynomialOf4Variables(_degree, fillWith);
}

PolynomialOf4Variables4::PolynomialOf4Variables4(int _degree, int _futureDegree, Complex fillWith)
{
    for(int i = 0; i < 4; ++i)
        subfunctions[i] = PolynomialOf4Variables(_degree, _futureDegree, fillWith);
}

PolynomialOf4Variables &PolynomialOf4Variables4::operator()(int index)
{
    if(index >= 4)
        throw new runtime_error("Index " + to_string(index) + " invalid for 4 dimensional function.");
    return subfunctions[index];
}

const PolynomialOf4Variables &PolynomialOf4Variables4::operator()(int index) const
{
    if(index >= 4)
        throw new runtime_error("Index " + to_string(index) + " invalid for 4 dimensional function.");
    return subfunctions[index];
}

std::string PolynomialOf4Variables4::toString(std::string var1, std::string var2, std::string var3, std::string var4) const
{
    string s = "{\n";
    for(int i = 0; i < 4; ++i)
        s += subfunctions[i].toString(var1, var2, var3, var4) + ",\n";

    s += "}";
    return s;
}
