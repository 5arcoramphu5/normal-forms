#include "PolynomialOf4Variables.h"

#include <iostream>
#include <sstream>
using namespace std;

PolynomialOf4Variables::PolynomialOf4Variables(int _degree, double fillWith) : PolynomialOf4Variables(_degree, _degree, fillWith) {}

PolynomialOf4Variables::PolynomialOf4Variables(int _degree, int _futureDegree, double fillWith) : degree(_degree), futureDegree(_futureDegree)
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

void PolynomialOf4Variables::extend(double fillWith)
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
    stringstream ss;
    for(int i = 0; i <= degree; ++i)
        for(int j = 0; j <= degree-i; ++j)
            for(int k = 0; k <= degree-i-j; ++k)
                for(int l = 0; l <= degree-i-j-k; ++l)
                {
                    if(i == 0 && j == 0 && k == 0 && l == 0 && coefficients[i][j][k][l] != 0)
                    {
                        ss << coefficients[i][j][k][l] << " + ";
                    }
                    else if(coefficients[i][j][k][l] != 0)
                    {
                        if(coefficients[i][j][k][l] != 1)
                        {
                            if(coefficients[i][j][k][l] == -1)
                                ss << "-";
                            else
                                ss << coefficients[i][j][k][l];
                        }
                        if(i > 0) ss << var1<<(i == 1 ? "" : "^"+to_string(i));
                        if(j > 0) ss << var2<<(j == 1 ? "" : "^"+to_string(j));
                        if(k > 0) ss << var3<<(k == 1 ? "" : "^"+to_string(k));
                        if(l > 0) ss << var4<<(l == 1 ? "" : "^"+to_string(l));

                        ss << " + ";
                    }
                }
    string str = ss.str();
    return str.substr(0, str.length()-3);
}

void PolynomialOf4Variables::set_coeff(int i, int j, int k, int l, double value)
{
    validate_indices(i, j, k, l);
    coefficients[i][j][k][l] = value;
}

double PolynomialOf4Variables::get_coeff(int i, int j, int k, int l) const
{
    validate_indices(i, j, k, l);
    return coefficients[i][j][k][l];
}

void PolynomialOf4Variables::validate_indices(int i, int j, int k, int l) const
{
    if(i+j+k+l > degree)
        throw runtime_error("Indices invalid for polynomial of degree "+to_string(degree)+".");
}