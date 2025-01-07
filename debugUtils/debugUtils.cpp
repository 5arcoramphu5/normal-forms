#include "debugUtils.h"
#include "../containers/Polynomial.cpp"

using namespace std;
using namespace capd;

// print in format compatible with Mathematica
// ostream& operator<<(ostream& os, const Complex& c)
// {
//     os << "(" << c.real() << " + " << c.imag() << "I)";
//     return os;
// }

template<ArithmeticType Coeff>
string toString(Polynomial<Coeff> polynomial, const string vars[])
{    
    const Multiindex zero(polynomial.dimension());
    string result = "\t";

    for(int i = 0; i < polynomial.imageDimension(); ++i)
    {
        stringstream ss;
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            Multiindex index(polynomial.dimension());
            index[0] = deg;

            do
            {
                Complex coeff = polynomial(i, index);
                if(coeff != (Coeff)0)
                {
                    if(index != zero)
                    {
                        if(coeff != Complex(1, 0))
                        {
                            if(coeff == Complex(-1, 0)) ss << "-";
                            else ss << coeff << " ";
                        }
                    }
                    else ss << coeff << " ";
                    
                    for(int j = 0; j < 4; ++j)
                        if(index[j] > 0) ss << vars[j] <<(index[j] == 1 ? "" : "^"+to_string(index[j])) << " ";

                    ss << " + ";
                }
            }while(index.hasNext());
        }

        string str = ss.str();
        if(str.length() == 0)
            str = "0";
        else 
            str = str.substr(0, str.length()-3);

        result += str + "\n\t";
    }

    return result;
}

template<ArithmeticType Coeff>
std::string toCoefficientString(Polynomial<Coeff> polynomial)
{
    stringstream result;
    result << "{\n";
    for(int i = 0; i < polynomial.imageDimension(); ++i)
    {
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            Multiindex index(polynomial.dimension());
            index[0] = deg;

            do
            {
                Complex coeff = polynomial(i, index);
                result << coeff << " ";
            }while(index.hasNext());

            result << endl;
        }
        if(i != polynomial.imageDimension()-1)
            result << "," << endl;
    }
    result << "}";

    return result.str();
}