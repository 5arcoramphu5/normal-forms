#include "debugUtils.h"

using namespace std;
using namespace capd;

string toString(CJet polynomial, string var1, string var2, string var3, string var4)
{    
    string vars[4] = {var1 , var2, var3, var4};
    const Multiindex zero({0, 0, 0, 0});
    string result = "\t";

    for(int i = 0; i < polynomial.dimension(); ++i)
    {
        stringstream ss;
        for(int deg = 0; deg <= polynomial.degree(); ++deg)
        {
            Multiindex index({deg, 0, 0, 0});
            do
            {
                Complex coeff = polynomial(i, index);
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
                    
                    for(int j = 0; j < 4; ++j)
                        if(index[j] > 0) ss << vars[j] <<(index[j] == 1 ? "" : "^"+to_string(index[j]));

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

void checkPseudoNormalCondition(const PseudoNormalForm &normalForm)
{
    // TODO
}
