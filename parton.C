#include "parton.h"
#include <string>
using namespace std;
Parton::Parton(double x, string type)
{
    _x=x;
    _type=type;

}

Parton::Parton()
{
    _x=0;
    _type="type";

}
