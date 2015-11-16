#ifndef PARTON_H
#define PARTON_H
#include <string>
#include <TRandom3.h>
using namespace std;
class Parton
{
private:
public:
    double _x;
    string _type;
    Parton(double x, string type);
    Parton();
};

#endif // PARTON_H
