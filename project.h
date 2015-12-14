#ifndef FILLRANDOM_HH
#define FILLRANDOM_HH
#include "parton.h"
#include <TLorentzVector.h>
void project();
void fillrandomUniform(Int_t ntimes);
Parton make_parton(double fmax);
double get_max();
inline double res(double E);
inline double P_trans(TLorentzVector v);
#endif
