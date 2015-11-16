#ifndef HISTS_H
#define HISTS_H
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
class hists
{
public:
    hists(TFile *file);
    TH2D *histoAUp;
    TH2D *histoUp;
    TH2D *histoGluon;
    TH2D *histoADown;
    TH2D *histoDown;
};

#endif // HISTS_H
