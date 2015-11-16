#include "hists.h"
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
hists::hists(TFile *file)
{
    file->cd();
    file->ls();
    histoAUp = (TH2D*) file->Get("histoAUp");
    histoUp = (TH2D*) file->Get("histoUp");
    histoGluon = (TH2D*) file->Get("histoGluon");
    histoADown = (TH2D*) file->Get("histoADown");
    histoDown = (TH2D*) file->Get("histoDown");
}
