#include "project.h"
#include <TCanvas.h>
#include <TSystem.h>
#include <TH1.h>
#include <TProfile.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TDirectory.h>
#include <iostream>
#include <math.h>

#include "parton.h"
TFile *file = new TFile("pdfHistograms.root","read");
file->cd();
// print what is in the file:
file->ls();

// get the histograms out. There are five.

TH2D *histoAUp = (TH2D*) file->Get("histoAUp");
TH2D *histoUp = (TH2D*) file->Get("histoUp");
TH2D *histoGluon = (TH2D*) file->Get("histoGluon");
TH2D *histoADown = (TH2D*) file->Get("histoADown");
TH2D *histoDown = (TH2D*) file->Get("histoDown");
TRandom3 random;

void project(){ // function with the same name as the .C file is automatically executed by root.

    Int_t ybin = histoDown->GetYaxis()->FindBin(90);
    TH1D* histo= histoDown->ProjectionX("Up",ybin,ybin+1);
    //histoGluon->FitSlicesX(0, ybin, ybin,0);
    //gDirectory->ls();
    //TH1D* histo = (TH1D*)gDirectory->Get("histoGluon_0");

    TCanvas *c1 = new TCanvas("c1","The Uniform random example");
    /* Step:
     * 1. Create 2 partons with 2 x
     * 2. Let partons create W/Z boson somehow
     * 3. Let W/Z boson decay
     */
    double _q=90.;
    //step 1: generate 2 x'es
    double part1=Parton(random);    

    c1->cd();//make sure you plot on the right canvas
    histoDown->Draw();
    histoDown->SetMinimum(0);
    c1->Update();

}
//cross section
double sigma(){
    //to be filled in
    return 1.0;

}

Parton make_parton(){
    double x=random.Uniform(0,1);
    Int_t ybin = histoDown->GetYaxis()->FindBin(90);
    Int_t xbin = histoDown->GetXaxis()->FindBin(-10*log10(x));
    doube fAUp =
}



//void fillrandomUniform(Int_t ntimes) {
//   //Fill a 1-D histogram from a uniform function
   
//    TCanvas *c1 = new TCanvas("c1","The Uniform random example");
////    TCanvas *c2 = new TCanvas("c2","the cross-check");

//    gBenchmark->Start("fillrandom");// this measures the time between now and the next gBenchmark of "fillrandom"

//    // create a histogram
//    TH1F *histo = new TH1F("histo","histogram of uniform distribution",100,0,1);
//    //TH3F *histo3d = new TH3F("histo3d","3D histogram for check of random generator",1000,0,1,1000,0,1,1000,0,1);
    
//   // create a random generator object. Advisable one to use is TRandom3. Other options are TRandom2 and TRandom1. DO NOT USE TRandom (so without a number) as it has bugs and is old
//   TRandom3 random;

//    Int_t iterator;
//    Float_t thisran,prevran,beforeran;
//    for(iterator=0; iterator<ntimes; iterator++){
//        // this is where you retrieve the random number
//        if(iterator%1000==0)
//            std::cout << iterator << "/" << ntimes << std::endl;
//        thisran = random.Uniform(0,1);
//        // and fill it into histograms
//        histo->Fill(thisran);
//      //  histo3d->Fill(thisran,prevran,beforeran);
//        // save these for next loop:
//        beforeran=prevran;
//        prevran=thisran;
        
//    }

//   // update the canvas so the plot shows
//    c1->cd();//make sure you plot on the right canvas
//    histo->Draw();
//    histo->SetMinimum(0);
//    c1->Update();
//   // c2->cd();// make sure you plot on the right canvas
//    //histo3d->Draw();
//    //c2->Update();

//   gBenchmark->Show("fillrandom"); // print out how long it took to run the program
//}
