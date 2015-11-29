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
#include <TLorentzVector.h>
#include <iostream>
#include <math.h>
#include "hists.h"
#include "parton.h"
using namespace std;

TFile *file = new TFile("pdfHistograms.root","read");
hists histos(file);
TRandom3 rng;
double q=90.;

void project(){ // function with the same name as the .C file is automatically executed by root.
    rng.SetSeed(0);
    //Int_t ybin = histos.histoDown->GetYaxis()->FindBin(90);
   // TH1D* histo= histos.histoDown->ProjectionX("Up",ybin,ybin+1);
    //histoGluon->FitSlicesX(0, ybin, ybin,0);
    //gDirectory->ls();
    //TH1D* histo = (TH1D*)gDirectory->Get("histoGluon_0");

    TCanvas *c1 = new TCanvas("c1","The Uniform random example");
    TH1D *histoZ = new TH1D("Z","histogram of uniform distribution",200,100,1000);
    TH1D *histoWp = new TH1D("W+","histogram of uniform distribution",200,100,1000);
    TH1D *histoWm = new TH1D("W-","histogram of uniform distribution",200,100,1000);
    /* Step:
     * 1. Create 2 partons with 2 x
     * 2. Let partons create W/Z boson somehow
     * 3. Let W/Z boson decay
     */
    double fmax=get_max();
    for(q=100;q<=1000;++q){
        for(int n=0;n<100000;++n){
            //step 1: generate 2 x'es
            Parton part1=make_parton(fmax);
            Parton part2=make_parton(fmax);
            //cout<<"part1 "<<part1._type<<" "<<part1._x<<endl;
            //cout<<"part2 "<<part2._type<<" "<<part2._x<<endl;
            TLorentzVector vec1;
            vec1.SetPx(q*part1._x);
            vec1.SetE(q*part1._x);
            TLorentzVector vec2;
            vec2.SetPx(-q*part2._x);
            vec2.SetE(q*part2._x);
            TLorentzVector vecW=vec1+vec2;
            //cout<<vecW.Px()<<" "<<vecW.E()<<endl;
            if (vecW.E()>=80.385){
                //correct types to make a Z
                if((part1._type=="Up" && part2._type=="AUp")||(part1._type=="AUp" && part2._type=="Up")||(part1._type=="Down" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Down")){
                    //enough mass to make a Z
                    if (vecW.E()>=91.1876){
                        //make Z
                        histoZ->Fill(q);
                        cout<<"Z "<<endl;
                    }
                    //correct types for W+
                }else if((part1._type=="Up" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Up")){
                    if (vecW.E()>=80.385){
                        //make W+
                        histoWp->Fill(q);
                        cout<<"Wp "<<endl;
                    }
                }else if((part1._type=="AUp" && part2._type=="Down")||(part1._type=="Down" && part2._type=="AUp")){
                    if (vecW.E()>=80.385){
                        //make W-
                        histoWm->Fill(q);
                        cout<<"Wm "<<endl;
                    }
                }
            }
        }
    }



    //if (part1._type ==& part2)
    /*TH1D *histo = new TH1D("Up_test","histogram of uniform distribution",200,0,1);
    int prod=0;
    while (prod<5000){
        Parton part=make_parton(fmax);
        if (part._type=="ADown"){
            histo->Fill(part._x);
            prod++;
            if(prod%100==0){
                cout<<prod<<endl;
            }
        }

    }*/
    c1->cd();//make sure you plot on the right canvas
    histoZ->Draw();
    histoZ->SetMinimum(0);
    c1->Update();

}
//cross section
double sigma(){
    //to be filled in
    return 1.0;
}

Parton make_parton(double fmax){
    double x=rng.Uniform(0,1);
    Int_t ybin = histos.histoDown->GetYaxis()->FindBin(q);
    Int_t xbin = histos.histoDown->GetXaxis()->FindBin(-10*log10(x));
    double fAUp = histos.histoAUp->GetBinContent(xbin,ybin);
    double fUp = histos.histoUp->GetBinContent(xbin,ybin);
    double fADown = histos.histoADown->GetBinContent(xbin,ybin);
    double fDown = histos.histoDown->GetBinContent(xbin,ybin);
    double fGluon = histos.histoGluon->GetBinContent(xbin,ybin);

    double f=rng.Uniform(0,fmax);
    if(f>=fAUp+fUp+fADown+fDown+fGluon){
        return make_parton(fmax);
    }else if(f>=fAUp+fUp+fADown+fDown){
        return Parton(x,"Gluon");
    }else if(f>=fAUp+fUp+fADown){
        return Parton(x,"Down");
    }else if(f>=fAUp+fUp){
        return Parton(x,"ADown");
    }else if(f>=fAUp){
        return Parton(x,"Up");
    }else{
        return Parton(x,"AUp");
    }
}

double get_max(){
    Int_t ybin = histos.histoDown->GetYaxis()->FindBin(q);
    Int_t xmax = histos.histoDown->GetXaxis()->GetNbins();
    double fmax=0;
    for(int xbin=0;xbin<=xmax;xbin++){
        double f= histos.histoAUp->GetBinContent(xbin,ybin)+histos.histoUp->GetBinContent(xbin,ybin)+histos.histoADown->GetBinContent(xbin,ybin)+histos.histoDown->GetBinContent(xbin,ybin)+histos.histoGluon->GetBinContent(xbin,ybin);
        if (f>fmax){
            fmax=f;
        }
    }
    cout<<fmax<<endl;
    return fmax;


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
