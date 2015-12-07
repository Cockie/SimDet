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
double zwidth=2.5;
double emass=0.000510998946;
double mumass=0.105658;

void project(){ // function with the same name as the .C file is automatically executed by root.
    rng.SetSeed(0);
    //Int_t ybin = histos.histoDown->GetYaxis()->FindBin(90);
   // TH1D* histo= histos.histoDown->ProjectionX("Up",ybin,ybin+1);
    //histoGluon->FitSlicesX(0, ybin, ybin,0);
    //gDirectory->ls();
    //TH1D* histo = (TH1D*)gDirectory->Get("histoGluon_0");

    TCanvas *c1 = new TCanvas("c1","electron");
    TCanvas *c2 = new TCanvas("c2","muon");
    TH1D *histoZ = new TH1D("Z","electron",30,-3,3);
    TH1D *histoZb = new TH1D("Z","muon",30,-3,3);
    TH1D *histoWp = new TH1D("W+","electron",30,-3,3);
    TH1D *histoWpb = new TH1D("W+","muon",30,-3,3);
    TH1D *histoWm = new TH1D("W-","electron",30,-3,3);
    TH1D *histoWmb = new TH1D("W-","muon",30,-3,3);
    /* Step:
     * 1. Create 2 partons with 2 x
     * 2. Let partons create W/Z boson somehow
     * 3. Let W/Z boson decay
     */
    double Etot=6500;
    double fmax=get_max();
    int wp=0;
    int wm=0;
    int z=0;
    for(q=90;q<=90;++q){
        for(int n=0;n<100000;++n){
            if(n%1000==0){
                cout<<"n="<<n<<endl;}
            //step 1: generate 2 x'es
            Parton part1=make_parton(fmax);
            Parton part2=make_parton(fmax);
            //cout<<"part1 "<<part1._type<<" "<<part1._x<<endl;
            //cout<<"part2 "<<part2._type<<" "<<part2._x<<endl;
            TLorentzVector vec1;
            vec1.SetPx(Etot*part1._x);
            vec1.SetE( Etot*part1._x );
            TLorentzVector vec2;
            vec2.SetPx(-Etot*part2._x);
            vec2.SetE( Etot*part2._x );
            TLorentzVector vecW=vec1+vec2;
            //cout<<vecW.Px()<<" "<<vecW.E()<<endl;
            if (vecW.E()>=75){
                //correct types to make a Z
                if((part1._type=="Up" && part2._type=="AUp")||(part1._type=="AUp" && part2._type=="Up")||(part1._type=="Down" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Down")){
                    //enough mass to make a Z
                    if (vecW.E()>=91.1876){
                        //make Z
                        ++z;
                        //get the boost vector
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        double x=acos(2*rng.Uniform(0,1)-1);
                        double E=vecW.E()/2;
                        double P=sqrt(E*E-emass*emass);
                        //double P=sqrt(E*E-0.1*0.1);
                        TLorentzVector veca(0, P*sin(x),P*cos(x),E );
                        TLorentzVector vecb(0, -P*sin(x),-P*cos(x),E );

                        veca.Boost(boost);
                        vecb.Boost(boost);
                        histoZ->Fill(veca.PseudoRapidity());
                        histoZ->Fill(vecb.PseudoRapidity());
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        P=sqrt(E*E-mumass*mumass);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );
                        vecb=TLorentzVector(0, -P*sin(x),-P*cos(x),E );

                        veca.Boost(boost);
                        vecb.Boost(boost);
                        histoZb->Fill(veca.PseudoRapidity());
                        histoZb->Fill(vecb.PseudoRapidity());


                    }
                    //correct types for W+
                }else if((part1._type=="Up" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Up")){
                    if (vecW.E()>=80.385){
                        //make W+
                        ++wp;
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        double x=acos(2*rng.Uniform(0,1)-1);
                        double E=(vecW.E()*vecW.E()+emass*emass)/(2*vecW.E());
                        double P=sqrt(E*E-emass*emass);
                        //double P=sqrt(E*E-0.1*0.1);
                        TLorentzVector veca(0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);;
                        histoWp->Fill(veca.PseudoRapidity());
                        /*//second electron
                        x=acos(2*rng.Uniform(0,1)-1);
                        veca=TLorentzVector (0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);;
                        histoWp->Fill(veca.PseudoRapidity());*/
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+mumass*mumass)/(2*vecW.E());
                        P=sqrt(E*E-mumass*mumass);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);
                        histoWpb->Fill(veca.PseudoRapidity());
                    }
                }else if((part1._type=="AUp" && part2._type=="Down")||(part1._type=="Down" && part2._type=="AUp")){
                    if (vecW.E()>=80.385){
                        //make W-
                        ++wm;
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        double x=acos(2*rng.Uniform(0,1)-1);
                        double E=(vecW.E()*vecW.E()+emass*emass)/(2*vecW.E());
                        double P=sqrt(E*E-emass*emass);
                        TLorentzVector veca(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);;
                        histoWm->Fill(veca.PseudoRapidity());
                        /*//second
                        x=acos(2*rng.Uniform(0,1)-1);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);;
                        histoWm->Fill(veca.PseudoRapidity());*/
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+mumass*mumass)/(2*vecW.E());
                        P=sqrt(E*E-mumass*mumass);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);
                        histoWmb->Fill(veca.PseudoRapidity());
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
    /*c1->cd();//make sure you plot on the right canvas
    histoWm->Draw();
    histoWm->SetMinimum(0);
    c1->Update();

    c2->cd();//make sure you plot on the right canvas
    histoWmb->Draw();
    histoWmb->SetMinimum(0);
    c2->Update();*/
    cout<<"Z="<<z<<" Wmin="<<wm<<" Wplus="<<wp<<endl;
    TH1D *ass = new TH1D("W-","electron",30,-3,3);
    TH1D *ass_sub = new TH1D("W-","electron",30,-3,3);
    TH1D *ass_sum = new TH1D("W-","electron",30,-3,3);
    ass_sub->Add(histoWm, histoWp,-1,2);
    ass_sum->Add(histoWm, histoWp,1,2);
    ass->Divide(ass_sub, ass_sum);
    c1->cd();//make sure you plot on the right canvas
    ass->Draw();
    ass->SetMinimum(0);
    c1->Update();
    TH1D *assb = new TH1D("W-","muon",30,-3,3);
    TH1D *assb_sub = new TH1D("W-","muon",30,-3,3);
    TH1D *assb_sum = new TH1D("W-","muon",30,-3,3);
    assb_sub->Add(histoWmb, histoWpb,-1,2);
    assb_sum->Add(histoWmb, histoWpb,1,2);
    assb->Divide(assb_sub, assb_sum);
    c2->cd();//make sure you plot on the right canvas
    assb->Draw();
    assb->SetMinimum(0);
    c2->Update();
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
