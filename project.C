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
double limit=3;
int previous_anti=0;
double lowPT=20;


void project(){ // function with the same name as the .C file is automatically executed by root.
    gBenchmark->Reset();
    gBenchmark->Start("fillrandom");
    rng.SetSeed(0);
    //Int_t ybin = histos.histoDown->GetYaxis()->FindBin(90);
   // TH1D* histo= histos.histoDown->ProjectionX("Up",ybin,ybin+1);
    //histoGluon->FitSlicesX(0, ybin, ybin,0);
    //gDirectory->ls();
    //TH1D* histo = (TH1D*)gDirectory->Get("histoGluon_0");
    TCanvas *c1 = new TCanvas("c1","Asymetry electron");
    TCanvas *c2 = new TCanvas("c2","Asymetry muon");
    TCanvas *cz = new TCanvas("cz","Z electron");
    TH1D *histoZ = new TH1D("Ze","Z electron",30,-limit,limit);
    TCanvas *czb = new TCanvas("czb","Z muon");
    TH1D *histoZb = new TH1D("Zm","Z muon",30,-limit,limit);
    TCanvas *cwp = new TCanvas("cwp","Wp electron");
    TH1D *histoWp = new TH1D("W+e","Wp electron",30,-limit,limit);
    TCanvas *cwpb = new TCanvas("cwpb","Wp muon");
    TH1D *histoWpb = new TH1D("W+m","Wp muon",30,-limit,limit);
    TCanvas *cwm = new TCanvas("cwm","Wm electron");
    TH1D *histoWm = new TH1D("W-e","Wm electron",30,-limit,limit);
    TCanvas *cwmb = new TCanvas("cwmb","Wm muon");
    TH1D *histoWmb = new TH1D("W-m","Wm muon",30,-limit,limit);
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
    double E=0;
    double P=0;
    double E1=0;
    double E2=0;
    double P1=0;
    double P2=0;
    double x=0;
    for(q=90;q<=90;++q){
        for(int n=0;n<10000;++n){
            if(n%100==0){
                cout<<"n="<<n<<endl;}
            //step 1: generate 2 x'es
            previous_anti=rng.Integer(2);
            //cout<<previous_anti<<endl;
            Parton part1=make_parton(fmax);
            //cout<<"part1 "<<part1._type<<" "<<part1._x<<endl;
            Parton part2=make_parton(fmax);
            //cout<<"part2 "<<part2._type<<" "<<part2._x<<endl;
            TLorentzVector vecW(0,0,Etot*part1._x-Etot*part2._x,Etot*part1._x+Etot*part2._x);
            //cout<<vecW.Px()<<" "<<vecW.E()<<endl;
            if (vecW.E()>=75){
                double M=sqrt(vecW.E()*vecW.E()-vecW.P()*vecW.P());
                //correct types to make a Z
                if((part1._type=="Up" && part2._type=="AUp")||(part1._type=="AUp" && part2._type=="Up")||(part1._type=="Down" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Down")){
                    //enough mass to make a Z
                    if (M>=91.1876){
                        //make Z
                        ++z;
                        //get the boost vector
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=vecW.E()/2;
                        E1=rng.Gaus(E, res(E));
                        P1=sqrt(E1*E1-emass*emass);
                        E2=rng.Gaus(E, res(E));
                        P2=sqrt(E2*E2-emass*emass);

                        TLorentzVector veca(0, P1*sin(x),P1*cos(x),E1 );
                        TLorentzVector vecb(0, -P2*sin(x),-P2*cos(x),E2 );

                        veca.Boost(boost);
                        vecb.Boost(boost);
                        if(P_trans(veca)>=lowPT){
                            histoZ->Fill(veca.PseudoRapidity());}
                        if(P_trans(vecb)>=lowPT){
                            histoZ->Fill(vecb.PseudoRapidity());}
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        E1=rng.Gaus(E, res(E));
                        P1=sqrt(E1*E1-mumass*mumass);
                        E2=rng.Gaus(E, res(E));
                        P2=sqrt(E2*E2-mumass*mumass);

                        veca=TLorentzVector(0, P1*sin(x),P1*cos(x),E1 );
                        vecb=TLorentzVector(0, -P2*sin(x),-P2*cos(x),E2 );

                        veca.Boost(boost);
                        vecb.Boost(boost);
                        if(P_trans(veca)>=lowPT){
                            histoZb->Fill(veca.PseudoRapidity());}
                        if(P_trans(vecb)>=lowPT){
                            histoZb->Fill(vecb.PseudoRapidity());}


                    }
                    //correct types for W+
                }else if((part1._type=="Up" && part2._type=="ADown")||(part1._type=="ADown" && part2._type=="Up")){
                    if (M>=80.385){
                        //make W+
                        ++wp;
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+emass*emass)/(2*vecW.E());
                        E=rng.Gaus(E, res(E));
                        P=sqrt(E*E-emass*emass);
                        TLorentzVector veca(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);
                        if(P_trans(veca)>=lowPT){
                            histoWp->Fill(veca.PseudoRapidity());}
                        /*//second electron
                        x=acos(2*rng.Uniform(0,1)-1);
                        veca=TLorentzVector (0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);;
                        histoWp->Fill(veca.PseudoPseudoRapidity());*/
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+mumass*mumass)/(2*vecW.E());
                        E=rng.Gaus(E, res(E));
                        P=sqrt(E*E-mumass*mumass);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);
                        if(P_trans(veca)>=lowPT){
                            histoWpb->Fill(veca.PseudoRapidity());}
                    }
                }else if((part1._type=="AUp" && part2._type=="Down")||(part1._type=="Down" && part2._type=="AUp")){
                    if (M>=80.385){
                        //make W-
                        ++wm;
                        TVector3 boost=vecW.BoostVector();

                        //boost system to com
                        vecW.Boost( -boost );

                        //decay into stuff
                        //electron
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+emass*emass)/(2*vecW.E());
                        E=rng.Gaus(E, res(E));
                        P=sqrt(E*E-emass*emass);
                        TLorentzVector veca(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);;
                        if(P_trans(veca)>=lowPT){
                            histoWm->Fill(veca.PseudoRapidity());}
                        /*//second
                        x=acos(2*rng.Uniform(0,1)-1);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );
                        veca.Boost(boost);;
                        histoWm->Fill(veca.PseudoPseudoRapidity());*/
                        //muon
                        x=acos(2*rng.Uniform(0,1)-1);
                        E=(vecW.E()*vecW.E()+mumass*mumass)/(2*vecW.E());
                        E=rng.Gaus(E, res(E));
                        P=sqrt(E*E-mumass*mumass);
                        veca=TLorentzVector(0, P*sin(x),P*cos(x),E );

                        veca.Boost(boost);
                        if(P_trans(veca)>=lowPT){
                            histoWmb->Fill(veca.PseudoRapidity());}
                    }
                }
            }
        }
    }



    //if (part1._type ==& part2)
    /*TH1D *histo = new TH1D("Up_test","histogram of uniform distribution",50,0,1);
    int prod=0;
    while (prod<10000){
        Parton part=make_parton(fmax);
        if (part._type=="Down"){
            histo->Fill(part._x);
            prod++;
            if(prod%100==0){
                cout<<prod<<endl;
            }
        }

    }
    c1->cd();//make sure you plot on the right canvas
    histo->Draw();
    histo->SetMinimum(0);
    c1->Update();
    */
    cwm->cd();//make sure you plot on the right canvas
    histoWm->Draw();
    histoWm->SetMinimum(0);
    cwm->Update();

    cwmb->cd();//make sure you plot on the right canvas
    histoWmb->Draw();
    histoWmb->SetMinimum(0);
    cwmb->Update();

    cz->cd();//make sure you plot on the right canvas
    histoZ->Draw();
    histoZ->SetMinimum(0);
    cz->Update();

    czb->cd();//make sure you plot on the right canvas
    histoZb->Draw();
    histoZb->SetMinimum(0);
    czb->Update();

    cwpb->cd();//make sure you plot on the right canvas
    histoWpb->Draw();
    histoWpb->SetMinimum(0);
    cwpb->Update();

    cwp->cd();//make sure you plot on the right canvas
    histoWp->Draw();
    histoWp->SetMinimum(0);
    cwp->Update();
    cout<<"Z="<<z<<" Wmin="<<wm<<" Wplus="<<wp<<endl;
    TH1D *ass = new TH1D("ass","electron",30,-limit,limit);
    TH1D *ass_sub = new TH1D("ass_sub","electron",30,-limit,limit);
    TH1D *ass_sum = new TH1D("ass_sum","electron",30,-limit,limit);
    ass_sub->Add(histoWm, histoWp,-1,1);
    ass_sum->Add(histoWm, histoWp,1,1);
    ass->Divide(ass_sub, ass_sum);
    c1->cd();//make sure you plot on the right canvas
    ass->Draw();
    ass->SetMinimum(0);
    c1->Update();
    TH1D *assb = new TH1D("assb","muon",30,-limit,limit);
    TH1D *assb_sub = new TH1D("assb_sub","muon",30,-limit,limit);
    TH1D *assb_sum = new TH1D("assb_sum","muon",30,-limit,limit);
    assb_sub->Add(histoWmb, histoWpb,-1,1);
    assb_sum->Add(histoWmb, histoWpb,1,1);
    assb->Divide(assb_sub, assb_sum);
    c2->cd();//make sure you plot on the right canvas
    assb->Draw();
    assb->SetMinimum(0);
    c2->Update();
    gBenchmark->Show("fillrandom");
}
//cross section
double sigma(){
    //to be filled in
    return 1.0;
}

inline double res(double E){
    return (0.028/sqrt(E)+0.12/E+0.3)*E;
}

inline double P_trans(TLorentzVector v){
    return sqrt( v.Px()*v.Px() +v.Py()*v.Py() );
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
    /*if(f>=fAUp+fUp+fADown+fDown+fGluon){
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
    }*/

    if(f>=fAUp+fUp+fADown+fDown+fGluon){
        //failed rejection
        return make_parton(fmax);
    }else if(f>=fAUp+fUp+fADown+fDown){
        //ignore gluons
        return make_parton(fmax);
    }else if(f>=fAUp+fUp+fADown){
        if(previous_anti==0){
            return make_parton(fmax);
        }else{
            previous_anti=0;
            return Parton(x,"Down");
        }
    }else if(f>=fAUp+fADown){
        if(previous_anti==0){
            return make_parton(fmax);
        }else{
            previous_anti=0;
            return Parton(x,"Up");
        }
    }else if(f>=fAUp){
        if(previous_anti==1){
            return make_parton(fmax);
        }else{
            previous_anti=1;
            return Parton(x,"ADown");
        }
    }else{
        if(previous_anti==1){
            return make_parton(fmax);
        }else{
            previous_anti=1;
            return Parton(x,"AUp");
        }
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
