{

  // example code on how to load histograms from a file and use them

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

  ///////////////////////////////////////////
  // the PDFs were collected from the PDF4LHC program. Specific description of this PDF:
    //  >>>>>> PDF description: <<<<<<
    //  CTEQ6.6                                                         
    //  Reference:                                                      
    //  P.M.Nadolsky, H.-L.Lai, Q.-H.Cao, J. Huston,                    
    //  J.Pumplin, D.Stump, W.-K.Tung and C.-P.Yuan.                    
    //  arXiv:0802.0007                                                 
    //  This set has 44 member PDFs.                                    
    //      mem=0 --> central value                                     
    //      mem=1-44 --> eigenvector sets                               
    //          use these for standard uncertainty calculations         
    //          all have alpha_s(mz) = 0.118                            

    
    //  as this is an example:  take an arbitrary value for x and q-squared and retrieve the xfx probabilities for gluons, up, anti-up, down and anti-down:
    
    // important: the x coordinate is in a different mapping and should be mapped as : xreal = pow(10,-0.1*xhisto). The y coordinate is q^2.

    // take some simple values like x = 90./8000.; (Assumption: Z boson at rest)
    // and q-squared (the scale of the process) = 90.;   

    // taking one histogram to get the coordinate system, as long as the axes are the same for all histograms this works (of course if they are different this should be done for every histogram)
    Float_t xhisto = -10*TMath::Log10(90./8000.);
    std::cout << "mapping " << 90./8000. << " to " << xhisto ;
    Int_t xbin = histoUp->GetXaxis()->FindBin(xhisto);
    std::cout << " which is in bin: " << xbin << std::endl;
    Int_t ybin = histoUp->GetYaxis()->FindBin(90.);

    // now retrieve the values from the various histograms.
    Double_t xfxUp = histoUp->GetBinContent(xbin,ybin);
    Double_t xfxAUp = histoAUp->GetBinContent(xbin,ybin);
    std::cout << "the xfx probability for Up and anti-Up quarks at these values is: " << xfxUp << " " << xfxAUp << std::endl;

    Double_t xfxDown = histoDown->GetBinContent(xbin,ybin);
    Double_t xfxADown = histoADown->GetBinContent(xbin,ybin);
    std::cout << "the xfx probability for Down and anti-Down quarks at these values is: " << xfxDown << " " << xfxADown << std::endl;



}
