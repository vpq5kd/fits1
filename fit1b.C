#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"

// chatGPT was used for root syntax/functions.
// fit1.C
// entries is the number of random samples filled into the histogram
void fit1(int entries=1000, bool save=false) {
   //Simple histogram fitting examples
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TFile *tf=0;
  if (save) tf=new TFile("histo.root","recreate");
  TH1F *meanHistChi2 = new TH1F("meanHistChi2","Mean #chi^{2} Fit",50,48,52);
  TH1F *meanHistNLL = new TH1F("meanHistNLL","Mean NLL fit",50,48,52);
  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock
  int ntrials = 1000;
  int point = 0; 
  for (int i = 0; i< ntrials; i++){
	 randomHist1->Reset();
 	 for (int i=0 ; i<entries ; i++){
   		 randomHist1->Fill(generator->Gaus(50,10)); // params: mean, sigma
	 }

 	 gStyle->SetOptFit(1111); // show reduced chi2, probability, and params
 	 randomHist1->Fit("gaus");  
 	 randomHist1->DrawCopy("e");  // "e" shows bin errors
 	 TF1 *fitfunc = randomHist1->GetFunction("gaus");
	 double meanChi2 = fitfunc->GetParameter(1);
	 meanHistChi2->Fill(meanChi2);

	 randomHist1->Fit("gaus","L");
 	 TF1 *fitfunc1 = randomHist1->GetFunction("gaus");
	 double meanNLL = fitfunc1->GetParameter(1);
	 meanHistNLL->Fill(meanNLL);
 }
 TCanvas *c1 = new TCanvas("c1","My Canvas", 800,600);
 c1->Divide(2,1);

 c1->cd(1);
 meanHistChi2->Draw();

 c1->cd(2);
 meanHistNLL->Draw();

  if (save) {
    tf->Write();
    tf->Close();
  }
  cout << "Use .q to exit root" << endl;
}
