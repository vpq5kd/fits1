#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"


// fit1.C
// entries is the number of random samples filled into the histogram
void fit1(int entries=1000, bool save=false) {
   //Simple histogram fitting examples
  gROOT->Reset();  // useful to reset ROOT to a cleaner state

  TFile *tf=0;
  if (save) tf=new TFile("histo.root","recreate");
  TH1F *meanHist = new TH1F("meanHist","Mean",50,48,52);
  TGraph *chi2Prob = new TGraph();
  TH1F *meanErrorHist = new TH1F("meanErrorHist","Mean Error",50,0,.6);
  TH1F *chi2hist = new TH1F("chi2hist","reduced #chi^{2}",50,0,3);
  TH1F *randomHist1 = new TH1F("randomHist1", "Random Histogram;x;frequency", 100, 0, 100);
  TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock
  int ntrials = 100;
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
	 double chi2 = fitfunc->GetChisquare();
	 double mean = fitfunc->GetParameter(1);
	 double meanError = fitfunc->GetParError(1);
	 double prob = TMath::Prob(fitfunc->GetChisquare(),fitfunc->GetNDF());
	 chi2hist->Fill(chi2/fitfunc->GetNDF());
	 meanHist->Fill(mean);
	 meanErrorHist->Fill(meanError);
	 chi2Prob->SetPoint(point,chi2/fitfunc->GetNDF(),prob);
	 point++;
 }
 TCanvas *c1 = new TCanvas("c1","My Canvas", 800,600);
 c1->Divide(2,2);

 c1->cd(1);
 chi2hist->Draw();

 c1->cd(2);
 meanHist->Draw();

 c1->cd(3);
 chi2Prob->SetTitle("#chi^{2} Probability vs Reduced #chi^{2};Reduced #chi^{2};P(#chi^{2})");
 chi2Prob->Draw("AP");

 c1->cd(4);
 meanErrorHist->Draw();
	
  if (save) {
    tf->Write();
    tf->Close();
  }
  cout << "Use .q to exit root" << endl;
}
