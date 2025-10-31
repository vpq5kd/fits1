#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include <TMath.h>


//chatGPT aided in creating this original nll function.
double ComputeNLL(TH1F* h, TF1* f) {
    double nll = 0.0;
    int nbins = h->GetNbinsX();

    for (int i = 1; i <= nbins; ++i) {
        double yi = h->GetBinContent(i);
        double xi = h->GetBinCenter(i);

        double lambda = f->Eval(xi) * h->GetBinWidth(i);

        if (lambda > 0) {
            nll += yi * log(lambda) - lambda - TMath::LnGamma(yi + 1.0);
        }
    }

    return -2.0 * nll;  // convention: multiply by -2
}

void exercise4(){
	TFile *f = new TFile("histo25.root");
	TH1F *data = (TH1F*)f->Get("randomHist1");
	data->Fit("gaus","L");
	TF1 *fitfunc = data->GetFunction("gaus");
	
	double bestMean = fitfunc->GetParameter(1);
	double bestErr = fitfunc->GetParError(1);

	int nPoints = 1000;

	TGraph *graph = new TGraph();
	for (int i = 0; i < nPoints; i++){
		double mean_test = bestMean - 4*bestErr + (8*bestErr * i)/(nPoints -1); // variable def from chatGPT
		fitfunc->FixParameter(1, mean_test);
		data->Fit(fitfunc, "L");
		double nll_val = fitfunc->GetChisquare(); // returns -2ln(l) according to chatGPT
		graph->SetPoint(i, mean_test, nll_val);
	}
	TCanvas *c1 = new TCanvas("c", "NLL and #chi^{2} Scans", 800, 600);
	graph->SetTitle("Test Mean vs NLL Value;Mean Value; -2lnL");
	graph->Draw("AP");
}
