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

void exercise3(){
	TFile *f = new TFile("histo25.root");
	TH1F *data = (TH1F*)f->Get("randomHist1");
	data->Fit("gaus","L");
	TF1 *fitfunc = data->GetFunction("gaus");
	
	double mean = fitfunc->GetParameter(0);
	double sigma = fitfunc->GetParameter(2);


	const int experiments = 1000;
  	TRandom2 *generator=new TRandom2(0);  // parameter == seed, 0->use clock
	double nll_for_data=ComputeNLL(data, fitfunc);

	TH1F *nllHist = new TH1F("NLL Distribution", "NLL Distribution from each Experiment", 100, nll_for_data *0.8, nll_for_data *1.2);
	//logic for this for loop developed in conjunction with ChatGPT
	for (int i = 0; i < experiments; i++){
		TH1F *experiment = (TH1F*)data->Clone(Form("experiment_%d",i));
		experiment->Reset();
		for (int i =1; i<=experiment->GetNbinsX();i++){
			double xi = experiment->GetBinCenter(i);
			double lambda = fitfunc->Eval(xi) * experiment->GetBinWidth(i);
			int yi = generator->Poisson(lambda);
			experiment->SetBinContent(i,yi);
		}

		double nll_experiment = ComputeNLL(experiment, fitfunc);
		nllHist->Fill(nll_experiment);
		
		delete experiment;

	}
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	nllHist->Draw();
	TLine *line = new TLine(nll_for_data, 0, nll_for_data, nllHist->GetMaximum());
	line->SetLineColor(kRed);
	line->Draw("same");

	//chatGPT taught me how to calculate a p value for this because I don't know a lot of stats lol
	int binData = nllHist->FindBin(nll_for_data);
	int nAbove = nllHist->Integral(binData, nllHist->GetNbinsX());
	double pVal = double(nAbove)/experiments;

	cout <<"p-value = " << pVal <<endl;
}
