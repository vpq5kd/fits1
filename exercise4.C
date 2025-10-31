#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom2.h"
#include "TROOT.h"
#include <TMath.h>
//I attempted to solve this as best I could but I'm not entirely sure that I did it correctly. I relied heavily on chatGPT to help me understand the assignment,and I did a lot of the coding in direct conjunction with chatGPT's help. I labeled specific parts that were chatGPTs idea. I am able to explain all of my code as well as my rational for how I interpreted the assignment and thus implemented it in my logic. 

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

	int nPoints = 50;

	TGraph *graph = new TGraph();
	for (int i = 0; i < nPoints; i++){
		double mean_test = bestMean - 4*bestErr + (8*bestErr * i)/(nPoints -1); // variable def from chatGPT
		fitfunc->FixParameter(1, mean_test);
		data->Fit(fitfunc, "L");
		double nll_val = ComputeNLL(data,fitfunc);
		graph->SetPoint(i, mean_test, nll_val);
		fitfunc->ReleaseParameter(1);
	}

	double nll_min = 1e99;
	for (int i = 0; i <graph->GetN(); i++){
		double x, y;
		graph->GetPoint(i,x,y);
		if (y < nll_min) nll_min = y;
	}
	for (int i = 0; i <graph->GetN(); i++){
		double x, y;
		graph->GetPoint(i, x, y);
		graph->SetPoint(i,x,y-nll_min);
	}

	TFile *f1 = new TFile("histo1k.root");
	TH1F *data1 = (TH1F*)f1->Get("randomHist1");
	data1->Fit("gaus");
	TF1 *fitfunc1 = data->GetFunction("gaus");
	
	double bestMean1 = fitfunc1->GetParameter(1);
	double bestErr1 = fitfunc1->GetParError(1);

	int nPoints1 = 50;

	TGraph *graphchi2 = new TGraph();
	for (int i = 0; i < nPoints1; i++){
		double mean_test = bestMean - 4*bestErr + (8*bestErr * i)/(nPoints -1); // variable def from chatGPT
		fitfunc1->FixParameter(1, mean_test);
		data->Fit(fitfunc1, "L");
		double nll_val = fitfunc1->GetChisquare();
		graphchi2->SetPoint(i, mean_test, nll_val);
		fitfunc1->ReleaseParameter(1);
	}
	nll_min = 1e99;
	for (int i = 0; i <graphchi2->GetN(); i++){
		double x, y;
		graphchi2->GetPoint(i,x,y);
		if (y < nll_min) nll_min = y;
	}
	for (int i = 0; i <graph->GetN(); i++){
		double x, y;
		graphchi2->GetPoint(i, x, y);
		graphchi2->SetPoint(i,x,y-nll_min);
	}

	TCanvas *c1 = new TCanvas("c", "NLL and #chi^{2} Scans", 800, 600);
	c1->Divide(1,2);

	c1->cd(1);
	graph->SetTitle("Test Mean vs NLL Value;Mean Value; -2lnL");
	graph->GetXaxis()->SetTitle("Mean Value");
	graph->GetYaxis()->SetTitle("-2lnL");
	graph->Draw("AL");
	
	c1->cd(2);
	graphchi2->SetTitle("Test Mean vs #chi^{2};Mean Value; #chi^{2}");
	graphchi2->GetXaxis()->SetTitle("Mean Value");
	graphchi2->GetYaxis()->SetTitle("#chi^{2}");
	graphchi2->Draw("AL");
	
	c1->SaveAs("result4.pdf");
}
