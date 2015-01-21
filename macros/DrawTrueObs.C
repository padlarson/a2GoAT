#include "TH1.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"

using namespace std;

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}


void DrawTrueObs()
{

	DrawEvTh();
	DrawPhysics();
}

	void DrawEvTh()
	{
	  	// Declaration of path to root file and scatterplots
	TString path = "/home/adlarson/MC/MyPhysics/";
	TFile *f = new TFile(path + "Phys_etap_pi0pi0eta2g.root"); 

	TH2F *protontrue 		= (TH2F*) f->Get("ThvE_p");
// TH2F *protontrue = (TH2F*) f->Get("ThvE_p");
	TH2F *gamma_etatrue 	= (TH2F*) f->Get("ThvE_eta_g");
	TH2F *gamma_piontrue 	= (TH2F*) f->Get("ThvE_pi0_g");

	set_plot_style();

	// Phase space

	TCanvas *c1 = new TCanvas("c1", "p eta' true observables", 200, 10, 1000, 700);
	c1->SetFillColor(16);
	pad1 = new TPad("pad1", "proton", 0.05, 0.5, 0.45, 0.95, 19);
	pad2 = new TPad("pad2", "#theta_{proton,CM} vs #theta_{#eta^{'},CM}", 0.5, 0.5, 0.95, 0.95, 18);
	pad3 = new TPad("pad3", "Energy vs #theta #gamma_{#eta}", 0.05, 0.05, 0.45, 0.45, 18);
	pad4 = new TPad("pad4", "Energy vs #theta #gamma_{#pi^{0}}", 0.5, 0.05, 0.95, 0.45, 18);

	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();

	pad1->cd();
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->GetFrame()->SetFillColor(9);
	pad1->GetFrame()->SetBorderMode(-1);
	pad1->GetFrame()->SetBorderSize(5);


	protontrue->SetLineColor(4); protontrue->SetLineWidth(3);
	protontrue->GetXaxis()->SetTitle("Proton Energy (GeV)"); protontrue->GetXaxis()->SetTitleSize(0.055); protontrue->GetXaxis()->SetTitleOffset(0.65);
	protontrue->GetYaxis()->SetTitle("#theta_{lab} (^{o})"); protontrue->GetYaxis()->SetTitleSize(0.055); protontrue->GetYaxis()->SetTitleOffset(0.8);
	protontrue->Draw("colz");

	pad3->cd();
	pad3->SetGridx();
	pad3->SetGridy();
	pad3->GetFrame()->SetFillColor(9);
	pad3->GetFrame()->SetBorderMode(-1);
	pad3->GetFrame()->SetBorderSize(5);
	gamma_etatrue->SetLineColor(4); gamma_etatrue->SetLineWidth(3);
	gamma_etatrue->GetXaxis()->SetTitle("Energy (GeV)"); gamma_etatrue->GetXaxis()->SetTitleSize(0.055); gamma_etatrue->GetXaxis()->SetTitleOffset(0.65);
	gamma_etatrue->GetYaxis()->SetTitle("#theta_{lab} (^{o})"); gamma_etatrue->GetYaxis()->SetTitleSize(0.055); gamma_etatrue->GetYaxis()->SetTitleOffset(0.8);
	gamma_etatrue->Draw("colz");

	pad4->cd();
	pad4->SetGridx();
	pad4->SetGridy();
	pad4->GetFrame()->SetFillColor(9);
	pad4->GetFrame()->SetBorderMode(-1);
	pad4->GetFrame()->SetBorderSize(5);
	gamma_piontrue->SetLineColor(4); gamma_piontrue->SetLineWidth(3);
	gamma_piontrue->GetXaxis()->SetTitle("Energy (GeV)"); gamma_piontrue->GetXaxis()->SetTitleSize(0.055); gamma_piontrue->GetXaxis()->SetTitleOffset(0.65);
	gamma_piontrue->GetYaxis()->SetTitle("#theta_{lab} (^{o})"); gamma_piontrue->GetYaxis()->SetTitleSize(0.055); gamma_piontrue->GetYaxis()->SetTitleOffset(0.8);
	gamma_piontrue->Draw("colz");

	c1->Print("/home/adlarson/Dropbox/postdoc Mainz/eta prime/Report/eta prime.tex/figure/MC/etapr6g_EvsTh.eps");

}

 
	void DrawPhysics()
	{
	// True physics 
	TCanvas *c2 = new TCanvas("c2", "Dalitz Plot", 200, 10, 1000, 700);
	c2->SetFillColor(18);
	pad1 = new TPad("pad1", "XvsY", 0.05, 0.5, 0.45, 0.95, 21);
	pad2 = new TPad("pad2", "#theta_{proton,CM} vs #theta_{#eta^{'},CM}", 0.5, 0.5, 0.95, 0.95, 21);
	pad3 = new TPad("pad3", "Energy vs #theta #gamma_{#eta}", 0.05, 0.05, 0.45, 0.45, 21);
	}



