#include<iomanip>
#include<sstream>
#include "/home/ilc/tianjp/bin/root_macros/style.C"

#define __electron__

void mH_IDR() {

  setILDStyle();
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);

#ifndef __electron__
  TFile *fL = new TFile("eps/HiggsMassNew_e2e2h_bb_eLpR.Large.root");
  TFile *fS = new TFile("eps/HiggsMassNew_e2e2h_bb_eLpR.Small.root");
#else
  TFile *fL = new TFile("eps/HiggsMassNew_e1e1h_bb_eLpR.Large.root");
  TFile *fS = new TFile("eps/HiggsMassNew_e1e1h_bb_eLpR.Small.root");
#endif  

  TH1D *hSL = (TH1D*)fL->Get("hS2");
  TH1D *hBL = (TH1D*)fL->Get("hB01");
  TH1D *hSS = (TH1D*)fS->Get("hS2");
  TH1D *hBS = (TH1D*)fS->Get("hB01");

  hSL->Add(hBL);
  hSS->Add(hBS);

#ifdef __electron__
  hSL->Rebin(2);
  hSS->Rebin(2);
#endif  
  
  hSS->SetLineStyle(kDashed);
  hSL->SetLineColor(4);  
  hSS->SetLineColor(2);  

  hSL->GetXaxis()->SetTitle("Higgs Mass [GeV]");
#ifndef __electron__
  hSL->GetYaxis()->SetTitle("Entries / 1.5 GeV");
#else
  hSL->GetYaxis()->SetTitle("Entries / 3 GeV");
#endif  

  Double_t ymax = hSL->GetMaximum() > hSS->GetMaximum() ? hSL->GetMaximum() : hSS->GetMaximum();
  hSL->SetMaximum(ymax*1.2);

  Double_t xmin=60, xmax=180;
  
#ifndef __electron__
  TString cvs="mH_e2e2h_bb_eLpR_IDR";
#else
  TString cvs="mH_e1e1h_bb_eLpR_IDR";
#endif  
  
  TCanvas *c = new TCanvas(cvs,cvs,600,600);
  hSL->Draw();    
  hSS->Draw("same");    
  
#if 1
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(hSL,"IDR-L","l");
  leg->AddEntry(hSS,"IDR-S","l");  
  leg->SetFillStyle(0);
  leg->Draw();

  TLatex *ltx = new TLatex(xmin+(xmax-xmin)/20,ymax*1.1,"ILD (preliminary)");
  ltx->SetTextFont(62);
  ltx->Draw();

#ifndef __electron__
  TLatex *ltx2 = new TLatex(xmax-(xmax-xmin)/4,ymax*0.6,"#mu#muH");
#else  
  TLatex *ltx2 = new TLatex(xmax-(xmax-xmin)/4,ymax*0.6,"eeH");  
#endif  
  ltx2->Draw();
  
  stringstream savename,savename2,savename3;
  savename << "eps/" << cvs << ".pdf" << ends;
  savename2 << "eps/" << cvs << ".eps" << ends;
  savename3 << "eps/" << cvs << ".C" << ends;
  c->SaveAs(savename.str().data());
  c->SaveAs(savename2.str().data());
  c->SaveAs(savename3.str().data());
#endif  
  
  
}
