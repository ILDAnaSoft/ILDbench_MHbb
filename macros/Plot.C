#include<iomanip>
#include<sstream>
#include "/home/ilc/tianjp/bin/root_macros/style.C"

#define __fit__

Int_t _nbin_fit = 5;

void PlotMH() {
  setILDStyle();

  //  TFile *f1 = new TFile("../root_merge/e2e2h_bb.eL.pR.root");
  TFile *f1 = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.eL.pR.mILD_l5_o1_v02.root");
  TNtupleD *n1 = f1->Get("hAnl");

  Int_t nbin=400;
  Double_t xmin=50,xmax=450;
  TH1D *h1 = new TH1D("h1","",nbin,xmin,xmax);
  TH1D *h2 = new TH1D("h2","",nbin,xmin,xmax);
  TH1D *h3 = new TH1D("h3","",nbin,xmin,xmax);

  n1->Draw("mhnew>>h1","nhbb==1");
  n1->Draw("mrecoil>>h2","nhbb==1");
  n1->Draw("mh>>h3","nhbb==1");

#if 0
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
#endif  

  h1->SetStats(0);
  h1->SetXTitle("Higgs Mass [GeV]");
  h1->SetYTitle("Entries / 1 GeV");
  h1->GetYaxis()->CenterTitle();

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(3);  

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);  

  //  Double_t ymax = h1->GetMaximum();
  //  h1->GetYaxis()->SetRangeUser(0.,ymax*1.3);

  TH1D *hh[3];
  hh[0] = h1;
  hh[1] = h2;
  hh[2] = h3;

  for (Int_t i=0;i<3;i++) {
    Double_t mx = hh[i]->GetMaximum();
    Int_t binmax1 = hh[i]->FindFirstBinAbove(mx-0.01);
    Int_t binmax2 = hh[i]->FindLastBinAbove(mx-0.01);
    Double_t vmax1 = hh[i]->GetBinCenter(binmax1);
    Double_t vmax2 = hh[i]->GetBinCenter(binmax2);
    Int_t binhmax1 = hh[i]->FindFirstBinAbove(mx/2-0.01);
    Int_t binhmax2 = hh[i]->FindLastBinAbove(mx/2-0.01);
    Double_t vhmax1 = hh[i]->GetBinCenter(binhmax1);
    Double_t vhmax2 = hh[i]->GetBinCenter(binhmax2);
    Double_t ntot = hh[i]->Integral();
    Double_t nfwhm = hh[i]->Integral(binhmax1,binhmax2);
    cerr << "maximum h" << i+1 << ": " << mx
	 << "; NBmax = " << binmax1 << " " << binmax2
	 << "; Vmax = " << vmax1 << " " << vmax2
	 << "; NBHmax = " << binhmax1 << " " << binhmax2
	 << "; VHmax = " << vhmax1 << " " << vhmax2
	 << "; wH = " << vhmax2 - vhmax1
	 << "; frac = " << nfwhm/ntot * 100
	 << endl;
  }

  
  Double_t ymax = h1->GetMaximum();
  h1->SetMaximum(ymax*1.2);
  
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");  

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry("h1","New Method");
  leg->AddEntry("h2","Recoil Mass");
  leg->AddEntry("h3","Direct M(JJ)");  
  leg->SetFillStyle(0);
  leg->Draw();

  TLatex *ltx = new TLatex(200.,0.5*ymax,"e^{+}+e^{-}#rightarrow#mu#muH, H#rightarrowb#bar{b} @ 500 GeV");
  ltx->Draw();

  TLatex *ltx2 = new TLatex(xmin+(xmax-xmin)/20,ymax*1.1,"ILD (preliminary)");
  ltx2->SetTextFont(62);
  ltx2->Draw();

  TLatex *ltx3 = new TLatex(300.,0.4*ymax,"IDR-L");
  ltx3->Draw();

  
  gPad->SaveAs("eps/HiggsMass_compare_mumuH_500_IDR-L.eps");
  gPad->SaveAs("eps/HiggsMass_compare_mumuH_500_IDR-L.pdf");
  gPad->SaveAs("eps/HiggsMass_compare_mumuH_500_IDR-L.C");

}

void PlotMH_SB() {
  setILDStyle();

  TFile *f1 = new TFile("../root_merge/e2e2h_bb.eL.pR.root");
  TFile *f2 = new TFile("../root_merge/4f_zz_sl.eL.pR.root");
  TNtupleD *n1 = f1->Get("hAnl");
  TNtupleD *n2 = f2->Get("hAnl");

  Int_t nbin=400;
  Double_t xmin=50,xmax=450;
  TH1D *h1 = new TH1D("h1","",nbin,xmin,xmax);
  TH1D *h2 = new TH1D("h2","",nbin,xmin,xmax);
  TH1D *h3 = new TH1D("h3","",nbin,xmin,xmax);
  TH1D *h4 = new TH1D("h4","",nbin,xmin,xmax);

  n1->Draw("mhnew>>h1");
  n2->Draw("mhnew>>h2");
  n1->Draw("mrecoil>>h3");
  n2->Draw("mrecoil>>h4");

#if 1
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
  h4->Scale(1./h4->Integral());
#endif  

  h1->SetStats(0);
  h1->SetXTitle("Higgs Mass [GeV]");
  h1->SetYTitle("Normalised to 1");
  h1->GetYaxis()->CenterTitle();

  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(2);  
  h4->SetLineColor(4);

  h3->SetLineStyle(kDashed);
  h4->SetLineStyle(kDashed);

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);
  h4->SetMarkerSize(0);    

  //  Double_t ymax = h1->GetMaximum();
  //  h1->GetYaxis()->SetRangeUser(0.,ymax*1.3);

  h1->Draw();
  h2->Draw("same");
  //  h3->Draw("same");  
  //  h4->Draw("same");  

  TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->AddEntry("h1","#mu#muH, H#rightarrowb#bar{b}");
  leg->AddEntry("h2","ZZ#rightarrowsemi-leptonic");
  leg->SetFillStyle(0);
  leg->Draw();

  Double_t ymax = h1->GetMaximum();
  TLatex *ltx = new TLatex(200.,0.5*ymax,"e^{+}+e^{-}#rightarrow#mu#muH, H#rightarrowb#bar{b} @ 500 GeV");
  //  ltx->Draw();

  gPad->SaveAs("eps/HiggsMass_compare_mumuH_zz_500.eps");

}

void PlotRes(Bool_t iLep = kFALSE) {
  setILDStyle();

  TFile *f1 = new TFile("../root_merge/e2e2h_bb.eL.pR.root");
  TNtupleD *n1 = f1->Get("hAnl");

  Int_t nbin=200;
  Double_t xmin=-1,xmax=1;
  TH1D *h1 = new TH1D("h1","",nbin,xmin,xmax);
  TH1D *h2 = new TH1D("h2","",nbin,xmin,xmax);
  TH1D *h3 = new TH1D("h3","",nbin,xmin,xmax);
  TH1D *h4 = new TH1D("h4","",nbin,xmin,xmax);
  TH1D *h5 = new TH1D("h5","",nbin,xmin,xmax);
  TH1D *h6 = new TH1D("h6","",nbin,xmin,xmax);
  TH1D *h7 = new TH1D("h7","",nbin,xmin,xmax);
  TH1D *h8 = new TH1D("h8","",nbin,xmin,xmax);

  n1->Draw("(cosj1-cosj1mc)/cosj1mc>>h1");
  n1->Draw("(phij1-phij1mc)/phij1mc>>h2");
  n1->Draw("(ej1-ej1mc)/ej1mc>>h3");
  n1->Draw("(cosj2-cosj2mc)/cosj2mc>>h4");
  n1->Draw("(phij2-phij2mc)/phij2mc>>h5");
  n1->Draw("(ej2-ej2mc)/ej2mc>>h6");
  if (iLep) {
    n1->Draw("(elep1-elep1mc)/elep1mc>>h7");
    n1->Draw("(elep2-elep2mc)/elep2mc>>h8");
  }

  h1->Add(h4);
  h2->Add(h5);
  h3->Add(h6);
  if (iLep) h7->Add(h8);
  
#if 1
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
  if (iLep)  h7->Scale(1./h7->Integral());
#endif  

  if (! iLep) {
  h2->SetStats(0);
  h2->SetXTitle("(rec.-truth)/truth");
  h2->SetYTitle("Normalised to 1");  
  h2->GetYaxis()->CenterTitle();
  }
  else {
  h7->SetStats(0);
  h7->SetXTitle("(rec.-truth)/truth");
  h7->SetYTitle("Normalised to 1");  
  h7->GetYaxis()->CenterTitle();
  }
  
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(3);  
  if (iLep) h7->SetLineColor(6);  

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);
  h3->SetMarkerSize(0);
  if (iLep) h7->SetMarkerSize(0);    

  //  Double_t ymax = h2->GetMaximum();
  //  h1->GetYaxis()->SetRangeUser(0.,ymax*1.2);

  if (iLep) {
    h7->Draw();
    h2->Draw("same");
  }
  else {
    h2->Draw();
  }
  h1->Draw("same");
  h3->Draw("same");  

  TLegend *leg = new TLegend(0.7,0.5,0.85,0.85);
  leg->AddEntry("h1","cos#theta_{j}");
  leg->AddEntry("h2","#phi_{j}");
  leg->AddEntry("h3","E_{j}");
  if (iLep) leg->AddEntry("h7","E_{#mu}");
  leg->SetFillStyle(0);
  leg->Draw();

  Double_t ymax = h1->GetMaximum();
  TLatex *ltx = new TLatex(200.,0.5*ymax,"e^{+}+e^{-}#rightarrow#mu#muH, H#rightarrowb#bar{b} @ 500 GeV");
  //  ltx->Draw();

  if (! iLep) {
  gPad->SaveAs("eps/Resolutions_compare_mumuH_500.eps");
  gPad->SetLogy();
  gPad->SaveAs("eps/Resolutions_compare_Logy_mumuH_500.eps");
  }
  else {
  gPad->SaveAs("eps/Resolutions_lep_compare_mumuH_500.eps");
  gPad->SetLogy();
  gPad->SaveAs("eps/Resolutions_lep_compare_Logy_mumuH_500.eps");
  }
}

void PlotPJ() {
  setILDStyle();

  TFile *f1 = new TFile("../root_merge/e2e2h_bb.eL.pR.root");
  TNtupleD *n1 = f1->Get("hAnl");

  Int_t nbin=200;
  Double_t xmin=-1,xmax=1;
  TH1D *h1 = new TH1D("h1","",nbin,xmin,xmax);
  TH1D *h2 = new TH1D("h2","",nbin,xmin,xmax);
  TH1D *h3 = new TH1D("h3","",nbin,xmin,xmax);
  TH1D *h4 = new TH1D("h4","",nbin,xmin,xmax);

  n1->Draw("(sqrt(pxnewj1**2+pynewj1**2+pznewj1**2)-sqrt(pxj1mc**2+pyj1mc**2+pzj1mc**2))/sqrt(pxj1mc**2+pyj1mc**2+pzj1mc**2)>>h1");
  n1->Draw("(sqrt(pxj1**2+pyj1**2+pzj1**2)-sqrt(pxj1mc**2+pyj1mc**2+pzj1mc**2))/sqrt(pxj1mc**2+pyj1mc**2+pzj1mc**2)>>h2");
  n1->Draw("(sqrt(pxnewj2**2+pynewj2**2+pznewj2**2)-sqrt(pxj2mc**2+pyj2mc**2+pzj2mc**2))/sqrt(pxj2mc**2+pyj2mc**2+pzj2mc**2)>>h3");
  n1->Draw("(sqrt(pxj2**2+pyj2**2+pzj2**2)-sqrt(pxj2mc**2+pyj2mc**2+pzj2mc**2))/sqrt(pxj2mc**2+pyj2mc**2+pzj2mc**2)>>h4");

  h1->Add(h3);
  h2->Add(h4);
#if 1
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
#endif  

  h1->SetStats(0);
  //  h1->SetXTitle("P_{j} [GeV]");
  h1->SetXTitle("(P_{j}-P_{truth})/P_{truth}");  
  h1->SetYTitle("Normalised to 1");  
  h1->GetYaxis()->CenterTitle();

  h1->SetLineColor(2);
  h2->SetLineColor(3);

  h1->SetMarkerSize(0);
  h2->SetMarkerSize(0);

  //  Double_t ymax = h1->GetMaximum();
  //  h1->GetYaxis()->SetRangeUser(0.,ymax*1.3);

  h1->Draw();
  h2->Draw("same");

  TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
  leg->AddEntry("h1","New Method");
  leg->AddEntry("h2","Direct Meas.");  
  leg->SetFillStyle(0);
  leg->Draw();

  Double_t ymax = h1->GetMaximum();
  TLatex *ltx = new TLatex(200.,0.5*ymax,"e^{+}+e^{-}#rightarrow#mu#muH, H#rightarrowb#bar{b} @ 500 GeV");
  //  ltx->Draw();

  gPad->SaveAs("eps/JetMomentum_compare_mumuH_500.eps");

}

void Comp(TString varname, TNtuple *n[3], Int_t nbin, Double_t xmin, Double_t xmax, TCut cut, TString cvs = "test", TString xtitle = "test") {
  setILDStyle();
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);

  TCanvas *c0 = new TCanvas("c0","",800,600);
  if (cvs == "test") cvs = varname;
  cerr << "cvs = " << cvs << endl;
  if (xtitle == "test") xtitle = varname;

  TH1D *h[3];
  Double_t ymax = -1;
  for (Int_t i=0;i<3;i++) {
    stringstream hname,drawvar;
    hname << "h" << i+1 << cvs << ends;
    h[i] = new TH1D(hname.str().data(),"",nbin,xmin,xmax);
    drawvar << varname << ">>" << hname.str().data() << ends;
    n[i]->Draw(drawvar.str().data(),cut);
    //    h[i]->SetLineColor(2+i*2);
    Double_t norm = h[i]->Integral();
    //    h[i]->Scale(1./norm);
    h[i]->SetXTitle(varname);
    //      h[i]->SetXTitle("z0 / #deltaz0");      
    if (h[i]->GetMaximum() > ymax) ymax = h[i]->GetMaximum();
  }

#if 1 // full width at the half maximum
  TH1D *hh[2];
  hh[0] = h[1];
  hh[1] = h[2];

  for (Int_t i=0;i<2;i++) {
    Double_t mx = hh[i]->GetMaximum();
    Int_t binmax1 = hh[i]->FindFirstBinAbove(mx-1.e-5);
    Int_t binmax2 = hh[i]->FindLastBinAbove(mx-1.e-5);
    Double_t vmax1 = hh[i]->GetBinCenter(binmax1);
    Double_t vmax2 = hh[i]->GetBinCenter(binmax2);
    Int_t binhmax1 = hh[i]->FindFirstBinAbove(mx/2);
    Int_t binhmax2 = hh[i]->FindLastBinAbove(mx/2);
    Double_t vhmax1 = hh[i]->GetBinCenter(binhmax1);
    Double_t vhmax2 = hh[i]->GetBinCenter(binhmax2);
    Double_t ntot = hh[i]->Integral();
    Double_t nfwhm = hh[i]->Integral(binhmax1,binhmax2);
    cerr << "maximum h" << i+1 << ": " << mx
	 << "; NBmax = " << binmax1 << " " << binmax2
	 << "; Vmax = " << vmax1 << " " << vmax2
	 << "; NBHmax = " << binhmax1 << " " << binhmax2
	 << "; VHmax = " << vhmax1 << " " << vhmax2
	 << "; wH = " << vhmax2 - vhmax1
	 << "; frac = " << nfwhm/ntot * 100
	 << endl;
#ifdef __fit__
    Int_t nbinfit = _nbin_fit;
    //    hh[i]->Fit("gaus","V","",vmax1-(xmax-xmin)/nbin*nbinfit,vmax2+(xmax-xmin)/nbin*nbinfit);
    hh[i]->Fit("gaus","","",vmax1-(xmax-xmin)/nbin*nbinfit,vmax2+(xmax-xmin)/nbin*nbinfit);    
    hh[i]->GetFunction("gaus")->SetLineColor(i==0? 4:2);
#endif  
  }
#endif


  h[0]->SetLineColor(3);
  h[1]->SetLineColor(4);
  h[2]->SetLineColor(2);
  h[2]->SetLineStyle(kDashed);


  h[1]->SetMaximum(ymax*1.2);
  h[1]->SetXTitle(xtitle);
  h[1]->SetYTitle("Normalized to 1");
  h[1]->GetYaxis()->CenterTitle();
  h[2]->SetMaximum(ymax*1.2);
  h[2]->SetXTitle(xtitle);
  h[2]->SetYTitle("Normalized to 1");
  //  h[2]->SetYTitle("Events");  
  h[2]->GetYaxis()->CenterTitle();
  
  //  stringstream cname;
  //  cname << "c" << varname << ends;
  //  TCanvas *c = new TCanvas(cname.str().data(),"",800,600);
  //  TCanvas *c = new TCanvas(cname.str().data(),cvs,800,600);
  //  TCanvas *c = new TCanvas(cvs,cvs,800,600);
  TCanvas *c = new TCanvas(cvs,cvs,600,600);
  //  TCanvas *c = new TCanvas(cvs,cvs,1024,768);        
  //  gPad->SetLogy();
  h[2]->Draw();    
  h[1]->Draw("same");
  //  h[0]->Draw("same");
  
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h[1],"IDR-L","l");
  leg->AddEntry(h[2],"IDR-S","l");  
  leg->SetFillStyle(0);
  leg->Draw();

  TLatex *ltx = new TLatex(xmin+(xmax-xmin)/20,ymax*1.1,"ILD (preliminary)");
  ltx->SetTextFont(62);
  ltx->Draw();
  
  stringstream savename,savename2,savename3;
  savename << "eps/" << cvs << ".pdf" << ends;
  savename2 << "eps/" << cvs << ".eps" << ends;
  savename3 << "eps/" << cvs << ".C" << ends;
  c->SaveAs(savename.str().data());
  c->SaveAs(savename2.str().data());
  c->SaveAs(savename3.str().data());
}

void CompTest(Int_t idet, Int_t ndraw, std::vector<TString> varnames, TNtuple *n[3], Int_t nbin, Double_t xmin, Double_t xmax, TCut cut, TString cvs, std::vector<TString> title, TString xtitle ) {
//void CompTest(Int_t idet, Int_t ndraw, TNtuple *n[3], Int_t nbin, Double_t xmin, Double_t xmax, TCut cut, TString cvs, TString xtitle ) {
  setILDStyle();
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);

  TCanvas *c0 = new TCanvas("c0","",800,600);
  return ;
}
void Comp(Int_t idet, Int_t ndraw, std::vector<TString> varnames, TNtuple *n[3], Int_t nbin, Double_t xmin, Double_t xmax, TCut cut, TString cvs, std::vector<TString> title, TString xtitle ) {
  setILDStyle();
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);

  TCanvas *c0 = new TCanvas("c0","",800,600);

  const Int_t ncomp = ndraw;
  TH1D *h[ncomp];
  Double_t ymax = -1;
  for (Int_t i=0;i<ncomp;i++) {
    stringstream hname,drawvar;
    hname << "h" << i+1 << cvs << ends;
    h[i] = new TH1D(hname.str().data(),"",nbin,xmin,xmax);
    drawvar << varnames.at(i) << ">>" << hname.str().data() << ends;
    n[idet]->Draw(drawvar.str().data(),cut);
    //    h[i]->SetLineColor(2+i*2);
    Double_t norm = h[i]->Integral();
    h[i]->Scale(1./norm);
    //    h[i]->SetXTitle(varname[i]);
    //      h[i]->SetXTitle("z0 / #deltaz0");      
    if (h[i]->GetMaximum() > ymax) ymax = h[i]->GetMaximum();
    h[i]->SetLineColor(i+1);
  }


  TCanvas *c = new TCanvas(cvs,cvs,800,600);
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);

  for (Int_t i=0;i<ncomp;i++) {
    h[i]->SetMaximum(ymax*1.2);
    h[i]->SetXTitle(xtitle);
    h[i]->SetYTitle("Normalized to 1");
    h[i]->GetYaxis()->CenterTitle();
    if (i==0) {
      h[i]->Draw();
    }
    else {
      h[i]->Draw("same");
    }
    leg->AddEntry(h[i], title.at(i), "l");
  }

  leg->SetFillStyle(0);
  leg->Draw();

  TLatex *ltx = new TLatex(xmin+(xmax-xmin)/20,ymax*1.1,"ILD (preliminary)");
  ltx->SetTextFont(62);
  ltx->Draw();
  
  stringstream savename,savename2,savename3;
  savename << "eps/" << cvs << ".pdf" << ends;
  savename2 << "eps/" << cvs << ".eps" << ends;
  savename3 << "eps/" << cvs << ".C" << ends;
  c->SaveAs(savename.str().data());
  c->SaveAs(savename2.str().data());
  c->SaveAs(savename3.str().data());
}

void fit_trkP(TNtuple *n) {

  setILDStyle();
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  Double_t xmin = -3e-4;
  Double_t xmax = 3e-4;
  Int_t nbin = 100;
  TH1D* h1 = new TH1D("h1","",nbin,xmin,xmax);

  TCanvas *c = new TCanvas("c","",600,600);
  n->Draw("(sqrt(pxlep1**2+pylep1**2)-sqrt(pxlep1mc**2+pylep1mc**2))/(pxlep1mc**2+pylep1mc**2)>>h1","nhbb==1&&elep1mc>100&&abs(coslep1mc)<0.3");
  //  n->Draw("(sqrt(pxlep1**2+pylep1**2)-sqrt(pxlep1mc**2+pylep1mc**2))/(pxlep1mc**2+pylep1mc**2)>>h1","nhbb==1&&sqrt(pxlep1mc**2+pylep1mc**2)>100&&abs(coslep1mc)<0.3");  

  TGaxis::SetMaxDigits(3);
  h1->SetXTitle("#delta(1/Pt) [GeV^{-1}]");  
  h1->SetYTitle("Events");  

  h1->Draw();
  h1->Fit("gaus");

  c->Update();
  //  c->SaveAs("eps/fit_trk_Pgt100_coslt0.3_Large.pdf");
  //  c->SaveAs("eps/fit_trk_Pgt100_coslt0.3_Small.pdf");  

}

void fit_JetAngle(TNtuple *n) {

  setILDStyle();
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  Double_t xmin = -0.2;
  Double_t xmax = 0.2;
  Int_t nbin = 100;
  TH1D* h1 = new TH1D("h1","",nbin,xmin,xmax);

  TCanvas *c = new TCanvas("c","",600,600);
  //  n->Draw("acos(cosj1)-acos(cosj1mc)>>h1","nhbb==1&&ej1mc>100&&abs(cosj1mc)<0.3");
  n->Draw("phij1-phij1mc>>h1","nhbb==1&&ej1mc>100&&abs(cosj1mc)<0.3");  

  TGaxis::SetMaxDigits(3);
  h1->SetXTitle("#theta_{j}^{Rec}-#theta^{MC} [radian]");  
  h1->SetYTitle("Events");  

  h1->Draw();
  h1->Fit("gaus");

  c->Update();
  //  c->SaveAs("eps/fit_trk_Pgt100_coslt0.3_Large.pdf");
  //  c->SaveAs("eps/fit_trk_Pgt100_coslt0.3_Small.pdf");
}

void Plot_JetAngle_E(TNtuple *n) {

  setILDStyle();
  gStyle->SetTitleOffset(1.1,"x");
#if 0
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadLeftMargin(0.14);
#endif  

#if 0
  Double_t xmin = 0;
  Double_t xmax = 200;
#else
  Double_t xmin = -1;
  Double_t xmax = 1;
#endif  
  Int_t nbinx = 100;
  Double_t ymin = -0.2;
  Double_t ymax = 0.2;
  Int_t nbiny = 100;
  TH2D* h1 = new TH2D("h1","",nbinx,xmin,xmax,nbiny,ymin,ymax);

  gStyle->SetPadRightMargin(0.14);
  //  gStyle->SetPadLeftMargin(0.14);

  TCanvas *c = new TCanvas("c","");
  //  n->Draw("acos(cosj1)-acos(cosj1mc):ej1mc>>h1","nhbb==1");
  n->Draw("acos(cosj1)-acos(cosj1mc):cosj1mc>>h1","nhbb==1");

  TGaxis::SetMaxDigits(3);
  //  h1->SetXTitle("E_{j}^{MC} [GeV]");
  h1->SetXTitle("cos#theta_{j}^{MC}");    
  h1->SetYTitle("#theta_{j}^{Rec}-#theta^{MC} [radian]");  

  //  gStyle->SetOptFit(1);
  h1->Draw("colz");

  c->Update();
  //  c->SaveAs("eps/JetAngle_theta_E_DBD.pdf");
  c->SaveAs("eps/JetAngle_theta_cos_DBD.pdf");  

}

void Plot() {
  //  PlotMH();
  //  PlotMH_SB();
  //  PlotRes();
  //  PlotRes(kTRUE);    
  //  PlotPJ();

  const Int_t ncomp = 3;
  TFile *f[ncomp];
  //  f[0] = new TFile("../root_merge.mILD_DBD/e2e2h.eL.pR.root");
  //  f[0] = new TFile("../root_merge.mILD_o1_v05/e2e2h.eL.pR.mILD_o1_v05.root");  
#if 0
  f[0] = new TFile("../root_merge.mILD_DBD/e2e2h_bb.eL.pR.root");    
  f[1] = new TFile("../root_merge.mILD_l5_o1_v02/e2e2h.eL.pR.mILD_l5_o1_v02.root");
  f[2] = new TFile("../root_merge.mILD_s5_o1_v02/e2e2h.eL.pR.mILD_s5_o1_v02.root");
#elif 1 // for final resolution plots in IDR note
  f[0] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.eL.pR.mILD_l5_o1_v02.root");
  f[1] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.eL.pR.mILD_l5_o1_v02.root");
  f[2] = new TFile("../../run_mc/root_merge.mILD_s5_o1_v02/e2e2h.eL.pR.mILD_s5_o1_v02.root");
  //  f[0] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e1e1h.eL.pR.mILD_l5_o1_v02.root");
  //  f[1] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e1e1h.eL.pR.mILD_l5_o1_v02.root");
  //  f[2] = new TFile("../../run_mc/root_merge.mILD_s5_o1_v02/e1e1h.eL.pR.mILD_s5_o1_v02.root");

#elif 0
  f[0] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.eR.pL.mILD_l5_o1_v02.root");
  f[1] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.eR.pL.mILD_l5_o1_v02.root");
  f[2] = new TFile("../../run_mc/root_merge.mILD_s5_o1_v02/e2e2h.eR.pL.mILD_s5_o1_v02.root");
#else // eLpR and eRpL merged, was used for resolution plot in IDR note V1
  f[0] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.mILD_l5_o1_v02.root");
  f[1] = new TFile("../../run_mc/root_merge.mILD_l5_o1_v02/e2e2h.mILD_l5_o1_v02.root");
  f[2] = new TFile("../../run_mc/root_merge.mILD_s5_o1_v02/e2e2h.mILD_s5_o1_v02.root");
#endif  

  TNtuple *n[ncomp],*ngen[ncomp];
  Double_t ntot[ncomp];
  for (Int_t i=0;i<ncomp;i++) {
    n[i] = (TNtuple*)f[i]->Get("hAnl");
    ngen[i] = (TNtuple*)f[i]->Get("hGen");
    ntot[i] = ngen[i]->GetEntries("nhbb==1");
  }


  TCut cut = "nhbb==1";
  //  TCut cut = "nhtt==1";
#if 0
  Comp("(cosj1-cosj1mc)/cosj1mc",n,100,-0.1,0.1,cut,"cosj1");
  Comp("(phij1-phij1mc)/phij1mc",n,100,-0.1,0.1,cut,"phij1");
  Comp("(ej1-ej1mc)/ej1mc",n,100,-0.3,0.3,cut,"ej1");  
  Comp("(elep1-elep1mc)/elep1mc",n,100,-0.01,0.01,cut,"elep1");

  Comp("mz",n,100,40,140,cut);  
  Comp("mh",n,100,60,160,cut);  
  Comp("mrecoil",n,100,100,200,cut);  
  Comp("mhnew",n,100,75,175,cut);
#endif  

  // for IDR
#if 0
  Comp("(elep1-elep1mc)/elep1mc",n,100,-0.02,0.02,cut,"elep1","(E_{#mu}^{Rec}-E^{MC}) / E^{MC}");
  Comp("mz",n,100,40,140,cut,"mz","M(#mu^{+}#mu^{-}) [GeV]");  
  Comp("acos(cosj1)-acos(cosj1mc)",n,100,-0.2,0.2,cut,"thetaj1","#theta_{j}^{Rec}-#theta_{q}^{MC} [radian]");
  Comp("acos(cosj2)-acos(cosj2mc)",n,100,-0.2,0.2,cut,"thetaj2","#theta_{j}^{Rec}-#theta_{q}^{MC} [radian]");  
  Comp("phij1-phij1mc",n,100,-0.2,0.2,cut,"phij1","#phi_{j}^{Rec}-#phi_{q}^{MC} [radian]");
  Comp("(ej1-ej1mc)/ej1mc",n,100,-0.5,0.5,cut,"ej1","(E_{j}^{Rec}-E^{MC}) / E^{MC}");  
  Comp("mh",n,100,60,160,cut,"mh","M(jj) [GeV]");  
  Comp("mrecoil",n,100,100,300,cut,"mrecoil","Recoil Mass [GeV]");  
  Comp("mhnew",n,100,75,175,cut,"mhnew","M(jj) [GeV]");
#endif  
  //  Comp("mhnew",n,80,75,175,cut,"mhnew","M(jj) [GeV]");
#if 0
  Comp("sqrt(ej1*ej1-pxj1*pxj1-pyj1*pyj1-pzj1*pzj1)",n,100,0,100,cut,"mj1","M(j) [GeV]");
  Comp("mj1",n,100,0,100,cut,"mj1","m_{j}^{Rec} [GeV]");
  Comp("sqrt(ej1mc*ej1mc-pxj1mc*pxj1mc-pyj1mc*pyj1mc-pzj1mc*pzj1mc)",n,100,0,100,cut,"mj1mc","m_{q}^{MC} [GeV]");    
  Comp("mj1h",n,100,0,100,cut,"mj1h","m_{j}^{MC} [GeV]");
  Comp("(mj1-mj1h)/mj1h",n,100,-1,1,cut,"mj1res","(m_{j}^{Rec}-m_{j}^{MC}) / m_{j}^{MC}");
  Comp("acos(cosj1)-acos(cosj1h)",n,100,-0.1,0.1,cut,"thetaj1h","#theta_{j}^{Rec}-#theta_{j}^{MC} [radian]");
  Comp("phij1-phij1h",n,100,-0.1,0.1,cut,"phij1h","#phi_{j}^{Rec}-#phi_{j}^{MC} [radian]");
#endif  
  //  Comp("(elep1-elep1mc)/elep1mc",n,100,-0.02,0.02,cut,"elep1","(E_{e}^{Rec}-E^{MC}) / E^{MC}");

  // finale  
  //  Comp("(elep1-elep1mc)/elep1mc",n,100,-0.02,0.02,cut,"elep1","(E_{#mu}^{Rec}-E^{MC}) / E^{MC}");
  //  Comp("phij1-phij1h",n,100,-0.1,0.1,cut,"phij1h","#phi_{j}^{Rec}-#phi_{j}^{MC} [radian]");
  //  Comp("acos(cosj1)-acos(cosj1h)",n,100,-0.1,0.1,cut,"thetaj1h","#theta_{j}^{Rec}-#theta_{j}^{MC} [radian]");

  //  Comp("acos(cosj1)-acos(cosj1h)",n,100,-0.05,0.05,cut+"abs(cosj1)<0.6&&ej1>150","thetaj1h_test","#theta_{j}^{Rec}-#theta_{j}^{MC} [radian]");
  //  Comp("phij1-phij1h",n,40,-0.03,0.03,cut+"abs(cosj1)<0.3","phij1h_test","#phi_{j}^{Rec}-#phi_{j}^{MC} [radian]");
  //  Comp("phij1-phij1h",n,40,-0.03,0.03,cut,"phij1h_test","#phi_{j}^{Rec}-#phi_{j}^{MC} [radian]");    
  //  Comp("mhnew",n,80,75,175,cut,"mhnew_test","M(jj) [GeV]");
  //  Comp("mhnew",n,80,75,175,cut,"mhnew","M(jj) [GeV]");  

  //  Comp("phij1-phij1mc",n,100,-0.1,0.1,cut,"phij1_test","#phi_{j}^{Rec}-#phi^{MC} [radian]");  
  

  
  //  Comp("(sqrt(pxlep1**2+pylep1**2)-sqrt(pxlep1mc**2+pylep1mc**2))/(pxlep1mc**2+pylep1mc**2)",n,100,-30e-5,30e-5,cut,"ptlep1","(Pt_{#mu^{-}}^{Rec}-Pt^{MC}) / (Pt^{MC})");  

  //  fit_trkP(n[0]);
  //  fit_trkP(n[1]);
  //  fit_trkP(n[2]);
  //  Plot_JetAngle_E(n[0]);
  //  fit_JetAngle(n[1]);
  
  
  //  Comp("(acos(cosj1)-acos(cosj1mc))/acos(cosj1mc)",n,100,-0.1,0.1,cut,"thetaj1");
  //  Comp("(phij1-phij1mc)/phij1mc",n,100,-0.1,0.1,cut,"phij1");

  
  //  n1->Draw("(cosj1-cosj1mc)/cosj1mc>>h1");
  //  n1->Draw("(phij1-phij1mc)/phij1mc>>h2");
  //  n1->Draw("(ej1-ej1mc)/ej1mc>>h3");

  //  PlotMH();

  //  Comp("mhnew",n,80,75,175,cut,"mhnew_test","M(jj) [GeV]");


  //  Comp("mhnew",n,80,75,175,cut,"mhnew_test","Higgs Mass [GeV]");
  //  Comp("mhnew",n,50,100,150,cut,"mhnew_test","Higgs Mass [GeV]");

  // fit
  _nbin_fit = 5;
  //  Comp("mhnew",n,50,100,150,cut,"mhnew_fit","Higgs Mass [GeV]");        
  _nbin_fit = 5;
  //  Comp("(elep1-elep1mc)/elep1mc",n,50,-0.02,0.02,cut,"elep1_fit","(E_{#mu}^{Rec}-E^{MC}) / E^{MC}");
  _nbin_fit = 3;
  //  Comp("phij1-phij1h",n,50,-0.05,0.05,cut,"phij1h_fit","#phi_{j}^{Rec}-#phi_{j}^{MC} [radian]");
  _nbin_fit = 3;
  //  Comp("acos(cosj1)-acos(cosj1h)",n,50,-0.05,0.05,cut,"thetaj1h_fit","#theta_{j}^{Rec}-#theta_{j}^{MC} [radian]");
  _nbin_fit = 5;
  //  Comp("mj1h",n,100,0,100,cut,"mj1h","m_{j}^{MC} [GeV]");
  //  Comp("(mj1-mj1h)/mj1h",n,50,-1,1,cut,"mj1res_fit","(m_{j}^{Rec}-m_{j}^{MC}) / m_{j}^{MC}");
  _nbin_fit = 4;
  Comp("mhjet",n,50,122,128,cut,"mhmc_jetmass_jetdir_compare_fit","Higgs Mass [GeV]");  
  
#if 0
  Comp("mhmuon",n,50,100,150,cut,"mhmc_muon_test","M(jj) [GeV]");  
  Comp("mhjet2",n,60,110,140,cut,"mhmc_jetdir_test","M(jj) [GeV]");  
  Comp("mhjet3",n,60,110,140,cut,"mhmc_jetmass_test","M(jj) [GeV]");  
  Comp("mhjet",n,50,120,130,cut,"mhmc_jetmass_jetdir_test","M(jj) [GeV]");  
  Comp("mhjet4",n,60,110,140,cut,"mhmc_jetdir_muon_test","M(jj) [GeV]");  
  Comp("mhjet5",n,60,110,140,cut,"mhmc_jetmass_muon_test","M(jj) [GeV]");
#endif  
  //  Comp("mhjet",n,50,120,130,cut,"mhmc_jetmass_jetdir_compare","Higgs Mass [GeV]");  

  
#if 0
  std::vector<TString> vars[10];
  vars->push_back("mhnew");
  for (Int_t i=0;i<1;i++) {
    TString var = vars->at(i);
    cerr << var  << endl;
  }
#endif  

#if 0
  std::vector<TString> vars, title;
  vars.push_back( "mhnew" );
  vars.push_back( "mhmuon" );
  vars.push_back( "mhjet2" );
  vars.push_back( "mhjet3" );
  vars.push_back( "mhjet4" );
  vars.push_back( "mhjet5" );
  vars.push_back( "mhjet" );

  title.push_back("Fully Reconstrcuted");
  title.push_back("Cheated Muon 4-momenta");
  title.push_back("Cheated Jet directions");
  title.push_back("Cheated Jet masses");
  title.push_back("Cheated Jet directions & Muon 4-momenta");
  title.push_back("Cheated Jet masses & Muon 4-momenta");
  title.push_back("Cheated Jet directions & Jet masses");
  for (Int_t i=0;i<4;i++) {
    TString var = vars.at(i);
    TString tit = title.at(i);
    cerr << var << " " << tit << endl;
  }
  Int_t nc = 7;
  Comp(1,nc,vars,n,60,110,140,cut,"mhmc_test",title,"Reconstructed Higgs Mass [GeV]"); // IDR-L
  //  Comp(2,2,vars,n,60,110,140,cut,"mhmc_test",title,"Reconstructed Higgs Mass [GeV]"); // IDR-S
#endif  

#if 0 // muon & jet direction
  std::vector<TString> vars, title;
  vars.push_back( "mhnew" );
  vars.push_back( "mhmuon" );
  vars.push_back( "mhjet2" );
  vars.push_back( "mhjet4" );

  title.push_back("Fully Reconstrcuted");
  title.push_back("MC for Muon 4-momenta");
  title.push_back("MC for Jet directions");
  title.push_back("MC for Both");
  for (Int_t i=0;i<4;i++) {
    TString var = vars.at(i);
    TString tit = title.at(i);
    cerr << var << " " << tit << endl;
  }
  Int_t nc = 4;
  Comp(1,nc,vars,n,120,100,160,cut,"mhmc_muon_jetdir_IDR-L",title,"Higgs Mass [GeV]"); // IDR-L
#endif  

#if 0 // muon & jet mass
  std::vector<TString> vars, title;
  vars.push_back( "mhnew" );
  vars.push_back( "mhmuon" );
  vars.push_back( "mhjet3" );
  vars.push_back( "mhjet5" );

  title.push_back("Fully Reconstrcuted");
  title.push_back("MC for Muon 4-momenta");
  title.push_back("MC for Jet masses");
  title.push_back("MC for Both");
  for (Int_t i=0;i<4;i++) {
    TString var = vars.at(i);
    TString tit = title.at(i);
    cerr << var << " " << tit << endl;
  }
  Int_t nc = 4;
  Comp(1,nc,vars,n,120,100,160,cut,"mhmc_muon_jetmass_IDR-L",title,"Higgs Mass [GeV]"); // IDR-L
#endif  

#if 0 // jet direction & jet mass
  std::vector<TString> vars, title;
  vars.push_back( "mhnew" );
  vars.push_back( "mhjet2" );
  vars.push_back( "mhjet3" );
  vars.push_back( "mhjet" );

  title.push_back("Fully Reconstrcuted");
  title.push_back("MC for Jet directions");
  title.push_back("MC for Jet masses");
  title.push_back("MC for Both");
  for (Int_t i=0;i<4;i++) {
    TString var = vars.at(i);
    TString tit = title.at(i);
    cerr << var << " " << tit << endl;
  }
  Int_t nc = 4;
  Comp(1,nc,vars,n,120,100,160,cut,"mhmc_jetdir_jetmass_IDR-L",title,"Higgs Mass [GeV]"); // IDR-L
#endif  
  
  
  //  CompTest(1,4,vars,n,60,110,140,cut,"mhmc_test",title,"Reconstructed Higgs Mass [GeV]");
  //  CompTest(1,4,n,60,110,140,cut,"mhmc_test","Reconstructed Higgs Mass [GeV]");    
}
