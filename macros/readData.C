#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "global.h"

using namespace std;

void read6f() {
  std::ifstream in;
  in.open("500_6f.list");
  Double_t ecm,xsec,xseclr,xsecrl;
  TString pname,ptype,pol;
  Int_t pid;
  while (1) {
    in >> ecm >> pid >> pname >> ptype >> pol >> xsec >> xseclr >> xsecrl;
    if (!in.good()) break;
    //    cerr << pname.Data() << " " << pol.Data() << " " << xsec << endl;
    pname6f.push_back(pname);
    pol6f.push_back(pol);
    xsec6f.push_back(xsec);    
  }
  in.close();

}

void readData(Int_t mode) {

  if (mode == 1) { // test data
#ifdef __IDRL__
#ifndef __electron__
    fS[ 1]   = new TFile("data_final/e2e2h.eL.pR.mILD_l5_o1_v02.root");
    fS[ 2]   = new TFile("data_final/e2e2h.eR.pL.mILD_l5_o1_v02.root");
#else
    fS[ 1]   = new TFile("data_final/e1e1h.eL.pL.mILD_l5_o1_v02.root");
    fS[ 2]   = new TFile("data_final/e1e1h.eL.pR.mILD_l5_o1_v02.root");
    fS[ 3]   = new TFile("data_final/e1e1h.eR.pL.mILD_l5_o1_v02.root");
    fS[ 4]   = new TFile("data_final/e1e1h.eR.pR.mILD_l5_o1_v02.root");
#endif    
    fB[ 1]   = new TFile("data_final/2f_z_bhabhag.eL.pL.mILD_l5_o1_v02.root");
    fB[ 2]   = new TFile("data_final/2f_z_bhabhag.eL.pR.mILD_l5_o1_v02.root");
    fB[ 3]   = new TFile("data_final/2f_z_bhabhag.eR.pL.mILD_l5_o1_v02.root");
    fB[ 4]   = new TFile("data_final/2f_z_bhabhag.eR.pR.mILD_l5_o1_v02.root");
    fB[ 5]   = new TFile("data_final/2f_z_l.eL.pR.mILD_l5_o1_v02.root");
    fB[ 6]   = new TFile("data_final/2f_z_l.eR.pL.mILD_l5_o1_v02.root");
    fB[ 7]   = new TFile("data_final/2f_z_h.eL.pR.mILD_l5_o1_v02.root");
    fB[ 8]   = new TFile("data_final/2f_z_h.eR.pL.mILD_l5_o1_v02.root");
    fB[ 9]   = new TFile("data_final/4f_zz_l.eL.pR.mILD_l5_o1_v02.root");
    fB[10]   = new TFile("data_final/4f_zz_l.eR.pL.mILD_l5_o1_v02.root");
    fB[11]   = new TFile("data_final/4f_zz_h.eL.pR.mILD_l5_o1_v02.root");
    fB[12]   = new TFile("data_final/4f_zz_h.eR.pL.mILD_l5_o1_v02.root");
    fB[13]   = new TFile("data_final/4f_zz_sl.eL.pR.mILD_l5_o1_v02.root");
    fB[14]   = new TFile("data_final/4f_zz_sl.eR.pL.mILD_l5_o1_v02.root");
    fB[15]   = new TFile("data_final/4f_ww_l.eL.pR.mILD_l5_o1_v02.root");
    fB[16]   = new TFile("data_final/4f_ww_l.eR.pL.mILD_l5_o1_v02.root");
    fB[17]   = new TFile("data_final/4f_ww_h.eL.pR.mILD_l5_o1_v02.root");
    fB[18]   = new TFile("data_final/4f_ww_h.eR.pL.mILD_l5_o1_v02.root");
    fB[19]   = new TFile("data_final/4f_ww_sl.eL.pR.mILD_l5_o1_v02.root");
    fB[20]   = new TFile("data_final/4f_ww_sl.eR.pL.mILD_l5_o1_v02.root");
    fB[21]   = new TFile("data_final/4f_zzorww_l.eL.pR.mILD_l5_o1_v02.root");
    fB[22]   = new TFile("data_final/4f_zzorww_l.eR.pL.mILD_l5_o1_v02.root");
    fB[23]   = new TFile("data_final/4f_zzorww_h.eL.pR.mILD_l5_o1_v02.root");
    fB[24]   = new TFile("data_final/4f_zzorww_h.eR.pL.mILD_l5_o1_v02.root");
    fB[25]   = new TFile("data_final/4f_sznu_l.eL.pR.mILD_l5_o1_v02.root");
    fB[26]   = new TFile("data_final/4f_sznu_l.eR.pL.mILD_l5_o1_v02.root");
    fB[27]   = new TFile("data_final/4f_sznu_sl.eL.pR.mILD_l5_o1_v02.root");
    fB[28]   = new TFile("data_final/4f_sznu_sl.eR.pL.mILD_l5_o1_v02.root");
    fB[29]   = new TFile("data_final/4f_sze_l.eL.pL.mILD_l5_o1_v02.root");
    fB[30]   = new TFile("data_final/4f_sze_l.eL.pR.mILD_l5_o1_v02.root");
    fB[31]   = new TFile("data_final/4f_sze_l.eR.pL.mILD_l5_o1_v02.root");
    fB[32]   = new TFile("data_final/4f_sze_l.eR.pR.mILD_l5_o1_v02.root");
    fB[33]   = new TFile("data_final/4f_sze_sl.eL.pL.mILD_l5_o1_v02.root");
    fB[34]   = new TFile("data_final/4f_sze_sl.eL.pR.mILD_l5_o1_v02.root");
    fB[35]   = new TFile("data_final/4f_sze_sl.eR.pL.mILD_l5_o1_v02.root");
    fB[36]   = new TFile("data_final/4f_sze_sl.eR.pR.mILD_l5_o1_v02.root");
    fB[37]   = new TFile("data_final/4f_sw_l.eL.pL.mILD_l5_o1_v02.root");
    fB[38]   = new TFile("data_final/4f_sw_l.eL.pR.mILD_l5_o1_v02.root");
    fB[39]   = new TFile("data_final/4f_sw_l.eR.pL.mILD_l5_o1_v02.root");
    fB[40]   = new TFile("data_final/4f_sw_l.eR.pR.mILD_l5_o1_v02.root");
    fB[41]   = new TFile("data_final/4f_sw_sl.eL.pL.mILD_l5_o1_v02.root");
    fB[42]   = new TFile("data_final/4f_sw_sl.eL.pR.mILD_l5_o1_v02.root");
    fB[43]   = new TFile("data_final/4f_sw_sl.eR.pL.mILD_l5_o1_v02.root");
    fB[44]   = new TFile("data_final/4f_sw_sl.eR.pR.mILD_l5_o1_v02.root");
    fB[45]   = new TFile("data_final/4f_szeorsw_l.eL.pL.mILD_l5_o1_v02.root");
    fB[46]   = new TFile("data_final/4f_szeorsw_l.eL.pR.mILD_l5_o1_v02.root");
    fB[47]   = new TFile("data_final/4f_szeorsw_l.eR.pL.mILD_l5_o1_v02.root");
    fB[48]   = new TFile("data_final/4f_szeorsw_l.eR.pR.mILD_l5_o1_v02.root");
#else
#ifndef __electron__
    fS[ 1]   = new TFile("data_final.small/e2e2h.eL.pR.mILD_s5_o1_v02.root");
    fS[ 2]   = new TFile("data_final.small/e2e2h.eR.pL.mILD_s5_o1_v02.root");
#else
    fS[ 1]   = new TFile("data_final.small/e1e1h.eL.pL.mILD_s5_o1_v02.root");
    fS[ 2]   = new TFile("data_final.small/e1e1h.eL.pR.mILD_s5_o1_v02.root");
    fS[ 3]   = new TFile("data_final.small/e1e1h.eR.pL.mILD_s5_o1_v02.root");
    fS[ 4]   = new TFile("data_final.small/e1e1h.eR.pR.mILD_s5_o1_v02.root");
#endif    

    fB[ 1]   = new TFile("data_final.small/2f_z_bhabhag.eL.pL.mILD_s5_o1_v02.root");
    fB[ 2]   = new TFile("data_final.small/2f_z_bhabhag.eL.pR.mILD_s5_o1_v02.root");
    fB[ 3]   = new TFile("data_final.small/2f_z_bhabhag.eR.pL.mILD_s5_o1_v02.root");
    fB[ 4]   = new TFile("data_final.small/2f_z_bhabhag.eR.pR.mILD_s5_o1_v02.root");
    fB[ 5]   = new TFile("data_final.small/2f_z_l.eL.pR.mILD_s5_o1_v02.root");
    fB[ 6]   = new TFile("data_final.small/2f_z_l.eR.pL.mILD_s5_o1_v02.root");
    fB[ 7]   = new TFile("data_final.small/2f_z_h.eL.pR.mILD_s5_o1_v02.root");
    fB[ 8]   = new TFile("data_final.small/2f_z_h.eR.pL.mILD_s5_o1_v02.root");
    fB[ 9]   = new TFile("data_final.small/4f_zz_l.eL.pR.mILD_s5_o1_v02.root");
    fB[10]   = new TFile("data_final.small/4f_zz_l.eR.pL.mILD_s5_o1_v02.root");
    fB[11]   = new TFile("data_final.small/4f_zz_h.eL.pR.mILD_s5_o1_v02.root");
    fB[12]   = new TFile("data_final.small/4f_zz_h.eR.pL.mILD_s5_o1_v02.root");
    fB[13]   = new TFile("data_final.small/4f_zz_sl.eL.pR.mILD_s5_o1_v02.root");
    fB[14]   = new TFile("data_final.small/4f_zz_sl.eR.pL.mILD_s5_o1_v02.root");
    fB[15]   = new TFile("data_final.small/4f_ww_l.eL.pR.mILD_s5_o1_v02.root");
    fB[16]   = new TFile("data_final.small/4f_ww_l.eR.pL.mILD_s5_o1_v02.root");
    fB[17]   = new TFile("data_final.small/4f_ww_h.eL.pR.mILD_s5_o1_v02.root");
    fB[18]   = new TFile("data_final.small/4f_ww_h.eR.pL.mILD_s5_o1_v02.root");
    fB[19]   = new TFile("data_final.small/4f_ww_sl.eL.pR.mILD_s5_o1_v02.root");
    fB[20]   = new TFile("data_final.small/4f_ww_sl.eR.pL.mILD_s5_o1_v02.root");
    fB[21]   = new TFile("data_final.small/4f_zzorww_l.eL.pR.mILD_s5_o1_v02.root");
    fB[22]   = new TFile("data_final.small/4f_zzorww_l.eR.pL.mILD_s5_o1_v02.root");
    fB[23]   = new TFile("data_final.small/4f_zzorww_h.eL.pR.mILD_s5_o1_v02.root");
    fB[24]   = new TFile("data_final.small/4f_zzorww_h.eR.pL.mILD_s5_o1_v02.root");
    fB[25]   = new TFile("data_final.small/4f_sznu_l.eL.pR.mILD_s5_o1_v02.root");
    fB[26]   = new TFile("data_final.small/4f_sznu_l.eR.pL.mILD_s5_o1_v02.root");
    fB[27]   = new TFile("data_final.small/4f_sznu_sl.eL.pR.mILD_s5_o1_v02.root");
    fB[28]   = new TFile("data_final.small/4f_sznu_sl.eR.pL.mILD_s5_o1_v02.root");
    fB[29]   = new TFile("data_final.small/4f_sze_l.eL.pL.mILD_s5_o1_v02.root");
    fB[30]   = new TFile("data_final.small/4f_sze_l.eL.pR.mILD_s5_o1_v02.root");
    fB[31]   = new TFile("data_final.small/4f_sze_l.eR.pL.mILD_s5_o1_v02.root");
    fB[32]   = new TFile("data_final.small/4f_sze_l.eR.pR.mILD_s5_o1_v02.root");
    fB[33]   = new TFile("data_final.small/4f_sze_sl.eL.pL.mILD_s5_o1_v02.root");
    fB[34]   = new TFile("data_final.small/4f_sze_sl.eL.pR.mILD_s5_o1_v02.root");
    fB[35]   = new TFile("data_final.small/4f_sze_sl.eR.pL.mILD_s5_o1_v02.root");
    fB[36]   = new TFile("data_final.small/4f_sze_sl.eR.pR.mILD_s5_o1_v02.root");
    fB[37]   = new TFile("data_final.small/4f_sw_l.eL.pL.mILD_s5_o1_v02.root");
    fB[38]   = new TFile("data_final.small/4f_sw_l.eL.pR.mILD_s5_o1_v02.root");
    fB[39]   = new TFile("data_final.small/4f_sw_l.eR.pL.mILD_s5_o1_v02.root");
    fB[40]   = new TFile("data_final.small/4f_sw_l.eR.pR.mILD_s5_o1_v02.root");
    fB[41]   = new TFile("data_final.small/4f_sw_sl.eL.pL.mILD_s5_o1_v02.root");
    fB[42]   = new TFile("data_final.small/4f_sw_sl.eL.pR.mILD_s5_o1_v02.root");
    fB[43]   = new TFile("data_final.small/4f_sw_sl.eR.pL.mILD_s5_o1_v02.root");
    fB[44]   = new TFile("data_final.small/4f_sw_sl.eR.pR.mILD_s5_o1_v02.root");
    fB[45]   = new TFile("data_final.small/4f_szeorsw_l.eL.pL.mILD_s5_o1_v02.root");
    fB[46]   = new TFile("data_final.small/4f_szeorsw_l.eL.pR.mILD_s5_o1_v02.root");
    fB[47]   = new TFile("data_final.small/4f_szeorsw_l.eR.pL.mILD_s5_o1_v02.root");
    fB[48]   = new TFile("data_final.small/4f_szeorsw_l.eR.pR.mILD_s5_o1_v02.root");
#endif    
  }
  else if (mode == 2) { // for training
  }

#if 0
  read6f();

  if (pname6f.size() != pol6f.size() || pname6f.size() != xsec6f.size()) {
    cerr << "something is inconsistent for 6f vectors!" << endl;
    return;
  }
  Int_t np6f = pname6f.size();
  for (Int_t i=0;i<np6f;i++) {
    stringstream file;
    file << "data_final/" << pname6f.at(i).Data() << "." << pol6f.at(i).Data()
	 << ".root" << ends;
    //    cerr << file.str().data() << endl;
    fB[49+i] = new TFile(file.str().data());
  }
#endif  

  for (Int_t i=1;i<nSignal;i++) {
    hAnl_S[i] = (TNtuple*)fS[i]->Get("hAnl");
    hGen_S[i] = (TNtuple*)fS[i]->Get("hGen");
    ngenS[i]  = hGen_S[i]->GetEntries();
  }

  for (Int_t i=1;i<nBackground;i++) {
    hAnl_B[i] = (TNtuple*)fB[i]->Get("hAnl");
    hGen_B[i] = (TNtuple*)fB[i]->Get("hGen");
    ngenB[i]  = hGen_B[i]->GetEntries();
  }

}
