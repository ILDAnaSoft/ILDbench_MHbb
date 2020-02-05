#include <iostream>
#include <sstream>
#include <iomanip>
#include "readPara.C"
#include "readData.C"
#include "readCut.C"
#include "global.h"

#include "TROOT.h"
#include "TNtuple.h"

#include "TMath.h"


using namespace std;

void cuts(Int_t iPol = -1,Bool_t iPrintAll = kFALSE) {

  if (iPol != 1 && iPol !=0 && iPol != -1) {
    cerr << "Please use a correct Polarization:  0,1,-1" << endl;
    return;
  }

  // get cuts
  readCut(myCut,iPol);

  // get ntuple data
  readData(1);

  // get x-section, # of generated events, weight, etc.
  readPara(1,iPol,inteLumi);
  // output the Nexp, Ngen, weight
  cerr << "Expected Number of Events" << " | " << "Generated Number of Events" << " | " << "Weight of Each Event" << endl;
  for (Int_t i=1;i<nSignal;i++) {
    cerr << "S" << setw(2) << i << ": " << setw(20) << nexpS[i] << setw(20) << ngenS[i] << setw(20) << weightS[i] << endl;
  }
  for (Int_t i=1;i<nBackground;i++) {
    cerr << "B" << setw(2) << i << ": " << setw(20) << nexpB[i] << setw(20) << ngenB[i] << setw(20) << weightB[i] << endl;
  }

  Int_t iZDecayMode = 13;
  // cut flow
  Double_t temp[nCutMax];
  Double_t temp0[nCutMax];
  Double_t nS[nSignal][nCutMax],nB[nBackground][nCutMax];
  Double_t nS0[nSignal][nCutMax];
  // signal
  for (Int_t i=1;i<nSignal;i++) {
    if (hAnl_S[i]) {
      TCut cut_temp = "";
      for (Int_t j=0;j<nCutMax;j++) {
	cut_temp += myCut[j];
	temp[j] = hAnl_S[i]->GetEntries(cut_temp);
	temp0[j] = hAnl_S[i]->GetEntries(cut_temp+"nhbb==1");
      }
    }
    else {
      for (Int_t k=0;k<nCutMax;k++) temp[k] = 0;
      for (Int_t k=0;k<nCutMax;k++) temp0[k] = 0;
    }
    for (Int_t j=0;j<nCutMax;j++) {
      nS[i][j] = weightS[i]*temp[j];
      nS0[i][j] = weightS[i]*temp0[j];
      cerr << "debug: " << nS0[i][j] << endl;
    }
  }
  // background
  for (Int_t i=1;i<nBackground;i++) {
    if (hAnl_B[i]) {
      TCut cut_temp = "";
      for (Int_t j=0;j<nCutMax;j++) {
	cut_temp += myCut[j];
	temp[j] = hAnl_B[i]->GetEntries(cut_temp);
      }
    }
    else {
      for (Int_t k=0;k<nCutMax;k++) temp[k] = 0;
    }
    for (Int_t j=0;j<nCutMax;j++) {
      nB[i][j] = weightB[i]*temp[j];
    }
  }
  if (iPrintAll) {
    cerr << "Expected" << " " << "selction" << endl;
    cerr << "------signal-------" << endl;
    for (Int_t i=1;i<nSignal;i++) {
      cerr << endl << "S" << i << ": " << nexpS[i] << " ";
      for (Int_t j=0;j<nCutMax;j++) {
	cerr << nS[i][j] << "  ";
      }
    }
    cerr << endl << "------background-------";
    for (Int_t i=1;i<nBackground;i++) {
      //      if (i>48) cerr << endl << pname6f.at(i-49).Data() << " " << pol6f.at(i-49).Data() << endl;
      cerr << endl << "B" << i << ": " << nexpB[i] << " ";
      for (Int_t j=0;j<nCutMax;j++) {
	cerr << nB[i][j] << "  ";
      }
    }
  }

  // group signal and background processes
  const Int_t nSSum = 1 + 1;
  const Int_t nBSum = 5 + 1;
  //# of expected
  Double_t nexpSS[nSSum];
  for (Int_t i=0;i<nSSum;i++) nexpSS[i] = 0.;
  Double_t nexpBS[nBSum];
  for (Int_t i=0;i<nBSum;i++) nexpBS[i] = 0.;
  // group
  for (Int_t i=1;i<=2;i++)    nexpSS[1] += nexpS[i];   // ax
  for (Int_t i=1;i<=6;i++)    nexpBS[1] += nexpB[i];   // 2f_l
  for (Int_t i=7;i<=8;i++)    nexpBS[2] += nexpB[i];   // 2f_h
  for (Int_t i=9;i<=10;i++)   nexpBS[3] += nexpB[i];   // 4f_zz_l
  for (Int_t i=15;i<=16;i++)  nexpBS[3] += nexpB[i];   // 4f_ww_l
  for (Int_t i=21;i<=22;i++)  nexpBS[3] += nexpB[i];   // 4f_zzorww_l
  for (Int_t i=25;i<=26;i++)  nexpBS[3] += nexpB[i];   // 4f_sznu_l
  for (Int_t i=29;i<=32;i++)  nexpBS[3] += nexpB[i];   // 4f_sze_l
  for (Int_t i=37;i<=40;i++)  nexpBS[3] += nexpB[i];   // 4f_sw_l
  for (Int_t i=45;i<=48;i++)  nexpBS[3] += nexpB[i];   // 4f_szeorsw_l
  for (Int_t i=13;i<=14;i++)  nexpBS[4] += nexpB[i];   // 4f_zz_sl
  for (Int_t i=19;i<=20;i++)  nexpBS[4] += nexpB[i];   // 4f_ww_sl
  for (Int_t i=27;i<=28;i++)  nexpBS[4] += nexpB[i];   // 4f_sznu_sl
  for (Int_t i=33;i<=36;i++)  nexpBS[4] += nexpB[i];   // 4f_sze_sl
  for (Int_t i=41;i<=44;i++)  nexpBS[4] += nexpB[i];   // 4f_sw_sl
  for (Int_t i=11;i<=12;i++)  nexpBS[5] += nexpB[i];   // 4f_zz_h
  for (Int_t i=17;i<=18;i++)  nexpBS[5] += nexpB[i];   // 4f_ww_h
  for (Int_t i=23;i<=24;i++)  nexpBS[5] += nexpB[i];   // 4f_zzorww_h
  //  for (Int_t i=49;i<nBackground;i++)  nexpBS[6] += nexpB[i];   // 6f
  Double_t nexpSG = 0.;
  for (Int_t i=1;i<nSSum;i++) {
    nexpSG += nexpSS[i];
  }
  Double_t nexpBG   =  0.;
  for (Int_t i=1;i<nBSum;i++) {
    nexpBG += nexpBS[i];
  }
  Double_t nexpSig  =  nexpSG/TMath::Sqrt(nexpSG+nexpBG);
  //# of generated
  Double_t ngenSS[nSSum];
  for (Int_t i=0;i<nSSum;i++) ngenSS[i] = 0.;
  Double_t ngenBS[nBSum];
  for (Int_t i=0;i<nBSum;i++) ngenBS[i] = 0.;
  // group
  for (Int_t i=1;i<=2;i++)    ngenSS[1] += ngenS[i];   // ax_inv
  for (Int_t i=1;i<=6;i++)    ngenBS[1] += ngenB[i];   // 2f_l
  for (Int_t i=7;i<=8;i++)    ngenBS[2] += ngenB[i];   // 2f_h
  for (Int_t i=9;i<=10;i++)   ngenBS[3] += ngenB[i];   // 4f_zz_l
  for (Int_t i=15;i<=16;i++)  ngenBS[3] += ngenB[i];   // 4f_ww_l
  for (Int_t i=21;i<=22;i++)  ngenBS[3] += ngenB[i];   // 4f_zzorww_l
  for (Int_t i=25;i<=26;i++)  ngenBS[3] += ngenB[i];   // 4f_sznu_l
  for (Int_t i=29;i<=32;i++)  ngenBS[3] += ngenB[i];   // 4f_sze_l
  for (Int_t i=37;i<=40;i++)  ngenBS[3] += ngenB[i];   // 4f_sw_l
  for (Int_t i=45;i<=48;i++)  ngenBS[3] += ngenB[i];   // 4f_szeorsw_l
  for (Int_t i=13;i<=14;i++)  ngenBS[4] += ngenB[i];   // 4f_zz_sl
  for (Int_t i=19;i<=20;i++)  ngenBS[4] += ngenB[i];   // 4f_ww_sl
  for (Int_t i=27;i<=28;i++)  ngenBS[4] += ngenB[i];   // 4f_sznu_sl
  for (Int_t i=33;i<=36;i++)  ngenBS[4] += ngenB[i];   // 4f_sze_sl
  for (Int_t i=41;i<=44;i++)  ngenBS[4] += ngenB[i];   // 4f_sw_sl
  for (Int_t i=11;i<=12;i++)  ngenBS[5] += ngenB[i];   // 4f_zz_h
  for (Int_t i=17;i<=18;i++)  ngenBS[5] += ngenB[i];   // 4f_ww_h
  for (Int_t i=23;i<=24;i++)  ngenBS[5] += ngenB[i];   // 4f_zzorww_h
  //  for (Int_t i=49;i<nBackground;i++)  ngenBS[6] += ngenB[i];   // 6f_zzorww_h  
  Double_t ngenSG = 0.;
  for (Int_t i=1;i<nSSum;i++) {
    ngenSG += ngenSS[i];
  }
  Double_t ngenBG   =  0.;
  for (Int_t i=1;i<nBSum;i++) {
    ngenBG += ngenBS[i];
  }
  //# after each cut
  Double_t nSS[nSSum][nCutMax];
  Double_t nSS0[nSSum][nCutMax];
  Double_t nBS[nBSum][nCutMax];
  Double_t nSG[nCutMax],nBG[nCutMax],nSig[nCutMax],nErrorBG[nCutMax];
  Double_t nSG0[nCutMax];
  for (Int_t j=0;j<nCutMax;j++) {
    nSG[j]  = 0.; 
    nBG[j]  = 0.; 
    for (Int_t i=0;i<nSSum;i++) nSS[i][j] = 0.;
    for (Int_t i=0;i<nSSum;i++) nSS0[i][j]= 0.;
    for (Int_t i=0;i<nBSum;i++) nBS[i][j] = 0.;
    for (Int_t i=1;i<=2;i++)    nSS[1][j] += nS[i][j];   // llh
    for (Int_t i=1;i<=2;i++)    nSS0[1][j]+= nS0[i][j];   // llh_bb
    for (Int_t i=1;i<=6;i++)    nBS[1][j] += nB[i][j];   // 2f_l
    for (Int_t i=7;i<=8;i++)    nBS[2][j] += nB[i][j];   // 2f_h
    for (Int_t i=9;i<=10;i++)   nBS[3][j] += nB[i][j];   // 4f_zz_l
    for (Int_t i=15;i<=16;i++)  nBS[3][j] += nB[i][j];   // 4f_ww_l
    for (Int_t i=21;i<=22;i++)  nBS[3][j] += nB[i][j];   // 4f_zzorww_l
    for (Int_t i=25;i<=26;i++)  nBS[3][j] += nB[i][j];   // 4f_sznu_l
    for (Int_t i=29;i<=32;i++)  nBS[3][j] += nB[i][j];   // 4f_sze_l
    for (Int_t i=37;i<=40;i++)  nBS[3][j] += nB[i][j];   // 4f_sw_l
    for (Int_t i=45;i<=48;i++)  nBS[3][j] += nB[i][j];   // 4f_szeorsw_l
    for (Int_t i=13;i<=14;i++)  nBS[4][j] += nB[i][j];   // 4f_zz_sl
    for (Int_t i=19;i<=20;i++)  nBS[4][j] += nB[i][j];   // 4f_ww_sl
    for (Int_t i=27;i<=28;i++)  nBS[4][j] += nB[i][j];   // 4f_sznu_sl
    for (Int_t i=33;i<=36;i++)  nBS[4][j] += nB[i][j];   // 4f_sze_sl
    for (Int_t i=41;i<=44;i++)  nBS[4][j] += nB[i][j];   // 4f_sw_sl
    for (Int_t i=11;i<=12;i++)  nBS[5][j] += nB[i][j];   // 4f_zz_h
    for (Int_t i=17;i<=18;i++)  nBS[5][j] += nB[i][j];   // 4f_ww_h
    for (Int_t i=23;i<=24;i++)  nBS[5][j] += nB[i][j];   // 4f_zzorww_h
    //    for (Int_t i=49;i<nBackground;i++)  nBS[6][j] += nB[i][j];   // 6f
    for (Int_t i=1;i<nSSum;i++) nSG[j]  += nSS[i][j];
    for (Int_t i=1;i<nSSum;i++) nSG0[j] += nSS0[i][j];
    for (Int_t i=1;i<nBSum;i++) nBG[j]  += nBS[i][j];
    nSig[j] = nSG[j]/TMath::Sqrt(nSG[j]+nBG[j]);
    for (Int_t i=1;i<nBackground;i++) nErrorBG[j] += nB[i][j]*weightB[i];
  }

  
  cerr << "----------Cuts----------" << endl;
  for (Int_t i=0;i<nCutMax;i++) {
    cerr << "Cut " << setw(2) << i << endl;
    myCut[i].Print();
  }
#if 1
  cerr << "Cuts :  " << setw(15) << "Backgrounds" << setw(15) << "Error_B" << setw(15) << "Signal" << setw(15) << "Significance" << endl;
  for (Int_t i=0;i<nCutMax;i++) {
    cerr << "Cut" << i << " :  " << setw(15) << nBG[i] << setw(15) << TMath::Sqrt(nErrorBG[i]) << setw(15) << nSG[i]  << setw(15) << nSig[i] << endl;
  }
#endif
  // reduction table
  if (iPol == 0) {
    cerr << "Polarization: (e-,e+) = (0,0)" << endl;
  }
  else if(iPol == -1) {
    cerr << "Polarization: (e-,e+) = (-0.8,+0.3)" << endl;
  }
  else if(iPol == +1) {
    cerr << "Polarization: (e-,e+) = (+0.8,-0.3)" << endl;
  }
  cerr << "---------------------------------------------------------------------------Reduction Table-------------------------------------------------------------------------------------" << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cerr << "Process          :  " << setw(12) << "2f_l" << setw(12) << "2f_h" << setw(12) << "4f_l" << setw(12) << "4f_sl" << setw(12) << "4f_h" << setw(12) << "BG" << setw(12) << "llh" << setw(12) << "Signal"  << setw(12) << "Signf" << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cerr << "Cross Section    :  ";
  for (Int_t i=1;i<nBSum;i++) {
    cerr << setw(12) << 1./inteLumi*nexpBS[i];
  } 
  cerr << setw(12) << 1./inteLumi*nexpBG << setw(12) <<  1./inteLumi*nexpSG << setw(12) <<  1./inteLumi*nexpSG*Brbb << setw(12) << " " << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cerr << "Generated        :  ";
  for (Int_t i=1;i<nBSum;i++) {
    cerr << setw(12) << ngenBS[i];
  }
  cerr << setw(12) << " " << setw(12) << ngenSG << setw(12) << ngenSG << setw(12) << " " << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cerr << "Expected         :  ";
  for (Int_t i=1;i<nBSum;i++) {
    cerr << setw(12) << nexpBS[i];
  }
  cerr << setw(12) << nexpBG << setw(12) << nexpSG << setw(12) << nexpSG*Brbb << setw(12) << nexpSig << endl;
  cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  for (Int_t i=0;i<nCutMax;i++) {
    cerr << "Cut" << i << "             :  ";
    for (Int_t j=1;j<nBSum;j++) {
      cerr << setw(12) << nBS[j][i];
    }
    cerr << setw(12) << nBG[i] << setw(12) << nSG[i] << setw(12) << nSG0[i] << setw(12) << nSig[i] << endl;
    cerr << "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  }

  // close root files
  for (Int_t i=1;i<nSignal;i++) {
    fS[i]->Close();
  }
  for (Int_t i=1;i<nBackground;i++) {
    fB[i]->Close();
  }

  return;

}

