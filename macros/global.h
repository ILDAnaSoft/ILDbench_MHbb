#ifndef __GLOBAL__
#define __GLOBAL__

//#define __electron__
//#define __IDRL__

#include "TROOT.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TString.h"


const Double_t inteLumi = 500.;

const Double_t kMass_H  = 125.;
//const Double_t BRinv = 0.1;
//const Double_t BRbb = 0.578;
const Double_t Brbb = 0.578, Brcc = 0.0268, Brtautau = 0.0637, Brmumu = 0.000221;  
const Double_t Brgg = 0.0856, Brgamgam = 0.0023, BrWW = 0.216, BrZZ = 0.0267;

#ifndef __electron__
const Int_t nSignal = 2 + 1;
#else
const Int_t nSignal = 4 + 1;
#endif
//const Int_t nBackground = 52 + 1;
const Int_t nBackground = 48 + 1;
//const Int_t nBackground = 48 + 1 + 188;
//const Int_t nBackground = 52 + 1 + 186;

const Int_t nCutMax = 10 + 1;
const Int_t nCutUsedMax = 8 + 1;

TCut myCut[nCutMax];
TFile *fS[nSignal], *fB[nBackground];
TNtuple *hAnl_S[nSignal],*hAnl_B[nBackground];
TNtuple *hGen_S[nSignal],*hGen_B[nBackground];

Double_t xsecS[nSignal],xsecB[nBackground];
Double_t nexpS[nSignal],nexpB[nBackground];
Double_t ngenS[nSignal],ngenB[nBackground];
Double_t weightS[nSignal],weightB[nBackground];


std::vector<TString> pname6f,pol6f;
std::vector<Double_t> xsec6f;


#endif
