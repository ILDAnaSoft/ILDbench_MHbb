// *****************************************************
// e+e- ------> ZH ------> (l+ l-) + X
// Processor for final selection
//                        ----Junping
// *****************************************************
#include "ZHll2JAnalysisProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>
#include "UTIL/PIDHandler.h"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "Utilities.h"
#include "JetAlgorithm.h"

//#define __ISR_subtract

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace mylib;
using namespace myjet;


ZHll2JAnalysisProcessor aZHll2JAnalysisProcessor ;


ZHll2JAnalysisProcessor::ZHll2JAnalysisProcessor() : Processor("ZHll2JAnalysisProcessor") {
  
  // modify processor description
  _description = "ZHll2JAnalysisProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticlesSkimmed") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "InputMCTruthLinkCollection" , 
			   "Name of the MCTruthLink collection"  ,
			   _colMCTL ,
			   std::string("RecoMCTruthLink") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPandoraPFOsCollection" , 
			   "Name of the PandoraPFOs collection"  ,
			   _colPFOs ,
			   std::string("PandoraPFOs") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputNewPFOsCollection",
			   "Name of the new PFOs collection after some pre-cuts",
			   _colNewPFOs,
			   std::string("NewPFOs_Uncluster") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputLeptonsCollection",
			   "Name of collection with the selected leptons",
			   _colLeptons,
			   std::string("Leptons") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputLeptonsOrigCollection",
			   "Name of collection with the selected leptons w/o FSR and BS recovery",
			   _colLeptonsOrig,
			   std::string("LeptonsOrig") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputJetsCollection",
			   "Name of the jets collection",
			   _colJets,
			   std::string("Durham_2Jets") );

  registerProcessorParameter("CenterOfMassEnergy",
			   "Center of mass energy"  ,
			   _ecm ,
			   float(500) ) ;
}

void ZHll2JAnalysisProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStatAnl = 0;
  
}

void ZHll2JAnalysisProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void ZHll2JAnalysisProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

  // some constants
  static const Double_t fEnergyCut = 0.05; // minimum energy of each pfo
  static const Double_t kMassW   =  80.4; // W mass by s.y
  static const Double_t kSigmaMw =  4.8; // W mass resolution
  static const Double_t kMassZ   = 91.187; // Z mass
  static const Double_t kSigmaMz =   6.0; // Z mass resolution
  static const Double_t kMassH   = 125.0; // H mass
  static const Double_t kSigmaMh =   7.2; // H mass resolution
  static const Double_t kMassT   = 174.0; // Top mass
  static const Double_t kSigmaMt =  20.0; // Top mass resolution
  static const Double_t kSigmaE  =   7.0; // Top mass resolution
  static const Double_t fMassZCut = 40.;     // mass cut for lepton pair from Z
  static const Double_t fMassHCut = 80.;     // mass cut for lepton pair from Z

  if (!hStatAnl) hStatAnl = new TH1D("hStatAnl", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  cerr << endl << "Hello, Analysis!" << " No: " << _nEvt << endl;

  static TNtupleD *hAnl = 0;
  if (!hAnl) {
    stringstream tupstr;
    tupstr << "nhbb:nhww:nhgg:nhtt:nhcc:nhzz:nhaa:nhmm" << ":"
	   << "elep1mc:elep2mc:coslep1mc:coslep2mc:ptzmc:coszmc:mzmc:acop12mc:mrecoilmc" << ":"
	   << "eisr1mc:eisr2mc" << ":"
	   << "elep1:elep2:coslep1:coslep2:ptz:cosz:mz:acop12:mrecoil" << ":"
	   << "elep1orig:elep2orig:coslep1orig:coslep2orig:ptzorig:coszorig:mzorig:acop12orig:mrecoilorig" << ":"
	   << "mvalep1:mvalep2:leptype:nphoton" << ":"
	   << "deltapt:deltaptorig:coslep12:coslep12orig:evis:cosmis:ephotonmax:econephoton" << ":"
	   << "mzcorr:mrecoilcorr:coszphoton:ifsr:iisr:ptvis:evischg:ptvischg:mmis" << ":"
	   << "pxlep1mc:pylep1mc:pzlep1mc:pxlep2mc:pylep2mc:pzlep2mc" << ":"
	   << "pxj1mc:pyj1mc:pzj1mc:ej1mc:pxj2mc:pyj2mc:pzj2mc:ej2mc:cosj1mc:phij1mc:cosj2mc:phij2mc" << ":"
	   << "pxnewj1mc:pynewj1mc:pznewj1mc:enewj1mc:pxnewj2mc:pynewj2mc:pznewj2mc:enewj2mc" << ":"
	   << "mhnewmc" << ":"
	   << "pxlep1:pylep1:pzlep1:pxlep2:pylep2:pzlep2" << ":"
	   << "pxj1:pyj1:pzj1:ej1:pxj2:pyj2:pzj2:ej2:cosj1:phij1:cosj2:phij2" << ":"
	   << "pxnewj1:pynewj1:pznewj1:enewj1:pxnewj2:pynewj2:pznewj2:enewj2" << ":"
	   << "mhnew"  << ":"
	   << "cosisr1mc:cosisr2mc:mh" << ":"
	   << "bmax1:bmax2:yminus:yplus:yminus4:yplus4:npfosc1:npfosc2" << ":"
	   << "mrecoilmc:cosj1h:phij1h:cosj2h:phij2h:mhcol:ej1h:ej2h" << ":"
	   << "mhjet:mhmuon:mhjet2:mj1h:mj2h:mj1:mj2:mhjet3:mhjet4:mhjet5" 
	   << ends;
    hAnl = new TNtupleD("hAnl","",tupstr.str().data());
  }

  // ------------------------------------------------
  // -- read out the MCParticles information
  // ------------------------------------------------
  LCCollection *colMC = evt->getCollection(_colMCP);
  // get the truth information
  Int_t nMCP = colMC->getNumberOfElements();
  TLorentzVector lortzLep1MC,lortzLep2MC,lortzZMC,lortzHMC,lortzISRMC,lortzHColMC;
  TLorentzVector lortzISR1MC, lortzISR2MC;
  TLorentzVector lortzJ1MC, lortzJ2MC;  
  LCCollectionVec *colMCH = new LCCollectionVec(LCIO::MCPARTICLE);
  std::vector<JJet> vjetsMCH;
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    Int_t pdg = mcPart->getPDG();
    Int_t nparents = mcPart->getParents().size();
    Int_t motherpdg = 0;
    if (nparents > 0) {
      MCParticle *mother = mcPart->getParents()[0];
      motherpdg = mother->getPDG();
    }
    Double_t energy = mcPart->getEnergy();
    TVector3 pv = TVector3(mcPart->getMomentum());
    Int_t ndaughters = mcPart->getDaughters().size();
    Int_t daughterpdg = 0;
    if (ndaughters > 0) {
      MCParticle *daughter = mcPart->getDaughters()[0];
      daughterpdg = daughter->getPDG();
    }
    TLorentzVector lortz = TLorentzVector(pv,energy);
    Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
    if ((pdg == 13 || pdg == 11) && motherpdg == 0 && ioverlay == 0) {
      lortzLep1MC = lortz;
      lortzZMC += lortz;
    }
    if ((pdg == -13 || pdg == -11) && motherpdg == 0 && ioverlay == 0) {
      lortzLep2MC = lortz;
      lortzZMC += lortz;
    }
    if (pdg == 25 && motherpdg == 0 && ioverlay == 0) {
      lortzHMC = lortz;
    }
    if (pdg > 0 && motherpdg == 25 && ioverlay == 0) {
      lortzJ1MC = lortz;
    }
    if (pdg < 0 && motherpdg == 25 && ioverlay == 0) {
      lortzJ2MC = lortz;
    }
    if (i == 0) {
      lortzISR1MC = lortz;
    }
    if (i == 1) {
      lortzISR2MC = lortz;
    }
    if (mcPart->getGeneratorStatus() == 1 && ioverlay == 0) {
      //      cerr << i << ": Origin = " << getOriginalPDG(mcPart,true) << endl;
#if 1
      if (getOriginalPDG(mcPart,true) == 25) {
	colMCH->addElement(mcPart);
	lortzHColMC += lortz;
	JJet jh(lortz);
	vjetsMCH.push_back(jh);
      }
#endif      
    }
  }
  // get Higgs decay modes
  std::vector<Int_t> nHDecay;
  nHDecay = getHiggsDecayModes(colMC);
  Double_t nHbb = nHDecay[0]; // H---> b b
  Double_t nHWW = nHDecay[1]; // H---> W W
  Double_t nHgg = nHDecay[2]; // H---> g g
  Double_t nHtt = nHDecay[3]; // H---> tau tau
  Double_t nHcc = nHDecay[4]; // H---> c c
  Double_t nHZZ = nHDecay[5]; // H---> Z Z
  Double_t nHaa = nHDecay[6]; // H---> gam gam
  Double_t nHmm = nHDecay[7]; // H---> mu mu

  TLorentzVector lortzEcm = getLorentzEcm(_ecm);
  TLorentzVector lortzRecoilMC = lortzEcm - lortzZMC;

  Double_t thetaJ1MC = lortzJ1MC.Theta();
  Double_t thetaJ2MC = lortzJ2MC.Theta();
  Double_t phiJ1MC = lortzJ1MC.Phi();  
  Double_t phiJ2MC = lortzJ2MC.Phi();  
  TVector2 pJ12NewMC = getMomentumNew(lortzRecoilMC.Px(),lortzRecoilMC.Py(),
				 thetaJ1MC,thetaJ2MC,phiJ1MC,phiJ2MC);
  Double_t pJ1NewMC = pJ12NewMC.X();
  TVector3 momentumJ1NewMC = TVector3(pJ1NewMC*TMath::Sin(thetaJ1MC)*TMath::Cos(phiJ1MC),
				      pJ1NewMC*TMath::Sin(thetaJ1MC)*TMath::Sin(phiJ1MC),
				      pJ1NewMC*TMath::Cos(thetaJ1MC));
  Double_t pJ2NewMC = pJ12NewMC.Y();
  TVector3 momentumJ2NewMC = TVector3(pJ2NewMC*TMath::Sin(thetaJ2MC)*TMath::Cos(phiJ2MC),
				      pJ2NewMC*TMath::Sin(thetaJ2MC)*TMath::Sin(phiJ2MC),
				      pJ2NewMC*TMath::Cos(thetaJ2MC));
  Double_t mJ1MC = lortzJ1MC.M();
  Double_t mJ2MC = lortzJ2MC.M();
  Double_t eJ1NewMC = TMath::Sqrt(momentumJ1NewMC.Mag()*momentumJ1NewMC.Mag()+
				  mJ1MC*mJ1MC);
  Double_t eJ2NewMC = TMath::Sqrt(momentumJ2NewMC.Mag()*momentumJ2NewMC.Mag()+
				  mJ2MC*mJ2MC);
  TLorentzVector lortzJ1NewMC = TLorentzVector(momentumJ1NewMC,eJ1NewMC);
  TLorentzVector lortzJ2NewMC = TLorentzVector(momentumJ2NewMC,eJ2NewMC);  
  TLorentzVector lortzHNewMC = lortzJ1NewMC + lortzJ2NewMC;

  // ------------------------------------------------
  // -- read out the PandoraPFOs information
  // ------------------------------------------------
  LCCollection *colPFO = evt->getCollection(_colPFOs);

  // ------------------------------------------------
  // -- read out the Thrust information
  // ------------------------------------------------
  LCCollection *colSelRecPart = evt->getCollection("SelectedReconstructedParticle");
  Double_t principleThrust = colSelRecPart->parameters().getFloatVal("principleThrustValue");
  Double_t majorThrust = colSelRecPart->parameters().getFloatVal("majorThrustValue");
  Double_t minorThrust = colSelRecPart->parameters().getFloatVal("minorThrustValue");
  FloatVec tAxis;
  FloatVec thrustAxis = colSelRecPart->parameters().getFloatVals("principleThrustAxis",tAxis);
  TVector3 principleAxis = TVector3(thrustAxis[0],thrustAxis[1],thrustAxis[2]);
  Double_t cosThrustAxis = principleAxis.CosTheta();

  // ------------------------------------------------
  // -- read out the NewPFOs information
  // ------------------------------------------------
  LCCollection *colNewPFO = evt->getCollection(_colNewPFOs);
  if (!colNewPFO) {
    cerr << "No NewPFOs Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "NewPFOs Collection Found" << ends;
  // get the visible energy
  TLorentzVector lortzVis(0.,0.,0.,0.);
  TLorentzVector lortzVisChg(0.,0.,0.,0.);
  Int_t nPFOs = colNewPFO->getNumberOfElements();
  Int_t nParticles = 0;
  ReconstructedParticle *photonMax = 0;
  Double_t energyPhotonMax = -1.;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colNewPFO->getElementAt(i));
    lortzVis += getLorentzVector(pfo);
    if (pfo->getEnergy() >= fEnergyCut) nParticles++;
    // find the photon with largest energy
    if (TMath::Abs(pfo->getCharge()) < 0.5) {
      Int_t leptonID = getLeptonID(pfo);
      if (leptonID == 22 && pfo->getEnergy() > energyPhotonMax) {
	energyPhotonMax = pfo->getEnergy();
	photonMax = pfo;
      }
    }
    else {
      lortzVisChg += getLorentzVector(pfo);
    }
  }

  // ------------------------------------------------
  // -- read out the MCTruthLink information
  // ------------------------------------------------
  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);
  
  // ------------------------------------------------
  // -- read out the Leptons information
  // ------------------------------------------------
  LCCollection *colLep = evt->getCollection(_colLeptons);
  Int_t nLeps = colLep->getNumberOfElements();
  if (!colLep || nLeps != 2) {
    cerr << "No Leptons Collection Found or nLeptons != 2 !" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "Leptons Collection Found" << ends;
  ReconstructedParticle *leptonMinus = dynamic_cast<ReconstructedParticle*>(colLep->getElementAt(0));
  ReconstructedParticle *leptonPlus  = dynamic_cast<ReconstructedParticle*>(colLep->getElementAt(1));
  ReconstructedParticle *lep_tmp = 0;
  if (leptonMinus->getCharge() > leptonPlus->getCharge()) {
    lep_tmp = leptonMinus;
    leptonMinus = leptonPlus;
    leptonPlus  = lep_tmp;
  }
  Double_t _mva_lep_minus = colLep->getParameters().getFloatVal("MVALepMinus");
  Double_t _mva_lep_plus  = colLep->getParameters().getFloatVal("MVALepPlus");
  Double_t _lep_type = colLep->getParameters().getIntVal("ISOLepType");
  Double_t nPhotonsRecovered = colLep->getParameters().getIntVal("NPhotons");

  // ------------------------------------------------
  // -- read out the LeptonsOrig information
  // ------------------------------------------------
  LCCollection *colLepOrig = evt->getCollection(_colLeptonsOrig);
  Int_t nLepsOrig = colLepOrig->getNumberOfElements();
  if (!colLepOrig || nLepsOrig != 2) {
    cerr << "No LeptonsOrig Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "LeptonsOrig Collection Found" << ends;
  ReconstructedParticle *leptonMinusOrig = dynamic_cast<ReconstructedParticle*>(colLepOrig->getElementAt(0));
  ReconstructedParticle *leptonPlusOrig  = dynamic_cast<ReconstructedParticle*>(colLepOrig->getElementAt(1));
  if (leptonMinusOrig->getCharge() > leptonPlusOrig->getCharge()) {
    lep_tmp = leptonMinusOrig;
    leptonMinusOrig = leptonPlusOrig;
    leptonPlusOrig  = lep_tmp;
  }

  // ------------------------------------------------
  // -- read out the Refined_2Jet information
  // ------------------------------------------------
  LCCollection *colJet = evt->getCollection(_colJets);
  if (!colJet) {
    cerr << "No Durham_2Jets Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "Durham_2Jets Collection Found" << ends;
  Int_t nJets = colJet->getNumberOfElements();
  if (nJets != 2) {
    cerr << "Number of Jets is not 2" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "2 Jets Found" << ends;

  ReconstructedParticle *jets[2];
  // flavor tagging information
  PIDHandler pidh (colJet);
  Int_t algo = pidh.getAlgorithmID("lcfiplus");
  Double_t FLV[2][11];
  Int_t nPFOsCJ1 = 0, nPFOsCJ2 = 0;
  for (Int_t i=0;i<nJets;i++) {
    jets[i] = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
    const ParticleID & jetID = pidh.getParticleID(jets[i], algo);
    FloatVec params = jetID.getParameters();
    FLV[i][0] = params[pidh.getParameterIndex(algo, "BTag")];
    FLV[i][1] = params[pidh.getParameterIndex(algo, "CTag")];
    FLV[i][2] = params[pidh.getParameterIndex(algo, "BCTag")];
    std::vector<lcio::ReconstructedParticle*> partVec = jets[i]->getParticles();
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPart=partVec.begin();iPart!=partVec.end();++iPart) {
      TVector3 momPart = TVector3((*iPart)->getMomentum());
      Double_t pTPart = momPart.Pt();
      if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 0) nPFOsCJ1++;
      if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 1) nPFOsCJ2++;
    }
  }
  Int_t algo_y = pidh.getAlgorithmID("yth");
  const ParticleID & ythID = pidh.getParticleID(jets[0], algo_y);
  FloatVec params_y = ythID.getParameters();
  Double_t yMinus = params_y[pidh.getParameterIndex(algo_y, "y12")];
  Double_t yPlus  = params_y[pidh.getParameterIndex(algo_y, "y23")];
  Double_t yMinus4= params_y[pidh.getParameterIndex(algo_y, "y34")];
  Double_t yPlus4 = params_y[pidh.getParameterIndex(algo_y, "y45")];
  Double_t bmax1 = FLV[0][0] > FLV[1][0] ? FLV[0][0] : FLV[1][0];
  Double_t bmax2 = FLV[0][0] > FLV[1][0] ? FLV[1][0] : FLV[0][0];

  // ------------------------------------------------
  // -- get the useful physical quantities and save them to ntuple
  // ------------------------------------------------
  TLorentzVector lortzLep1 = getLorentzVector(leptonMinus);
  TLorentzVector lortzLep2 = getLorentzVector(leptonPlus);
  TLorentzVector lortzLepOrig1 = getLorentzVector(leptonMinusOrig);
  TLorentzVector lortzLepOrig2 = getLorentzVector(leptonPlusOrig);
  TLorentzVector lortzZ = lortzLep1 + lortzLep2;
  TLorentzVector lortzZOrig = lortzLepOrig1 + lortzLepOrig2;
  // get the Pt balance between Z and Photon
  Double_t deltaPtbal = lortzZ.Pt();
  Double_t deltaPtbalOrig = lortzZOrig.Pt();
  Double_t coneEnergyPhoton = -1.;
  TLorentzVector lortzPhotonMax(0.,0.,0.,0.);
  Double_t cosThetaZPhoton = -1.;
  if (energyPhotonMax > 10) {
    lortzPhotonMax = getLorentzVector(photonMax);
    deltaPtbal -= lortzPhotonMax.Pt();
    deltaPtbalOrig -= lortzPhotonMax.Pt();
    coneEnergyPhoton = getConeEnergy(photonMax,colNewPFO,0.98);
    cosThetaZPhoton = getCosTheta(lortzZ,lortzPhotonMax);
  }
  TLorentzVector lortzMis = lortzEcm - lortzVis - lortzZ;

  // recover the remained FSR/ISR if a energetic and isolated photon is found
  TLorentzVector lortzZCorr = lortzZ + lortzPhotonMax;
  Double_t massZCorr = lortzZ.M();
  Double_t recoilMassCorr = getRecoilMass(lortzEcm,lortzZ);
  Int_t iFSRCorr = 0, iISRCorr = 0;
  if (energyPhotonMax > 10 && coneEnergyPhoton/energyPhotonMax < 0.1 && cosThetaZPhoton > 0.) {
    if (TMath::Abs(lortzZCorr.M()-kMassZ) < TMath::Abs(lortzZ.M()-kMassZ)) {
      iFSRCorr = 1;
      massZCorr = lortzZCorr.M();
    }
  }
  if (energyPhotonMax > 10 && coneEnergyPhoton/energyPhotonMax < 0.1 ) {
    Double_t recoilMass = getRecoilMass(lortzEcm,lortzZ);
    Double_t recoilMassNew = getRecoilMass(lortzEcm,lortzZCorr);
    if (TMath::Abs(recoilMassNew-kMassH) < TMath::Abs(recoilMass-kMassH)) {
      iISRCorr = 1;
      recoilMassCorr = recoilMassNew;
    }
  }
  // jets information
  ReconstructedParticle *jet1 = jets[0];
  ReconstructedParticle *jet2 = jets[1];
  TVector3 momentum_j1 = TVector3(jet1->getMomentum());
  TVector3 momentum_j2 = TVector3(jet2->getMomentum());
  Double_t e_j1  = jet1->getEnergy();
  Double_t e_j2  = jet2->getEnergy();
  TLorentzVector lortzJ1 = TLorentzVector(momentum_j1,e_j1);
  TLorentzVector lortzJ2 = TLorentzVector(momentum_j2,e_j2);
  TLorentzVector lortzH = lortzJ1 + lortzJ2;
  // match
  Double_t cosj1_1 = getCosTheta(lortzJ1,lortzJ1MC);
  Double_t cosj1_2 = getCosTheta(lortzJ1,lortzJ2MC);
  if (cosj1_1 < cosj1_2) {
    TLorentzVector lortz_tmp = lortzJ1;
    lortzJ1 = lortzJ2;
    lortzJ2 = lortz_tmp;
  }
  // apply new method
  TLorentzVector lortzRecoil = lortzEcm - lortzZ;
#ifdef __ISR_subtract
  LCCollection *colPhoton = evt->getCollection("photons");
  if (colPhoton) {
    for (Int_t i=0;i<colPhoton->getNumberOfElements();i++) {
      ReconstructedParticle *ph = dynamic_cast<ReconstructedParticle*>(colPhoton->getElementAt(i));
      lortzRecoil -= getLorentzVector(ph);
    }
  }
#endif  
  Double_t thetaJ1 = lortzJ1.Theta();
  Double_t thetaJ2 = lortzJ2.Theta();
  Double_t phiJ1 = lortzJ1.Phi();  
  Double_t phiJ2 = lortzJ2.Phi();  
  TVector2 pJ12New = getMomentumNew(lortzRecoil.Px(),lortzRecoil.Py(),
				 thetaJ1,thetaJ2,phiJ1,phiJ2);
  Double_t pJ1New = pJ12New.X();
  TVector3 momentumJ1New = TVector3(pJ1New*TMath::Sin(thetaJ1)*TMath::Cos(phiJ1),
				      pJ1New*TMath::Sin(thetaJ1)*TMath::Sin(phiJ1),
				      pJ1New*TMath::Cos(thetaJ1));
  Double_t pJ2New = pJ12New.Y();
  TVector3 momentumJ2New = TVector3(pJ2New*TMath::Sin(thetaJ2)*TMath::Cos(phiJ2),
				      pJ2New*TMath::Sin(thetaJ2)*TMath::Sin(phiJ2),
				      pJ2New*TMath::Cos(thetaJ2));
  Double_t mJ1 = lortzJ1.M();
  Double_t mJ2 = lortzJ2.M();
  Double_t eJ1New = TMath::Sqrt(momentumJ1New.Mag()*momentumJ1New.Mag()+
				  mJ1*mJ1);
  Double_t eJ2New = TMath::Sqrt(momentumJ2New.Mag()*momentumJ2New.Mag()+
				  mJ2*mJ2);
  TLorentzVector lortzJ1New = TLorentzVector(momentumJ1New,eJ1New);
  TLorentzVector lortzJ2New = TLorentzVector(momentumJ2New,eJ2New);  
  TLorentzVector lortzHNew = lortzJ1New + lortzJ2New;

  // for cheating study
  // get two jets from MC Particles from Higgs decay
  //  JJets jetsMCH(colMCH);
  JJets jetsMCH(vjetsMCH,1);  
  //  jetsMCH.SetAlgorithm(1); // 1=Durham fixedN
  jetsMCH.DoClustering(2);
    //    double yMinus = jets.GetYMinus(); // current smallest Yij
    //    double yPlus  = jets.GetYPlus();  // last smallest Yij
  LCCollection *colJetMCH = jetsMCH.GetJetsCol(); // get output jets collection
  ReconstructedParticle *jetsh[2];
  for (Int_t i=0;i<2;i++) {
    jetsh[i] = dynamic_cast<ReconstructedParticle*>(colJetMCH->getElementAt(i));
  }
  // jets information
  ReconstructedParticle *jet1h = jetsh[0];
  ReconstructedParticle *jet2h = jetsh[1];
  TVector3 momentum_j1h = TVector3(jet1h->getMomentum());
  TVector3 momentum_j2h = TVector3(jet2h->getMomentum());
  Double_t e_j1h  = jet1h->getEnergy();
  Double_t e_j2h  = jet2h->getEnergy();
  TLorentzVector lortzJ1h = TLorentzVector(momentum_j1h,e_j1h);
  TLorentzVector lortzJ2h = TLorentzVector(momentum_j2h,e_j2h);
  TLorentzVector lortzHh = lortzJ1h + lortzJ2h;
  // match
  Double_t cosj1h_1 = getCosTheta(lortzJ1h,lortzJ1MC);
  Double_t cosj1h_2 = getCosTheta(lortzJ1h,lortzJ2MC);
  if (cosj1h_1 < cosj1h_2) {
    TLorentzVector lortz_tmp = lortzJ1h;
    lortzJ1h = lortzJ2h;
    lortzJ2h = lortz_tmp;
  }

  // cheat jets direction
  Double_t mH_jet  = getHMass(lortzJ1h,lortzJ2h,lortzRecoil);   // cheat jet mass & jet direction
  Double_t mH_muon = getHMass(lortzJ1,lortzJ2,lortzRecoilMC);   // cheat muon 4-momentum
  Double_t mH_jet2 = getHMass(lortzJ1h,lortzJ2h,lortzRecoil,lortzJ1,lortzJ2); // cheat jet direction
  Double_t mH_jet3 = getHMass(lortzJ1,lortzJ2,lortzRecoil,lortzJ1h,lortzJ2h); // cheat jet mass
  Double_t mH_jet4 = getHMass(lortzJ1h,lortzJ2h,lortzRecoilMC,lortzJ1,lortzJ2); // cheat jet direction & muon 4-momentum
  Double_t mH_jet5 = getHMass(lortzJ1,lortzJ2,lortzRecoilMC,lortzJ1h,lortzJ2h); // cheat jet mass & muon 4-momentum
  

  Double_t data[200];
  data[ 0] = nHbb;
  data[ 1] = nHWW;
  data[ 2] = nHgg;
  data[ 3] = nHtt;
  data[ 4] = nHcc;
  data[ 5] = nHZZ;
  data[ 6] = nHaa;
  data[ 7] = nHmm;
  data[ 8] = lortzLep1MC.E();
  data[ 9] = lortzLep2MC.E();
  data[10] = lortzLep1MC.CosTheta();
  data[11] = lortzLep2MC.CosTheta();
  data[12] = lortzZMC.Pt();
  data[13] = lortzZMC.CosTheta();
  data[14] = lortzZMC.M();
  data[15] = getAcoPlanarity(lortzLep1MC,lortzLep2MC);
  data[16] = getRecoilMass(lortzEcm,lortzZMC);
  data[17] = lortzISR1MC.E();
  data[18] = lortzISR2MC.E();
  data[19] = lortzLep1.E();
  data[20] = lortzLep2.E();
  data[21] = lortzLep1.CosTheta();
  data[22] = lortzLep2.CosTheta();
  data[23] = lortzZ.Pt();
  data[24] = lortzZ.CosTheta();
  data[25] = lortzZ.M();
  data[26] = getAcoPlanarity(lortzLep1,lortzLep2);
  data[27] = getRecoilMass(lortzEcm,lortzZ);
  data[28] = lortzLepOrig1.E();
  data[29] = lortzLepOrig2.E();
  data[30] = lortzLepOrig1.CosTheta();
  data[31] = lortzLepOrig2.CosTheta();
  data[32] = lortzZOrig.Pt();
  data[33] = lortzZOrig.CosTheta();
  data[34] = lortzZOrig.M();
  data[35] = getAcoPlanarity(lortzLepOrig1,lortzLepOrig2);
  data[36] = getRecoilMass(lortzEcm,lortzZOrig);
  data[37] = _mva_lep_minus;
  data[38] = _mva_lep_plus;
  data[39] = _lep_type;
  data[40] = nPhotonsRecovered;
  data[41] = deltaPtbal;
  data[42] = deltaPtbalOrig;
  data[43] = getCosTheta(leptonMinus,leptonPlus);
  data[44] = getCosTheta(leptonMinusOrig,leptonPlusOrig);
  data[45] = lortzVis.E();
  data[46] = lortzMis.CosTheta();
  data[47] = energyPhotonMax;
  data[48] = coneEnergyPhoton;
  data[49] = massZCorr;
  data[50] = recoilMassCorr;
  data[51] = cosThetaZPhoton;
  data[52] = iFSRCorr;
  data[53] = iISRCorr;
  data[54] = lortzVis.Pt();
  data[55] = lortzVisChg.E();
  data[56] = lortzVisChg.Pt();
  data[57] = (lortzEcm-lortzZ-lortzVis).M();
  data[ 58]= lortzLep1MC.Px();  
  data[ 59]= lortzLep1MC.Py();  
  data[ 60]= lortzLep1MC.Pz();  
  data[ 61]= lortzLep2MC.Px();  
  data[ 62]= lortzLep2MC.Py();  
  data[ 63]= lortzLep2MC.Pz();  
  data[ 64]= lortzJ1MC.Px();  
  data[ 65]= lortzJ1MC.Py();  
  data[ 66]= lortzJ1MC.Pz();  
  data[ 67]= lortzJ1MC.E();  
  data[ 68]= lortzJ2MC.Px();  
  data[ 69]= lortzJ2MC.Py();  
  data[ 70]= lortzJ2MC.Pz();  
  data[ 71]= lortzJ2MC.E();  
  data[ 72]= lortzJ1MC.CosTheta();  
  data[ 73]= lortzJ1MC.Phi();  
  data[ 74]= lortzJ2MC.CosTheta();  
  data[ 75]= lortzJ2MC.Phi();  
  data[ 76]= lortzJ1NewMC.Px();  
  data[ 77]= lortzJ1NewMC.Py();  
  data[ 78]= lortzJ1NewMC.Pz();  
  data[ 79]= lortzJ1NewMC.E();  
  data[ 80]= lortzJ2NewMC.Px();  
  data[ 81]= lortzJ2NewMC.Py();  
  data[ 82]= lortzJ2NewMC.Pz();  
  data[ 83]= lortzJ2NewMC.E();
  data[ 84]= lortzHNewMC.M();    
  data[ 85]= lortzLep1.Px();  
  data[ 86]= lortzLep1.Py();  
  data[ 87]= lortzLep1.Pz();  
  data[ 88]= lortzLep2.Px();  
  data[ 89]= lortzLep2.Py();  
  data[ 90]= lortzLep2.Pz();  
  data[ 91]= lortzJ1.Px();  
  data[ 92]= lortzJ1.Py();  
  data[ 93]= lortzJ1.Pz();  
  data[ 94]= lortzJ1.E();  
  data[ 95]= lortzJ2.Px();  
  data[ 96]= lortzJ2.Py();  
  data[ 97]= lortzJ2.Pz();  
  data[ 98]= lortzJ2.E();  
  data[ 99]= lortzJ1.CosTheta();  
  data[100]= lortzJ1.Phi();  
  data[101]= lortzJ2.CosTheta();  
  data[102]= lortzJ2.Phi();  
  data[103]= lortzJ1New.Px();  
  data[104]= lortzJ1New.Py();  
  data[105]= lortzJ1New.Pz();  
  data[106]= lortzJ1New.E();  
  data[107]= lortzJ2New.Px();  
  data[108]= lortzJ2New.Py();  
  data[109]= lortzJ2New.Pz();  
  data[110]= lortzJ2New.E();
  data[111]= lortzHNew.M();    
  data[112]= lortzISR1MC.CosTheta();
  data[113]= lortzISR2MC.CosTheta();
  data[114]= lortzH.M();
  data[115]= bmax1;
  data[116]= bmax2;
  data[117]= yMinus;
  data[118]= yPlus;
  data[119]= yMinus4;
  data[120]= yPlus4;
  data[121]= nPFOsCJ1;
  data[122]= nPFOsCJ2;
  data[123]= lortzRecoilMC.M();
  data[124]= lortzJ1h.CosTheta();  
  data[125]= lortzJ1h.Phi();  
  data[126]= lortzJ2h.CosTheta();  
  data[127]= lortzJ2h.Phi();  
  data[128]= lortzHColMC.M();  
  data[129]= lortzJ1h.E();  
  data[130]= lortzJ2h.E();  
  data[131]= mH_jet;
  data[132]= mH_muon;
  data[133]= mH_jet2;
  data[134]= lortzJ1h.M();  
  data[135]= lortzJ2h.M();  
  data[136]= lortzJ1.M();  
  data[137]= lortzJ2.M();  
  data[138]= mH_jet3;
  data[139]= mH_jet4;
  data[140]= mH_jet5;
  hAnl->Fill(data);

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  //  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  //		       << "   in run:  " << evt->getRunNumber() 
  //		       << std::endl ;

  //  _nEvt ++ ;

  last->cd();
}



void ZHll2JAnalysisProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHll2JAnalysisProcessor::end(){ 

  cerr << "ZHll2JAnalysisProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  //  cerr << endl;
  cerr << "  =============" << endl;
  cerr << "   Cut Summary " << endl;
  cerr << "  =============" << endl;
  cerr << "   ll+4 Jet    " << endl;
  cerr << "  =============" << endl;
  cerr << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStatAnl->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;
  
}

Double_t ZHll2JAnalysisProcessor::getHMass(TLorentzVector lortzJ1, TLorentzVector lortzJ2, TLorentzVector lortzRecoil) {
  Double_t thetaJ1 = lortzJ1.Theta();
  Double_t thetaJ2 = lortzJ2.Theta();
  Double_t phiJ1 = lortzJ1.Phi();  
  Double_t phiJ2 = lortzJ2.Phi();  
  TVector2 pJ12New = getMomentumNew(lortzRecoil.Px(),lortzRecoil.Py(),
				 thetaJ1,thetaJ2,phiJ1,phiJ2);
  Double_t pJ1New = pJ12New.X();
  TVector3 momentumJ1New = TVector3(pJ1New*TMath::Sin(thetaJ1)*TMath::Cos(phiJ1),
				      pJ1New*TMath::Sin(thetaJ1)*TMath::Sin(phiJ1),
				      pJ1New*TMath::Cos(thetaJ1));
  Double_t pJ2New = pJ12New.Y();
  TVector3 momentumJ2New = TVector3(pJ2New*TMath::Sin(thetaJ2)*TMath::Cos(phiJ2),
				      pJ2New*TMath::Sin(thetaJ2)*TMath::Sin(phiJ2),
				      pJ2New*TMath::Cos(thetaJ2));
  Double_t mJ1 = lortzJ1.M();
  Double_t mJ2 = lortzJ2.M();
  Double_t eJ1New = TMath::Sqrt(momentumJ1New.Mag()*momentumJ1New.Mag()+
				  mJ1*mJ1);
  Double_t eJ2New = TMath::Sqrt(momentumJ2New.Mag()*momentumJ2New.Mag()+
				  mJ2*mJ2);
  TLorentzVector lortzJ1New = TLorentzVector(momentumJ1New,eJ1New);
  TLorentzVector lortzJ2New = TLorentzVector(momentumJ2New,eJ2New);  
  TLorentzVector lortzHNew = lortzJ1New + lortzJ2New;

  return lortzHNew.M();
}

Double_t ZHll2JAnalysisProcessor::getHMass(TLorentzVector lortzJ1, TLorentzVector lortzJ2, TLorentzVector lortzRecoil,TLorentzVector lortzJ1m, TLorentzVector lortzJ2m) {
  Double_t thetaJ1 = lortzJ1.Theta();
  Double_t thetaJ2 = lortzJ2.Theta();
  Double_t phiJ1 = lortzJ1.Phi();  
  Double_t phiJ2 = lortzJ2.Phi();  
  TVector2 pJ12New = getMomentumNew(lortzRecoil.Px(),lortzRecoil.Py(),
				 thetaJ1,thetaJ2,phiJ1,phiJ2);
  Double_t pJ1New = pJ12New.X();
  TVector3 momentumJ1New = TVector3(pJ1New*TMath::Sin(thetaJ1)*TMath::Cos(phiJ1),
				      pJ1New*TMath::Sin(thetaJ1)*TMath::Sin(phiJ1),
				      pJ1New*TMath::Cos(thetaJ1));
  Double_t pJ2New = pJ12New.Y();
  TVector3 momentumJ2New = TVector3(pJ2New*TMath::Sin(thetaJ2)*TMath::Cos(phiJ2),
				      pJ2New*TMath::Sin(thetaJ2)*TMath::Sin(phiJ2),
				      pJ2New*TMath::Cos(thetaJ2));
  Double_t mJ1 = lortzJ1m.M();
  Double_t mJ2 = lortzJ2m.M();
  Double_t eJ1New = TMath::Sqrt(momentumJ1New.Mag()*momentumJ1New.Mag()+
				  mJ1*mJ1);
  Double_t eJ2New = TMath::Sqrt(momentumJ2New.Mag()*momentumJ2New.Mag()+
				  mJ2*mJ2);
  TLorentzVector lortzJ1New = TLorentzVector(momentumJ1New,eJ1New);
  TLorentzVector lortzJ2New = TLorentzVector(momentumJ2New,eJ2New);  
  TLorentzVector lortzHNew = lortzJ1New + lortzJ2New;

  return lortzHNew.M();
}
