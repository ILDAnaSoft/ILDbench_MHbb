// *****************************************************
// e+e- ------> ZH ------> (l+ l-) + X
// Processor for leptons selection
//                        ----Junping
// *****************************************************
#include "ZHll2JLeptonSelectionMVAProcessor.h"
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
#include <EVENT/Vertex.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "TMVA/Reader.h"

#include "Utilities.h"

//#define __ISR_subtract

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace TMVA;
using namespace mylib;

ZHll2JLeptonSelectionMVAProcessor aZHll2JLeptonSelectionMVAProcessor ;


ZHll2JLeptonSelectionMVAProcessor::ZHll2JLeptonSelectionMVAProcessor() : Processor("ZHll2JLeptonSelectionMVAProcessor") {
  
  // modify processor description
  _description = "ZHll2JLeptonSelectionMVAProcessor does whatever it does ..." ;
  

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

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputNewPFOsCollection",
			    "Name of the new PFOs collection without isolated lepton ",
			    _colNewPFOs,
			    std::string("newPandoraPFOs") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputLeptonsCollection",
			    "Name of collection with the selected leptons",
			    _colLeptons,
			    std::string("leptons") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputLeptonsOrigCollection",
			    "Name of collection with the selected leptons w/o FSR and BS recovery",
			    _colLeptonsOrig,
			    std::string("leptonsOrig") );

  registerProcessorParameter("DirOfISOElectronWeights",
			   "Directory of Weights for the Isolated Electron MVA Classification"  ,
			   _isolated_electron_weights ,
			   std::string("isolated_electron_weights") ) ;

  registerProcessorParameter("DirOfISOMuonWeights",
			   "Directory of Weights for the Isolated Muon MVA Classification"  ,
			   _isolated_muon_weights ,
			   std::string("isolated_muon_weights") ) ;

  registerProcessorParameter("CutOnTheISOElectronMVA",
			   "Cut on the mva output of isolated electron selection"  ,
			   _mvacut_electron ,
			   float(0.5) ) ;

  registerProcessorParameter("CutOnTheISOMuonMVA",
			   "Cut on the mva output of isolated muon selection"  ,
			   _mvacut_muon ,
			   float(0.7) ) ;

  registerProcessorParameter("MethodForISOLeptonSelection",
			   "Which mehtod for isolated lepton selection: 1--MVA; 2-traditonal cone energy"  ,
			   _method_iso ,
			   int(1) ) ;
}

void ZHll2JLeptonSelectionMVAProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStat = 0;

  //  TFile *outRootFile = new TFile("output.root","RECREATE");

  for (Int_t i=0;i<2;i++) {
    TMVA::Reader *reader = new TMVA::Reader( "Color:Silent" );    
    // add variables
    if (i == 0) {  // electron
      reader->AddVariable( "coneec",          &_coneec );
      reader->AddVariable( "coneen",          &_coneen );
      reader->AddVariable( "momentum",        &_momentum );
      reader->AddVariable( "coslarcon",       &_coslarcon );
      reader->AddVariable( "energyratio",     &_energyratio );
      reader->AddVariable( "ratioecal",       &_ratioecal );
      reader->AddVariable( "ratiototcal",     &_ratiototcal );
      reader->AddVariable( "nsigd0",          &_nsigd0 );
      reader->AddVariable( "nsigz0",          &_nsigz0 );
    }
    else {
      reader->AddVariable( "coneec",          &_coneec );
      reader->AddVariable( "coneen",          &_coneen );
      reader->AddVariable( "momentum",        &_momentum );
      reader->AddVariable( "coslarcon",       &_coslarcon );
      reader->AddVariable( "energyratio",     &_energyratio );
      //      reader->AddVariable( "yokeenergy",      &_yokeenergy );
      reader->AddVariable( "nsigd0",          &_nsigd0 );
      reader->AddVariable( "nsigz0",          &_nsigz0 );
      reader->AddVariable( "totalcalenergy",  &_totalcalenergy );
    }
    
    // book the reader (method, weights)
    TString dir    = _isolated_electron_weights;
    if (i == 1) dir = _isolated_muon_weights;
    TString prefix = "TMVAClassification";
    TString methodName = "MLP method";
    TString weightfile = dir + "/" + prefix + "_" + "MLP.weights.xml";
    reader->BookMVA( methodName, weightfile ); 
    _readers.push_back(reader);
  }
  
}

void ZHll2JLeptonSelectionMVAProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void ZHll2JLeptonSelectionMVAProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

#if 1
  Double_t fEtrackCut = -1.;        // lower edge of each PFO energy 
#else
  Double_t fEtrackCut = 0.05;        // lower edge of each PFO energy 
#endif
  Double_t fElectronCut1 = 0.5;      // lower edge of totalCalEnergy/momentum
  Double_t fElectronCut2 = 1.3;      // upper edge of totalCalEnergy/momentum
  Double_t fElectronCut3 = 0.9;     // lower edge of ecalEnergy/totalCalEnergy
  Double_t fMuonCut1 = 0.3;      // upper edge of totalCalEnergy/momentum
  //  Double_t fMuonCut2 = 0.5;      // upper edge of ecalEnergy/totalCalEnergy
  //  Double_t fMuonCut3 = 1.2;      // lower edge of yoke energy
  Double_t fMuonCut3 = -999;      // lower edge of yoke energy  
  Double_t fCosConeCut = 0.98;   // the angle of cone around the direction of pfo
  Double_t fCosLargeConeCut = 0.95; // angel of large cone around the pfo
  Double_t fCosFSRCut = 0.99;   // the angle of BS and FSR around the direction of charged lepton
  Double_t kMassZ = 91.187;     // Z mass
  Double_t fMassZCut = 40.;     // mass cut for lepton pair from Z
  Double_t fEpsilon = 1.E-10;
  // Fisher coefficients
  Double_t c0Electron = 12.2; 
  Double_t c1Electron = 0.87; 
  Double_t c0Muon = 12.6; 
  Double_t c1Muon = 4.62; 
  Double_t fPElectronCut = 5.;
  Double_t fPMuonCut = 5.;

  // cut table
  if (!hStat) hStat = new TH1D("hStat", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  streamlog_out(DEBUG) << "Hello, MVA Lepton Pair Selection!" << endl;

  static TNtupleD *hGen = 0;
  if (!hGen) {
    stringstream tupstr_gen;
    tupstr_gen << "nhbb:nhww:nhgg:nhtt:nhcc:nhzz:nhaa:nhmm" << ":"
	       << "elep1:elep2:coslep1:coslep2:ptz:cosz:mz:acop12:mrecoil" << ":"
	       << "eisr1:eisr2" << ":"
	       << "pxlep1:pylep1:pzlep1:pxlep2:pylep2:pzlep2" << ":"
	       << "pxj1:pyj1:pzj1:ej1:pxj2:pyj2:pzj2:ej2:cosj1:phij1:cosj2:phij2" << ":"
	       << "pxnewj1:pynewj1:pznewj1:enewj1:pxnewj2:pynewj2:pznewj2:enewj2" << ":"
	       << "mhnew:cosisr1:cosisr2"
	       << ends;
    hGen = new TNtupleD("hGen","",tupstr_gen.str().data());
  }
  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr_gen;
    tupstr_gen << "nhbb:elep1mc:elep2mc:coslep1mc:coslep2mc:ptzmc:coszmc:mzmc:acop12mc:mrecoilmc" << ":"
	       << "eisr1mc:eisr2mc" << ":"
	       << "elep1:elep2:coslep1:coslep2:ptz:cosz:mz:acop12:mrecoil" << ":"
	       << "elep1orig:elep2orig:coslep1orig:coslep2orig:ptzorig:coszorig:mzorig:acop12orig:mrecoilorig" << ":"
	       << "mvalep1:mvalep2:leptype:nphoton"
	       << ends;
    hEvt = new TNtupleD("hEvt","",tupstr_gen.str().data());
  }
  static TNtupleD *hVtx = 0;
  if (!hVtx) {
    stringstream tupstr_gen;
    tupstr_gen << "zip:dip:ntrks:zipmc:xipmc:yipmc:zipvf:xipvf:yipvf"
	       << ends;
    hVtx = new TNtupleD("hVtx","",tupstr_gen.str().data());
  }
  // -- Get the MCTruth Linker --
  //  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  //  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

  // -- Read out MC information --  
  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nMCP = colMC->getNumberOfElements();
  TLorentzVector lortzLep1MC,lortzLep2MC,lortzZMC,lortzHMC,lortzISRMC;
  TLorentzVector lortzISR1MC, lortzISR2MC;
  TLorentzVector lortzJ1MC, lortzJ2MC;  
  Double_t z_IPMC=999.,x_IPMC=999.,y_IPMC=999.;
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
#if 1
      TVector3 vip = TVector3(mcPart->getVertex());
      z_IPMC = vip[2];
      x_IPMC = vip[0];
      y_IPMC = vip[1];
#endif      
    }
    if (i == 1) {
      lortzISR2MC = lortz;
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

  const Double_t fEcm = 500.;
  TLorentzVector lortzEcm = getLorentzEcm(fEcm);
  TLorentzVector lortzRecoilMC = lortzEcm - lortzZMC;
#ifdef __ISR_subtract
  lortzRecoilMC -= lortzISR1MC;
  lortzRecoilMC -= lortzISR2MC;
#endif  

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
  
  Double_t data_gen[50];
  data_gen[ 0]= nHbb;
  data_gen[ 1]= nHWW;
  data_gen[ 2]= nHgg;
  data_gen[ 3]= nHtt;
  data_gen[ 4]= nHcc;
  data_gen[ 5]= nHZZ;
  data_gen[ 6]= nHaa;
  data_gen[ 7]= nHmm;
  data_gen[ 8]= lortzLep1MC.E();
  data_gen[ 9]= lortzLep2MC.E();
  data_gen[10]= lortzLep1MC.CosTheta();
  data_gen[11]= lortzLep2MC.CosTheta();
  data_gen[12]= lortzZMC.Pt();
  data_gen[13]= lortzZMC.CosTheta();
  data_gen[14]= lortzZMC.M();
  data_gen[15]= getAcoPlanarity(lortzLep1MC,lortzLep2MC);
  data_gen[16]= getRecoilMass(lortzEcm,lortzZMC);
  data_gen[17]= lortzISR1MC.E();
  data_gen[18]= lortzISR2MC.E();
  data_gen[19]= lortzLep1MC.Px();  
  data_gen[20]= lortzLep1MC.Py();  
  data_gen[21]= lortzLep1MC.Pz();  
  data_gen[22]= lortzLep2MC.Px();  
  data_gen[23]= lortzLep2MC.Py();  
  data_gen[24]= lortzLep2MC.Pz();  
  data_gen[25]= lortzJ1MC.Px();  
  data_gen[26]= lortzJ1MC.Py();  
  data_gen[27]= lortzJ1MC.Pz();  
  data_gen[28]= lortzJ1MC.E();  
  data_gen[29]= lortzJ2MC.Px();  
  data_gen[30]= lortzJ2MC.Py();  
  data_gen[31]= lortzJ2MC.Pz();  
  data_gen[32]= lortzJ2MC.E();  
  data_gen[33]= lortzJ1MC.CosTheta();  
  data_gen[34]= lortzJ1MC.Phi();  
  data_gen[35]= lortzJ2MC.CosTheta();  
  data_gen[36]= lortzJ2MC.Phi();  
  data_gen[37]= lortzJ1NewMC.Px();  
  data_gen[38]= lortzJ1NewMC.Py();  
  data_gen[39]= lortzJ1NewMC.Pz();  
  data_gen[40]= lortzJ1NewMC.E();  
  data_gen[41]= lortzJ2NewMC.Px();  
  data_gen[42]= lortzJ2NewMC.Py();  
  data_gen[43]= lortzJ2NewMC.Pz();  
  data_gen[44]= lortzJ2NewMC.E();
  data_gen[45]= lortzHNewMC.M();    
  data_gen[46]= lortzISR1MC.CosTheta();
  data_gen[47]= lortzISR2MC.CosTheta();
  hGen->Fill(data_gen);

  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
    std::cerr << "No PFO Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "MCParticle and PandoraPFOs Collections found!" << ends;
  Int_t nPFOs = colPFO->getNumberOfElements();
  //  cerr << "Number of PFOs: " << nPFOs << endl;
  LCCollectionVec *pNewPFOsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *pLeptonsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *pLeptonsOrigCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  pNewPFOsCollection->setSubset(true);
  pLeptonsOrigCollection->setSubset(true);
  std::vector<lcio::ReconstructedParticle*> newPFOs;
  std::vector<lcio::ReconstructedParticle*> leptons;
  std::vector<lcio::ReconstructedParticle*> electrons;
  std::vector<lcio::ReconstructedParticle*> muons;

  // primary vertex from MCTruth: z_IPMC

  // primary vertex from VertexFinder (LCFIPlus)
  LCCollection *colPVtx = evt->getCollection("PrimaryVertex");
  Int_t npvtx = colPVtx->getNumberOfElements();
  Vertex *pvtx = dynamic_cast<Vertex*>(colPVtx->getElementAt(0));
  TVector3 v_pvtx = TVector3(pvtx->getPosition());
  Double_t z_pvtx = v_pvtx[2];
  Double_t x_pvtx = v_pvtx[0];
  Double_t y_pvtx = v_pvtx[1];

  
  // primary vertex from simple algorithm
  // d = 0
  // z = pt weighted average of candidate tracks from primary IP
  Double_t z_sum = 0., w_sum = 0.;
  Int_t nTrksIP = 0;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    Double_t energy = recPart->getEnergy();
    Double_t charge = recPart->getCharge();
    TVector3 momentum = TVector3(recPart->getMomentum());
    Double_t pt = momentum.Pt();
    TrackVec tckvec = recPart->getTracks();
    Int_t ntracks = tckvec.size();
    Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
    if (ntracks > 0) {
      d0 = tckvec[0]->getD0();
      z0 = tckvec[0]->getZ0();
      deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
      deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
      nsigd0 = d0/deltad0;
      nsigz0 = z0/deltaz0;
      if (TMath::Abs(nsigd0) > 10) continue; // candidate track: |nsigd0| < 10
      Double_t wpt = pt/deltaz0/deltaz0;
      w_sum += wpt;
      z_sum += wpt*z0;
      nTrksIP++;
    }
  }
  Double_t z_IP = 0., d_IP = 0.;
  if (nTrksIP > 0) z_IP = z_sum/w_sum;
  Double_t data_vtx[10];
  data_vtx[ 0] = z_IP;
  data_vtx[ 1] = d_IP;
  data_vtx[ 2] = nTrksIP;
  data_vtx[ 3] = z_IPMC;
  data_vtx[ 4] = x_IPMC;
  data_vtx[ 5] = y_IPMC;
  data_vtx[ 6] = z_pvtx;  
  data_vtx[ 7] = x_pvtx;  
  data_vtx[ 8] = y_pvtx;  
  hVtx->Fill(data_vtx);

  // loop all the PFOs
  const Int_t nLepMax = 20;
  ReconstructedParticle *leps[nLepMax];
  Double_t mva_outs[nLepMax];
  Int_t nLeps = 0;
  float _mva_lep_minus = -1.;
  float _mva_lep_plus = -1.;
  int _lep_type = 0;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    Double_t energy = recPart->getEnergy();
    Double_t charge = recPart->getCharge();
    TrackVec tckvec = recPart->getTracks();
    Int_t ntracks = tckvec.size();
    Double_t d0=0.,z0=0.,deltad0=0.,deltaz0=0.,nsigd0=0.,nsigz0=0.;
    if (ntracks > 0) {
      d0 = tckvec[0]->getD0();
      z0 = tckvec[0]->getZ0();
      //      z0 -= z_IP;
      //      z0 -= z_IPMC;
      z0 -= z_pvtx;
      //      d0 -= d_IP;
      deltad0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[0]);
      deltaz0 = TMath::Sqrt(tckvec[0]->getCovMatrix()[9]);
      nsigd0 = d0/deltad0;
      nsigz0 = z0/deltaz0;
    }
    //    Double_t r0 = TMath::Sqrt(d0*d0+z0*z0);
    //    Double_t nsigr0 = TMath::Sqrt(nsigd0*nsigd0+nsigz0*nsigz0);
    if (energy > fEtrackCut) {
      newPFOs.push_back(recPart);
      Double_t ecalEnergy = 0;
      Double_t hcalEnergy = 0;
      Double_t yokeEnergy = 0;
      Double_t totalCalEnergy = 0;
      Int_t nHits = 0;
      std::vector<lcio::Cluster*> clusters = recPart->getClusters();
      for (std::vector<lcio::Cluster*>::const_iterator iCluster=clusters.begin();iCluster!=clusters.end();++iCluster) {
	ecalEnergy += (*iCluster)->getSubdetectorEnergies()[0];
	hcalEnergy += (*iCluster)->getSubdetectorEnergies()[1];
	yokeEnergy += (*iCluster)->getSubdetectorEnergies()[2];
	ecalEnergy += (*iCluster)->getSubdetectorEnergies()[3];
	hcalEnergy += (*iCluster)->getSubdetectorEnergies()[4];
	CalorimeterHitVec calHits = (*iCluster)->getCalorimeterHits();
	//	nHits += (*iCluster)->getCalorimeterHits().size();
	nHits = calHits.size();
      }
      totalCalEnergy = ecalEnergy + hcalEnergy;
      TVector3 momentum = TVector3(recPart->getMomentum());
      Double_t momentumMagnitude = momentum.Mag();
      //      Double_t cosTheta = momentum.CosTheta();
      //get cone information
      Bool_t woFSR = kTRUE;
      Double_t coneEnergy0[3] = {0.,0.,0.};
      Double_t pFSR[4] = {0.,0.,0.,0.};
      Double_t pLargeCone[4]  = {0.,0.,0.,0.};
      Int_t nConePhoton = 0;
      getConeEnergy(recPart,colPFO,fCosConeCut,woFSR,coneEnergy0,pFSR,fCosLargeConeCut,pLargeCone,nConePhoton);
      //      Double_t coneEnergy = coneEnergy0[0];
      Double_t coneEN     = coneEnergy0[1];
      Double_t coneEC     = coneEnergy0[2];
      TLorentzVector lortzFSR = TLorentzVector(pFSR[0],pFSR[1],pFSR[2],pFSR[3]);
      TLorentzVector lortzLargeCone = TLorentzVector(pLargeCone[0],pLargeCone[1],pLargeCone[2],pLargeCone[3]);
      TVector3 momentumLargeCone = lortzLargeCone.Vect();
      Double_t cosThetaWithLargeCone = 1.;
      if (momentumLargeCone.Mag() > 0.0000001) {
	cosThetaWithLargeCone = momentum.Dot(momentumLargeCone)/momentumMagnitude/momentumLargeCone.Mag();
      }
      Double_t energyRatioWithLargeCone = energy/(energy+lortzLargeCone.E());
      Double_t ratioECal = 0., ratioTotalCal = 0.;
      if (ecalEnergy > 0.) ratioECal = ecalEnergy/totalCalEnergy;
      ratioTotalCal = totalCalEnergy/momentumMagnitude;
      //      Int_t nConePFOs    = conePFOs.size();
      //      Int_t nConePFOs    = 0;
      //      Int_t nConeCharged = 0;
      //      Int_t nConeNeutral = 0;
      // evaluate the neural-net output of isolated-lepton classfication
      Double_t mva_electron = -1.,mva_muon = -1.;
      _coneec      = coneEC;
      _coneen      = coneEN;
      _momentum    = momentumMagnitude;
      _coslarcon   = cosThetaWithLargeCone;
      _energyratio = energyRatioWithLargeCone;
      _ratioecal   = ratioECal;
      _ratiototcal = ratioTotalCal;
      _nsigd0      = nsigd0;
      _nsigz0      = nsigz0;
      _yokeenergy  = yokeEnergy;
      _totalcalenergy = totalCalEnergy;
      if (charge != 0 && 
	  totalCalEnergy/momentumMagnitude > fElectronCut1 && totalCalEnergy/momentumMagnitude < fElectronCut2 &&
	  ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && 
	  (momentumMagnitude > fPElectronCut)) {
	if (TMath::Abs(nsigd0) < 50 && TMath::Abs(nsigz0) < 50) {   // contraint to primary vertex
	  mva_electron = _readers[0]->EvaluateMVA( "MLP method"           );
	  leps[nLeps] = recPart;
	  mva_outs[nLeps] = mva_electron;
	  nLeps++;
	}
      }
      if (charge != 0 && 
	  totalCalEnergy/momentumMagnitude < fMuonCut1 && 
	  yokeEnergy > fMuonCut3 && 
	  (momentumMagnitude > fPMuonCut)) {
	if (TMath::Abs(nsigd0) < 20 && TMath::Abs(nsigz0) < 20) {
	  mva_muon = _readers[1]->EvaluateMVA( "MLP method"           );
	  leps[nLeps] = recPart;
	  mva_outs[nLeps] = mva_muon;
	  nLeps++;
	}
      }
      // select the leptons
      if (_method_iso == 1) {  // MVA
	if (mva_electron > _mvacut_electron) {
	  leptons.push_back(recPart);
	  electrons.push_back(recPart);
	}
	if (mva_muon > _mvacut_muon) {
	  leptons.push_back(recPart);
	  muons.push_back(recPart);
	}
      }
      else if (_method_iso == 2) {  // Econe.vs.P
	if (charge != 0 && 
	    totalCalEnergy/momentumMagnitude > fElectronCut1 && totalCalEnergy/momentumMagnitude < fElectronCut2 &&
	    ecalEnergy/(totalCalEnergy + fEpsilon) > fElectronCut3 && 
	    (momentumMagnitude > c0Electron + c1Electron*coneEC)) {
	  if (TMath::Abs(nsigd0) > 50 || TMath::Abs(nsigz0) > 5) continue;   // contraint to primary vertex
	  leptons.push_back(recPart);
	  electrons.push_back(recPart);
	}
	if (charge != 0 && 
	    totalCalEnergy/momentumMagnitude < fMuonCut1 && 
	    yokeEnergy > fMuonCut3 && 
	    (momentumMagnitude > c0Muon + c1Muon*coneEC)) {
	  if (TMath::Abs(nsigd0) > 5 || TMath::Abs(nsigz0) > 5) continue;  // contraint to primary vertex
	  leptons.push_back(recPart);
	  muons.push_back(recPart);
	}
      }
      else {
	std::cerr << "Not proper method specified for isolated lepton selection!" << std::endl;
	return;
      }
    }
  }
  Int_t nelectrons = electrons.size();
  Int_t nmuons = muons.size();
  if (nelectrons < 2 && nmuons < 2) throw marlin::SkipEventException(this);
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "nLeptons >= 2" << ends;

  // find the lepton pair nearest to the Z mass
  std::vector<lcio::ReconstructedParticle*> electronZ;
  std::vector<lcio::ReconstructedParticle*> muonZ;
  std::vector<lcio::ReconstructedParticle*> leptonZ;
  std::vector<lcio::ReconstructedParticle*> photonSplit;
  Double_t deltaMassZ_M = fMassZCut + 20; // for muon
  Double_t deltaMassZ_E = fMassZCut + 40; // for electron
  Double_t massZ=0.;
  if (electrons.size() > 1) {
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iElectron=electrons.begin();iElectron<electrons.end()-1;iElectron++) {
      for (std::vector<lcio::ReconstructedParticle*>::const_iterator jElectron=iElectron+1;jElectron<electrons.end();jElectron++) {
	if ((*iElectron)->getCharge() != (*jElectron)->getCharge()) {
	  Double_t mass = getInvariantMass((*iElectron),(*jElectron));
	  if (TMath::Abs(mass-kMassZ) < deltaMassZ_E) {
	    deltaMassZ_E = TMath::Abs(mass-kMassZ);
	    massZ = mass;
	    electronZ.clear();
	    electronZ.push_back(*iElectron);
	    electronZ.push_back(*jElectron);
	  }
	}
      }
    }
  }
  if (muons.size() > 1) {
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iMuon=muons.begin();iMuon<muons.end()-1;iMuon++) {
      for (std::vector<lcio::ReconstructedParticle*>::const_iterator jMuon=iMuon+1;jMuon<muons.end();jMuon++) {
	if ((*iMuon)->getCharge() != (*jMuon)->getCharge()) {
	  Double_t mass = getInvariantMass((*iMuon),(*jMuon));
	  if (TMath::Abs(mass-kMassZ) < deltaMassZ_M) {
	    deltaMassZ_M = TMath::Abs(mass-kMassZ);
	    massZ = mass;
	    muonZ.clear();
	    muonZ.push_back(*iMuon);
	    muonZ.push_back(*jMuon);
	  }
	}
      }
    }
  }

  ReconstructedParticle* leptonMinus = 0;
  ReconstructedParticle* leptonPlus  = 0;
  if (muonZ.size() == 0 && electronZ.size() == 0) {
    cerr << "no lepton pair candidate found!" << endl;
    throw marlin::SkipEventException(this);
  }
  else if (muonZ.size() > 0 && deltaMassZ_M <= deltaMassZ_E) {
    _lep_type = 13;
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iMuon=muonZ.begin();iMuon<muonZ.end();iMuon++) {
      if ((*iMuon)->getCharge() < 0) {
	leptonMinus = (*iMuon);
      }
      else {
	leptonPlus = (*iMuon);
      }
    }
  }
  else {
    _lep_type = 11;
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iElectron=electronZ.begin();iElectron<electronZ.end();iElectron++) {
      if ((*iElectron)->getCharge() < 0) {
	leptonMinus = (*iElectron);
      }
      else {
	leptonPlus = (*iElectron);
      }
    }
  }
  for (Int_t i=0;i<nLeps;i++) {
    if (leps[i] == leptonMinus) _mva_lep_minus = mva_outs[i];
    if (leps[i] == leptonPlus)  _mva_lep_plus  = mva_outs[i];
  }

  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "|M(l+l-) - M(Z)| < " << fMassZCut+20 << " ("  << fMassZCut+40 << ")" << ends;

  // a little bit tighter cut on sum of mva_out of two leptons
  if (_lep_type == 11 && _mva_lep_minus + _mva_lep_plus < 1.1) throw marlin::SkipEventException(this);
  if (_lep_type == 13 && _mva_lep_minus + _mva_lep_plus < 1.5) throw marlin::SkipEventException(this);
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "MVA(l+) + MVA(l-) > 1.5 (1.1)" << ends;

  pLeptonsOrigCollection->addElement(leptonMinus);
  pLeptonsOrigCollection->addElement(leptonPlus);

  // recovery of FSR and BS
  std::vector<lcio::ReconstructedParticle*> photons;
  ReconstructedParticleImpl * recoLeptonMinus = new ReconstructedParticleImpl();
  doPhotonRecovery(leptonMinus,colPFO,recoLeptonMinus,fCosFSRCut,_lep_type,photons);
  ReconstructedParticleImpl * recoLeptonPlus = new ReconstructedParticleImpl();
  doPhotonRecovery(leptonPlus,colPFO,recoLeptonPlus,fCosFSRCut,_lep_type,photons);

  // create the lepton collection
  pLeptonsCollection->addElement(recoLeptonMinus);
  pLeptonsCollection->addElement(recoLeptonPlus);
  pLeptonsCollection->parameters().setValue( "MVALepMinus", _mva_lep_minus );
  pLeptonsCollection->parameters().setValue( "MVALepPlus", _mva_lep_plus );
  pLeptonsCollection->parameters().setValue( "ISOLepType", _lep_type );
  Int_t nPhotonsRecovered = photons.size();
  pLeptonsCollection->parameters().setValue( "NPhotons", nPhotonsRecovered );

  // save the other PFOs to a new collection
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=newPFOs.begin();iObj<newPFOs.end();++iObj) {
    if ((*iObj) == leptonMinus || (*iObj) == leptonPlus) continue;
    if (isFoundInVector((*iObj),photons)) continue;
    pNewPFOsCollection->addElement(*iObj);
  }
  evt->addCollection(pNewPFOsCollection,_colNewPFOs.c_str());
  evt->addCollection(pLeptonsCollection,_colLeptons.c_str());
  evt->addCollection(pLeptonsOrigCollection,_colLeptonsOrig.c_str());

  photons.clear();

  TLorentzVector lortzLep1 = getLorentzVector(recoLeptonMinus);
  TLorentzVector lortzLep2 = getLorentzVector(recoLeptonPlus);
  TLorentzVector lortzLepOrig1 = getLorentzVector(leptonMinus);
  TLorentzVector lortzLepOrig2 = getLorentzVector(leptonPlus);
  TLorentzVector lortzZ = lortzLep1 + lortzLep2;
  TLorentzVector lortzZOrig = lortzLepOrig1 + lortzLepOrig2;

  Double_t data_evt[50];
  data_evt[ 0] = nHbb;
  data_evt[ 1] = lortzLep1MC.E();
  data_evt[ 2] = lortzLep2MC.E();
  data_evt[ 3] = lortzLep1MC.CosTheta();
  data_evt[ 4] = lortzLep2MC.CosTheta();
  data_evt[ 5] = lortzZMC.Pt();
  data_evt[ 6] = lortzZMC.CosTheta();
  data_evt[ 7] = lortzZMC.M();
  data_evt[ 8] = getAcoPlanarity(lortzLep1MC,lortzLep2MC);
  data_evt[ 9] = getRecoilMass(lortzEcm,lortzZMC);
  data_evt[10] = lortzISR1MC.E();
  data_evt[11] = lortzISR2MC.E();
  data_evt[12] = lortzLep1.E();
  data_evt[13] = lortzLep2.E();
  data_evt[14] = lortzLep1.CosTheta();
  data_evt[15] = lortzLep2.CosTheta();
  data_evt[16] = lortzZ.Pt();
  data_evt[17] = lortzZ.CosTheta();
  data_evt[18] = lortzZ.M();
  data_evt[19] = getAcoPlanarity(lortzLep1,lortzLep2);
  data_evt[20] = getRecoilMass(lortzEcm,lortzZ);
  data_evt[21] = lortzLepOrig1.E();
  data_evt[22] = lortzLepOrig2.E();
  data_evt[23] = lortzLepOrig1.CosTheta();
  data_evt[24] = lortzLepOrig2.CosTheta();
  data_evt[25] = lortzZOrig.Pt();
  data_evt[26] = lortzZOrig.CosTheta();
  data_evt[27] = lortzZOrig.M();
  data_evt[28] = getAcoPlanarity(lortzLepOrig1,lortzLepOrig2);
  data_evt[29] = getRecoilMass(lortzEcm,lortzZOrig);
  data_evt[30] = _mva_lep_minus;
  data_evt[31] = _mva_lep_plus;
  data_evt[32] = _lep_type;
  data_evt[33] = nPhotonsRecovered;
  hEvt->Fill(data_evt);



  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  last->cd();
}



void ZHll2JLeptonSelectionMVAProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ZHll2JLeptonSelectionMVAProcessor::end(){ 

  for (std::vector<TMVA::Reader*>::const_iterator ireader=_readers.begin();ireader!=_readers.end();++ireader) {
    delete *ireader;
  }

  cerr << "ZHll2JLeptonSelectionMVAProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  
  cerr << "  =============" << endl;
  cerr << "   Cut Summary " << endl;
  cerr << "  =============" << endl;
  cerr << "   ll+X    " << endl;
  cerr << "  =============" << endl;
  cerr << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;

}
