// *****************************************************
// Processor for Overlay Removal By MC
//                        ----Junping
// *****************************************************
#include "OverlayRemovalByMCProcessor.h"

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

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

OverlayRemovalByMCProcessor aOverlayRemovalByMCProcessor ;


OverlayRemovalByMCProcessor::OverlayRemovalByMCProcessor() : Processor("OverlayRemovalByMCProcessor") {
  
  // modify processor description
  _description = "OverlayRemovalByMCProcessor does whatever it does ..." ;
  

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
			    "OutputPFOsWithoutOverlayCollection",
			    "Name of the PFOs collection without overlaid gam-gam to hadron background",
			    _colPFOsWithoutOverlay,
			    std::string("PFOsWithoutOverlay") );

}

void OverlayRemovalByMCProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

}

void OverlayRemovalByMCProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void OverlayRemovalByMCProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

  // -- Read out MC information --  
  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }

  // -- Get the MCTruth Linker --
  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);

  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
    std::cerr << "No PFO Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }

  LCCollectionVec *pPFOsWithoutOverlayCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  pPFOsWithoutOverlayCollection->setSubset(true);
  Int_t nPFOs = colPFO->getNumberOfElements();
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    LCObjectVec vecMCTL = navMCTL->getRelatedToObjects(pfo);
    //    Double_t weight = 1.;
    //    Int_t nMCTL = 1;
    //    MCParticle *mc = getLinkedMCParticle(pfo,colMCTL,weight,nMCTL);
    Int_t nMCTL = vecMCTL.size();
    Bool_t iOverlay = kFALSE;
    if (vecMCTL.size() > 0) {
      Int_t iMCTL = -1;
      Double_t maxEnergyMCTL = -1.;
      for (Int_t im=0;im<nMCTL;im++) {
	MCParticle *mcLink = dynamic_cast<MCParticle *>(vecMCTL[im]);
	if (mcLink->getEnergy() > maxEnergyMCTL) {
	  maxEnergyMCTL = mcLink->getEnergy();
	  iMCTL = im;
	}
      }
      MCParticle *mcPart = dynamic_cast<MCParticle *>(vecMCTL[iMCTL]);
      iOverlay = mcPart->isOverlay();
      if (!iOverlay) {
	pPFOsWithoutOverlayCollection->addElement(pfo);
      }
    }
  }
  // add new collections
  evt->addCollection(pPFOsWithoutOverlayCollection,_colPFOsWithoutOverlay.c_str());

  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  //  _nEvt ++ ;

}



void OverlayRemovalByMCProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void OverlayRemovalByMCProcessor::end(){ 

  cerr << "OverlayRemovalByMCProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  
}
