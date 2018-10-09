#ifndef ZHll2JLeptonSelectionMVAProcessor_h
#define ZHll2JLeptonSelectionMVAProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TH1D.h"

#include "TMVA/Reader.h"


/**  processor for selecting the two isolated leptons from Z decay, with BS and FSR recovered.
 * @version $Id: ZHll2JLeptonSelectionMVAProcessor.h,v 1.0 2013-10-31 12:57:39 Junping Exp $ 
 */

class ZHll2JLeptonSelectionMVAProcessor : public marlin::Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ZHll2JLeptonSelectionMVAProcessor ; }
  
  
  ZHll2JLeptonSelectionMVAProcessor() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
 protected:

  /** Input collection name.
   */
  std::string _colMCP ;
  std::string _colMCTL ;
  std::string _colPFOs ;
  std::string _colNewPFOs ;
  std::string _colLeptons ;
  std::string _colLeptonsOrig ;
  std::string _isolated_electron_weights ;
  std::string _isolated_muon_weights ;

  std::vector<TMVA::Reader*> _readers;
  Float_t _coneec, _coneen, _momentum, _coslarcon, _energyratio;
  Float_t _ratioecal, _ratiototcal, _nsigd0, _nsigz0, _yokeenergy, _totalcalenergy;

  float _mvacut_electron, _mvacut_muon;
  int _method_iso;

  int _nRun ;
  int _nEvt ;

  TH1D *hStat;
  std::stringstream gCutName[20];

} ;

#endif



