#ifndef JETALGORITHM_H
#define JETALGORITHM_H

// *******************************************************
// some classes for jet clustering
//                ----junping.tian
// *******************************************************

#include <math.h>

#include "TLorentzVector.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h" 
#include "EVENT/LCCollection.h"

#include "TVector2.h"

namespace myjet{

  // unary predicate to be used with std::remove_if
  bool true_if_value_minus_999(double value);

  // --------------------------------------------------------------------
  class JJet : public TLorentzVector{
  public:
    JJet(double px,double py,double pz,double e);
    JJet(const TLorentzVector j);
    JJet(EVENT::ReconstructedParticle *j);
    virtual ~JJet();

    inline int GetFlavor() const {return _fFlavor;};
    inline int GetAlgorithm() {return _fA;};
    double GetYij(const JJet &j);
    void Merge(const JJet &j);
    inline std::vector<EVENT::ReconstructedParticle*> GetParticles() const {return _fParticles;};
    void AddParticle(EVENT::ReconstructedParticle* p);
    IMPL::ReconstructedParticleImpl * GetRPJet() const;
    inline int GetN() const {return _fParticles.size();};
    TVector2 GetPull();

    inline void SetFlavor(int i) {_fFlavor = i;};
    inline void SetAlgorithm(int i) {_fA = i;};
    
  private:
    std::vector<EVENT::ReconstructedParticle*> _fParticles;
    int _fFlavor,_fA;

    //    ClassDef(JJet, 1)
  };

  // --------------------------------------------------------------------
  class JJets : public TLorentzVector{
  public:
    JJets();
    JJets(std::vector<JJet> js, int alg = 1);
    JJets(EVENT::LCCollection* colPFOs, int alg = 1);
    virtual ~JJets();

    void SetAlgorithm(int i);

    inline double GetY() {return _fY;};
    inline double GetYMinus() {return _fY/E()/E();};
    inline double GetYPlus() {return _fYPlus/E()/E();};
    inline int GetN() {return _fJets.size();};
    inline std::vector<JJet> GetJets() {return _fJets;};
    inline int GetAlgorithm() {return _fA;};

    void Initialize();
    void Process();
    void UpdateAll();
    void UpdateMinimum();
    void Add(JJet j);
    void DoClusteringN(int ncut);
    void DoClusteringY(double ycut);
    void DoClusteringY0(double ycut);
    void DoClustering(double cut);
    void DoClusteringNN(int ncut1, int ncut2, int alg1, int alg2);
    IMPL::LCCollectionVec* GetJetsCol();
    
  private:
    double _fY,_fYPlus;
    std::vector<JJet> _fJets;
    std::vector<double> _fYs;
    int _fI,_fJ;
    int _fA; // jet algorithm -- 1: Durham fixed Njet; 2: Durham fixed Ycut
        
    //    ClassDef(JJets, 1)
  };
}

#endif
