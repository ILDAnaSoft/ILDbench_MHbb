// *******************************************************
// some classes for jet clustering
//                ----junping.tian
// *******************************************************
#include <iostream>

#include "JetAlgorithm.h"

#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCCollectionVec.h" 
#include "EVENT/LCIO.h"
#include "EVENT/LCCollection.h"

#include "TLorentzVector.h"
#include "TVector2.h"

//ClassImp(myjet::JJet)
//ClassImp(myjet::JJets)

using namespace std;

using namespace EVENT;
using namespace IMPL;

namespace myjet{

  // --------------------------
  //  JJet
  // --------------------------
  JJet::JJet(double px,double py,double pz,double e)
    : TLorentzVector(px,py,pz,e),_fParticles(),_fFlavor(0),_fA(1)
  {
  }
  JJet::JJet(const TLorentzVector j)
    : TLorentzVector(j),_fParticles(),_fFlavor(0),_fA(1)
  {
  }
  JJet::JJet(ReconstructedParticle *j)
    : TLorentzVector(j->getMomentum(),j->getEnergy()),
      _fParticles(),_fFlavor(0),_fA(1)
  {
    for (int i=0;i<j->getParticles().size();i++) {
      ReconstructedParticle* p = dynamic_cast<ReconstructedParticle*>(j->getParticles().at(i));
      AddParticle(p);
    }
  }
  JJet::~JJet()
  {
  }
  // --------------------------
  double JJet::GetYij(const JJet &j)
  {
    //  Calculate Yij, the distance to another jet
    //
    TVector3 v1 = (TVector3)Vect();
    TVector3 v2 = (TVector3)j.Vect();
    double costheta = v1.Dot(v2)/P()/j.P();
    TVector3 v12 = v1 + v2;
    double e12 = E() + j.E();
    double q2 = e12*e12 - v12.Mag2();
    double emin = j.E() > E() ? E() : j.E();
    double x = emin/e12;
    //    double z = E() > j.E() ? E()/e12 : j.E()/e12;
    double z = 1.-x;
    //
    int alg = GetAlgorithm();
    // Durham algorithm
    if (alg == 1 || alg == 2) {
      return 2.*emin*emin*(1.-costheta);
    }
    // JADE_E algorithm
    if (alg == 3) {
      return q2;
    }
    // gluon splitting algorithm
    if (alg == 4) {
      return q2*(1.-z)/(1.+z*z);
    }
    // JADE algorithm
    if (alg == 5) {
      return 2.*E()*j.E()*(1.-costheta);
    }
    // gluon splitting algorithm - Durham variant
    if (alg == 6) {
      return 2.*emin*emin*(1.-costheta)*(1.-z)/(1.+z*z);
    }
  }
  // --------------------------
  void JJet::Merge(const JJet &j)
  {
    //  Merge with another jet
    //
    double e = E() + j.E();
    double x = X() + j.X();
    double y = Y() + j.Y();
    double z = Z() + j.Z();
    SetPxPyPzE(x,y,z,e);    
    std::vector<ReconstructedParticle*> jps = j.GetParticles();
    int npj = jps.size();
    for (int jj=0;jj<npj;jj++) _fParticles.push_back(jps[jj]);
  }
  // --------------------------
  void JJet::AddParticle(ReconstructedParticle* p)
  {
    //  add particle to jet
    //
    ReconstructedParticleVec ps = p->getParticles();
    int np = ps.size();
    if (np == 0) {
      _fParticles.push_back(p);
    }
    else {
      for (int i=0;i<np;i++) _fParticles.push_back(ps[i]);
    }
  }
  // --------------------------
  ReconstructedParticleImpl * JJet::GetRPJet() const
  {
    //  return the jet as reconstructed particle for use in marlin processor
    //
    ReconstructedParticleImpl* j = new ReconstructedParticleImpl;
    double momentum[3] = {Px(),Py(),Pz()};
    j->setMomentum(momentum);
    j->setEnergy(E());
    for( int i=0;i<GetN();i++){
      ReconstructedParticle* ji = dynamic_cast<ReconstructedParticle*>(_fParticles[i]);
      j->addParticle(ji);
    }
    return j;    
  }
  // --------------------------
  TVector2 JJet::GetPull()
  {
    //  return the jet pull as defined in arxiv:1001.5027
    //
    TVector2 pull(0.,0.);
    TVector2 r_j = TVector2(Rapidity(),Phi());
    for (int i=0;i<GetN();i++)    {
      ReconstructedParticle* ji = dynamic_cast<ReconstructedParticle*>(_fParticles[i]);
      TLorentzVector li = TLorentzVector(ji->getMomentum(),ji->getEnergy());
      TVector2 r_i = TVector2(li.Rapidity(),li.Phi());
      pull += (r_i-r_j)*(r_i-r_j).Mod()*li.Pt()/Pt();
    }
    return pull;
  }

  
  // --------------------------
  //  JJets
  // --------------------------
  JJets::JJets()
    : TLorentzVector(),_fY(-1),_fYPlus(-1),_fJets(),_fYs(),_fI(0),_fJ(0),_fA(1)
  {
  }
  JJets::JJets(std::vector<JJet> js, int alg)
    : TLorentzVector(),_fY(-1),_fYPlus(-1),_fJets(js),_fYs(),_fI(0),_fJ(0),_fA(alg)
  {
    SetAlgorithm(alg);
    Initialize();
  }
  JJets::JJets(LCCollection* colPFOs, int alg)
    : TLorentzVector(),_fY(-1),_fYPlus(-1),_fJets(),_fYs(),_fI(0),_fJ(0),_fA(alg)
  {
    for (int i=0;i<colPFOs->getNumberOfElements();i++) {
      ReconstructedParticle* p = dynamic_cast<ReconstructedParticle*>(colPFOs->getElementAt(i));
      TLorentzVector lortz = TLorentzVector(p->getMomentum(),p->getEnergy());
      JJet j(lortz);
      j.AddParticle(p);
      j.SetAlgorithm(_fA);
      Add(j);
      //      _fJets.push_back(j);
    }
    //    Initialize();
  }
  JJets::~JJets()
  {
  }
  // --------------------------
  void JJets::Initialize()
  {
    //  initialize when given all input particles/jets
    //  notice the ordering of Yij
    //
    TLorentzVector lortz = _fJets[0];
    _fY = -1.;
    for (int j=1;j<GetN();j++) {
      lortz += _fJets[j];
      for (int i=0;i<j;i++) {
	const double y = _fJets[i].GetYij(_fJets[j]);	  
	_fYs.push_back(y);
	if (y < _fY || _fY < 0) {
	  _fY = y;
	  _fI = i;
	  _fJ = j;
	}
      }
    }
    SetPxPyPzE(lortz.Px(),lortz.Py(),lortz.Pz(),lortz.E());

  }
  // --------------------------
  void JJets::SetAlgorithm(int i)
  {
    _fA = i;
    for (int j=0;j<GetN();j++) _fJets[j].SetAlgorithm(i);
  }
  // --------------------------
  void JJets::UpdateMinimum()
  {
    //  find current closest two particles/jets 
    //
    _fYPlus = _fY;
    _fY = -1.;
    for (int i=0;i<GetN()-1;i++) {
      for (int j=i+1;j<GetN();j++) {
	double y = _fYs[j*(j-1)/2+i];
	if (y < _fY || _fY < 0) {
	  _fY = y;
	  _fI = i;
	  _fJ = j;
	}
      }
    }
  }
  // --------------------------
  void JJets::UpdateAll()
  {
    //  usually only for purpose of using mixed clusterig algorithms
    //  re-calculate all Yij for current particles/jets
    //  find current closest two particles/jets 
    //
    _fYPlus = _fY;
    _fY = -1.;
    for (int i=0;i<GetN()-1;i++) {
      for (int j=i+1;j<GetN();j++) {
	double y = _fJets[i].GetYij(_fJets[j]);
	_fYs[j*(j-1)/2+i] = y;
	if (y < _fY || _fY < 0) {
	  _fY = y;
	  _fI = i;
	  _fJ = j;
	}
      }
    }
  }
  // --------------------------
  void JJets::Add(JJet j)
  {
    //  add a new particle/jet
    //
    _fJets.push_back(j);
    SetPxPyPzE(Px()+j.Px(),Py()+j.Py(),Pz()+j.Pz(),E()+j.E());    
    for (int i=0;i<GetN()-1;i++) {
      double y = _fJets[i].GetYij(j);
      _fYs.push_back(y);
      if (y < _fY || _fY < 0) {
	_fY = y;
	_fI = i;
	_fJ = GetN()-1;
      }
    }
  }
  // --------------------------
  void JJets::DoClusteringN(int ncut)
  {
    //  do fixed n jet clustering
    //
    while (GetN()>=2 && GetN() > ncut) Process();
  }
  // --------------------------
  void JJets::DoClusteringY(double ycut)
  {
    //  do fixed y cut clustering
    //
    while (GetN()>=2 && GetY() < ycut*E()*E()) Process();
  }
  // --------------------------
  void JJets::DoClusteringY0(double ycut)
  {
    //  do fixed y cut clustering
    //
    while (GetN()>=2 && GetY() < ycut) Process();
  }
  // --------------------------
  void JJets::DoClustering(double cut)
  {
    if (_fA==1) DoClusteringN((int) cut);
    if (_fA==2) DoClusteringY(cut);
  }
  // --------------------------
  void JJets::DoClusteringNN(int ncut1, int ncut2, int alg1, int alg2)
  {
    //  do two steps mini-jet clustering
    //
    if (alg1 != GetAlgorithm()) {
      SetAlgorithm(alg1);
      UpdateAll();
    }
    while (GetN()>=2 && GetN() > ncut1) Process();
    if (alg2 != alg1) {
      SetAlgorithm(alg2);
      UpdateAll();
    }
    while (GetN()>=2 && GetN() > ncut2) Process();
  }
  // --------------------------
  void JJets::Process()
  {
    //  iterate operations to merge the two closest particle/jet
    //  such as: add four momentum, remove later one (j) from jet vector,
    //           remove Yij values related to jet j (make sure later first),
    //           update Yij values related to jet i,
    //           and find new closest particle/jet
    //
    _fJets[_fI].Merge(_fJets[_fJ]);

#if 0
    for (int i=GetN()-1;i>_fJ;i--) _fYs.erase(_fYs.begin()+i*(i-1)/2+_fJ);
    for (int i=_fJ-1;i>=0;i--) _fYs.erase(_fYs.begin()+_fJ*(_fJ-1)/2+i);
#else    
    const double fdrop_val = -999.0;
    for (int i=GetN()-1; i>_fJ; i--)
      _fYs.at(i*(i-1)/2 + _fJ) = fdrop_val;

    for (int i=_fJ-1; i>=0; i--)
      _fYs.at(_fJ*(_fJ-1)/2 + i) = fdrop_val;

    std::vector<double>::iterator it = _fYs.end();
    it = std::remove_if(_fYs.begin(), _fYs.end(), true_if_value_minus_999);
    if (it != _fYs.end()) // we need to resize _fYs
      _fYs.resize(it - _fYs.begin());
#endif    

    _fJets.erase(_fJets.begin()+_fJ);

    for (int j=_fI+1;j<GetN();j++) {
      _fYs[j*(j-1)/2+_fI] = _fJets[_fI].GetYij(_fJets[j]);
    }
    for (int j=0;j<_fI;j++) {
      _fYs[_fI*(_fI-1)/2+j] = _fJets[_fI].GetYij(_fJets[j]);
    }
    UpdateMinimum();
  }
  // --------------------------
  LCCollectionVec* JJets::GetJetsCol()
  {
    //  return directly jets collection for convenient use in marlin processors
    //
    LCCollectionVec* jetsCol= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    for (int i=0;i<GetN();i++) {
      ReconstructedParticle* j = dynamic_cast<ReconstructedParticle*>(_fJets[i].GetRPJet());
      jetsCol->addElement(j);
    }
    jetsCol->parameters().setValue( "YMinus", (float)GetYMinus() );
    jetsCol->parameters().setValue( "YPlus" , (float)GetYPlus()  );
    return jetsCol;
  }

  //=========================================================================
  bool true_if_value_minus_999(double value)
  {
    return (value == -999.0);
  }
  
} // namespace myjet
