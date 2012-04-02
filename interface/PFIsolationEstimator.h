//--------------------------------------------------------------------------------------------------
// $Id $
//
// PFIsolationEstimator 
//
// Helper Class for calculating PFIsolation for Photons & Electron onthe fly. This class takes 
//        PF Particle collection and the reconstructed vertex collection as input.
//
// Authors: Vasundhara Chetluru
//--------------------------------------------------------------------------------------------------


/// --> NOTE if you want to use this class as standalone without the CMSSW part 
///  you need to uncomment the below line and compile normally with scramv1 b 
///  Then you need just to load it in your root macro the lib with the correct path, eg:
///  gSystem->Load("/data/benedet/CMSSW_5_2_2/lib/slc5_amd64_gcc462/pluginEGammaEGammaAnalysisTools.so");

//#define STANDALONE   // <---- this line

#ifndef PFIsolationEstimator_H
#define PFIsolationEstimator_H

#ifndef STANDALONE
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#endif
#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

using namespace std;

class PFIsolationEstimator{
 public:
  PFIsolationEstimator();
  ~PFIsolationEstimator(); 
  
  enum VetoType {
    kElectron = -1,      // MVA for non-triggering electrons          
    kPhoton  =  1            // MVA for triggering electrons
  };
  
  void     initialize( Bool_t  bApplyVeto, int iParticleType);
  void     initializeElectronIsolation( Bool_t  bApplyVeto );
  void     initializePhotonIsolation( Bool_t  bApplyVeto );

  Bool_t   isInitialized() const { return fisInitialized; }
  

  /*  float fElectronIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  float fPhotonIsolation(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);

  float* fElectronIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  float* fPhotonIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  */

  float fGetIsolation(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl,const reco::VertexCollection& vertices );
  vector<float >  fGetIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices );
  
 
  int chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const;

  void setConeSize(float fValue = 0.4){ fConeSize = fValue;};

  void setParticleType(int iValue){iParticleType = iValue;};

  //Veto booleans
  void setApplyVeto(Bool_t bValue = kTRUE){  bApplyVeto = bValue;};
  void setDeltaRVetoBarrel(Bool_t bValue = kTRUE){  bDeltaRVetoBarrel = bValue;};
  void setDeltaRVetoEndcap(Bool_t bValue = kTRUE){  bDeltaRVetoEndcap = bValue;};
  void setRectangleVetoBarrel(Bool_t bValue = kTRUE){  bRectangleVetoBarrel = bValue;};
  void setRectangleVetoEndcap(Bool_t bValue = kTRUE){  bRectangleVetoEndcap = bValue;};

  //Veto Values
  void setDeltaRVetoBarrelPhotons(float fValue = -1.){fDeltaRVetoBarrelPhotons=fValue;};
  void setDeltaRVetoBarrelNeutrals(float fValue = -1.){fDeltaRVetoBarrelNeutrals=fValue;};
  void setDeltaRVetoBarrelCharged(float fValue = -1.){fDeltaRVetoBarrelPhotons=fValue;};
  void setDeltaRVetoEndcapPhotons(float fValue = -1.){fDeltaRVetoEndcapPhotons=fValue;};
  void setDeltaRVetoEndcapNeutrals(float fValue = -1.){fDeltaRVetoEndcapNeutrals=fValue;};
  void setDeltaRVetoEndcapCharged(float fValue = -1.){fDeltaRVetoEndcapPhotons=fValue;};

  
  void setRectangleDeltaPhiVetoBarrelPhotons(float fValue = -1.){fRectangleDeltaPhiVetoBarrelPhotons=fValue;};
  void setRectangleDeltaPhiVetoBarrelNeutrals(float fValue = -1.){fRectangleDeltaPhiVetoBarrelNeutrals=fValue;};
  void setRectangleDeltaPhiVetoBarrelCharged(float fValue = -1.){fRectangleDeltaPhiVetoBarrelPhotons=fValue;};
  void setRectangleDeltaPhiVetoEndcapPhotons(float fValue = -1.){fRectangleDeltaPhiVetoEndcapPhotons=fValue;};
  void setRectangleDeltaPhiVetoEndcapNeutrals(float fValue = -1.){fRectangleDeltaPhiVetoEndcapNeutrals=fValue;};
  void setRectangleDeltaPhiVetoEndcapCharged(float fValue = -1.){fRectangleDeltaPhiVetoEndcapPhotons=fValue;};
  

  void setRectangleDeltaEtaVetoBarrelPhotons(float fValue = -1.){fRectangleDeltaEtaVetoBarrelPhotons=fValue;};
  void setRectangleDeltaEtaVetoBarrelNeutrals(float fValue = -1.){fRectangleDeltaEtaVetoBarrelNeutrals=fValue;};
  void setRectangleDeltaEtaVetoBarrelCharged(float fValue = -1.){fRectangleDeltaEtaVetoBarrelPhotons=fValue;};
  void setRectangleDeltaEtaVetoEndcapPhotons(float fValue = -1.){fRectangleDeltaEtaVetoEndcapPhotons=fValue;};
  void setRectangleDeltaEtaVetoEndcapNeutrals(float fValue = -1.){fRectangleDeltaEtaVetoEndcapNeutrals=fValue;};
  void setRectangleDeltaEtaVetoEndcapCharged(float fValue = -1.){fRectangleDeltaEtaVetoEndcapPhotons=fValue;};

  float getIsolationPhoton(){ return fIsolationPhoton; };
  float getIsolationNeutral(){ return fIsolationNeutral; };
  float getIsolationCharged(){ return fIsolationCharged; };
  float getIsolationChargedAll(){ return fIsolationChargedAll; };

  vector<float >  getIsolationInRingsPhoton(){ return fIsolationInRingsPhoton; };
  vector<float >  getIsolationInRingsNeutral(){ return fIsolationInRingsNeutral; };
  vector<float >  getIsolationInRingsCharged(){ return fIsolationInRingsCharged; };
  vector<float >  getIsolationInRingsChargedAll(){ return fIsolationInRingsChargedAll; };


  void SetNumbersOfRings(int iValue = 1){iNumberOfRings = iValue;};
  void setRingSize(float fValue = 0.4){fRingSize = fValue;};

#ifndef STANDALONE
 
#endif
  
 
 private:
 

  int                     iParticleType;

  Bool_t                    fisInitialized;
  float                   fIsolation;
  float                   fIsolationPhoton;
  float                   fIsolationNeutral;
  float                   fIsolationCharged;
  float                   fIsolationChargedAll;
  
  vector<float >          fIsolationInRings;
  vector<float >          fIsolationInRingsPhoton;
  vector<float >          fIsolationInRingsNeutral;
  vector<float >          fIsolationInRingsCharged;  
  vector<float >          fIsolationInRingsChargedAll;

  Bool_t                    checkClosestZVertex;
  float                   fConeSize;
  Bool_t                    bApplyVeto;
  
  Bool_t                    bDeltaRVetoBarrel; 
  Bool_t                    bDeltaRVetoEndcap; 
  
  Bool_t                    bRectangleVetoBarrel; 
  Bool_t                    bRectangleVetoEndcap; 
  
  float                   fDeltaRVetoBarrelPhotons; 
  float                   fDeltaRVetoBarrelNeutrals;
  float                   fDeltaRVetoBarrelCharged;

  float                   fDeltaRVetoEndcapPhotons; 
  float                   fDeltaRVetoEndcapNeutrals;
  float                   fDeltaRVetoEndcapCharged;  

  float                   fRectangleDeltaPhiVetoBarrelPhotons; 
  float                   fRectangleDeltaPhiVetoBarrelNeutrals;
  float                   fRectangleDeltaPhiVetoBarrelCharged;

  float                   fRectangleDeltaPhiVetoEndcapPhotons; 
  float                   fRectangleDeltaPhiVetoEndcapNeutrals;
  float                   fRectangleDeltaPhiVetoEndcapCharged;
  
  float                   fRectangleDeltaEtaVetoBarrelPhotons; 
  float                   fRectangleDeltaEtaVetoBarrelNeutrals;
  float                   fRectangleDeltaEtaVetoBarrelCharged;

  float                   fRectangleDeltaEtaVetoEndcapPhotons; 
  float                   fRectangleDeltaEtaVetoEndcapNeutrals;
  float                   fRectangleDeltaEtaVetoEndcapCharged;

  int                     iNumberOfRings;
  float                   fRingSize;

};

#endif
