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

class PFIsolationEstimator{
 public:
  PFIsolationEstimator();
  ~PFIsolationEstimator(); 
  
  enum VetoType {
    kElectron = -1,      // MVA for non-triggering electrons          
    kPhoton  =  1            // MVA for triggering electrons
  };
  
  void     initialize( Bool_t  bApplyVeto, Int_t iParticleType );
  
  Bool_t   isInitialized() const { return fisInitialized; }
  

  /*  Float_t fElectronIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  Float_t fPhotonIsolation(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);

  Float_t* fElectronIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  Float_t* fPhotonIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  */

  Float_t fGetIsolation(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  Float_t* fGetIsolationInRings(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  
 
  int chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const;

  void setConeSize(Float_t fValue = 0.4){ fConeSize = fValue;};

  //Veto booleans
  void setApplyVeto(Bool_t bValue = kTRUE){  bApplyVeto = bValue;};
  void setDeltaRVetoBarrel(Bool_t bValue = kTRUE){  bDeltaRVetoBarrel = bValue;};
  void setDeltaRVetoEndcap(Bool_t bValue = kTRUE){  bDeltaRVetoEndcap = bValue;};
  void setRectangleVetoBarrel(Bool_t bValue = kTRUE){  bRectangleVetoBarrel = bValue;};
  void setRectangleVetoEndcap(Bool_t bValue = kTRUE){  bRectangleVetoEndcap = bValue;};

  //Veto Values
  void setDeltaRVetoBarrelPhotons(Float_t fValue = -1.){fDeltaRVetoBarrelPhotons=fValue;};
  void setDeltaRVetoBarrelNeutrals(Float_t fValue = -1.){fDeltaRVetoBarrelNeutrals=fValue;};
  void setDeltaRVetoBarrelCharged(Float_t fValue = -1.){fDeltaRVetoBarrelPhotons=fValue;};
  void setDeltaRVetoEndcapPhotons(Float_t fValue = -1.){fDeltaRVetoEndcapPhotons=fValue;};
  void setDeltaRVetoEndcapNeutrals(Float_t fValue = -1.){fDeltaRVetoEndcapNeutrals=fValue;};
  void setDeltaRVetoEndcapCharged(Float_t fValue = -1.){fDeltaRVetoEndcapPhotons=fValue;};

  
  void setRectangleDeltaPhiVetoBarrelPhotons(Float_t fValue = -1.){fRectangleDeltaPhiVetoBarrelPhotons=fValue;};
  void setRectangleDeltaPhiVetoBarrelNeutrals(Float_t fValue = -1.){fRectangleDeltaPhiVetoBarrelNeutrals=fValue;};
  void setRectangleDeltaPhiVetoBarrelCharged(Float_t fValue = -1.){fRectangleDeltaPhiVetoBarrelPhotons=fValue;};
  void setRectangleDeltaPhiVetoEndcapPhotons(Float_t fValue = -1.){fRectangleDeltaPhiVetoEndcapPhotons=fValue;};
  void setRectangleDeltaPhiVetoEndcapNeutrals(Float_t fValue = -1.){fRectangleDeltaPhiVetoEndcapNeutrals=fValue;};
  void setRectangleDeltaPhiVetoEndcapCharged(Float_t fValue = -1.){fRectangleDeltaPhiVetoEndcapPhotons=fValue;};
  

  void setRectangleDeltaEtaVetoBarrelPhotons(Float_t fValue = -1.){fRectangleDeltaEtaVetoBarrelPhotons=fValue;};
  void setRectangleDeltaEtaVetoBarrelNeutrals(Float_t fValue = -1.){fRectangleDeltaEtaVetoBarrelNeutrals=fValue;};
  void setRectangleDeltaEtaVetoBarrelCharged(Float_t fValue = -1.){fRectangleDeltaEtaVetoBarrelPhotons=fValue;};
  void setRectangleDeltaEtaVetoEndcapPhotons(Float_t fValue = -1.){fRectangleDeltaEtaVetoEndcapPhotons=fValue;};
  void setRectangleDeltaEtaVetoEndcapNeutrals(Float_t fValue = -1.){fRectangleDeltaEtaVetoEndcapNeutrals=fValue;};
  void setRectangleDeltaEtaVetoEndcapCharged(Float_t fValue = -1.){fRectangleDeltaEtaVetoEndcapPhotons=fValue;};



#ifndef STANDALONE
 
#endif
  
 
 private:
 

  Int_t                     iParticleType;

  Bool_t                    fisInitialized;
  Float_t                   fIsolation;
  Float_t                   fIsolationPhoton;
  Float_t                   fIsolationNeutral;
  Float_t                   fIsolationCharged;
  Float_t                   fIsolationChargedAll;
  
  Float_t                   fIsolationInRings[10];
  
  
  Bool_t                    checkClosestZVertex;
  Float_t                   fConeSize;
  Bool_t                    bApplyVeto;
  
  Bool_t                    bDeltaRVetoBarrel; 
  Bool_t                    bDeltaRVetoEndcap; 
  
  Bool_t                    bRectangleVetoBarrel; 
  Bool_t                    bRectangleVetoEndcap; 
  
  Float_t                   fDeltaRVetoBarrelPhotons; 
  Float_t                   fDeltaRVetoBarrelNeutrals;
  Float_t                   fDeltaRVetoBarrelCharged;

  Float_t                   fDeltaRVetoEndcapPhotons; 
  Float_t                   fDeltaRVetoEndcapNeutrals;
  Float_t                   fDeltaRVetoEndcapCharged;  

  Float_t                   fRectangleDeltaPhiVetoBarrelPhotons; 
  Float_t                   fRectangleDeltaPhiVetoBarrelNeutrals;
  Float_t                   fRectangleDeltaPhiVetoBarrelCharged;

  Float_t                   fRectangleDeltaPhiVetoEndcapPhotons; 
  Float_t                   fRectangleDeltaPhiVetoEndcapNeutrals;
  Float_t                   fRectangleDeltaPhiVetoEndcapCharged;
  
  Float_t                   fRectangleDeltaEtaVetoBarrelPhotons; 
  Float_t                   fRectangleDeltaEtaVetoBarrelNeutrals;
  Float_t                   fRectangleDeltaEtaVetoBarrelCharged;

  Float_t                   fRectangleDeltaEtaVetoEndcapPhotons; 
  Float_t                   fRectangleDeltaEtaVetoEndcapNeutrals;
  Float_t                   fRectangleDeltaEtaVetoEndcapCharged;

};

#endif
