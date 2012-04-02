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
  
  void     initialize( );
  
  Bool_t   isInitialized() const { return fisInitialized; }
  

  Float_t fElectronIsolation(const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  Float_t fPhotonIsolation(const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);

  Float_t* fElectronIsolationInRings(const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  Float_t* fPhotonIsolationInRings(const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >);
  

#ifndef STANDALONE
 
#endif
  
 
 private:
 
  Bool_t                    fisInitialized;
  Float_t                   fIsolation;
  Float_t                   fIsolationInRings[10];
  

};

#endif
