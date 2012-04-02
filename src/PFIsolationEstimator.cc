#include <TFile.h>
#include "../interface/PFIsolationEstimator.h"
#include <cmath>
using namespace std;

#ifndef STANDALONE
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
using namespace reco;
#endif





//--------------------------------------------------------------------------------------------------
PFIsolationEstimator::PFIsolationEstimator() :
fisInitialized(kFALSE)
{
  // Constructor.
}



//--------------------------------------------------------------------------------------------------
PFIsolationEstimator::~PFIsolationEstimator()
{
 
}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initialize( ) {
  
  fisInitialized = kTRUE;

}



//--------------------------------------------------------------------------------------------------
Float_t PFIsolationEstimator::fElectronIsolation(const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection > ) {
  
  return fIsolation;
}
 
//--------------------------------------------------------------------------------------------------
Float_t PFIsolationEstimator::fPhotonIsolation( const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >) {
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
Float_t* PFIsolationEstimator::fElectronIsolationInRings( const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >) {
  
  return fIsolationInRings;
}
 
//--------------------------------------------------------------------------------------------------
Float_t* PFIsolationEstimator::fPhotonIsolationInRings( const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection >) {
  
  return fIsolationInRings;
}
