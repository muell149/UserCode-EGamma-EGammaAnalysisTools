#include <TFile.h>
#include "../interface/PFIsolationEstimator.h"
#include <cmath>
#include "DataFormats/Math/interface/deltaR.h"
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
void PFIsolationEstimator::initialize( Bool_t  bApplyVeto, int iParticleType ) {

  setParticleType(iParticleType);

  //By default check for an option vertex association
  checkClosestZVertex = kTRUE;

  //Default cone size is 0.4 for electrons and photons
  setConeSize(0.4);
  
  //Apply vetoes
  setApplyVeto(bApplyVeto);
  
  setDeltaRVetoBarrelPhotons();
  setDeltaRVetoBarrelNeutrals();
  setDeltaRVetoBarrelCharged();
  setDeltaRVetoEndcapPhotons();
  setDeltaRVetoEndcapNeutrals();
  setDeltaRVetoEndcapCharged();

  
  setRectangleDeltaPhiVetoBarrelPhotons();
  setRectangleDeltaPhiVetoBarrelNeutrals();
  setRectangleDeltaPhiVetoBarrelCharged();
  setRectangleDeltaPhiVetoEndcapPhotons();
  setRectangleDeltaPhiVetoEndcapNeutrals();
  setRectangleDeltaPhiVetoEndcapCharged();
  

  setRectangleDeltaEtaVetoBarrelPhotons();
  setRectangleDeltaEtaVetoBarrelNeutrals();
  setRectangleDeltaEtaVetoBarrelCharged();
  setRectangleDeltaEtaVetoEndcapPhotons();
  setRectangleDeltaEtaVetoEndcapNeutrals();
  setRectangleDeltaEtaVetoEndcapCharged();


  if(bApplyVeto && iParticleType==kElectron){
    //Setup veto conditions for electrons
    setDeltaRVetoBarrel(kFALSE);
    setDeltaRVetoEndcap(kTRUE);
    setRectangleVetoBarrel(kFALSE);
    setRectangleVetoEndcap(kFALSE);
    
    //Current recommended default value for the electrons
    setDeltaRVetoEndcapPhotons(0.08);
    setDeltaRVetoEndcapCharged(0.015);


  }else{
    //Setup veto conditions for photons
    setDeltaRVetoBarrel(kTRUE);
    setDeltaRVetoEndcap(kTRUE);
    setRectangleVetoBarrel(kTRUE);
    setRectangleVetoEndcap(kTRUE);
  }
  

  fisInitialized = kTRUE;

}


void PFIsolationEstimator::initializeElectronIsolation( Bool_t  bApplyVeto ){
  initialize(bApplyVeto,kElectron);
}

void PFIsolationEstimator::initializePhotonIsolation( Bool_t  bApplyVeto ){
  initialize(bApplyVeto,kPhoton);
}

/*

//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fElectronIsolation(const reco::PFCandidate * pfCandidate,const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices) 
{
  
  fIsolation=0.0;

  
  
  return fIsolation;
}
 
//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fPhotonIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl,const reco::VertexCollection& vertices) {
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
float* PFIsolationEstimator::fElectronIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices) {
  
  return fIsolationInRings;
}
 
//--------------------------------------------------------------------------------------------------
float* PFIsolationEstimator::fPhotonIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices) {
  
  return fIsolationInRings;
}
*/

 
//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fGetIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices) {
  setRingSize(fConeSize);
  SetNumbersOfRings(1);

  fGetIsolationInRings( pfCandidate, pfParticlesColl, vertices);

  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
vector<float >  PFIsolationEstimator::fGetIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, const reco::VertexCollection& vertices) {

  
  
  return fIsolationInRings;
}


int PFIsolationEstimator::chargedHadronVertex( const reco::VertexCollection& vertices, const reco::PFCandidate& pfcand ) const {

  //code copied from Florian's PFNoPU class


  reco::TrackBaseRef trackBaseRef( pfcand.trackRef() );

  size_t  iVertex = 0;
  unsigned index=0;
  unsigned nFoundVertex = 0;
  typedef reco::VertexCollection::const_iterator IV;
  typedef reco::Vertex::trackRef_iterator IT;
  float bestweight=0;
  for(IV iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

    const reco::Vertex& vtx = *iv;

    // loop on tracks in vertices
    for(IT iTrack=vtx.tracks_begin();
        iTrack!=vtx.tracks_end(); ++iTrack) {

      const reco::TrackBaseRef& baseRef = *iTrack;

      // one of the tracks in the vertex is the same as 
      // the track considered in the function
      if(baseRef == trackBaseRef ) {
        float w = vtx.trackWeight(baseRef);
        //select the vertex for which the track has the highest weight
        if (w > bestweight){
          bestweight=w;
          iVertex=index;
          nFoundVertex++;
        }
      }
    }
  }

  if (nFoundVertex>0){
    if (nFoundVertex!=1)
      edm::LogWarning("TrackOnTwoVertex")<<"a track is shared by at least two verteces. Used to be an assert";
    return iVertex;
  }
  // no vertex found with this track. 

  // optional: as a secondary solution, associate the closest vertex in z
  if ( checkClosestZVertex ) {

    double dzmin = 10000;
    double ztrack = pfcand.vertex().z();
    bool foundVertex = false;
    index = 0;
    for(IV iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

      double dz = fabs(ztrack - iv->z());
      if(dz<dzmin) {
        dzmin = dz;
        iVertex = index;
        foundVertex = true;
      }
    }

    if( foundVertex ) 
      return iVertex;  
  
  }


  return -1 ;
}


