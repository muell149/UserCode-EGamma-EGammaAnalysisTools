#include <TFile.h>
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"
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

#endif

using namespace reco;



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
  //  setConeSize(0.4);
  
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


//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeElectronIsolation( Bool_t  bApplyVeto,float  fConeSize){
  initialize(bApplyVeto,kElectron);
  initializeRings(1, fConeSize);

}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializePhotonIsolation( Bool_t  bApplyVeto, float fConeSize ){
  initialize(bApplyVeto,kPhoton);
  initializeRings(1, fConeSize);
}


//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeElectronIsolationInRings( Bool_t  bApplyVeto, int iNumberOfRings, float fRingSize ){
  initialize(bApplyVeto,kElectron);
  initializeRings(iNumberOfRings, fRingSize);
}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializePhotonIsolationInRings( Bool_t  bApplyVeto, int iNumberOfRings, float fRingSize  ){
  initialize(bApplyVeto,kPhoton);
  initializeRings(iNumberOfRings, fRingSize);
}


//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initializeRings(int iNumberOfRings, float fRingSize){
 
  setRingSize(fRingSize);
  setNumbersOfRings(iNumberOfRings);
 
  fIsolationInRings.clear();
  for(int isoBin =0;isoBin<iNumberOfRings;isoBin++){
    float fTemp = 0.0;
    fIsolationInRings.push_back(fTemp);
    
    float fTempPhoton = 0.0;
    fIsolationInRingsPhoton.push_back(fTempPhoton);

    float fTempNeutral = 0.0;
    fIsolationInRingsNeutral.push_back(fTempNeutral);

    float fTempCharged = 0.0;
    fIsolationInRingsCharged.push_back(fTempCharged);

    float fTempChargedAll = 0.0;
    fIsolationInRingsChargedAll.push_back(fTempChargedAll);

  }

  fConeSize = fRingSize * (float)iNumberOfRings;

}
  
 
//--------------------------------------------------------------------------------------------------
float PFIsolationEstimator::fGetIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection > vertices) {
 
  fGetIsolationInRings( pfCandidate, pfParticlesColl, vertices);
  fIsolation = fIsolationInRings[0];
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
vector<float >  PFIsolationEstimator::fGetIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, edm::Handle< reco::VertexCollection > vertices) {

  int isoBin;

  
  for(isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]=0.;
    fIsolationInRingsPhoton[isoBin]=0.;
    fIsolationInRingsNeutral[isoBin]=0.;
    fIsolationInRingsCharged[isoBin]=0.;
    fIsolationInRingsChargedAll[isoBin]=0.;
  }
  
 
  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {

    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 

    if(pfParticle.pdgId()==22){
      
      if(isPhotonParticleVetoed(  pfCandidate, &pfParticle)>=0.){
	isoBin = (int)(fDeltaR*10.);
	isoBin = 0;
	fIsolationInRingsPhoton[isoBin]  = fIsolationInRingsPhoton[isoBin] + pfParticle.pt();
      }
    }else if(abs(pfParticle.pdgId())==130){
        
      if(isNeutralParticleVetoed(  pfCandidate, &pfParticle)>=0.){
	
	isoBin = (int)(fDeltaR*10.);
	isoBin = 0;
	//	cout<<fDeltaR<<" : "<<isoBin<<endl;
	fIsolationInRingsNeutral[isoBin]  = fIsolationInRingsNeutral[isoBin] + pfParticle.pt();
      }

      //    }else if(abs(pfParticle.pdgId()) == 11 ||abs(pfParticle.pdgId()) == 13 || abs(pfParticle.pdgId()) == 211){
    }else if(abs(pfParticle.pdgId()) == 211){
      if(isChargedParticleVetoed(  pfCandidate, &pfParticle, vertices)>=0.){
	//isoBin = (int)(fDeltaR*10.);
	isoBin = 0;
	fIsolationInRingsCharged[isoBin]  = fIsolationInRingsCharged[isoBin] + pfParticle.pt();
      }
      
    }
  }

 
  for(int isoBin =0;isoBin<iNumberOfRings;isoBin++){
    fIsolationInRings[isoBin]= fIsolationInRingsPhoton[isoBin]+ fIsolationInRingsNeutral[isoBin] +  fIsolationInRingsCharged[isoBin];
  }

  return fIsolationInRings;
}



//--------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isPhotonParticleVetoed(const reco::PFCandidate* pfcand , const reco::PFCandidate* pfIsoCand ){
  
  
  fDeltaR = deltaR(pfcand->eta(),pfcand->phi(),pfIsoCand->eta(),pfIsoCand->phi()); 



  if(fDeltaR > fConeSize)
    return -999.;
  
  if(!bApplyVeto)
    return fDeltaR;
  

  fDeltaPhi = abs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = abs(pfcand->eta()-pfIsoCand->eta()); 

  if(abs(pfcand->eta())<1.44442){
    if(bDeltaRVetoBarrel){
      if(fDeltaR < fDeltaRVetoBarrelPhotons)
        return -999.;
    }
    
    if(bRectangleVetoBarrel){
      if(fDeltaEta < fRectangleDeltaEtaVetoBarrelPhotons && fDeltaPhi < fRectangleDeltaPhiVetoBarrelPhotons){
	return -999.;
      }
    }
  }else{
    if(bDeltaRVetoEndcap){
      if(fDeltaR < fDeltaRVetoEndcapPhotons)
	return -999.;
    }
    if(bRectangleVetoEndcap){
      if(fDeltaEta < fRectangleDeltaEtaVetoEndcapPhotons && fDeltaPhi < fRectangleDeltaPhiVetoEndcapPhotons){
	 return -999.;
      }
    }
  }

  return fDeltaR;
}

//--------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isNeutralParticleVetoed(const reco::PFCandidate* pfcand , const reco::PFCandidate* pfIsoCand ){

  fDeltaR = deltaR(pfcand->eta(),pfcand->phi(),pfIsoCand->eta(),pfIsoCand->phi()); 
  
  if(fDeltaR > fConeSize)
    return -999;
  
  if(!bApplyVeto)
    return fDeltaR;


  fDeltaPhi = abs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = abs(pfcand->eta()-pfIsoCand->eta()); 

  if(abs(pfcand->eta())<1.44442){
    if(!bDeltaRVetoBarrel&&!bRectangleVetoBarrel){
      return fDeltaR;
    }
    
    if(bDeltaRVetoBarrel){
	if(fDeltaR < fDeltaRVetoBarrelNeutrals)
	  return -999.;
      }
      if(bRectangleVetoBarrel){
	if(fDeltaEta < fRectangleDeltaEtaVetoBarrelNeutrals && fDeltaPhi < fRectangleDeltaPhiVetoBarrelNeutrals){
	    return -999.;
	}
      }
      
    }else{
     if(!bDeltaRVetoEndcap&&!bRectangleVetoEndcap){
       return fDeltaR;
     }
      if(bDeltaRVetoEndcap){
	if(fDeltaR < fDeltaRVetoEndcapNeutrals)
	  return -999.;
      }
      if(bRectangleVetoEndcap){
	if(fDeltaEta < fRectangleDeltaEtaVetoEndcapNeutrals && fDeltaPhi < fRectangleDeltaPhiVetoEndcapNeutrals){
	  return -999.;
	}
      }
  }

  return fDeltaR;
}


//--------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isChargedParticleVetoed(const reco::PFCandidate* pfcand , const reco::PFCandidate* pfIsoCand, edm::Handle< reco::VertexCollection >  vertices  ){
  
  
  VertexRef vtx = chargedHadronVertex(vertices,  *pfIsoCand );
  if(vtx.isNull())
    return -999.;
  
  if(iParticleType==kElectron){
    // math::XYZPoint pfParticleVtx(pfIsoCand->momentum());
    
    
    float dz = fabs(pfIsoCand->vz() -pfcand->vz());
    if (dz > 1.)
      return -999.;
    
    // double dxy = ( -(pfIsoCand->vx() - pfcand->vx())*pfIsoCand->py() + (pfIsoCand->vy() - pfcand->vy())*pfIsoCand->px()) / pfIsoCand->pt();
     double dxy = ( -(vtx->x() - pfcand->vx())*pfIsoCand->py() + (vtx->y() - pfcand->vy())*pfIsoCand->px()) / pfIsoCand->pt();

    if(fabs(dxy) > 0.1)
      return -999.;

    fDeltaR = deltaR(pfIsoCand->eta(),pfIsoCand->phi(),pfcand->eta(),pfcand->phi()); 

    
  }else{
    /*
    math::XYZVector photonWrtVertex(pfcand->superClusterRef()->x() - vtx->x(),
				    pfcand->superClusterRef()->y() - vtx->y(),
				    pfcand->superClusterRef()->x() - vtx->z());
    */

    float dz = fabs(pfIsoCand->vz() -pfcand->vz());
    if (dz > 1.)
      return -999.;
    
    // double dxy = ( -(pfIsoCand->vx() - pfcand->vx())*pfIsoCand->py() + (pfIsoCand->vy() - pfcand->vy())*pfIsoCand->px()) / pfIsoCand->pt();

    double dxy = ( -(pfIsoCand->vx() - pfcand->vx())*pfIsoCand->py() + (pfIsoCand->vy() - pfcand->vy())*pfIsoCand->px()) / pfIsoCand->pt();
    if(fabs(dxy) > 0.1)
      return -999.;

    fDeltaR = deltaR(pfIsoCand->eta(),pfIsoCand->phi(),pfcand->eta(),pfcand->phi()); 

  }

    
  if(fDeltaR > fConeSize)
    return -999.;
  
  if(!bApplyVeto)
    return fDeltaR;
  
  fDeltaPhi = fabs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = fabs(pfcand->eta()-pfIsoCand->eta()); 
  
  
  if(abs(pfcand->eta())<1.44442){
    if(!bDeltaRVetoBarrel&&!bRectangleVetoBarrel){
      return fDeltaR;
    }
    
    if(bDeltaRVetoBarrel){
	if(fDeltaR < fDeltaRVetoBarrelCharged)
	  return -999.;
      }
      if(bRectangleVetoBarrel){
	if(fDeltaEta < fRectangleDeltaEtaVetoBarrelCharged && fDeltaPhi < fRectangleDeltaPhiVetoBarrelCharged){
	    return -999.;
	}
      }
      
    }else{
     if(!bDeltaRVetoEndcap&&!bRectangleVetoEndcap){
       return fDeltaR;
     }
      if(bDeltaRVetoEndcap){
	if(fDeltaR < fDeltaRVetoEndcapCharged)
	  return -999.;
      }
      if(bRectangleVetoEndcap){
	if(fDeltaEta < fRectangleDeltaEtaVetoEndcapCharged && fDeltaPhi < fRectangleDeltaPhiVetoEndcapCharged){
	  return -999.;
	}
      }
  }
		   
  

  return fDeltaR;
}


//--------------------------------------------------------------------------------------------------
 VertexRef  PFIsolationEstimator::chargedHadronVertex(  edm::Handle< reco::VertexCollection > verticesColl, const reco::PFCandidate& pfcand ){

  //code copied from Florian's PFNoPU class
    
  reco::TrackBaseRef trackBaseRef( pfcand.trackRef() );

  size_t  iVertex = 0;
  unsigned index=0;
  unsigned nFoundVertex = 0;

  float bestweight=0;
  
  const reco::VertexCollection& vertices = *(verticesColl.product());

  for( reco::VertexCollection::const_iterator iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {
    
    const reco::Vertex& vtx = *iv;
    
    // loop on tracks in vertices
    for(reco::Vertex::trackRef_iterator iTrack=vtx.tracks_begin();iTrack!=vtx.tracks_end(); ++iTrack) {

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
    return  VertexRef( verticesColl, iVertex);
  }
  // no vertex found with this track. 

  // optional: as a secondary solution, associate the closest vertex in z
  if ( checkClosestZVertex ) {

    double dzmin = 10000.;
    double ztrack = pfcand.vertex().z();
    bool foundVertex = false;
    index = 0;
    for( reco::VertexCollection::const_iterator  iv=vertices.begin(); iv!=vertices.end(); ++iv, ++index) {

      double dz = fabs(ztrack - iv->z());
      if(dz<dzmin) {
        dzmin = dz;
        iVertex = index;
        foundVertex = true;
      }
    }

    if( foundVertex ) 
      return  VertexRef( verticesColl, iVertex);  
  
  }
   
  return  VertexRef( );
}


