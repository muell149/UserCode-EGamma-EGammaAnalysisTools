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
  TFile *fOutput = new TFile("IsolationMaps.root","recreate");
  
  hDeltaEtaDeltaPhiPhoton[0]->Write();
  hDeltaEtaDeltaPhiPhotonPtWeighted[0]->Write();
  hDeltaEtaDeltaPhiNeutral[0]->Write();
  hDeltaEtaDeltaPhiNeutralPtWeighted[0]->Write();
  hDeltaEtaDeltaPhiCharged[0]->Write();
  hDeltaEtaDeltaPhiChargedPtWeighted[0]->Write();


  hDeltaRPhoton[0]->Write();
  hDeltaRPhotonPtWeighted[0]->Write();
  hDeltaRNeutral[0]->Write();
  hDeltaRNeutralPtWeighted[0]->Write();
  hDeltaRCharged[0]->Write();
  hDeltaRChargedPtWeighted[0]->Write();

  hDeltaEtaDeltaPhiPhoton[1]->Write();
  hDeltaEtaDeltaPhiPhotonPtWeighted[1]->Write();
  hDeltaEtaDeltaPhiNeutral[1]->Write();
  hDeltaEtaDeltaPhiNeutralPtWeighted[1]->Write();
  hDeltaEtaDeltaPhiCharged[1]->Write();
  hDeltaEtaDeltaPhiChargedPtWeighted[1]->Write();


  hDeltaRPhoton[1]->Write();
  hDeltaRPhotonPtWeighted[1]->Write();
  hDeltaRNeutral[1]->Write();
  hDeltaRNeutralPtWeighted[1]->Write();
  hDeltaRCharged[1]->Write();
  hDeltaRChargedPtWeighted[1]->Write();

  fOutput->Close();

}

//--------------------------------------------------------------------------------------------------
void PFIsolationEstimator::initialize( Bool_t  bApplyVeto, int iParticleType ) {

  setParticleType(iParticleType);

  //By default check for an option vertex association
  checkClosestZVertex = kTRUE;
  
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
    setDeltaRVetoBarrel(kFALSE);
    setDeltaRVetoEndcap(kFALSE);
    setRectangleVetoBarrel(kFALSE);
    setRectangleVetoEndcap(kFALSE);

  }
 
  hDeltaEtaDeltaPhiPhoton[0] = new TH2F("hDeltaEtaDeltaPhiPhoton","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiPhotonPtWeighted[0] = new TH2F("hDeltaEtaDeltaPhiPhotonPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiNeutral[0] = new TH2F("hDeltaEtaDeltaPhiNeutral","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiNeutralPtWeighted[0] = new TH2F("hDeltaEtaDeltaPhiNeutralPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiCharged[0] = new TH2F("hDeltaEtaDeltaPhiCharged","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiChargedPtWeighted[0] = new TH2F("hDeltaEtaDeltaPhiChargedPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);


  hDeltaRPhoton[0] = new TH1F("hDeltaRPhoton","",40,0,1.0);
  hDeltaRPhotonPtWeighted[0] = new TH1F("hDeltaRPhotonPtWeighted","",40,0,1.0);
  hDeltaRNeutral[0] = new TH1F("hDeltaRNeutral","",40,0,1.0);
  hDeltaRNeutralPtWeighted[0] = new TH1F("hDeltaRNeutralPtWeighted","",40,0,1.0);
  hDeltaRCharged[0] = new TH1F("hDeltaRCharged","",40,0,1.0);
  hDeltaRChargedPtWeighted[0] = new TH1F("hDeltaRChargedPtWeighted","",40,0,1.0);


  hDeltaEtaDeltaPhiPhoton[1] = new TH2F("hEEDeltaEtaDeltaPhiPhoton","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiPhotonPtWeighted[1] = new TH2F("hEEDeltaEtaDeltaPhiPhotonPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiNeutral[1] = new TH2F("hEEDeltaEtaDeltaPhiNeutral","",80,-1.0,1.0,80,-1.0,1.0);
  hDeltaEtaDeltaPhiNeutralPtWeighted[1] = new TH2F("hEEDeltaEtaDeltaPhiNeutralPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiCharged[1] = new TH2F("hEEDeltaEtaDeltaPhiCharged","",100,-1.0,1.0,100,-1.0,1.0);
  hDeltaEtaDeltaPhiChargedPtWeighted[1] = new TH2F("hEEDeltaEtaDeltaPhiChargedPtWeighted","",100,-1.0,1.0,100,-1.0,1.0);


  hDeltaRPhoton[1] = new TH1F("hEEDeltaRPhoton","",40,0,1.0);
  hDeltaRPhotonPtWeighted[1] = new TH1F("hEEDeltaRPhotonPtWeighted","",40,0,1.0);
  hDeltaRNeutral[1] = new TH1F("hEEDeltaRNeutral","",40,0,1.0);
  hDeltaRNeutralPtWeighted[1] = new TH1F("hEEDeltaRNeutralPtWeighted","",40,0,1.0);
  hDeltaRCharged[1] = new TH1F("hEEDeltaRCharged","",40,0,1.0);
  hDeltaRChargedPtWeighted[1] = new TH1F("hEEDeltaRChargedPtWeighted","",40,0,1.0);

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
float PFIsolationEstimator::fGetIsolation(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices,bool kDoPlots) {
 
  fGetIsolationInRings( pfCandidate, pfParticlesColl, vtx, vertices);
  fIsolation = fIsolationInRings[0];
  
  return fIsolation;
}


//--------------------------------------------------------------------------------------------------
vector<float >  PFIsolationEstimator::fGetIsolationInRings(const reco::PFCandidate * pfCandidate, const reco::PFCandidateCollection* pfParticlesColl, reco::Vertex& vtx, edm::Handle< reco::VertexCollection > vertices,bool kDoPlots) {

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

    if(&pfParticle==(pfCandidate))
      continue;

    if(pfParticle.pdgId()==22){
      
      if(isPhotonParticleVetoed(  pfCandidate, &pfParticle)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsPhoton[isoBin]  = fIsolationInRingsPhoton[isoBin] + pfParticle.pt();
      }

      if(kDoPlots){
	 if(abs(pfParticle.eta())<1.44442){
	   hDeltaEtaDeltaPhiPhotonPtWeighted[0]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRPhotonPtWeighted[0]->Fill(fDeltaR,pfParticle.pt());
	    hDeltaEtaDeltaPhiPhoton[0]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRPhoton[0]->Fill(fDeltaR);
	 }else{
	   hDeltaEtaDeltaPhiPhotonPtWeighted[1]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRPhotonPtWeighted[1]->Fill(fDeltaR,pfParticle.pt());
	   hDeltaEtaDeltaPhiPhoton[1]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRPhoton[1]->Fill(fDeltaR);
	 }
  
      }
      
    }else if(abs(pfParticle.pdgId())==130){
        
      if(isNeutralParticleVetoed(  pfCandidate, &pfParticle)>=0.){
       	isoBin = (int)(fDeltaR/fRingSize);
	cout<<isoBin<<endl;
	fIsolationInRingsNeutral[isoBin]  = fIsolationInRingsNeutral[isoBin] + pfParticle.pt();
      }
      if(kDoPlots){
	 if(abs(pfParticle.eta())<1.44442){
	   hDeltaEtaDeltaPhiNeutralPtWeighted[0]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRNeutralPtWeighted[0]->Fill(fDeltaR,pfParticle.pt());
	   hDeltaEtaDeltaPhiNeutral[0]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRNeutral[0]->Fill(fDeltaR);
	 }else{
	   hDeltaEtaDeltaPhiNeutralPtWeighted[1]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRNeutralPtWeighted[1]->Fill(fDeltaR,pfParticle.pt());
	   hDeltaEtaDeltaPhiNeutral[1]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRNeutral[1]->Fill(fDeltaR);
	 }
  
      }

      //}else if(abs(pfParticle.pdgId()) == 11 ||abs(pfParticle.pdgId()) == 13 || abs(pfParticle.pdgId()) == 211){
    }else if(abs(pfParticle.pdgId()) == 211){
      if(isChargedParticleVetoed(  pfCandidate, &pfParticle, vtx, vertices)>=0.){
	isoBin = (int)(fDeltaR/fRingSize);
	fIsolationInRingsCharged[isoBin]  = fIsolationInRingsCharged[isoBin] + pfParticle.pt();
      }
      if(kDoPlots){
	 if(abs(pfParticle.eta())<1.44442){
	   hDeltaEtaDeltaPhiChargedPtWeighted[0]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRChargedPtWeighted[0]->Fill(fDeltaR,pfParticle.pt());
	   hDeltaEtaDeltaPhiCharged[0]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRCharged[0]->Fill(fDeltaR);
	 }else{
	   hDeltaEtaDeltaPhiChargedPtWeighted[1]->Fill(fDeltaPhi,fDeltaEta,pfParticle.pt());
	   hDeltaRChargedPtWeighted[1]->Fill(fDeltaR,pfParticle.pt());
	   hDeltaEtaDeltaPhiCharged[1]->Fill(fDeltaPhi,fDeltaEta);
	   hDeltaRCharged[1]->Fill(fDeltaR);
	 }
  
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
  
  fDeltaPhi = abs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = abs(pfcand->eta()-pfIsoCand->eta()); 

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
  
  fDeltaPhi = abs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = abs(pfcand->eta()-pfIsoCand->eta()); 

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


//----------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isChargedParticleVetoed(const reco::PFCandidate* pfcand , const reco::PFCandidate* pfIsoCand, edm::Handle< reco::VertexCollection >  vertices  ){
  //need code to handle special conditions
  
  return -999;
}

//-----------------------------------------------------------------------------------------------------
float  PFIsolationEstimator::isChargedParticleVetoed(const reco::PFCandidate* pfcand , const reco::PFCandidate* pfIsoCand, reco::Vertex& vtxMain, edm::Handle< reco::VertexCollection >  vertices  ){
  

  
  
  VertexRef vtx = chargedHadronVertex(vertices,  *pfIsoCand );
  if(vtx.isNull())
    return -999.;
  
  float fVtxMainX = vtxMain.x();
  float fVtxMainY = vtxMain.y();
  float fVtxMainZ = vtxMain.z();

  if(iParticleType==kPhoton){

    
    //this piece of code does not use the chargedhadronvertex function. 
    /*  float dz = fabs(pfIsoCand->vz() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
    double dxy = ( -(pfIsoCand->vx() - fVtxMainX)*pfIsoCand->py() + (pfIsoCand->vy() - fVtxMainY)*pfIsoCand->px()) / pfIsoCand->pt();

    if(fabs(dxy) > 0.1)
      return -999.;
    */
    
    float dz = fabs(vtx->z() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
   
    double dxy = ( -(vtx->x() - fVtxMainX)*pfIsoCand->py() + (vtx->y() - fVtxMainY)*pfIsoCand->px()) / pfIsoCand->pt();

    if(fabs(dxy) > 0.1)
      return -999.;
    
  }else{
  

    float dz = fabs(vtx->z() - fVtxMainZ);
    if (dz > 1.)
      return -999.;
    
    double dxy = ( -(vtx->x() - pfcand->vx())*pfIsoCand->py() + (vtx->y() - pfcand->vy())*pfIsoCand->px()) / pfIsoCand->pt();
    if(fabs(dxy) > 0.1)
      return -999.;
  }
    
  fDeltaR = deltaR(pfIsoCand->eta(),pfIsoCand->phi(),pfcand->eta(),pfcand->phi()); 

  if(fDeltaR > fConeSize)
    return -999.;

  fDeltaPhi = abs(deltaPhi(pfcand->phi(),pfIsoCand->phi())); 
  fDeltaEta = abs(pfcand->eta()-pfIsoCand->eta()); 
  
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


