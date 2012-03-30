//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronMVAEstimator
//
// Helper Class for applying MVA electron ID selection
//
// Authors: D.Benedetti, E.DiMaro, S.Xie
//--------------------------------------------------------------------------------------------------


/// --> NOTE if you want to use this class as standalone without the CMSSW part 
///  you need to uncomment the below line and compile normally with scramv1 b 
///  Then you need just to load it in your root macro the lib with the correct path, eg:
///  gSystem->Load("/data/benedet/CMSSW_5_2_2/lib/slc5_amd64_gcc462/pluginEGammaEGammaAnalysisTools.so");

//#define STANDALONE   // <---- this line

#ifndef ElectronMVAEstimator_H
#define ElectronMVAEstimator_H

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

class ElectronMVAEstimator{
 public:
  ElectronMVAEstimator();
  ~ElectronMVAEstimator(); 
  
  enum MVAType {
    kTrig = 0,      // MVA for non-triggering electrons          
    kNonTrig              // MVA for triggering electrons
  };
  
  void     initialize(std::string methodName,
		      std::string weightsfile , 
		      ElectronMVAEstimator::MVAType type );
  
  Bool_t   isInitialized() const { return fisInitialized; }
  
  void SetPrintMVADebug(bool b) { fPrintMVADebug = b; }
  
  void bindVariables();

#ifndef STANDALONE
  Double_t mvaValue(const reco::GsfElectron *ele, 
		    const reco::Vertex vertex, 
		    const TransientTrackBuilder *transientTrackBuilder,
		    EcalClusterLazyTools myEcalCluster,
		    bool printDebug = kFALSE);
#endif
  
  Double_t mvaValue(Double_t fbrem, 
		    Double_t deta,
		    Double_t dphi,
		    Double_t see,
		    Double_t etawidth,
		    Double_t phiwidth,
		    Double_t HoE,
		    Double_t EoP,
		    Double_t e1x5e5x5,
		    Double_t EoPout,
		    Double_t eleEoPout,
		    Double_t detacalo,
		    Double_t kfchi2,
		    Int_t    kfhits,
		    Double_t spp,
		    Double_t IoEmIoP,
		    Int_t    nbrems,
		    Double_t R9,
		    Double_t dphicalo,
		    Double_t gsfchi2,
		    Double_t PreShowerOverRaw,
		    Double_t d0,
		    Double_t ip3d,
		    Double_t eta,
		    Double_t pt,
		    Bool_t printDebug = kFALSE );
 
  Double_t mvaValue(Double_t fbrem, 
		    Double_t deta,
		    Double_t dphi,
		    Double_t see,
		    Double_t etawidth,
		    Double_t phiwidth,
		    Double_t HoE,
		    Double_t EoP,
		    Double_t e1x5e5x5,
		    Double_t EoPout,
		    Double_t eleEoPout,
		    Double_t detacalo,
		    Double_t kfchi2,
		    Int_t    kfhits,
		    Double_t spp,
		    Double_t IoEmIoP,
		    Int_t    nbrems,
		    Double_t R9,
		    Double_t dphicalo,
		    Double_t gsfchi2,
		    Double_t PreShowerOverRaw,
		    Double_t eta,
		    Double_t pt,
		    Bool_t printDebug = kFALSE );
 

 
 private:
  TMVA::Reader             *fTMVAReader;
  std::string               fMethodname;
  MVAType                   fMVAType;
  
  Bool_t                    fisInitialized;
  Bool_t                    fPrintMVADebug;
  Float_t                   fMVAVar_fbrem; 
  Float_t                   fMVAVar_deta;
  Float_t                   fMVAVar_dphi;
  Float_t                   fMVAVar_see;
  Float_t                   fMVAVar_etawidth;
  Float_t                   fMVAVar_phiwidth;
  Float_t                   fMVAVar_HoE;
  Float_t                   fMVAVar_EoP;
  Float_t                   fMVAVar_e1x5e5x5;
  Float_t                   fMVAVar_EoPout;
  Float_t                   fMVAVar_eleEoPout;
  Float_t                   fMVAVar_detacalo;
  Float_t                   fMVAVar_kfchi2;
  Float_t                   fMVAVar_kfhits;
  Float_t                   fMVAVar_spp;
  Float_t                   fMVAVar_IoEmIoP;
  Float_t                   fMVAVar_nbrems;
  Float_t                   fMVAVar_R9;
  Float_t                   fMVAVar_dphicalo;
  Float_t                   fMVAVar_gsfchi2;
  Float_t                   fMVAVar_PreShowerOverRaw;
  Float_t                   fMVAVar_d0;
  Float_t                   fMVAVar_ip3d;
  Float_t                   fMVAVar_eta;
  Float_t                   fMVAVar_pt;
  Int_t                     fMVAVar_matchConv;
  

};

#endif
