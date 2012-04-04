#include <TFile.h>
#include "../interface/ElectronMVAEstimator.h"
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
ElectronMVAEstimator::ElectronMVAEstimator() :
fMethodname("BDTG method"),
fisInitialized(kFALSE)
{
  // Constructor.
  fTMVAReader = 0;
}



//--------------------------------------------------------------------------------------------------
ElectronMVAEstimator::~ElectronMVAEstimator()
{
    if (fTMVAReader) delete fTMVAReader;
}

//--------------------------------------------------------------------------------------------------
void ElectronMVAEstimator::initialize( std::string methodName,
				       std::string weightsfile,
				       ElectronMVAEstimator::MVAType type) {
  
  fisInitialized = kTRUE;
  fMVAType = type;

  fMethodname = methodName;
    
  if (fTMVAReader) delete fTMVAReader;
  
  fTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
  fTMVAReader->SetVerbose(kTRUE);
  
  if (type == kTrig) {
    // Pure tracking variables
    fTMVAReader->AddVariable("fbrem",           &fMVAVar_fbrem);
    fTMVAReader->AddVariable("kfchi2",          &fMVAVar_kfchi2);
    fTMVAReader->AddVariable("kfhits",          &fMVAVar_kfhits);
    fTMVAReader->AddVariable("gsfchi2",         &fMVAVar_gsfchi2);

    // Geometrical matchings
    fTMVAReader->AddVariable("deta",            &fMVAVar_deta);
    fTMVAReader->AddVariable("dphi",            &fMVAVar_dphi);
    fTMVAReader->AddVariable("detacalo",        &fMVAVar_detacalo);
    // fTMVAReader->AddVariable("dphicalo",        &fMVAVar_dphicalo);   // Pruned but save in your ntuple. 
    
    // Pure ECAL -> shower shapes
    fTMVAReader->AddVariable("see",             &fMVAVar_see);
    fTMVAReader->AddVariable("spp",             &fMVAVar_spp);
    fTMVAReader->AddVariable("etawidth",        &fMVAVar_etawidth);
    fTMVAReader->AddVariable("phiwidth",        &fMVAVar_phiwidth);
    fTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_e1x5e5x5);
    fTMVAReader->AddVariable("R9",              &fMVAVar_R9);
    // fTMVAReader->AddVariable("nbrems",          &fMVAVar_nbrems); // Pruned but save in your ntuple. 

    // Energy matching
    fTMVAReader->AddVariable("HoE",             &fMVAVar_HoE);
    fTMVAReader->AddVariable("EoP",             &fMVAVar_EoP); 
    fTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_IoEmIoP);
    fTMVAReader->AddVariable("eleEoPout",       &fMVAVar_eleEoPout);
    //  fTMVAReader->AddVariable("EoPout",          &fMVAVar_EoPout); // Pruned but save in your ntuple.    
    fTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);

    // IP
    fTMVAReader->AddVariable("d0",              &fMVAVar_d0);
    fTMVAReader->AddVariable("ip3d",            &fMVAVar_ip3d);
    
    fTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
    fTMVAReader->AddSpectator("pt",             &fMVAVar_pt);
  }
  
  if (type == kNonTrig) {
    // Pure tracking variables
    fTMVAReader->AddVariable("fbrem",           &fMVAVar_fbrem);
    fTMVAReader->AddVariable("kfchi2",          &fMVAVar_kfchi2);
    fTMVAReader->AddVariable("kfhits",          &fMVAVar_kfhits);
    fTMVAReader->AddVariable("gsfchi2",         &fMVAVar_gsfchi2);

    // Geometrical matchings
    fTMVAReader->AddVariable("deta",            &fMVAVar_deta);
    fTMVAReader->AddVariable("dphi",            &fMVAVar_dphi);
    fTMVAReader->AddVariable("detacalo",        &fMVAVar_detacalo);
    // fTMVAReader->AddVariable("dphicalo",        &fMVAVar_dphicalo);   // Pruned but save in your ntuple. 
    
    // Pure ECAL -> shower shapes
    fTMVAReader->AddVariable("see",             &fMVAVar_see);
    fTMVAReader->AddVariable("spp",             &fMVAVar_spp);
    fTMVAReader->AddVariable("etawidth",        &fMVAVar_etawidth);
    fTMVAReader->AddVariable("phiwidth",        &fMVAVar_phiwidth);
    fTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_e1x5e5x5);
    fTMVAReader->AddVariable("R9",              &fMVAVar_R9);
    // fTMVAReader->AddVariable("nbrems",          &fMVAVar_nbrems); // Pruned but save in your ntuple. 

    // Energy matching
    fTMVAReader->AddVariable("HoE",             &fMVAVar_HoE);
    fTMVAReader->AddVariable("EoP",             &fMVAVar_EoP); 
    fTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_IoEmIoP);
    fTMVAReader->AddVariable("eleEoPout",       &fMVAVar_eleEoPout);
    //  fTMVAReader->AddVariable("EoPout",          &fMVAVar_EoPout); // Pruned but save in your ntuple. 
    fTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);
    
    fTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
    fTMVAReader->AddSpectator("pt",             &fMVAVar_pt);
  }
  
  fTMVAReader->BookMVA(fMethodname , weightsfile);
  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << weightsfile << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t ElectronMVAEstimator::mvaValue(Double_t fbrem, 
					Double_t kfchi2,
					Int_t    kfhits,
					Double_t gsfchi2,
					Double_t deta,
					Double_t dphi,
					Double_t detacalo,
					//Double_t dphicalo,
					Double_t see,
					Double_t spp,
					Double_t etawidth,
					Double_t phiwidth,
					Double_t e1x5e5x5,
					Double_t R9,
					//Int_t    nbrems,
					Double_t HoE,
					Double_t EoP,
					Double_t IoEmIoP,
					Double_t eleEoPout,
					Double_t PreShowerOverRaw,
					//Double_t EoPout,
					Double_t d0,
					Double_t ip3d,
					Double_t eta,
					Double_t pt,
					Bool_t printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: ElectronMVAEstimator not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_fbrem           = fbrem; 
  fMVAVar_kfchi2          = kfchi2;
  fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
  fMVAVar_gsfchi2         = gsfchi2;

  fMVAVar_deta            = deta;
  fMVAVar_dphi            = dphi;
  fMVAVar_detacalo        = detacalo;
  // fMVAVar_dphicalo        = dphicalo;


  fMVAVar_see             = see;
  fMVAVar_spp             = spp;
  fMVAVar_etawidth        = etawidth;
  fMVAVar_phiwidth        = phiwidth;
  fMVAVar_e1x5e5x5        = e1x5e5x5;
  fMVAVar_R9              = R9;
  //fMVAVar_nbrems          = float(nbrems);   // BTD does not support int variables


  fMVAVar_HoE             = HoE;
  fMVAVar_EoP             = EoP;
  fMVAVar_IoEmIoP         = IoEmIoP;
  fMVAVar_eleEoPout       = eleEoPout;
  fMVAVar_PreShowerOverRaw= PreShowerOverRaw;
  //fMVAVar_EoPout          = EoPout; 

  fMVAVar_d0              = d0;
  fMVAVar_ip3d            = ip3d;
  fMVAVar_eta             = eta;
  fMVAVar_pt              = pt;


  bindVariables();
  Double_t mva = -9999;  
  mva = fTMVAReader->EvaluateMVA(fMethodname);



  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << endl;
    cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
      // << " dphicalo " << fMVAVar_dphicalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
	 << " R9 " << fMVAVar_R9  
      // << " mynbrems " << fMVAVar_nbrems  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
      //<< " EoPout " << fMVAVar_EoPout  
	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
	 << " d0 " << fMVAVar_d0  
	 << " ip3d " << fMVAVar_ip3d  
	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }


  return mva;
}
//--------------------------------------------------------------------------------------------------
Double_t ElectronMVAEstimator::mvaValue(Double_t fbrem, 
					Double_t kfchi2,
					Int_t    kfhits,
					Double_t gsfchi2,
					Double_t deta,
					Double_t dphi,
					Double_t detacalo,
					//Double_t dphicalo,
					Double_t see,
					Double_t spp,
					Double_t etawidth,
					Double_t phiwidth,
					Double_t e1x5e5x5,
					Double_t R9,
					//Int_t    nbrems,
					Double_t HoE,
					Double_t EoP,
					Double_t IoEmIoP,
					Double_t eleEoPout,
					Double_t PreShowerOverRaw,
					//Double_t EoPout,
					Double_t eta,
					Double_t pt,
					Bool_t printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: ElectronMVAEstimator not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_fbrem           = fbrem; 
  fMVAVar_kfchi2          = kfchi2;
  fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
  fMVAVar_gsfchi2         = gsfchi2;

  fMVAVar_deta            = deta;
  fMVAVar_dphi            = dphi;
  fMVAVar_detacalo        = detacalo;
  // fMVAVar_dphicalo        = dphicalo;


  fMVAVar_see             = see;
  fMVAVar_spp             = spp;
  fMVAVar_etawidth        = etawidth;
  fMVAVar_phiwidth        = phiwidth;
  fMVAVar_e1x5e5x5        = e1x5e5x5;
  fMVAVar_R9              = R9;
  //fMVAVar_nbrems          = float(nbrems);   // BTD does not support int variables


  fMVAVar_HoE             = HoE;
  fMVAVar_EoP             = EoP;
  fMVAVar_IoEmIoP         = IoEmIoP;
  fMVAVar_eleEoPout       = eleEoPout;
  fMVAVar_PreShowerOverRaw= PreShowerOverRaw;
  //fMVAVar_EoPout          = EoPout; 

  fMVAVar_eta             = eta;
  fMVAVar_pt              = pt;


  bindVariables();
  Double_t mva = -9999;  
  mva = fTMVAReader->EvaluateMVA(fMethodname);



  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << endl;
    cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
      // << " dphicalo " << fMVAVar_dphicalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
	 << " R9 " << fMVAVar_R9  
      // << " mynbrems " << fMVAVar_nbrems  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
      //<< " EoPout " << fMVAVar_EoPout  
	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }


  return mva;
}


//--------------------------------------------------------------------------------------------------
#ifndef STANDALONE
Double_t ElectronMVAEstimator::mvaValue(const reco::GsfElectron& ele, 
					const reco::Vertex& vertex, 
					const TransientTrackBuilder& transientTrackBuilder,					
					EcalClusterLazyTools myEcalCluster,
					bool printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: ElectronMVAEstimator not properly initialized.\n"; 
    return -9999;
  }
  
  bool validKF= false; 
  reco::TrackRef myTrackRef = ele.closestCtfTrackRef();
  validKF = (myTrackRef.isAvailable());
  validKF = (myTrackRef.isNonnull());  

  // Pure tracking variables
  fMVAVar_fbrem           =  ele.fbrem();
  fMVAVar_kfchi2          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
  fMVAVar_kfhits          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; 
  //  fMVAVar_kfhits          =  (validKF) ? myTrackRef->numberOfValidHits() : -1. ;   //  save also this in your ntuple as possible alternative
  fMVAVar_gsfchi2         =  ele.gsfTrack()->normalizedChi2();  

  
  // Geometrical matchings
  fMVAVar_deta            =  ele.deltaEtaSuperClusterTrackAtVtx();
  fMVAVar_dphi            =  ele.deltaPhiSuperClusterTrackAtVtx();
  fMVAVar_detacalo        =  ele.deltaEtaSeedClusterTrackAtCalo();
  // fMVAVar_dphicalo        =  ele.deltaPhiSeedClusterTrackAtCalo();   //  save also this in your ntuple 


  // Pure ECAL -> shower shapes
  fMVAVar_see             =  ele.sigmaIetaIeta();    //EleSigmaIEtaIEta
  std::vector<float> vCov = myEcalCluster.localCovariances(*(ele.superCluster()->seed())) ;
  if (!isnan(vCov[2])) fMVAVar_spp = sqrt (vCov[2]);   //EleSigmaIPhiIPhi
  else fMVAVar_spp = 0.;    
  // fMVAVar_sigmaIEtaIPhi = vCov[1];  //  save also this in your ntuple 

  fMVAVar_etawidth        =  ele.superCluster()->etaWidth();
  fMVAVar_phiwidth        =  ele.superCluster()->phiWidth();
  fMVAVar_e1x5e5x5        =  (ele.e5x5()) !=0. ? 1.-(ele.e1x5()/ele.e5x5()) : -1. ;
  fMVAVar_R9              =  myEcalCluster.e3x3(*(ele.superCluster()->seed())) / ele.superCluster()->rawEnergy();
  //fMVAVar_nbrems          =  fabs(ele.numberOfBrems());    //  save also this in your ntuple 

  // Energy matching
  fMVAVar_HoE             =  ele.hadronicOverEm();
  fMVAVar_EoP             =  ele.eSuperClusterOverP();
  fMVAVar_IoEmIoP         =  (1.0/(ele.superCluster()->energy())) - (1.0 / ele.p());  // in the future to be changed with ele.gsfTrack()->p()
  fMVAVar_eleEoPout       =  ele.eEleClusterOverPout();
  fMVAVar_PreShowerOverRaw=  ele.superCluster()->preshowerEnergy() / ele.superCluster()->rawEnergy();
  // fMVAVar_EoPout          =  ele.eSeedClusterOverPout();     //  save also this in your ntuple 


  // Spectators
  fMVAVar_eta             =  ele.superCluster()->eta();         
  fMVAVar_pt              =  ele.pt();                          

 

  // for triggering electrons get the impact parameteres
  if(fMVAType == kTrig) {
    //d0
    if (ele.gsfTrack().isNonnull()) {
      fMVAVar_d0 = (-1.0)*ele.gsfTrack()->dxy(vertex.position()); 
    } else if (ele.closestCtfTrackRef().isNonnull()) {
      fMVAVar_d0 = (-1.0)*ele.closestCtfTrackRef()->dxy(vertex.position()); 
    } else {
      fMVAVar_d0 = -9999.0;
    }
    
    //default values for IP3D
    fMVAVar_ip3d = -999.0; 
    // fMVAVar_ip3dSig = 0.0;
    if (ele.gsfTrack().isNonnull()) {
      const double gsfsign   = ( (-ele.gsfTrack()->dxy(vertex.position()))   >=0 ) ? 1. : -1.;
      
      const reco::TransientTrack &tt = transientTrackBuilder.build(ele.gsfTrack()); 
      const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,vertex);
      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	//double ip3derr = ip3dpv.second.error();  
	fMVAVar_ip3d = ip3d; 
	// fMVAVar_ip3dSig = ip3d/ip3derr;
      }
    }
  }
  

  // evaluate
  bindVariables();
  Double_t mva = -9999;  
  mva = fTMVAReader->EvaluateMVA(fMethodname);



  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << " fMVAType " << fMVAType << endl;
    cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
      // << " dphicalo " << fMVAVar_dphicalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_e1x5e5x5  
	 << " R9 " << fMVAVar_R9  
      // << " mynbrems " << fMVAVar_nbrems  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
      // << " EoPout " << fMVAVar_EoPout  
	 << " d0 " << fMVAVar_d0  
	 << " ip3d " << fMVAVar_ip3d  
	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }



  return mva;
}
#endif
void ElectronMVAEstimator::bindVariables() {

  // this binding is needed for variables that sometime diverge. 


  if(fMVAVar_fbrem < -1.)
    fMVAVar_fbrem = -1.;	
  
  fMVAVar_deta = fabs(fMVAVar_deta);
  if(fMVAVar_deta > 0.06)
    fMVAVar_deta = 0.06;
  
  
  fMVAVar_dphi = fabs(fMVAVar_dphi);
  if(fMVAVar_dphi > 0.6)
    fMVAVar_dphi = 0.6;
  
  
//   if(fMVAVar_EoPout > 20.)
//     fMVAVar_EoPout = 20.;
  
  if(fMVAVar_EoP > 20.)
    fMVAVar_EoP = 20.;
  
  if(fMVAVar_eleEoPout > 20.)
    fMVAVar_eleEoPout = 20.;
  
  
  fMVAVar_detacalo = fabs(fMVAVar_detacalo);
  if(fMVAVar_detacalo > 0.2)
    fMVAVar_detacalo = 0.2;
  
  
//   fMVAVar_dphicalo = fabs(fMVAVar_dphicalo);
//   if(fMVAVar_dphicalo > 0.4)
//     fMVAVar_dphicalo = 0.4;
  
  
  if(fMVAVar_e1x5e5x5 < -1.)
    fMVAVar_e1x5e5x5 = -1;
  
  if(fMVAVar_e1x5e5x5 > 2.)
    fMVAVar_e1x5e5x5 = 2.; 
  
  
  
  if(fMVAVar_R9 > 5)
    fMVAVar_R9 = 5;
  
  if(fMVAVar_gsfchi2 > 200.)
    fMVAVar_gsfchi2 = 200;
  
  
  if(fMVAVar_kfchi2 > 10.)
    fMVAVar_kfchi2 = 10.;
  
  
  // Needed for a bug in CMSSW_420, fixed in more recent CMSSW versions
  if(std::isnan(fMVAVar_spp))
    fMVAVar_spp = 0.;	
  
  
  return;
}








