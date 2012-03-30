// -*- C++ -*-
//
// Package:    ElectronAnalyzer
// Class:      ElectronAnalyzer
// 
/**\class ElectronAnalyzer

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Daniele Benedetti



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "RecoParticleFlow/PFProducer/interface/Utils.h"

#include "../interface/ElectronMVAEstimator.h"



#include <vector>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
//
// class decleration
//

using namespace edm;
using namespace reco;
using namespace std;
class ElectronAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ElectronAnalyzer(const edm::ParameterSet&);
      ~ElectronAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  
  ParameterSet conf_;

  //  ElectronMVAEstimator *fMVASiDanV2;
  ElectronMVAEstimator* fMVASiDanV2;

  unsigned int ev;
      // ----------member data ---------------------------

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig):
  conf_(iConfig)

{
   // = new ElectronMVAEstimator();
 //  fMVASiDanV2->initialize("BDTCat_BDTG_SiDanV2",
// 			  "/mnt/data2/OutputBatchTMVA/EmanueleV12/weights/DanieleMVA_BDTCat_BDTG_SiDanV2.weights.xml",
// 			  ElectronMVAEstimator::kNonTrig);
  
  edm::Service<TFileService> fs;

  //  h_mva_ele  = fs->make<TH1F>("h_mva_ele"," ",50,-1.1,1.1);

}


ElectronAnalyzer::~ElectronAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  InputTag egammaLabel(string("gsfElectrons"));
  Handle<GsfElectronCollection> theEGammaCollection;
  iEvent.getByLabel(egammaLabel,theEGammaCollection);
  const GsfElectronCollection theEGamma = *(theEGammaCollection.product());


  InputTag  MCTruthCollection(string("generator"));
  edm::Handle<edm::HepMCProduct> pMCTruth;
  iEvent.getByLabel(MCTruthCollection,pMCTruth);
  const HepMC::GenEvent* genEvent = pMCTruth->GetEvent();



  bool debug = true;


  ev++;

  if(debug)
    cout << "************************* New Event:: " << ev << " *************************" << endl;

  // Validation from generator events 

  for(HepMC::GenEvent::particle_const_iterator cP = genEvent->particles_begin(); 
      cP != genEvent->particles_end(); cP++ ) {

    float etamc= (*cP)->momentum().eta();
    float phimc= (*cP)->momentum().phi();
    float ptmc = (*cP)->momentum().perp();


    if(abs((*cP)->pdg_id())==11 && 
       (*cP)->status()==1       &&
       ptmc > 5.                && 
       fabs(etamc) < 2.5 ){
 
      
      for (uint j=0; j<theEGamma.size();j++) {
	float etareco = theEGamma[j].eta();
	float phireco = theEGamma[j].phi();
	float deta = etamc - etareco;
	float dphi = Utils::mpi_pi(phimc - phireco);
	float dR = sqrt(deta*deta + dphi*dphi);

	if(dR < 0.1) {
	  if(debug)
	    cout << " niente " << endl;

	} 
      } // End Loop on RECO electrons
    } // End if MC electrons selection
  } //End Loop Generator Particles 

}
// ------------ method called once each job just before starting event loop  ------------
void 
ElectronAnalyzer::beginJob(const edm::EventSetup&)
{

  ev = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronAnalyzer::endJob() {
  cout << " endJob:: #events " << ev << endl;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ElectronAnalyzer);
