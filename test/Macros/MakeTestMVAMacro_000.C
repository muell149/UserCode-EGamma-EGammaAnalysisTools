#define MakeTestMVAMacro_000_cxx
#include "MakeTestMVAMacro_000.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../../interface/ElectronMVAEstimator.h"


void MakeTestMVAMacro_000::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L MakeTestMVAMacro_000.C
//      Root > MakeTestMVAMacro_000 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


  if(sig_)
    cout << " I am processing the signal ... input file " << myfile_
	 << " ...  output file name " << myoutputfile_ << endl;
  else
    cout << " I am processing the background ... input file " << myfile_
	 << "... output file name " << myoutputfile_ << endl;
 
 

  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;
  cout << " * " << endl;

  TFile *file= new TFile(myoutputfile_.c_str(),"Recreate");
  TH1D *h_norm_lowbarrel = new  TH1D("h_norm_lowbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_norm_highbarrel = new  TH1D("h_norm_highbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_norm_endcap = new  TH1D("h_norm_endcap"," BDT Cat ",1,1.,2.);
  
  TH1D *h_lowet_norm_lowbarrel = new  TH1D("h_lowet_norm_lowbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_lowet_norm_highbarrel = new  TH1D("h_lowet_norm_highbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_lowet_norm_endcap = new  TH1D("h_lowet_norm_endcap"," BDT Cat ",1,1.,2.);

  TH1D *h_midet_norm_lowbarrel = new  TH1D("h_midet_norm_lowbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_midet_norm_highbarrel = new  TH1D("h_midet_norm_highbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_midet_norm_endcap = new  TH1D("h_midet_norm_endcap"," BDT Cat ",1,1.,2.);

  TH1D *h_highet_norm_lowbarrel = new  TH1D("h_highet_norm_lowbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_highet_norm_highbarrel = new  TH1D("h_highet_norm_highbarrel"," BDT Cat ",1,1.,2.);
  TH1D *h_highet_norm_endcap = new  TH1D("h_highet_norm_endcap"," BDT Cat ",1,1.,2.);


  TH1D *h_mvaNonTrigV0_lowbarrel = new  TH1D("h_mvaNonTrigV0_lowbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_mvaNonTrigV0_highbarrel = new  TH1D("h_mvaNonTrigV0_highbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_mvaNonTrigV0_endcap = new  TH1D("h_mvaNonTrigV0_endcap"," BDT ",400.,-1.1,1.1);

  TH1D *h_lowet_mvaNonTrigV0_lowbarrel = new  TH1D("h_lowet_mvaNonTrigV0_lowbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_lowet_mvaNonTrigV0_highbarrel = new  TH1D("h_lowet_mvaNonTrigV0_highbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_lowet_mvaNonTrigV0_endcap = new  TH1D("h_lowet_mvaNonTrigV0_endcap"," BDT ",400.,-1.1,1.1);

  TH1D *h_midet_mvaNonTrigV0_lowbarrel = new  TH1D("h_midet_mvaNonTrigV0_lowbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_midet_mvaNonTrigV0_highbarrel = new  TH1D("h_midet_mvaNonTrigV0_highbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_midet_mvaNonTrigV0_endcap = new  TH1D("h_midet_mvaNonTrigV0_endcap"," BDT ",400.,-1.1,1.1);

  TH1D *h_highet_mvaNonTrigV0_lowbarrel = new  TH1D("h_highet_mvaNonTrigV0_lowbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_highet_mvaNonTrigV0_highbarrel = new  TH1D("h_highet_mvaNonTrigV0_highbarrel"," BDT ",400.,-1.1,1.1);
  TH1D *h_highet_mvaNonTrigV0_endcap = new  TH1D("h_highet_mvaNonTrigV0_endcap"," BDT ",400.,-1.1,1.1);





  ElectronMVAEstimator *myMVANonTrigV0 = new ElectronMVAEstimator();
  myMVANonTrigV0->initialize("BDTCat_BDTG_NonTrigV0",
			   "/afs/cern.ch/cms/data/CMSSW/RecoEgamma/ElectronIdentification/data/Electrons_BDTGCat_NonTrigV0.weights.xml", 
			   ElectronMVAEstimator::kNonTrig);
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      if ((jentry < 10000 && jentry%1000 == 0) || 
	  jentry%10000 == 0 ) 
	cout<<"... entry "<< jentry << endl;
      
      if(matchConv == 1)
	continue;
      
    

      if(myMVANonTrigV0 == 0)
	cout << " MVA not init!!! " <<endl;
      
      bool print = false;
      
      double myMVAOut = myMVANonTrigV0->mvaValue(fbrem, 
						 kfchi2,
						 kfhits,
						 gsfchi2,
						 deta,
						 dphi,
						 detacalo,
						 // dphicalo,
						 see,
						 spp,
						 etawidth,
						 phiwidth,
						 e1x5e5x5,
						 R9,
						 //nbrems,
						 HoE,
						 EoP,
						 IoEmIoP,
						 eleEoPout,
						 PreShowerOverRaw,
						 // EoPout,
						 eta,
						 pt,
						 print);
      
      if(print) 
	cout << " MVA Standalone Class " << myMVAOut << endl;



      if(fabs(eta) <= 0.8) {
	h_norm_lowbarrel->Fill(1.);
	h_mvaNonTrigV0_lowbarrel->Fill(myMVAOut);
	
	if(pt < 10.) {
	  h_lowet_norm_lowbarrel->Fill(1.);
	  h_lowet_mvaNonTrigV0_lowbarrel->Fill(myMVAOut);
	}
	if(pt >= 10 && pt < 20.) {
	  h_midet_norm_lowbarrel->Fill(1.);
	  h_midet_mvaNonTrigV0_lowbarrel->Fill(myMVAOut);
	}
	else if(pt >= 20.) {
	  h_highet_norm_lowbarrel->Fill(1.);
	  h_highet_mvaNonTrigV0_lowbarrel->Fill(myMVAOut);
	}
	
      }
      else if( fabs(eta) > 0.8 && fabs(eta) <= 1.485) {
	h_norm_highbarrel->Fill(1.);
	h_mvaNonTrigV0_highbarrel->Fill(myMVAOut);

	if(pt < 10.) {
	  h_lowet_norm_highbarrel->Fill(1.);
	  h_lowet_mvaNonTrigV0_highbarrel->Fill(myMVAOut);
	}
	if(pt >= 10 && pt < 20.) {
	  h_midet_norm_highbarrel->Fill(1.);
	  h_midet_mvaNonTrigV0_highbarrel->Fill(myMVAOut);
	}
	else if(pt >= 20.) {
	  h_highet_norm_highbarrel->Fill(1.);
	  h_highet_mvaNonTrigV0_highbarrel->Fill(myMVAOut);
	}
	
      }
      else if(  fabs(eta) > 1.485 ){	
	h_norm_endcap->Fill(1.);
	h_mvaNonTrigV0_endcap->Fill(myMVAOut);

	if(pt < 10.) {
	  h_lowet_norm_endcap->Fill(1.);
	  h_lowet_mvaNonTrigV0_endcap->Fill(myMVAOut);
	}
	if(pt >= 10 && pt < 20.) {
	  h_midet_norm_endcap->Fill(1.);
	  h_midet_mvaNonTrigV0_endcap->Fill(myMVAOut);
	}
	else if(pt >= 20.) {
	  h_highet_norm_endcap->Fill(1.);
	  h_highet_mvaNonTrigV0_endcap->Fill(myMVAOut);
	}
      }


   }

   file->cd();
   file->Write();
   file->Close();

}
