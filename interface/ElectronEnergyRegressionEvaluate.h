//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronEnergyRegressionEvaluate
//
// Helper Class for applying electron energy regression calculation
//
// Authors: A.Takeda, S.Xie
//--------------------------------------------------------------------------------------------------


/// --> NOTE if you want to use this class as standalone without the CMSSW part 
///  you need to uncomment the below line and compile normally with scramv1 b 
///  Then you need just to load it in your root macro the lib with the correct path, eg:
///  gSystem->Load("/data/benedet/CMSSW_5_2_2/lib/slc5_amd64_gcc462/pluginEGammaEGammaAnalysisTools.so");

//#define STANDALONE   // <---- this line


#ifndef ELECTRONENERGYREGRESSIONEVALUATE_H
#define ELECTRONENERGYREGRESSIONEVALUATE_H

#include "TFile.h"
#include "TTree.h"

// For applying regression
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#ifndef STANDALONE
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#endif


class ElectronEnergyRegressionEvaluate{
  public:
    ElectronEnergyRegressionEvaluate();
    ~ElectronEnergyRegressionEvaluate();

    enum ElectronEnergyRegressionType {
      kNoTrkVar,
      kWithTrkVar,
      kNoTrkVarTwoPtBins,
      kWithTrkVarTwoPtBins
    };

    void initialize(std::string weightsFile,
                    ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType type);

    bool isInitialized() const {return fIsInitialized;}
                
#ifndef STANDALONE
    double calculateRegressionEnergy(const reco::GsfElectron *ele, 
                                     EcalClusterLazyTools &myEcalCluster, 
                                     const edm::EventSetup &setup,
                                     double rho, double nvertices, 
                                     bool printDebug = false);
#endif

    // Evaluates regression without tracker variables
    double regressionValueNoTrkVar(
      double SCRawEnergy,
      double scEta,
      double scPhi,
      double R9,
      double etawidth,
      double phiwidth,
      double NClusters,
      double HoE,
      double rho,
      double vertices,
      double EtaSeed,
      double PhiSeed,
      double ESeed,
      double E3x3Seed,
      double E5x5Seed,
      double see,
      double spp,
      double sep,
      double EMaxSeed,
      double E2ndSeed,
      double ETopSeed,
      double EBottomSeed,
      double ELeftSeed,
      double ERightSeed,
      double E2x5MaxSeed,
      double E2x5TopSeed,
      double E2x5BottomSeed,
      double E2x5LeftSeed,
      double E2x5RightSeed,
      double pt,
      double IEtaSeed,
      double IPhiSeed,
      double EtaCrySeed,
      double PhiCrySeed,
      double PreShowerOverRaw,
      bool printDebug = false);

    // Evaluates regression using tracker variables
    double regressionValueWithTrkVar(				
      double SCRawEnergy,
      double scEta,
      double scPhi,
      double R9,
      double etawidth,
      double phiwidth,
      double NClusters,
      double HoE,
      double rho,
      double vertices,
      double EtaSeed,
      double PhiSeed,
      double ESeed,
      double E3x3Seed,
      double E5x5Seed,
      double see,
      double spp,
      double sep,
      double EMaxSeed,
      double E2ndSeed,
      double ETopSeed,
      double EBottomSeed,
      double ELeftSeed,
      double ERightSeed,
      double E2x5MaxSeed,
      double E2x5TopSeed,
      double E2x5BottomSeed,
      double E2x5LeftSeed,
      double E2x5RightSeed,
      double pt,
      double GsfTrackPIn,
      double fbrem,
      double Charge,
      double EoP,
      double IEtaSeed,
      double IPhiSeed,
      double EtaCrySeed,
      double PhiCrySeed,
      double PreShowerOverRaw,
      bool printDebug = false );

  private:
    bool fIsInitialized;
    ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType fVersionType;
    GBRForest *forest_eb;		// Pointer to the GBRForest for barrel
    GBRForest *forest_ee;		// Pointer to the GBRForest for endcap
    GBRForest *forest_lowPt_eb;	// Pointer to the GBRForest for barrel, low Pt
    GBRForest *forest_lowPt_ee;	// Pointer to the GBRForest for endcap, low Pt
    GBRForest *forest_highPt_eb;	// Pointer to the GBRForest for barrel, high Pt
    GBRForest *forest_highPt_ee;	// Pointer to the GBRForest for endcap, high Pt


};

#endif
