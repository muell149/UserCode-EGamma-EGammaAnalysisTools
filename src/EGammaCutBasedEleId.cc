#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include <algorithm>

#ifndef STANDALONEID

bool EgammaCutBasedEleId::PassWP(WorkingPoint workingPoint,
    const reco::GsfElectronRef &ele,
    const edm::Handle<reco::ConversionCollection> &conversions,
    const math::XYZPoint &beamspot,
    const reco::Vertex &vertex,
    const IsoDepositVals &isoVals,
    const double &rho)
{
    unsigned int mask = TestWP(workingPoint, ele, conversions, beamspot, vertex, isoVals, rho);
    if ((mask & PassAll) == PassAll) return true;
    return false;
}

unsigned int EgammaCutBasedEleId::TestWP(WorkingPoint workingPoint,
    const reco::GsfElectronRef &ele,
    const edm::Handle<reco::ConversionCollection> &conversions,
    const math::XYZPoint &beamspot,
    const reco::Vertex & vertex,
    const IsoDepositVals &isoVals,
    const double &rho)
{

    // get the ID variables from the electron object
    bool isEB           = ele->isEB() ? true : false;
    float pt            = ele->pt();
    float eta           = ele->superCluster()->eta();
    float dEtaIn        = ele->deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn        = ele->deltaPhiSuperClusterTrackAtVtx();
    float sigmaIEtaIEta = ele->sigmaIetaIeta();
    float hoe           = ele->hadronicOverEm();
    float ooemoop       = (1.0/ele->superCluster()->energy() - 1.0/ele->p());
    float d0vtx         = ele->gsfTrack()->dxy(vertex.position());
    float dzvtx         = ele->gsfTrack()->dz(vertex.position());
    float iso_nh        = (*(isoVals[3]))[ele];
    float iso_ch        = (*(isoVals[0]))[ele];
    float iso_em        = (*(isoVals[1]))[ele];
    float vtxFit        = ConversionTools::hasMatchedConversion(*ele, conversions, beamspot);
    float mHits         = ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits(); 

    // get the mask value
    unsigned int mask = EgammaCutBasedEleId::TestWP(workingPoint, isEB, pt, eta, dEtaIn, dPhiIn,
        sigmaIEtaIEta, hoe, ooemoop, d0vtx, dzvtx, iso_nh, iso_ch, iso_em, vtxFit, mHits, rho);

    // return the mask value
    return mask;

}

#endif

unsigned int EgammaCutBasedEleId::PassWP(WorkingPoint workingPoint, const bool isEB, const float pt, const float eta,
    const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe,
    const float ooemoop, const float d0vtx, const float dzvtx, const float iso_nh, const float iso_ch, const float iso_em, 
    const bool vtxFit, const unsigned int mHits, const double rho)
{
    unsigned int mask = EgammaCutBasedEleId::TestWP(workingPoint, isEB, pt, eta, dEtaIn, dPhiIn,
        sigmaIEtaIEta, hoe, ooemoop, d0vtx, dzvtx, iso_nh, iso_ch, iso_em, vtxFit, mHits, rho);
    if ((mask & PassAll) == PassAll) return true;
    return false;
}

unsigned int EgammaCutBasedEleId::TestWP(WorkingPoint workingPoint, const bool isEB, const float pt, const float eta,
    const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe, 
    const float ooemoop, const float d0vtx, const float dzvtx, const float iso_nh, const float iso_ch, const float iso_em, 
    const bool vtxFit, const unsigned int mHits, const double rho)
{

    unsigned int mask = 0;
    float cut_dEtaIn[2]         = {999.9, 999.9};
    float cut_dPhiIn[2]         = {999.9, 999.9};
    float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
    float cut_hoe[2]            = {999.9, 999.9};
    float cut_ooemoop[2]        = {999.9, 999.9};
    float cut_d0vtx[2]          = {999.9, 999.9};
    float cut_dzvtx[2]          = {999.9, 999.9};
    float cut_iso[2]            = {999.9, 999.9};
    bool cut_vtxFit[2]          = {false, false};
    unsigned int cut_mHits[2]   = {999, 999};

    if (workingPoint == EgammaCutBasedEleId::VETO) {
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.010;
        cut_dPhiIn[0]        = 0.800; cut_dPhiIn[1]        = 0.700;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.150; cut_hoe[1]           = 999.9;
        cut_ooemoop[0]       = 999.9; cut_ooemoop[1]       = 999.9;
        cut_d0vtx[0]         = 0.040; cut_d0vtx[1]         = 0.040;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = false; cut_vtxFit[1]        = false;
        cut_mHits[0]         = 999  ; cut_mHits[1]         = 999;
        cut_iso[0]           = 0.150; cut_iso[1]           = 0.150;
    } 
    else if (workingPoint == EgammaCutBasedEleId::LOOSE) {
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.009;
        cut_dPhiIn[0]        = 0.015; cut_dPhiIn[1]        = 0.010;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    } 
    else if (workingPoint == EgammaCutBasedEleId::MEDIUM) {
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.007;
        cut_dPhiIn[0]        = 0.006; cut_dPhiIn[1]        = 0.003;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    } 
    else if (workingPoint == EgammaCutBasedEleId::TIGHT) {
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.005;
        cut_dPhiIn[0]        = 0.003; cut_dPhiIn[1]        = 0.002;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 0    ; cut_mHits[1]         = 0;
        if (pt >= 20.0) {
            cut_iso[0] = 0.100; cut_iso[1] = 0.100;
        }
        else {
            cut_iso[0] = 0.100; cut_iso[1] = 0.070;
        }
    } 
    else {
        std::cout << "[EgammaCutBasedEleId::TestWP] Undefined working point" << std::endl;
    }

    // choose cut if barrel or endcap
    unsigned int idx = isEB ? 0 : 1;
    float etaAbs = fabs(eta);

    // effective area for isolation
    float AEff = 0.18;
    if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.19;
    if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.21;
    if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.38;
    if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.61;
    if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.73;
    if (etaAbs > 2.4) AEff = 0.78;

    // apply to neutrals
    double rhoPrime = std::max(rho, 0.0);
    double iso_n = std::max(iso_nh + iso_em - rhoPrime * AEff, 0.0);

    // compute final isolation
    double iso = (iso_n + iso_ch) / pt;

    // test cuts
    if (fabs(dEtaIn) < cut_dEtaIn[idx])             mask |= DETAIN;
    if (fabs(dPhiIn) < cut_dPhiIn[idx])             mask |= DPHIIN; 
    if (sigmaIEtaIEta < cut_sigmaIEtaIEta[idx])     mask |= SIGMAIETAIETA;
    if (hoe < cut_hoe[idx])                         mask |= HOE;
    if (fabs(ooemoop) < cut_ooemoop[idx])           mask |= OOEMOOP;
    if (fabs(d0vtx) < cut_d0vtx[idx])               mask |= D0VTX;
    if (fabs(dzvtx) < cut_dzvtx[idx])               mask |= DZVTX;
    if (cut_vtxFit[idx] && vtxFit)                  mask |= VTXFIT;
    if (mHits <= cut_mHits[idx])                    mask |= MHITS;
    if (iso < cut_iso[idx])                         mask |= ISO;

    // return the mask
    return mask;

}

