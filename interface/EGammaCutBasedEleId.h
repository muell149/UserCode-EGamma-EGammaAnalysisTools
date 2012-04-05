
//
// to use this code outside of CMSSW
// set this definition
//

//#define STANDALONEID

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include <vector>

namespace EgammaCutBasedEleId {

//
// typedefs
//

typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > >   IsoDepositMaps;
typedef std::vector< edm::Handle< edm::ValueMap<double> > >             IsoDepositVals;

//
// defined ID working points
//

enum WorkingPoint {
    VETO,
    LOOSE,
    MEDIUM,
    TIGHT
};

//
// cuts used within working points
//

enum CutType {
    DETAIN          = (1<<0),
    DPHIIN          = (1<<1),
    SIGMAIETAIETA   = (1<<2),
    HOE             = (1<<3),
    OOEMOOP         = (1<<4),
    D0VTX           = (1<<5),
    DZVTX           = (1<<6),
    ISO             = (1<<7),
    VTXFIT          = (1<<8),
    MHITS           = (1<<9),
    EOPFBREM        = (1<<10)
};

//
// all possible cuts pass
//

static const unsigned int PassAll = DETAIN | DPHIIN | SIGMAIETAIETA | HOE | OOEMOOP | D0VTX | DZVTX | ISO | VTXFIT | MHITS;

//
// CMSSW interface
//

#ifndef STANDALONEID

bool PassWP(const WorkingPoint workingPoint,
    const reco::GsfElectronRef &ele,
    const edm::Handle<reco::ConversionCollection> &conversions,
    const math::XYZPoint &beamspot,
    const reco::Vertex & vertex,
    const IsoDepositVals &isoVals,
    const double &rho);

unsigned int TestWP(const WorkingPoint workingPoint,
    const reco::GsfElectronRef &ele,
    const edm::Handle<reco::ConversionCollection> &conversions,
    const math::XYZPoint &beamspot,
    const reco::Vertex & vertex,
    const IsoDepositVals &isoVals,
    const double &rho);

#endif

//
// implementation of working points
// used by CMSSW interface, does not 
// itself depend on CMSSW code
//

unsigned int PassWP(WorkingPoint workingPoint, const bool isEB, const float pt, const float eta,
    const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe,
    const float ooemoop, const float d0vtx, const float dzvtx, const float iso_nh, const float iso_ch, const float iso_em, 
    const bool vtxFit, const unsigned int mHits, const double rho);

unsigned int TestWP(WorkingPoint workingPoint, const bool isEB, const float pt, const float eta,
    const float dEtaIn, const float dPhiIn, const float sigmaIEtaIEta, const float hoe,
    const float ooemoop, const float d0vtx, const float dzvtx, const float iso_nh, const float iso_ch, const float iso_em, 
    const bool vtxFit, const unsigned int mHits, const double rho);

}

