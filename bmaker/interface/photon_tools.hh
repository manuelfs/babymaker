// photon_tools: Functions related to RA2/b photons
// Copied from RA2/b's TreeMaker
// https://github.com/TreeMaker/TreeMaker/blob/74848c2a9fd16acaaf8d711842c6e2dc1d5bd56b/Utils/src/PhotonIDisoProducer.cc

#ifndef H_PHOTON_TOOLS
#define H_PHOTON_TOOLS

// System include files
#include <vector>

// FW physics include files
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class photon_tools{
public:

  ///////////////// PHOTON CUTS ///////////////////////
  const float PhotonPtCut  = 100.0;

  std::vector<double> etaHigh;
  std::vector<double> etaLow;
  std::vector< std::vector<double> > effA;
  enum isoType {pfCh,pfNu,pfGam,coneCh,coneNu,coneGam,NisoTypes};

  bool idPhoton(const pat::Photon &photon, edm::Handle<std::vector<pat::Electron> > &electrons,
		edm::Handle<reco::ConversionCollection> &conversions,
		edm::Handle<reco::BeamSpot> &beamspot, double rho);
  bool electronMatch(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &electrons,
		     const edm::Handle<reco::ConversionCollection> &conversions, const math::XYZPoint &beamspot,
		     double lxyMin=2.0, double probMin=1e-6, unsigned int nHitsBeforeVtxMax=0) ;
  double rhoCorrectedIso(isoType isoType_, double isoVal, double eta, double rho);
  void addEffA(double etaLow_,double etaHigh_,double effA_pfCh_,double effA_pfNu_,double effA_pfGa_);

  photon_tools();
  ~photon_tools();
};



#endif
