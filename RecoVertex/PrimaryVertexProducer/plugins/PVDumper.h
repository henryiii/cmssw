// -*- C++ -*-
//
// Package:    PrimaryVertexProducer
// Class:      PrimaryVertexProducer
//
/**\class PrimaryVertexProducer PrimaryVertexProducer.cc
 RecoVertex/PrimaryVertexProducer/src/PrimaryVertexProducer.cc

 Description: steers tracker primary vertex reconstruction and storage

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Pascal Vanlaer
//         Created:  Tue Feb 28 11:06:34 CET 2006
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <fstream>
#include <iostream>

//#include
//"RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZT_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include <TFile.h>
#include <TTree.h>
#include <algorithm>

struct Writer;

template <class T> struct BranchVector : public std::vector<T> { BranchVector(Writer *parent, TString name); };

struct Writer {
  using Vec = BranchVector<double>;

  TTree *tree;
  std::map<TString, Vec *> known_branches;

  /// Truth PVs and SVs
  // PV locations, randomly generated from LHCb expected Gaussian
  Vec pvr_x{this, "pvr_x"};
  Vec pvr_y{this, "pvr_y"};
  Vec pvr_z{this, "pvr_z"};
  // SV locations, randomly generated from LHCb expected Gaussian
  Vec svr_x{this, "svr_x"};
  Vec svr_y{this, "svr_y"};
  Vec svr_z{this, "svr_z"};
  // The ID of the owning PV for the SV
  Vec svr_pvr{this, "svr_pvr"};

  /// Hit information
  // Leave empty
  Vec prt_pid{this, "prt_pid"};
  Vec prt_hits{this, "prt_hits"};

  /// Truth tracks
  // Momentum (direction) of particle
  Vec prt_px{this, "prt_px"};
  Vec prt_py{this, "prt_py"};
  Vec prt_pz{this, "prt_pz"};
  // Leave 0
  Vec prt_e{this, "prt_e"};
  // Location of particle
  Vec prt_x{this, "prt_x"};
  Vec prt_y{this, "prt_y"};
  Vec prt_z{this, "prt_z"};
  // The ID of the owning PV
  Vec prt_pvr{this, "prt_pvr"};
  // Number of prompt tracks in event
  Vec ntrks_prompt{this, "ntrks_prompt"};

  /// Reconstructed tracks
  Vec recon_x{this, "recon_x"};
  Vec recon_y{this, "recon_y"};
  Vec recon_z{this, "recon_z"};
  Vec recon_tx{this, "recon_tx"};
  Vec recon_ty{this, "recon_ty"};
  Vec recon_chi2{this, "recon_chi2"};
  Vec recon_pocax{this, "recon_pocax"};
  Vec recon_pocay{this, "recon_pocay"};
  Vec recon_pocaz{this, "recon_pocaz"};
  Vec recon_sigmapocaxy{this, "recon_sigmapocaxy"};
  Vec recon_errz0{this, "recon_errz0"};

  Writer(TTree *tree) : tree(tree) {}

  inline void clear() {
    for (auto &[_, value] : known_branches)
      value->clear();
  }
};

//
// class declaration
//

class PVDumper : public edm::one::EDProducer<> {
public:
  PVDumper(const edm::ParameterSet &);
  ~PVDumper() override;

  void produce(edm::Event &, const edm::EventSetup &) override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  // access to config
  edm::ParameterSet config() const { return theConfig; }

private:
  // ----------member data ---------------------------

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
  TrackFilterForPVFindingBase *theTrackFilter;
  edm::ParameterSet theConfig;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
  edm::EDGetTokenT<edm::ValueMap<float>> trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float>> trkTimeResosToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> vecPileupSummaryInfoToken_;
  edm::EDGetTokenT<TrackingParticleCollection> vec_TrackingParticle_Token_;
  edm::EDGetTokenT<TrackingVertexCollection> vec_VertexParticle_Token_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> trkTokenView;

  std::ofstream f_trk, f_vtx, f_truth_vtx, f_truth_trk;

  TFile myfile;
  TTree mytree;
  Writer writer;
  ULong64_t evt_num = 0;
};

template <class T> BranchVector<T>::BranchVector(Writer *parent, TString name) : std::vector<T>() {
  parent->tree->Branch(name, static_cast<std::vector<T> *>(this));
  parent->known_branches.emplace(name, this);
}
