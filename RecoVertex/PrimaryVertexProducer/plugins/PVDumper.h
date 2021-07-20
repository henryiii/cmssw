// -*- C++ -*-
//
// Package:    PrimaryVertexProducer
// Class:      PrimaryVertexProducer
//
/**\class PrimaryVertexProducer PrimaryVertexProducer.cc RecoVertex/PrimaryVertexProducer/src/PrimaryVertexProducer.cc

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
#include <iostream>
#include <fstream>

//#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZT_vect.h"

#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include <algorithm>
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

//
// class declaration
//

class PVDumper : public edm::one::EDProducer<> {
public:
  PVDumper(const edm::ParameterSet&);
  ~PVDumper() override;

  void produce(edm::Event&, const edm::EventSetup&) override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // access to config
  edm::ParameterSet config() const { return theConfig; }

private:
  // ----------member data ---------------------------

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
  TrackFilterForPVFindingBase* theTrackFilter;
  edm::ParameterSet theConfig;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > vecPileupSummaryInfoToken_;
  edm::EDGetTokenT<TrackingParticleCollection> vec_TrackingParticle_Token_;
  edm::EDGetTokenT<TrackingVertexCollection> vec_VertexParticle_Token_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> trkTokenView;

  std::ofstream f_trk, f_vtx, f_truth_vtx, f_truth_trk;
};

