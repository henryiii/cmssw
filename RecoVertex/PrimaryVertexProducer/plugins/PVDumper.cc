#include "PVDumper.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"
#include <fstream>
#include <iostream>

typedef TrackingParticleRefVector::iterator tp_iterator;

PVDumper::PVDumper(const edm::ParameterSet &conf)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf),
      myfile("output_tree.root", "RECREATE"), mytree("trks", "Tracks and such"), writer(&mytree) {

  mytree.Branch("evt_num", &evt_num, "l");

  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  vtxToken = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
  trkTokenView = consumes<edm::View<reco::Track>>(conf.getParameter<edm::InputTag>("TrackLabel"));

  // Open files
  f_trk.open("tracks_info.txt", std::ios::out);
  f_vtx.open("vertex_info.txt", std::ios::out);
  f_truth_vtx.open("truth_vertex_info.txt", std::ios::out);
  f_truth_trk.open("truth_track_info.txt", std::ios::out);

  vec_TrackingParticle_Token_ = consumes<TrackingParticleCollection>(edm::InputTag("mix", "MergedTrackTruth"));
  vec_VertexParticle_Token_ = consumes<TrackingVertexCollection>(edm::InputTag("mix", "MergedTrackTruth"));
  associatorToken_ = consumes<reco::TrackToTrackingParticleAssociator>(edm::InputTag("quickTrackAssociatorByHits"));

  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  produces<int>();
}

PVDumper::~PVDumper() {}

void PVDumper::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {

  edm::Handle<TrackingParticleCollection> TruthTrackContainer;
  edm::Handle<TrackingVertexCollection> TruthVertexContainer;
  iEvent.getByToken(vec_TrackingParticle_Token_, TruthTrackContainer);
  iEvent.getByToken(vec_VertexParticle_Token_, TruthVertexContainer);
  const TrackingParticleCollection *tPC = TruthTrackContainer.product();
  const TrackingVertexCollection *tVC = TruthVertexContainer.product();

  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  iEvent.getByToken(associatorToken_, theAssociator);

  evt_num = iEvent.id().event();

  unsigned int i = 0;
  std::map<const TrackingParticle *, unsigned int> trktru_map;
  for (TrackingParticleCollection::const_iterator t = tPC->begin(); t != tPC->end(); ++t) {

    if (t->eventId().bunchCrossing() != 0)
      continue;

    /// Truth Track info
    writer.prt_x.push_back(t->vx());
    writer.prt_y.push_back(t->vy());
    writer.prt_z.push_back(t->vz());
    writer.prt_px.push_back(t->px());
    writer.prt_py.push_back(t->py());
    writer.prt_pz.push_back(t->pz());
    writer.prt_e.push_back(0); // TODO: Remove?

    writer.prt_pvr.push_back(0);
    // TODO: missing pvr_prt

    f_truth_trk << " TruthTrack ";
    f_truth_trk << iEvent.luminosityBlock() << " " << iEvent.id().event() << " " << i << " "; //<< &(*t) ;
    // f_truth_trk << t -> eventId().bunchCrossing() << " ";
    f_truth_trk << " " << t->charge() << " " << t->px() << " " << t->py() << " " << t->pz() << " ";
    f_truth_trk << t->vx() << " " << t->vy() << " " << t->vz() << " dud ";

    // It appears that the GEN information is only available for the primary
    // vertex not the pileup ones so I am removing it
    // if ( t->genParticle_begin() != t->genParticle_end() ) {
    //  const reco::GenParticle *gp = &(**(t->genParticle_begin()));
    //
    //  f_truth_trk << (*(t->genParticle_begin()))->momentum().rho() << " "
    //		  << (*(t->genParticle_begin()))->vertex().x() << " ";
    //  const reco::GenParticle *mother=gp;
    //  // I want the vertex of the first decay - so the daughter of the top particle
    //  while ( (mother->mother() != nullptr) &&
    //  (mother->mother()->mother() != nullptr) ) mother=(const
    //  reco::GenParticle *)mother->mother(); f_truth_trk <<
    //  mother->vertex().x() << " ";
    //
    //}
    // else{
    //  f_truth_trk << "0 0 0";
    //}

    f_truth_trk << std::endl;
    trktru_map[&(*t)] = i;
    i++;
  }

  i = 0;
  for (TrackingVertexCollection::const_iterator v = tVC->begin(); v != tVC->end(); ++v) {

    if (v->eventId().bunchCrossing() != 0)
      continue;

    // Vertex info
    writer.pvr_x.push_back(v->position().x());
    writer.pvr_y.push_back(v->position().y());
    writer.pvr_z.push_back(v->position().z());
    writer.ntrks_prompt.push_back(v->daughterTracks().size());

    f_truth_vtx << " TruthVertex ";
    f_truth_vtx << iEvent.luminosityBlock() << " " << iEvent.id().event() << " " << i << " "; // << &(*v) << " ";
    // f_truth_vtx << v -> eventId().bunchCrossing() << " ";
    f_truth_vtx << v->position().x() << " " << v->position().y() << " " << v->position().z() << " "
                << v->daughterTracks().size() << " ";

    for (tp_iterator iTP = v->daughterTracks_begin(); iTP != v->daughterTracks_end(); ++iTP) {
      f_truth_vtx << trktru_map[&(*(*iTP))] << " "; // << &(*(*iTP)) << " " ;
    }

    f_truth_vtx << std::endl;
    i++;
  }

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken, tks);

  edm::Handle<edm::View<reco::Track>> trackCollectionV;
  iEvent.getByToken(trkTokenView, trackCollectionV);

  reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackCollectionV, TruthTrackContainer);

  // get the BeamSpot, it will always be needed, even when not used as a
  // constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  // interface RECO tracks to vertex reconstruction
  const auto &theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;
  t_tks = (*theB).build(tks, beamSpot);

  // select tracks
  std::vector<reco::TransientTrack> &&seltks = theTrackFilter->select(t_tks);

  int i_val = 0;
  //    f_trk << "TrkInfo1 " << &i << " " << i.pt() << std::endl;
  // }
  // for ( auto & i : seltks) {
  //  const auto &my_trk= i;//.track();

  std::map<const reco::Track *, unsigned int> trk_map;
  for (auto &i : (*tks)) {
    int is_good = 0;
    const auto &my_trk = i;
    trk_map[&my_trk] = i_val;
    for (auto &j : seltks) {
      const auto &my_trkbase = j.trackBaseRef();
      if (&(*my_trkbase) == &i) {
        is_good = true;
        break;
      }
    }

    const math::XYZPoint &inner_position = my_trk.innerPosition();
    const math::XYZVector &inner_mom = my_trk.innerMomentum();
    const math::XYZVector &unit_mom = inner_mom.unit();

    double errd0 = my_trk.d0Error();
    double errz0 = my_trk.dzError();

    // Reco tracks
    writer.recon_x.push_back(inner_position.x());
    writer.recon_y.push_back(inner_position.y());
    writer.recon_z.push_back(inner_position.z());
    writer.recon_tx.push_back(unit_mom.x());
    writer.recon_ty.push_back(unit_mom.y());
    writer.recon_chi2.push_back(my_trk.chi2());
    // TODO: Use actual POCA
    writer.recon_pocax.push_back(inner_position.x());
    writer.recon_pocay.push_back(inner_position.y());
    writer.recon_pocaz.push_back(inner_position.z());
    writer.recon_sigmapocaxy.push_back(errd0);
    writer.recon_errz0.push_back(errz0);

    f_trk << "TrkInfo " << iEvent.luminosityBlock() << " " << iEvent.id().event() << " ";
    f_trk << i_val << " "; //<< " " << &my_trk << " "
    f_trk << is_good
          //<< &my_trk << " " << &(*(my_trkbase))
          << " " << my_trk.pt() << " "; // << my_trkbase->pt() << " ";

    for (unsigned ij = 0; ij < 5; ij++)
      f_trk << my_trk.parameter(ij) << " ";
    const reco::Track::CovarianceMatrix innerStateCovariance = my_trk.innerStateCovariance();

    // f_trk << "[ ";
    for (int ii = 0; ii < 5; ii++) {
      // f_trk << "[ ";
      for (int ij = 0; ij < 5; ij++) {
        f_trk << innerStateCovariance(ii, ij);
        if (ij != 4)
          f_trk << ", ";
      }
      // f_trk << " ]";
      if (ii != 4)
        f_trk << ", ";
      // else f_trk << " ";
    }
    // f_trk << " ] ";

    edm::RefToBase<reco::Track> trackref(trackCollectionV, i_val);
    reco::RecoToSimCollection::const_iterator f = recSimColl.find(trackref);

    if (f != recSimColl.end()) {
      TrackingParticleRef tp = f->val.front().first;
      int tpVal = -1;
      if (trktru_map.find(&(*tp)) != trktru_map.end()) {
        tpVal = trktru_map[&(*tp)];
      }
      f_trk << tpVal << " ";
    }

    f_trk << std::endl;
    i_val++;
  }

  edm::Handle<reco::VertexCollection> vColl;
  iEvent.getByToken(vtxToken, vColl);

  int i_vtx = 0;
  for (reco::VertexCollection::const_iterator v = vColl->begin(); v != vColl->end(); ++v) {
    f_vtx << "VtxInfo "
          << " " << iEvent.luminosityBlock() << " " << iEvent.id().event() << " " << i_vtx << " ";
    f_vtx << v->position().x() << " " << v->position().y() << " " << v->position().z() << " ";
    f_vtx << v->xError() << " " << v->yError() << " " << v->zError() << " ";
    f_vtx << v->tracksSize() << " ";
    for (auto iv = v->tracks_begin(); iv != v->tracks_end(); iv++) {
      f_vtx << trk_map[&(*(*iv))] << " ";
      // f_vtx << &(*(*iv)) << " ";
      // f_vtx << (*iv)->pt() << " ";
    }
    f_vtx << std::endl;
    i_vtx++;
  }

  auto result = std::make_unique<int>(1);
  iEvent.put(std::move(result));
  mytree.Fill();
  myfile.Write();
  writer.clear();
}

void PVDumper::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // offlinePrimaryVertices
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<std::string>("label", "");
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("minNdof", 0.0);
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 0.0);
      temp1.push_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 2.0);
      temp1.push_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    TrackFilterForPVFinding::fillPSetDescription(psd0);
    psd0.add<int>("numTracksThreshold", 0); // HI only
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel",
                          edm::InputTag("dummy_default")); // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel",
                          edm::InputTag("dummy_default")); // 4D only

  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZT_vect::fillPSetDescription(psd1);
      psd0.add<edm::ParameterSetDescription>("TkDAClusParameters", psd1);

      edm::ParameterSetDescription psd2;
      GapClusterizerInZ::fillPSetDescription(psd2);
      psd0.add<edm::ParameterSetDescription>("TkGapClusParameters", psd2);
    }
    psd0.add<std::string>("algorithm", "DA_vect");
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});

  descriptions.add("pvDumper", desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(PVDumper);
