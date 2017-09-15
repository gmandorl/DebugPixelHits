// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTrackerEvent.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/MeasurementDet/interface/TempMeasurements.h"
#include "MagneticField/Engine/interface/MagneticField.h"

//#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
//#include "CalibFormats/SiPixelObjects/interface/SiPixelQuality.h"
//#include "CalibTracker/Records/interface/SiPixelQualityRcd.h"
#include "RecoTracker/MeasurementDet/src/TkMeasurementDetSet.h"
#include "../interface/BadComponents.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h" 
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit2DPosConstraint.h"


#include "DebugPixelHits/DebugPixelHits/interface/VarsOfTrack.h"


class DebugPixelHits_TTbar : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit DebugPixelHits_TTbar(const edm::ParameterSet&);
        ~DebugPixelHits_TTbar();

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) ;

        // ----------member data ---------------------------
//         const edm::EDGetTokenT<edm::View<reco::Candidate>> pairs_;   
        edm::EDGetTokenT<edm::View<reco::Candidate>> candidates_;  
        const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterLabel_;
        const edm::EDGetTokenT<MeasurementTrackerEvent> tracker_;
        const edm::EDGetTokenT<LumiScalersCollection> lumiScaler_;
        const edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_; 
        /// Track Transformer
        TrackTransformer refitter_;
        const std::string propagatorOpposite_;
        const std::string propagatorAlong_;
        //std::string pixelQualityLabel_;
        const std::string estimatorName_;
        const double rescaleError_;
        

        edm::ESHandle<TrackerTopology> theTrkTopo;
        edm::ESHandle<Propagator> thePropagatorOpposite;
        edm::ESHandle<Propagator> thePropagatorAlong;
        edm::ESHandle<Chi2MeasurementEstimatorBase> theEstimator;
        edm::ESHandle<TrajectoryStateUpdator>         updator_;        

        //edm::ESHandle<SiPixelQuality> thePixelQuality;
        BadComponents theBadComponents;

        TTree *tree_;
        int run_, lumi_, bx_, instLumi_, npv_; uint64_t event_; 
        float pair_mass_, track_pt_, track_eta_, track_phi_;
        float pv_ndof_, pv_z0_, tag_dz_, tag_dxy_, pv_xyErr_, pv_zErr_, track_dz_, track_dxy_, track_dzErr_, track_dxyErr_;
        int source_det_, source_layer_;
        int detid_, roc_, hitFound_, hitOnTrack_, detIsActive_, maybeBadROC_, trackHasHit_, trackHasLostHit_, hitInRandomWindow_;
        float track_global_phi_, track_global_z_, track_local_x_, track_local_y_, track_exp_sizeX_, track_exp_sizeY_, track_exp_charge_;
        float hit_global_phi_, hit_global_z_, hit_local_x_, hit_local_y_, hit_sizeX_, hit_sizeY_, hit_firstpixel_x_, hit_firstpixel_y_, hit_chi2_, hit_charge_, hitInRandomWindowDistance_;
        
        int cluster_center_x_, cluster_center_y_;
        const static int nHit = 147;//21*7
        
        int cluster_charge_in_hits_[nHit];
        int cluster_column_ON_[nHit];
        
        int cluster_chargeBroken_in_hits_[nHit];
        int cluster_columnBroken_ON_[nHit];
        
        int cluster_x_inModule_[nHit];
        int cluster_y_inModule_[nHit];

        int cluster_hits_x_[nHit], cluster_hits_y_[nHit];
        
        float trackfromPV_pt_,trackfromPV_eta_,trackfromPV_phi_,trackfromPV_global_phi_,trackfromPV_global_z_,trackfromPV_local_x_,trackfromPV_local_y_;
        int trackfromPV_detIsActive_, trackfromPV_roc_, trackfromPV_detid_;
        float trackfromPV_alpha_,trackfromPV_beta_;
        float track_alpha_, track_beta_;
        
        float track_global_x_, track_global_y_, trackfromPV_global_x_, trackfromPV_global_y_;
        
        
        const static int module_x = 160;//21*7
        const static int module_y = 416;//21*7
        
        bool column1_has_hit_[module_y]={ false };
        bool column1_status_[module_y]={ false };
        
        bool column2_has_hit_[module_y]={ false };
        bool column2_status_[module_y]={ false };
        
        std::vector<int> All_hits_charge;
        std::vector<int> All_hits_Px;
        std::vector<int> All_hits_Py;
        
        int debug_;
        
        std::vector<VarsOfTrack> VarsTrack_PXB2;
        std::vector<VarsOfTrack> VarsTrack_wPV;
        
        float trackfromPV_localPixel_x_,trackfromPV_localPixel_y_, track_localPixel_x_,track_localPixel_y_; 
        float hit_localPixel_x_,hit_localPixel_y_, cluster_localPixel_x_,cluster_localPixel_y_;
        float trackfromPV_local_Dx_,trackfromPV_local_Dy_,track_local_Dx_,track_local_Dy_;
        
        bool broken_cluster_ = false;
        
        
};

//
// constructors and destructor
//
DebugPixelHits_TTbar::DebugPixelHits_TTbar(const edm::ParameterSet& iConfig):
    //pairs_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pairs"))),
    //pixelClusterLabel_(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("pixelClusters"))),
    tracker_(consumes<MeasurementTrackerEvent>(iConfig.getParameter<edm::InputTag>("tracker"))),
    lumiScaler_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"))),
    vertices_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
    refitter_(iConfig),
    propagatorOpposite_(iConfig.getParameter<std::string>("PropagatorOpposite")),
    propagatorAlong_(iConfig.getParameter<std::string>("PropagatorAlong")),
    estimatorName_(iConfig.getParameter<std::string>("Chi2MeasurementEstimator")),
    rescaleError_(iConfig.getParameter<double>("rescaleError")),
    debug_(iConfig.getUntrackedParameter<int>("debug",0))
    //pixelQualityLabel_(iConfig.getParameter<std::string>("SiPixelQuality"))
{
    if (iConfig.existsAs<std::string>("badComponentsFile")) {
        theBadComponents.init(iConfig.getParameter<std::string>("badComponentsFile"));
    }

    usesResource("TFileService");
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("tree","tree");
    tree_->Branch("run",  &run_,  "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("bx", &bx_, "bx/i");
    tree_->Branch("instLumi", &instLumi_, "instLumi/i");
    tree_->Branch("npv", &npv_, "npv/i");
    tree_->Branch("pair_mass", &pair_mass_, "pair_mass/F");
    tree_->Branch("track_pt", &track_pt_, "track_pt/F");
    tree_->Branch("track_eta", &track_eta_, "track_eta/F");
    tree_->Branch("track_phi", &track_phi_, "track_phi/F");
    tree_->Branch("track_global_phi", &track_global_phi_, "track_global_phi/F");
    tree_->Branch("track_global_z", &track_global_z_, "track_global_z/F");
    tree_->Branch("track_local_x", &track_local_x_, "track_local_x/F");
    tree_->Branch("track_local_y", &track_local_y_, "track_local_y/F");
    tree_->Branch("track_exp_sizeX", &track_exp_sizeX_, "track_exp_sizeX/F");
    tree_->Branch("track_exp_sizeY", &track_exp_sizeY_, "track_exp_sizeY/F");
    tree_->Branch("track_exp_charge", &track_exp_charge_, "track_exp_charge/F");
    tree_->Branch("hit_global_phi", &hit_global_phi_, "hit_global_phi/F");
    tree_->Branch("hit_global_z", &hit_global_z_, "hit_global_z/F");
    tree_->Branch("hit_local_x", &hit_local_x_, "hit_local_x/F");
    tree_->Branch("hit_local_y", &hit_local_y_, "hit_local_y/F");
    tree_->Branch("hit_firstpixel_x", &hit_firstpixel_x_, "hit_firstpixel_x/F");
    tree_->Branch("hit_firstpixel_y", &hit_firstpixel_y_, "hit_firstpixel_y/F");
    tree_->Branch("hit_sizeX", &hit_sizeX_, "hit_sizeX/F");
    tree_->Branch("hit_sizeY", &hit_sizeY_, "hit_sizeY/F");
    tree_->Branch("hit_chi2", &hit_chi2_, "hit_chi2/F");
    tree_->Branch("hit_charge", &hit_charge_, "hit_charge/F");

    tree_->Branch("pv_ndof", &pv_ndof_, "pv_ndof/F");
    tree_->Branch("pv_z0", &pv_z0_, "pv_z0/F");
    tree_->Branch("tag_dz", &tag_dz_, "tag_dz/F");
    tree_->Branch("tag_dxy", &tag_dxy_, "tag_dxy/F");
    tree_->Branch("pv_xyErr", &pv_xyErr_, "pv_xyErr/F");
    tree_->Branch("pv_zErr", &pv_zErr_, "pv_zErr/F");
    tree_->Branch("track_dz", &track_dz_, "track_dz/F");
    tree_->Branch("track_dxy", &track_dxy_, "track_dxy/F");
    tree_->Branch("track_dzErr", &track_dzErr_, "track_dzErr/F");
    tree_->Branch("track_dxyErr", &track_dxyErr_, "track_dxyErr/F");

    tree_->Branch("detid", &detid_, "detid/I");
    tree_->Branch("roc", &roc_, "roc/I");
    tree_->Branch("source_det", &source_det_, "source_det/I");
    tree_->Branch("source_layer", &source_layer_, "source_layer/I");
    tree_->Branch("hitFound", &hitFound_, "hitFound/I");
    tree_->Branch("hitOnTrack", &hitOnTrack_, "hitOnTrack/I");
    tree_->Branch("hitInRandomWindow", &hitInRandomWindow_, "hitInRandomWindow/I");
    tree_->Branch("hitInRandomWindowDistance", &hitInRandomWindowDistance_, "hitInRandomWindowDistance/F");
    tree_->Branch("detIsActive", &detIsActive_, "detIsActive/I");
    tree_->Branch("maybeBadROC", &maybeBadROC_, "maybeBadROC/I");
    tree_->Branch("trackHasHit", &trackHasHit_, "trackHasHit/I");
    tree_->Branch("trackHasLostHit", &trackHasLostHit_, "trackHasLostHit/I");
    
    //nuove vars
    
    tree_->Branch("track_alpha", &track_alpha_, "track_alpha/F");
    tree_->Branch("track_beta", &track_beta_, "track_beta/F");
    

    tree_->Branch("cluster_center_x", &cluster_center_x_, "cluster_center_x/I");
    tree_->Branch("cluster_center_y", &cluster_center_y_, "cluster_center_y/I");
    //tree_->Branch("cluster_hits_x", cluster_hits_x_, "cluster_hits_x/I");  
    //tree_->Branch("cluster_hits_y", cluster_hits_y_, "cluster_hits_y/I");
    tree_->Branch("cluster_charge_in_hits", &cluster_charge_in_hits_, "cluster_charge_in_hits[147]/I");
    tree_->Branch("cluster_column_ON", &cluster_column_ON_, "cluster_column_ON[147]/I");
    tree_->Branch("cluster_chargeBroken_in_hits", &cluster_chargeBroken_in_hits_, "cluster_chargeBroken_in_hits[147]/I");
    tree_->Branch("cluster_columnBroken_ON", &cluster_columnBroken_ON_, "cluster_columnBroken_ON[147]/I");
//     tree_->Branch("cluster_x_inModule", &cluster_x_inModule_, "cluster_x_inModule[147]/I");
//     tree_->Branch("cluster_y_inModule", &cluster_y_inModule_, "cluster_y_inModule[147]/I");
        
    //nuove variabili

    tree_->Branch("trackfromPV_global_phi", &trackfromPV_global_phi_, "trackfromPV_global_phi/F");
    tree_->Branch("trackfromPV_global_z", &trackfromPV_global_z_, "trackfromPV_global_z/F");
    tree_->Branch("trackfromPV_local_x", &trackfromPV_local_x_, "trackfromPV_local_x/F");
    tree_->Branch("trackfromPV_local_y", &trackfromPV_local_y_, "trackfromPV_local_y/F");
    tree_->Branch("trackfromPV_detid", &trackfromPV_detid_, "trackfromPV_detid/I");
    tree_->Branch("trackfromPV_roc", &trackfromPV_roc_, "trackfromPV_roc/I");
    tree_->Branch("trackfromPV_detIsActive", &trackfromPV_detIsActive_, "trackfromPV_detIsActive/I");
    tree_->Branch("trackfromPV_alpha", &trackfromPV_alpha_, "trackfromPV_alpha/F");
    tree_->Branch("trackfromPV_beta", &trackfromPV_beta_, "trackfromPV_beta/F");
    
    
    tree_->Branch("column1_has_hit",&column1_has_hit_,"column1_has_hit[416]/O");
    tree_->Branch("column1_status",&column1_status_,"column1_status[416]/O");
        
    tree_->Branch("column2_has_hit",&column2_has_hit_,"column2_has_hit[416]/O");
    tree_->Branch("column2_status",&column2_status_,"column2_status[416]/O");
        
    tree_->Branch("All_hits_charge", &All_hits_charge);
    tree_->Branch("All_hits_Px", &All_hits_Px);
    tree_->Branch("All_hits_Py", &All_hits_Py);
    
    tree_->Branch("trackfromPV_localPixel_x", &trackfromPV_localPixel_x_, "trackfromPV_localPixel_x/F");
    tree_->Branch("trackfromPV_localPixel_y", &trackfromPV_localPixel_y_, "trackfromPV_localPixel_y/F");
    tree_->Branch("track_localPixel_x", &track_localPixel_x_, "track_localPixel_x/F");
    tree_->Branch("track_localPixel_y", &track_localPixel_y_, "track_localPixel_y/F");
    tree_->Branch("hit_localPixel_x", &hit_localPixel_x_, "hit_localPixel_x/F");
    tree_->Branch("hit_localPixel_y", &hit_localPixel_y_, "hit_localPixel_y/F");
    
    tree_->Branch("cluster_localPixel_x", &cluster_localPixel_x_, "cluster_localPixel_x/F");
    tree_->Branch("cluster_localPixel_y", &cluster_localPixel_y_, "cluster_localPixel_y/F");
    tree_->Branch("trackfromPV_local_Dy", &trackfromPV_local_Dy_, "trackfromPV_local_Dy/F");
    tree_->Branch("trackfromPV_local_Dx", &trackfromPV_local_Dx_, "trackfromPV_local_Dx/F");
    tree_->Branch("track_local_Dx", &track_local_Dx_, "track_local_Dx/F");
    tree_->Branch("track_local_Dy", &track_local_Dy_, "track_local_Dy/F");
    
    tree_->Branch("broken_cluster", &broken_cluster_, "broken_cluster/O");
    
    candidates_ = consumes<edm::View<reco::Candidate>>(edm::InputTag("particleFlow"));
   
}


DebugPixelHits_TTbar::~DebugPixelHits_TTbar()
{
}

    void
DebugPixelHits_TTbar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    run_  = iEvent.id().run();
    lumi_ = iEvent.id().luminosityBlock();
    bx_   = iEvent.bunchCrossing();

    edm::Handle<LumiScalersCollection> lumiScaler;
    iEvent.getByToken(lumiScaler_, lumiScaler);
    instLumi_ = (!lumiScaler->empty()) ? lumiScaler->front().instantLumi() : -99;

    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(vertices_, vertices);
    npv_ = vertices->size();    
    
    edm::ESHandle<MagneticField> theMFH;
    iSetup.get<IdealMagneticFieldRecord>().get("", theMFH);
    const MagneticField *theMF = theMFH.product();
    
    // read input
//     Handle<View<reco::Candidate> > pairs;
//     iEvent.getByToken(pairs_, pairs);
    
    Handle<View<reco::Candidate> > candidates;
    iEvent.getByToken(candidates_, candidates);

//     if (pairs->empty()) return;

    iSetup.get<TrackerTopologyRcd>().get(theTrkTopo);
    refitter_.setServices(iSetup);

    iSetup.get<TrackingComponentsRecord>().get(estimatorName_, theEstimator);
    iSetup.get<TrackingComponentsRecord>().get(propagatorOpposite_, thePropagatorOpposite);
    iSetup.get<TrackingComponentsRecord>().get(propagatorAlong_, thePropagatorAlong);
    iSetup.get<TrackingComponentsRecord>().get("KFUpdator",updator_);
    
    //iSetup.get<SiPixelQualityRcd>().get(pixelQualityLabel_, thePixelQuality);

    //edm::Handle<edmNew::DetSetVector<SiPixelCluster> > pixelC; 
    //iEvent.getByToken(pixelClusterLabel_, pixelC); 

    Handle<MeasurementTrackerEvent> tracker;
    iEvent.getByToken(tracker_, tracker);
    //const PxMeasurementConditionSet & pixelConds = tracker->measurementTracker().pixelDetConditions();

    std::vector<int> badROCs;    
    
    
    
//     for (const reco::Candidate & pair : *pairs) {
    for (const reco::Candidate & candidate : *candidates) {
//         pair_mass_ = pair.mass();
        pair_mass_ = 0;
        
        if (candidate.bestTrack() == 0 ) continue;
        if (candidate.pt() < 2) continue;
//         const reco::Muon &tag = dynamic_cast<const reco::Muon &>(*pair.daughter(0)->masterClone());
//         if (tag.innerTrack().isNull()) continue;
//         const reco::Muon &mu = dynamic_cast<const reco::Muon &>(*pair.daughter(1)->masterClone());
//         if (mu.innerTrack().isNull()) continue;

        const reco::Track & mutk = *candidate.bestTrack(); //questa non è una mutk ma una traccia qualunque--> cambiare nome
        const reco::Vertex *bestPV = &vertices->front();
        
//         float bestDZ = std::abs(tag.innerTrack()->dz(bestPV->position()));
        
        float bestDZ = std::abs(mutk.dz(bestPV->position()));
        for (const reco::Vertex &pv : *vertices) {
//             float thisDZ = std::abs(tag.innerTrack()->dz(pv.position()));
            float thisDZ = std::abs(mutk.dz(pv.position()));
            if (thisDZ < bestDZ) { bestDZ = thisDZ; bestPV = &pv; }            
        }
        pv_ndof_ = bestPV->ndof();
        pv_z0_ = bestPV->z();
        pv_xyErr_ = std::hypot(bestPV->xError(), bestPV->yError());
        pv_zErr_ = bestPV->zError();
        tag_dz_ = -1;//tag.innerTrack()->dz(bestPV->position());
        tag_dxy_ = -1;//tag.innerTrack()->dxy(bestPV->position());
//         const reco::Track & mutk = *mu.innerTrack();
        track_pt_ = mutk.pt();
        track_eta_ = mutk.eta();
        track_phi_ = mutk.phi();
        track_dz_ = mutk.dz(bestPV->position());
        track_dxy_ = mutk.dxy(bestPV->position());
        track_dzErr_ = mutk.dzError();
        track_dxyErr_ = mutk.dxyError();
        
        if (track_dxy_/track_dxyErr_>2||track_dz_/track_dzErr_>2) continue;
        
        trackHasHit_ = mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
        trackHasLostHit_ = mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS);
        bool PXB1Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
        bool PXB2Valid =  mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2);
        if (debug_) printf("\n\nMuon with pt %f eta %f, found PXB hits: %d, lost PXB inner hits along: %d, lost PXB inner hits before: %d, PXB1 hit: %d, PXB2 hit: %d, PXB3 hit: %d, PXF1 hit: %d, PXF2 hit: %d\n",
//                             mu.pt(),  mu.eta(),
                            candidate.pt(),  candidate.eta(),
                            mutk.hitPattern().numberOfValidPixelBarrelHits(), 
                            mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::TRACK_HITS), mutk.hitPattern().numberOfLostPixelBarrelHits(reco::HitPattern::MISSING_INNER_HITS), 
                            PXB1Valid, PXB2Valid, mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3),
                            mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1), 
                            mutk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2));
        int nhits = mutk.recHitsSize();
        std::vector<Trajectory> traj  = refitter_.transform(mutk);
        if (traj.size() != 1) continue; 
        std::vector<const SiPixelCluster *> PXB1Clusters;
        if (PXB1Valid) {
            for (int i = 0; i < nhits; ++i) {
                const TrackingRecHit *hit = &* mutk.recHit(i);
                if (hit->isValid() && hit->geographicalId().subdetId() == PixelSubdetector::SubDetector::PixelBarrel && 
                        theTrkTopo->layer(hit->geographicalId()) == 1) { 
                    const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
                    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                    auto clustref = pixhit->cluster();
                    if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                    PXB1Clusters.push_back(&*clustref);
                }
            }
            if (PXB1Clusters.empty()) std::cout << "WARNING: did not find Cluster for PXB1 Hit" << std::endl;
        }
        
                
        TrajectoryStateOnSurface tsosPXB2;
        for (const auto &tm : traj.front().measurements()) {
            if (tm.recHit().get() && tm.recHitR().isValid()) {
                DetId where = tm.recHitR().geographicalId();
                source_det_ = where.subdetId();
                source_layer_ = theTrkTopo->layer(where);
                if (source_det_ != PixelSubdetector::SubDetector::PixelBarrel ||  source_layer_ != 1) {
                    tsosPXB2 = tm.updatedState().isValid() ? tm.updatedState() : tm.backwardPredictedState();
                    if (debug_) printf("starting state on det %d, layer %d, r = %5.1f, z = %+6.1f\n", 
                                        source_det_, source_layer_, tsosPXB2.globalPosition().perp(), tsosPXB2.globalPosition().z());
                    break;
                }
            }
        }
        if (!tsosPXB2.isValid()) std::cout << "WARNING: did not find state for PXB2 Hit" << std::endl;
        if (!tsosPXB2.isValid()) continue; // for now
        tsosPXB2.rescaleError(rescaleError_);
               
        const GeometricSearchTracker * gst = tracker->geometricSearchTracker();
        //const auto & PXBLs = gst->pixelBarrelLayers();
        //for (const auto * PXBLayer : PXBLs) { std::cout << "PXB Layer with radius = " << PXBLayer->specificSurface().radius() << std::endl; }
        const auto *pxbLayer1 = gst->pixelBarrelLayers().front();
        auto compDets = pxbLayer1->compatibleDets(tsosPXB2, *thePropagatorOpposite, *theEstimator);

       
        VarsTrack_PXB2.clear();
        for (const auto & detAndState : compDets) {
            VarsOfTrack trackFromLayer2; 
            bool found = false;
            if (debug_) printf("Compatible module %u\n", detAndState.first->geographicalId().rawId()); 
            detid_ = detAndState.first->geographicalId().rawId();
            const auto &mdet = tracker->idToDet(detAndState.first->geographicalId());
            detIsActive_ = mdet.isActive();
            if (debug_) printf("   module center: rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", mdet.position().perp(),  mdet.position().z(), mdet.position().phi().value());
            if (debug_) printf("   track state:   rho = % 6.2f  z = %+6.3f  phi = %+6.3f\n", detAndState.second.globalPosition().perp(), detAndState.second.globalPosition().z(), detAndState.second.globalPosition().phi().value());
            const auto & tkpos = detAndState.second.localPosition();
            const auto & tkerr = detAndState.second.localError().positionError();
            if (debug_) printf("              local x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f \n",  tkpos.x(), std::sqrt(tkerr.xx()),  tkpos.y(), std::sqrt(tkerr.yy()));
            LocalVector localDir = detAndState.second.localMomentum()/detAndState.second.localMomentum().mag();
            float chargeCorr = std::abs(localDir.z());
            track_global_phi_ = detAndState.second.globalPosition().phi();
            track_global_z_   = detAndState.second.globalPosition().z();
            track_global_x_   = detAndState.second.globalPosition().x();
            track_global_y_   = detAndState.second.globalPosition().y();
            track_local_x_ = tkpos.x();
            track_local_y_ = tkpos.y();
            roc_ = (tkpos.x() > 0)*8 + std::min(std::max<int>(0,std::floor(tkpos.y()+3.24)),7);
            track_exp_sizeX_ = 1.5f;
            track_exp_sizeY_ = std::max<float>(1.f, std::hypot(1.3f, 1.9f*mutk.pz()/mutk.pt()));
            track_exp_charge_ = std::sqrt(1.08f + std::pow(mutk.pz()/mutk.pt(),2)) * 26000;
            theBadComponents.fetch(iEvent.id().run(), iEvent.id().luminosityBlock(),  detAndState.first->geographicalId().rawId(), badROCs);
            maybeBadROC_ = false;
            track_local_Dx_=std::sqrt(tkerr.xx());
            track_local_Dy_=std::sqrt(tkerr.yy());
            
                        
//                 posizione locale
            LocalPoint xytrack(track_local_x_, track_local_y_);
            track_localPixel_x_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack).first;
            track_localPixel_y_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack).second;
            
            for (int badROC: badROCs) {
                float rocx = ((badROC/8)-0.5)*82*0.0100; // one ROC is 80 pixels, each 100um wide, but there are big pixels
                float rocy = ((badROC%8)-3.5)*54*0.0150; // one ROC is 52 pixels, each 150um wide, but there are big pixels
                if (debug_) printf("      bad ROC local x = %+6.3f            y = %+6.3f    (gio) \n",  rocx, rocy);
                //MeasurementPoint rocixiy((badROC/8)*80+40, (badROC%8)*52 + 26);
                //LocalPoint rocxy = (dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().localPosition(MeasurementPoint(rocixiy));
                //printf("      bad ROC local x = %+6.3f            y = %+6.3f    (top) \n",  rocxy.x(), rocxy.y());
 
                
                if (std::abs(rocx - tkpos.x()) < (40*0.0100 + 2*0.0100 + 5*std::sqrt(tkerr.xx()))) maybeBadROC_ = true; // 1 half-ROC + 2 pixels + 5 sigma
                if (std::abs(rocy - tkpos.y()) < (26*0.0150 + 2*0.0150 + 5*std::sqrt(tkerr.yy()))) maybeBadROC_ = true; // 1 half-ROC + 2 pixels + 5 sigma
            }
            //auto rechits = mdet.recHits(detAndState.second);
            //for (const auto & hit : rechits) {
            std::vector<std::pair<float, TrackingRecHit::ConstRecHitPointer>> rechitsAndChi2 ;
            for (const auto & hitAndChi2 : mdet.fastMeasurements(detAndState.second, tsosPXB2, *thePropagatorOpposite, *theEstimator)) {
                if (hitAndChi2.recHit()->isValid()) rechitsAndChi2.emplace_back(hitAndChi2.estimate(), hitAndChi2.recHit());                
                if (hitAndChi2.recHit()->isValid()) std::cout << "rechitsAndChi2 FILLING"<<std::endl;
            } 
            const auto & allhits = mdet.recHits(detAndState.second);
            if (rechitsAndChi2.empty()) {
                for (const auto & hit : allhits) { 
                    float distance = std::max(std::abs(hit->localPosition().x()-tkpos.x()), std::abs(hit->localPosition().y()-tkpos.y()));
                    if (hit->isValid() && distance < 0.2) {
                         rechitsAndChi2.emplace_back(99 + distance, hit);
                    }
                }
            }
            std::sort(rechitsAndChi2.begin(), rechitsAndChi2.end());
            // now we make the matching around a fake position, in the same TBM; go inwards by one ROC if one is not in the innermost ROC
            bool innermostROC = std::abs(tkpos.y()) < 0.8;
            float fakey = tkpos.y() + std::copysign(tkpos.y(), 0.8)*(innermostROC ? +1 : -1); 
            hitInRandomWindow_ = false; hitInRandomWindowDistance_ = 0.2;
            for (const auto & hit : allhits) { 
                float distance = std::max(std::abs(hit->localPosition().x()-tkpos.x()), std::abs(hit->localPosition().y()-fakey));
                if (hit->isValid() && distance < hitInRandomWindowDistance_) {
                    hitInRandomWindow_ = true; hitInRandomWindowDistance_ = distance;
                }
            }
            
            trackFromLayer2.set_found(hitFound_);trackFromLayer2.set_detid(detid_);trackFromLayer2.set_detIsActive(detIsActive_);
            trackFromLayer2.set_roc(roc_);
            trackFromLayer2.set_chargeCorr(chargeCorr);
            trackFromLayer2.set_track_global_phi(track_global_phi_);trackFromLayer2.set_track_global_z(track_global_z_);
            trackFromLayer2.set_track_local_x(track_local_x_);trackFromLayer2.set_track_local_y(track_local_y_);
            trackFromLayer2.set_track_local_Dx(track_local_Dx_);trackFromLayer2.set_track_local_Dy(track_local_Dy_);
            trackFromLayer2.set_track_localPixel_x(track_localPixel_x_);trackFromLayer2.set_track_localPixel_y(track_localPixel_y_);
            trackFromLayer2.set_track_exp_sizeX(track_exp_sizeX_);trackFromLayer2.set_track_exp_sizeY(track_exp_sizeY_);
            trackFromLayer2.set_track_exp_charge(track_exp_charge_);        
            trackFromLayer2.set_source_det(source_det_); trackFromLayer2.set_source_layer(source_layer_);        
            trackFromLayer2.set_maybeBadROC(maybeBadROC_); trackFromLayer2.set_trackHasHit(trackHasHit_); 
            trackFromLayer2.set_trackHasLostHit(trackHasLostHit_); 
            trackFromLayer2.set_alpha(acos(track_global_z_/std::hypot(track_global_x_, track_global_z_))); 
            trackFromLayer2.set_beta(acos(track_global_z_/std::hypot(track_global_y_, track_global_z_))); 
            
            for (const auto & hitAndChi2 : rechitsAndChi2) {
                const auto & hit = hitAndChi2.second;
                const auto &hitpos = hit->localPosition();
                const auto &hiterr = hit->localPositionError();
                const auto * pixhit = dynamic_cast<const SiPixelRecHit*>(&*hit);    if (!pixhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                auto clustref = pixhit->cluster();                                  if (clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                auto hitgpos = detAndState.first->toGlobal(hitpos);
                hit_global_phi_ = hitgpos.phi();
                hit_global_z_   = hitgpos.z();
                hit_local_x_ = hitpos.x();
                hit_local_y_ = hitpos.y();
                
                //                 posizione locale
                LocalPoint xytrack(hit_local_x_, hit_local_y_);
                hit_localPixel_x_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack).first;
                hit_localPixel_y_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack).second;
                
                cluster_localPixel_x_=clustref->x();
                cluster_localPixel_y_=clustref->y();
                
                
                cluster_center_x_ = (int) round(track_localPixel_x_/*clustref->x()*/);
                cluster_center_y_ = (int) round(track_localPixel_y_/*clustref->y()*/);
                hit_firstpixel_x_ = clustref->minPixelRow();
                hit_firstpixel_y_ = clustref->minPixelCol();
                hit_chi2_ = hitAndChi2.first;
                hit_charge_ = clustref->charge();
                hit_sizeX_ = clustref->sizeX();
                hit_sizeY_ = clustref->sizeY();
                hitFound_ = true;
                hitOnTrack_ = (std::find(PXB1Clusters.begin(), PXB1Clusters.end(), clustref.get()) != PXB1Clusters.end());

//                 for (const auto &P : clustref->pixels()) {
//                     if(abs(P.x-cluster_center_x_) < 4 &&  abs(P.y-cluster_center_y_) < 11)
//                     cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.adc;
//                     
//                     printf("      casella = %d \n", (P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2));
//                     printf("      pixel at x = %3d - %3d  \t  y = %3d - %3d \t adc: %5d\n", P.x, cluster_center_x_, P.y, cluster_center_y_, P.adc);
//                     printf("      adc1 = %3d  \t ad 2 = %3d \n", cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)], P.adc);
// 
//                 }
//                 for (const auto & all_hit : allhits) {
                
                for (unsigned int c=0; c< module_y; c++){
                column1_has_hit_[c]=false;
                column1_status_[c]=false; 
                column2_has_hit_[c]=false;
                column2_status_[c]=false;
                }
                
                for (int n=0; n<nHit ;++n) cluster_charge_in_hits_[n]=0;
                for (int n=0; n<nHit ;++n) cluster_column_ON_[n]=0;
                for (int n=0; n<nHit ;++n) cluster_chargeBroken_in_hits_[n]=0;
                for (int n=0; n<nHit ;++n) cluster_columnBroken_ON_[n]=0;
//                 for (int n=0; n<nHit ;++n) cluster_x_inModule_[n]=0;
//                 for (int n=0; n<nHit ;++n) cluster_y_inModule_[n]=0;
                
                
                std::vector<int> cluster_double_columns;
                
                for (const auto &P : clustref->pixels()) { 
                    
                    if(abs(P.x-cluster_center_x_) < 4 &&  abs(P.y-cluster_center_y_) < 11 && P.adc>0)
                    {
                        float double_column=P.y/2*1000+P.x/80;
                        std::cout<<P.x<<" "<<P.y<<" "<<double_column<<std::endl;
                        if ( std::find(cluster_double_columns.begin(), cluster_double_columns.end(), double_column) == cluster_double_columns.end() ){
                            std::cout<<P.x<<" "<<P.y<<" "<<double_column<< "  if not in vec already " <<std::endl;
                            cluster_double_columns.push_back(double_column);
                            
                        }
                            
                    }
                    
                }
                
                for (unsigned int ss=0; ss<cluster_double_columns.size(); ++ss) std::cout<<cluster_double_columns[ss]<< "col in columns"<<std::endl;
                
                
                
                int column_to_remove=-1;
                if (cluster_double_columns.size()>1) {
                    column_to_remove= std::rand() % cluster_double_columns.size(); 
                    std::cout<<column_to_remove<< "  column to remove " <<std::endl;
                    
                }
                
                
                std::cout<< "All hits"<<std::endl; 
                for (const auto & all_hit : allhits) {
                   //                 for (const auto & rec_hit : rechitsAndChi2) {
//                     auto & all_hit =rec_hit.second;

                    const auto * pixAllhit = dynamic_cast<const SiPixelRecHit*>(&*all_hit);    if (!pixAllhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                    auto all_clustref = pixAllhit->cluster();                                  if (all_clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");
                    for (const auto &P : all_clustref->pixels()) {
                             std::cout << "colums "<<P.adc << "  " << P.x << "  " << P.y << std::endl;
                             if (P.adc>0 && !column1_has_hit_[P.y] && P.x<80){
                                 std::cout << column1_has_hit_[P.y] << "set value for column " << P.y << std::endl;
                                 column1_has_hit_[P.y]=true;
                                 std::cout << column1_has_hit_[P.y] << "set value for column " << P.y << std::endl;
                            }
                            if (P.adc>0 && !column2_has_hit_[P.y] && P.x>=80){
                                 std::cout << column2_has_hit_[P.y] << "set value for column " << P.y << std::endl;
                                 column2_has_hit_[P.y]=true;
                                 std::cout << column2_has_hit_[P.y] << "set value for column " << P.y << std::endl;
                            }
                             
                            
                    }
                }
                
                for (unsigned int c=0; c< module_y; c++){
                std::cout <<column1_has_hit_[c]<<" " ;
                std::cout <<column2_has_hit_[c]<<" " ;
                if(c%2==0){
                    if(column1_has_hit_[c] || column1_has_hit_[c+1]){
                        column1_status_[c]=  true;    
                        column1_status_[c+1]=  true;  
                    }
                    
                }
                if(c%2==0){
                    if(column2_has_hit_[c] || column2_has_hit_[c+1]){
                        column2_status_[c]=  true;    
                        column2_status_[c+1]=  true;  
                    }
                    
                }
                }
                std::cout <<" "<<std::endl;  
                
                for (unsigned int c=0; c< module_y; c++){
                std::cout <<column1_status_[c]<<" " ;
                }
                std::cout <<" "<<std::endl; 
                
                for (unsigned int c=0; c< module_y; c++){
                std::cout <<column2_status_[c]<<" " ;
                }
                std::cout <<" "<<std::endl; 
                
                
                for (unsigned int c=0; c< module_y; c++){
                    
//                     std::cout<<c<<" "<<cluster_center_y_<<std::endl;
//                     std::cout<<abs((int)c-cluster_center_y_)<<std::endl;
                if(abs((int)c-cluster_center_y_)<11) {
//                     std::cout<<c<<" "<<cluster_center_y_<<std::endl;
                    std::cout <<" has hit: " <<column1_has_hit_[c] << " status :"  << column1_status_[c]<<" "<<c<<std::endl; 
                    std::cout <<" has hit: " <<column2_has_hit_[c] << " status :"  << column2_status_[c]<<" "<<c<<std::endl; 
                    for (int x_size=-3; x_size<4; x_size++) 
                        {   
                            std::cout <<" xxx: " <<x_size+cluster_center_x_ <<" , " << c <<std::endl; 
                            printf("      casella status = %d \n", (x_size) + 7*(c-cluster_center_y_) + ((int) 147/2));
                            if ((x_size+cluster_center_x_)<80) cluster_column_ON_[(x_size) + 7*(c-cluster_center_y_) + ((int) 147/2)]=column1_status_[c];
                            if ((x_size+cluster_center_x_)>=80) cluster_column_ON_[(x_size) + 7*(c-cluster_center_y_) + ((int) 147/2)]=column2_status_[c];
                            
                            
                            cluster_columnBroken_ON_[(x_size) + 7*(c-cluster_center_y_) + ((int) 147/2)]=cluster_column_ON_[(x_size) + 7*(c-cluster_center_y_) + ((int) 147/2)];
                            
                            if (column_to_remove>=0){
                                if (((int)(c)/2)==(cluster_double_columns[column_to_remove]/1000) && ((x_size+cluster_center_x_)/80)==(cluster_double_columns[column_to_remove]%1000)) 
                                    cluster_columnBroken_ON_[(x_size) + 7*(c-cluster_center_y_) + ((int) 147/2)]=0;
                                
                            }
                            
                            
                            //riempire status
                        }
                    }
                }

                All_hits_charge.clear();
                All_hits_Px.clear();
                All_hits_Py.clear();
                for (const auto & all_hit : allhits) {
                   //                 for (const auto & rec_hit : rechitsAndChi2) {
//                     auto & all_hit =rec_hit.second;

                    const auto * pixAllhit = dynamic_cast<const SiPixelRecHit*>(&*all_hit);    if (!pixAllhit) throw cms::Exception("CorruptData", "Valid PXB1 hit that is not a SiPixelRecHit");
                    auto all_clustref = pixAllhit->cluster();                                  if (all_clustref.isNull()) throw cms::Exception("CorruptData", "Valid PXB1 SiPixelRecHit with null cluster ref");             
                    
                    
                    for (const auto &P : all_clustref->pixels()) { 
                        All_hits_charge.push_back(P.adc);
                        All_hits_Px.push_back(P.x);
                        All_hits_Py.push_back(P.y);
                        
                        if(abs(P.x-cluster_center_x_) < 4 &&  abs(P.y-cluster_center_y_) < 11)
                        {
                         cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.adc;
                         cluster_chargeBroken_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.adc;
                         
                         if (column_to_remove>=0){
                             if ((P.y/2)==(cluster_double_columns[column_to_remove]/1000) && (P.x/80)==(cluster_double_columns[column_to_remove]%1000)) 
                                cluster_chargeBroken_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = 0; 
                             
                        }
                         
                        
//                          cluster_x_inModule_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.x;
//                          cluster_y_inModule_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.y;
                        printf("      casella = %d \n", (P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2));
                        printf("      pixel at x = %3d - %3d  \t  y = %3d - %3d \t adc: %5d\n", P.x, cluster_center_x_, P.y, cluster_center_y_, P.adc);
                        printf("      adc1 = %3d  \t ad 2 = %3d \n", cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)], P.adc);
                        printf("      active = %3d  \t ad 2 = %3d ad 2 = %3d \n ", cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)], column1_status_[P.y], column2_status_[P.y]);
                        }


                    }
                    
                    
                }                
                
                
                
                trackFromLayer2.set_hit_global_phi(hit_global_phi_);trackFromLayer2.set_hit_global_z(hit_global_z_);
                trackFromLayer2.set_hit_local_x(hit_local_x_);trackFromLayer2.set_hit_local_y(hit_local_y_);
                trackFromLayer2.set_hit_localPixel_x(hit_localPixel_x_);trackFromLayer2.set_hit_localPixel_y(hit_localPixel_y_);
                trackFromLayer2.set_cluster_localPixel_x(cluster_localPixel_x_);trackFromLayer2.set_cluster_localPixel_y(cluster_localPixel_y_);
                trackFromLayer2.set_cluster_center_x(cluster_center_x_);trackFromLayer2.set_cluster_center_y(cluster_center_y_);
                trackFromLayer2.set_hit_firstpixel_x(hit_firstpixel_x_);trackFromLayer2.set_hit_firstpixel_y(hit_firstpixel_y_);
                trackFromLayer2.set_hit_chi2(hit_chi2_); trackFromLayer2.set_hit_charge(hit_charge_);
                trackFromLayer2.set_hit_sizeX(hit_sizeX_);trackFromLayer2.set_hit_sizeY(hit_sizeY_);
                trackFromLayer2.set_hitFound(hitFound_); trackFromLayer2.set_hitOnTrack(hitOnTrack_);
                trackFromLayer2.set_hitInRandomWindow(hitInRandomWindow_); trackFromLayer2.set_hitInRandomWindowDistance(hitInRandomWindowDistance_);
                trackFromLayer2.set_cluster_charge_in_hits(cluster_charge_in_hits_);
                trackFromLayer2.set_cluster_col_status(cluster_column_ON_);
                trackFromLayer2.set_cluster_chargeBroken_in_hits(cluster_chargeBroken_in_hits_);
                trackFromLayer2.set_cluster_colBroken_status(cluster_columnBroken_ON_);
                trackFromLayer2.set_column1_has_hit(column1_has_hit_);
                trackFromLayer2.set_column1_status(column1_status_);       
                trackFromLayer2.set_column2_has_hit(column2_has_hit_);
                trackFromLayer2.set_column2_status(column2_status_);   

                trackFromLayer2.set_All_hits_charge(All_hits_charge);
                trackFromLayer2.set_All_hits_Px(All_hits_Px);
                trackFromLayer2.set_All_hits_Py(All_hits_Py);
//                 trackFromLayer2.set_cluster_xs(cluster_x_inModule_);
//                 trackFromLayer2.set_cluster_ys(cluster_y_inModule_);
                
                if (!found) {VarsTrack_PXB2.push_back(trackFromLayer2);found = true; }
                    
//                 if (!found) { tree_->Fill(); found = true; } Qui fill vecchio
                
                if (debug_) printf("             rechit x = %+6.3f +- %5.3f   y = %+6.3f +- %5.3f     dx = %+6.3f +- %5.3f   dy = %+6.3f +- %5.3f   chi2 = %7.2f    ontrack = %c  ",  
                                    hitpos.x(), std::sqrt(hiterr.xx()), hitpos.y(), std::sqrt(hiterr.yy()),
                                    hitpos.x()-tkpos.x(), std::sqrt(hiterr.xx()+tkerr.xx()), hitpos.y()-tkpos.y(), std::sqrt(hiterr.yy()+tkerr.yy()), hitAndChi2.first,
                                    (hitOnTrack_ ? 'Y' : 'n'));
                if (debug_) printf(" size: %3d (%2d x %2d), raw charge %8d   corr charge %8.1f\n", clustref->size(), clustref->sizeX(), clustref->sizeY(), clustref->charge(), clustref->charge()*chargeCorr);
//                 if (debug_) {
//                     for (const auto &P : clustref->pixels()) {
//                         printf("      pixel at x = %3d y = %3d \n", P.x, P.y);
//                         cluster_hits_x_.push_back(P.x-cluster_center_x_);
//                         cluster_hits_y_.push_back(P.y-cluster_center_y_);
//                         if(abs(P.x-cluster_center_x_) < 4 &&  abs(P.y-cluster_center_y_) < 11)
//                         cluster_charge_in_hits_[(P.x-cluster_center_x_) + 7*(P.y-cluster_center_y_) + ((int) 147/2)] = P.adc;
                        
//                     }
//                 }
            }
            //another old fill for not found case
            if (!found) {
                hitFound_ = false;
                trackFromLayer2.set_hitFound(hitFound_);
                VarsTrack_PXB2.push_back(trackFromLayer2);
//                 tree_->Fill();
             }
        }
        
        
        
        ///parte nuova
        
        std::cout << "analytical propagator" <<std::endl;        
        AnalyticalPropagator myPropagator(theMF,anyDirection);// no material surface inside 1st hit
        
        std::cout << "TransverseImpactPointExtrapolator" <<std::endl; 
        TransverseImpactPointExtrapolator extrapolator(myPropagator);
        
        
        std::cout << "fts" <<std::endl; 
        //questa fts si userebbe senza propagare dal layer2, ma non va bene fts---> tsosPXB2
//         GlobalPoint tkpos(mutk.vertex().x(),mutk.vertex().y(),mutk.vertex().z());
//         GlobalVector tkmom(mutk.momentum().x(),mutk.momentum().y(),mutk.momentum().z());
//         FreeTrajectoryState fts(tkpos,tkmom,mutk.charge(),theMF);
        /////         
        
        std::cout << "global Point PV" <<std::endl; 
        GlobalPoint pos(bestPV->x(), bestPV->y(), bestPV->z());
        std::cout << "global Error PV" <<std::endl;
        GlobalError err(bestPV->covariance(0,0), 
                        bestPV->covariance(0,1),bestPV->covariance(1,1),
                        bestPV->covariance(0,2),bestPV->covariance(1,2),bestPV->covariance(2,2));
        
        std::cout << "Extrapolation" <<std::endl; 
        std::cout<< tsosPXB2.isValid()<< " is valid "<<std::endl; 
        TrajectoryStateOnSurface tsosAtVtx = extrapolator.extrapolate(*tsosPXB2.freeState(),pos);
        std::cout << "Extrapolated" <<std::endl; 
        const Surface& surfAtVtx = tsosAtVtx.surface();
        LocalPoint vtxPosL = surfAtVtx.toLocal(pos);
        LocalError vtxErrL = ErrorFrameTransformer().transform(err,surfAtVtx);
        TransientTrackingRecHit::RecHitPointer vtxhit = TRecHit2DPosConstraint::build(vtxPosL,vtxErrL,&surfAtVtx);
        
        //update 
        std::cout << "update" <<std::endl; 
        TrajectoryStateOnSurface updated = updator_->update(tsosAtVtx, *vtxhit);
        
        //ora prendere "updated" e propagare su L1       
        //propagazione su PXB1 di updated
        auto compDets_updated = pxbLayer1->compatibleDets(updated, *thePropagatorAlong, *theEstimator);
        
        VarsTrack_wPV.clear();
        for (const auto & detAndState : compDets_updated) {
//             bool found = false;
            int detid2 = detAndState.first->geographicalId().rawId();
            const auto &mdet = tracker->idToDet(detAndState.first->geographicalId());
            int detIsActive2 = mdet.isActive();
            const auto & tkpos = detAndState.second.localPosition();
            const auto & tkerr = detAndState.second.localError().positionError();
//             LocalVector localDir = detAndState.second.localMomentum()/detAndState.second.localMomentum().mag();
            
            float track_local_x2 = tkpos.x();
            float track_local_y2 = tkpos.y();
            
            //                 posizione locale
            LocalPoint xytrack2(track_local_x2, track_local_y2);
            trackfromPV_localPixel_x_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack2).first;
            trackfromPV_localPixel_y_=(dynamic_cast<const PixelGeomDetUnit*>(detAndState.first))->specificTopology().pixel(xytrack2).second;
            
            std::cout << "propagation part"<< std::endl;
            std::cout << track_local_x2 << " " << track_local_y2 << " "<< detid2 << " "<< detIsActive2<< " "<< tkerr.xx() << " "<< std::endl;
            std::cout << track_local_x_ << " " << track_local_y_ << " "<< detid_ << " "<< detIsActive_<< " "<< std::endl;
            
            std::cout <<track_global_phi_ <<" "<< detAndState.second.globalPosition().phi()<<std::endl;
            std::cout <<track_global_z_   <<" "<< detAndState.second.globalPosition().z()<<std::endl;
            
            trackfromPV_global_phi_=detAndState.second.globalPosition().phi();
            trackfromPV_global_z_=detAndState.second.globalPosition().z();
            trackfromPV_global_x_=detAndState.second.globalPosition().x();
            trackfromPV_global_y_=detAndState.second.globalPosition().y();
            
            
            
            VarsOfTrack trackFromPV; 
            
            trackFromPV.set_found(-1);
            trackFromPV.set_detid(detid2); trackFromPV.set_detIsActive(detIsActive2);
            trackFromPV.set_roc((tkpos.x() > 0)*8 + std::min(std::max<int>(0,std::floor(tkpos.y()+3.24)),7));
            trackFromPV.set_track_global_phi(trackfromPV_global_phi_);
            trackFromPV.set_track_global_z(trackfromPV_global_z_);
            trackFromPV.set_track_local_x(track_local_x2);
            trackFromPV.set_track_local_y(track_local_y2);
            trackFromPV.set_track_local_Dx(std::sqrt(tkerr.xx()));
            trackFromPV.set_track_local_Dy(std::sqrt(tkerr.yy()));
            trackFromPV.set_track_localPixel_x(trackfromPV_localPixel_x_);trackFromPV.set_track_localPixel_y(trackfromPV_localPixel_y_);
            trackFromPV.set_alpha(acos(trackfromPV_global_z_/std::hypot(trackfromPV_global_x_, trackfromPV_global_z_))); 
            trackFromPV.set_beta(acos(trackfromPV_global_z_/std::hypot(trackfromPV_global_y_, trackfromPV_global_z_))); 
            
            VarsTrack_wPV.push_back(trackFromPV);
            
            
        }
    //     fill tree multiple times here
    for (unsigned int i=0; i<VarsTrack_PXB2.size(); i++){

        std::cout<< VarsTrack_PXB2.size() << " " << VarsTrack_wPV.size()<< std::endl;
        hitFound_=VarsTrack_PXB2[i].hitFound; 
        detid_=VarsTrack_PXB2[i].detid; detIsActive_=VarsTrack_PXB2[i].detIsActive; roc_=VarsTrack_PXB2[i].roc; 
        track_global_phi_=VarsTrack_PXB2[i].track_global_phi; track_global_z_=VarsTrack_PXB2[i].track_global_z;
        track_local_x_=VarsTrack_PXB2[i].track_local_x; track_local_y_=VarsTrack_PXB2[i].track_local_y;
        track_local_Dx_=VarsTrack_PXB2[i].track_local_Dx; track_local_Dy_=VarsTrack_PXB2[i].track_local_Dy;
        track_localPixel_x_=VarsTrack_PXB2[i].track_localPixel_x; track_localPixel_y_=VarsTrack_PXB2[i].track_localPixel_y;
        track_exp_sizeX_=VarsTrack_PXB2[i].track_exp_sizeX; track_exp_sizeY_=VarsTrack_PXB2[i].track_exp_sizeY; track_exp_charge_=VarsTrack_PXB2[i].track_exp_charge;		
        hit_global_phi_=VarsTrack_PXB2[i].hit_global_phi; hit_global_z_=VarsTrack_PXB2[i].hit_global_z;
        hit_local_x_=VarsTrack_PXB2[i].hit_local_x; hit_local_y_=VarsTrack_PXB2[i].hit_local_y;
        cluster_center_x_=VarsTrack_PXB2[i].cluster_center_x; cluster_center_y_=VarsTrack_PXB2[i].cluster_center_y;
        hit_firstpixel_x_=VarsTrack_PXB2[i].hit_firstpixel_x; hit_firstpixel_y_ =VarsTrack_PXB2[i].hit_firstpixel_y;
        hit_chi2_=VarsTrack_PXB2[i].hit_chi2; hit_charge_=VarsTrack_PXB2[i].hit_charge;
        hit_sizeX_=VarsTrack_PXB2[i].hit_sizeX; hit_sizeY_=VarsTrack_PXB2[i].hit_sizeY; hitFound_=VarsTrack_PXB2[i].hitFound; hitOnTrack_=VarsTrack_PXB2[i].hitOnTrack;
        hitInRandomWindow_=VarsTrack_PXB2[i].hitInRandomWindow;  hitInRandomWindowDistance_=VarsTrack_PXB2[i].hitInRandomWindowDistance;
        hit_localPixel_x_=VarsTrack_PXB2[i].hit_localPixel_x;hit_localPixel_y_=VarsTrack_PXB2[i].hit_localPixel_y;
        cluster_localPixel_x_=VarsTrack_PXB2[i].cluster_localPixel_x;cluster_localPixel_y_=VarsTrack_PXB2[i].cluster_localPixel_y;
        source_det_=VarsTrack_PXB2[i].source_det; source_layer_=VarsTrack_PXB2[i].source_layer;
        maybeBadROC_=VarsTrack_PXB2[i].maybeBadROC; trackHasHit_=VarsTrack_PXB2[i].trackHasHit; trackHasLostHit_=VarsTrack_PXB2[i].trackHasLostHit; 
        track_alpha_=VarsTrack_PXB2[i].alpha; track_beta_=VarsTrack_PXB2[i].beta;
        for (unsigned int k=0; k< nHit; k++){
//             std::cout<<VarsTrack_PXB2[i].cluster_charge_in_hits[k]<<std::endl;
            cluster_charge_in_hits_[k]=VarsTrack_PXB2[i].cluster_charge_in_hits[k];
            cluster_column_ON_[k]=VarsTrack_PXB2[i].cluster_col_status[k];  
            cluster_chargeBroken_in_hits_[k]=VarsTrack_PXB2[i].cluster_chargeBroken_in_hits[k];
            cluster_columnBroken_ON_[k]=VarsTrack_PXB2[i].cluster_colBroken_status[k];  
//             cluster_x_inModule_[k]=VarsTrack_PXB2[i].cluster_xs[k];
//             cluster_y_inModule_[k]=VarsTrack_PXB2[i].cluster_ys[k];
        }
        
        broken_cluster_=!(std::equal(cluster_charge_in_hits_, cluster_charge_in_hits_+147, cluster_chargeBroken_in_hits_)); 
        
       
        for (unsigned int k=0; k< 416; k++)
        {
//             std::cout<<k<<" qui stampaaa stato prima " <<column1_has_hit_[k]<<std::endl;
//             std::cout<<k<<" qui stampaaa stato " <<VarsTrack_PXB2[i].column1_has_hit[k]<<std::endl;
            column1_has_hit_[k]=VarsTrack_PXB2[i].column1_has_hit[k];
//             std::cout<<k<<" qui stampaaa stato dopo " <<column1_has_hit_[k]<<std::endl;
            column1_status_[k]=VarsTrack_PXB2[i].column1_status[k];
            column2_has_hit_[k]=VarsTrack_PXB2[i].column2_has_hit[k];
            column2_status_[k]=VarsTrack_PXB2[i].column2_status[k];
            
        }
        All_hits_charge.clear();
        All_hits_Px.clear();
        All_hits_Py.clear();
        for (unsigned int k=0; k< VarsTrack_PXB2[i].All_hits_charge.size(); k++)
        {
            std::cout<<VarsTrack_PXB2[i].All_hits_charge[k]<<" qui stampaaa carica " <<std::endl;
            All_hits_charge.push_back(VarsTrack_PXB2[i].All_hits_charge[k]);
            All_hits_Px.push_back(VarsTrack_PXB2[i].All_hits_Px[k]);
            All_hits_Py.push_back(VarsTrack_PXB2[i].All_hits_Py[k]);
        }
        
       
        if (VarsTrack_wPV.size()){
            
            float min_distance=1000;
            int min_dist_index=-1;
            
            for (unsigned int j=0; j<VarsTrack_wPV.size(); j++){
                
                if (std::hypot(VarsTrack_wPV[j].track_local_x-track_local_x_,VarsTrack_wPV[j].track_local_y-track_local_y_ )<min_distance){
                    min_distance=std::hypot(VarsTrack_wPV[j].track_local_x-track_local_x_,VarsTrack_wPV[j].track_local_y-track_local_y_ );
                    min_dist_index=j;
                }
                }
                
                trackfromPV_global_phi_=VarsTrack_wPV[min_dist_index].track_global_phi; trackfromPV_global_z_=VarsTrack_wPV[min_dist_index].track_global_z;
                trackfromPV_local_x_=VarsTrack_wPV[min_dist_index].track_local_x; trackfromPV_local_y_=VarsTrack_wPV[min_dist_index].track_local_y;
                trackfromPV_local_Dx_=VarsTrack_wPV[min_dist_index].track_local_Dx; trackfromPV_local_Dy_=VarsTrack_wPV[min_dist_index].track_local_Dy;
                trackfromPV_localPixel_x_=VarsTrack_wPV[min_dist_index].track_localPixel_x; trackfromPV_localPixel_y_=VarsTrack_wPV[min_dist_index].track_localPixel_y;
                trackfromPV_detid_=VarsTrack_wPV[min_dist_index].detid; trackfromPV_roc_=VarsTrack_wPV[min_dist_index].roc; trackfromPV_detIsActive_=VarsTrack_wPV[min_dist_index].detIsActive;
                trackfromPV_alpha_=VarsTrack_wPV[min_dist_index].alpha; trackfromPV_beta_=VarsTrack_wPV[min_dist_index].beta;
                tree_->Fill();
                
                //per ora stiamo considerando punti incontrati dalla trackfromPV < punti incontrati dalla traccia dal BPIX2
            
        }
        else{
            trackfromPV_global_phi_=-10; trackfromPV_global_z_=-10; trackfromPV_local_x_=-10; trackfromPV_local_y_=-10;
            tree_->Fill();
        } 
    }    
        
    } 
    

    
    if (debug_ > 0) debug_--;
}


//define this as a plug-in
DEFINE_FWK_MODULE(DebugPixelHits_TTbar);
//  edm::ESHandle<MagneticField> theMF;
