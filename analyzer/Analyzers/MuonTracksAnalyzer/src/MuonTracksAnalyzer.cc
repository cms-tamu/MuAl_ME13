// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
//#include "TrackingTools/PatternTools/interface/Trajectory.h"
//#include "TrackingTools/GeomPropagators/interface/Propagator.h"

//
// class declaration
//

class MuonTracksAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MuonTracksAnalyzer(const edm::ParameterSet&);
      ~MuonTracksAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  edm::InputTag m_tracksTag;
  edm::InputTag m_muonsTag;
  double m_minTrackPt;
  double m_minTrackEta;
  double m_maxTrackEta;
  int m_minTrackerHits;
  int m_minDTHits;
  int m_minCSCHits;
  int m_DTWheel;
  int m_DTStation;
  int m_DTSector;
  int m_CSCEndcap;
  int m_CSCStation;
  int m_CSCRing;
  int m_CSCChamber;
  
  int nLayersHit;
  int nInnerHitValid;
  
  edm::Service<TFileService> fs;
  
	TTree * mual_ttree;
	Float_t h_localx;
	Float_t h_localy;
	Float_t h_localz;

	Float_t h_pt;
	Float_t h_eta;
	Float_t h_pz;

    
    Float_t v_hitx[6], v_hity[6], v_hitz[6];
    Float_t v_resx[6], v_resy[6];

    UChar_t t_endcap;
    UChar_t t_station;
    UChar_t t_ring;
    UChar_t t_chamber;

    UInt_t h_nlayers;
    UInt_t h_validTrackerHits;
    
    Int_t h_charge;

};

MuonTracksAnalyzer::MuonTracksAnalyzer(const edm::ParameterSet& iConfig)
  : m_tracksTag(      iConfig.getParameter<edm::InputTag>("tracksTag"))
  , m_muonsTag(       iConfig.getParameter<edm::InputTag>("muonsTag"))
  , m_minTrackPt(     iConfig.getParameter<double>("minTrackPt"))
  , m_minTrackEta(    iConfig.getParameter<double>("minTrackEta"))
  , m_maxTrackEta(    iConfig.getParameter<double>("maxTrackEta"))
  , m_minTrackerHits( iConfig.getParameter<int>("minTrackerHits"))
  , m_minDTHits(      iConfig.getParameter<int>("minDTHits"))
  , m_minCSCHits(     iConfig.getParameter<int>("minCSCHits"))
  , m_DTWheel(        iConfig.getParameter<int>("dtWheel"))
  , m_DTStation(      iConfig.getParameter<int>("dtStation"))
  , m_DTSector(       iConfig.getParameter<int>("dtSector"))
  , m_CSCEndcap(      iConfig.getParameter<int>("cscEndcap"))
  , m_CSCStation(     iConfig.getParameter<int>("cscStation"))
  , m_CSCRing(        iConfig.getParameter<int>("cscRing"))
  , m_CSCChamber(     iConfig.getParameter<int>("cscChamber"))
{
	
	mual_ttree = fs->make<TTree>("mual_ttree", "mual_ttree");
    
    mual_ttree->Branch("localx",&h_localx,"localx/F");
    mual_ttree->Branch("localy",&h_localy,"localy/F");
    mual_ttree->Branch("localz",&h_localz,"localz/F");
    mual_ttree->Branch("pz",&h_pz,"h_pz/F");
    mual_ttree->Branch("pt",&h_pt,"h_pt/F");
    mual_ttree->Branch("eta",&h_eta,"h_eta/F");
    mual_ttree->Branch("nlayers",&h_nlayers,"nlayers/i");
    mual_ttree->Branch("validTrackerHits",&h_validTrackerHits,"validTrackerHits/i");
    mual_ttree->Branch("charge",&h_charge,"charge/i");   


    char buff[20], buff2[20];
    for(int i = 0; i < 6; i++) {
        sprintf(buff, "lay%i_x", i+1);
        sprintf(buff2, "lay%i_x/F", i+1);
        std::cout << buff << " " << buff2 << std::endl;
        mual_ttree->Branch(buff,&v_hitx[i],buff2);
        
        sprintf(buff, "lay%i_y", i+1);
        sprintf(buff2, "lay%i_y/F", i+1);
        mual_ttree->Branch(buff,&v_hity[i],buff2);
        
        sprintf(buff, "lay%i_z", i+1);
        sprintf(buff2, "lay%i_z/F", i+1);
        mual_ttree->Branch(buff,&v_hitz[i],buff2);
        
        sprintf(buff, "lay%i_res_x", i+1);
        sprintf(buff2, "lay%i_res_x/F", i+1);
        mual_ttree->Branch(buff,&v_resx[i],buff2);
        
        sprintf(buff, "lay%i_res_y", i+1);
        sprintf(buff2, "lay%i_res_y/F", i+1);
        mual_ttree->Branch(buff,&v_resy[i],buff2);
    }

    mual_ttree->Branch("endcap",&t_endcap,"endcap/b");
    mual_ttree->Branch("station",&t_station,"station/b");
    mual_ttree->Branch("ring",&t_ring,"ring/b");
    mual_ttree->Branch("chamber",&t_chamber,"chamber/b");
    
}


MuonTracksAnalyzer::~MuonTracksAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


// ------------ method called for each event  ------------
void MuonTracksAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
//  edm::Handle<reco::TrackCollection> trackColl;
//  iEvent.getByLabel(m_tracksTag, trackColl);
  
  edm::Handle<reco::MuonCollection> recoMuonCollection;
  iEvent.getByLabel(m_muonsTag, recoMuonCollection);
  
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByLabel("offlineBeamSpot", beamspot);
  
  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(globalGeometry);

  edm::ESHandle<Propagator> propagator;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagator);

  edm::ESHandle<TrackerGeometry> trackerGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometry);

  edm::ESHandle<CSCGeometry> cscGeometry;
  iSetup.get<MuonGeometryRecord>().get(cscGeometry);

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

    std::vector<const reco::Muon*> recoMuonsSelected;
    std::cout << "recoMuonCollection size: " << recoMuonCollection->size() << std::endl;
    if ( recoMuonCollection.isValid() ) {
	  for (reco::MuonCollection::const_iterator muon = recoMuonCollection->begin();  muon != recoMuonCollection->end();  ++muon) {
		  if ( muon->isGlobalMuon() && muon->isStandAloneMuon() ) {
			  
			  // m_recoMu_pt = muon->pt();
			  // m_recoMu_eta = muon->eta();
			  // m_recoMu_phi = muon->phi();
			  h_pt = muon->pt();
			  h_eta = muon->eta();
			  h_pz = muon->pz();
			  
			  // m_tree_RecoMuons->Fill();
			  
			  if (    muon->globalTrack()->normalizedChi2()   < 10
			       && muon->innerTrack()->numberOfValidHits() > 10
			       // && muon->numberOfMatchedStations() > 1
//				       && fabs( muon->innerTrack()->dxy( beamspot->position() ) ) < 0.2) {
			       && fabs( muon->innerTrack()->dxy( beamspot->position() ) < 0.2 ) ) {
                  
                  recoMuonsSelected.push_back(&*muon);
			  }
		  }
	  }
  }
    
  //  std::cout << "Start loop over recoMuonsSelected" << std::endl;
  //std::cout << "recoMuonsSelected collection size: " << recoMuonsSelected.size() << std::endl;
  for ( unsigned i = 0; i < recoMuonsSelected.size(); ++i ) {
    const reco::Track* innerTrack = recoMuonsSelected[i]->innerTrack().get();
    
    int nInnerHit = 0;
    int nInnerHitValid = 0;
    for ( trackingRecHit_iterator innerHit = innerTrack->recHitsBegin(); innerHit != innerTrack->recHitsEnd();  ++innerHit ) {
      //std::cout << "Inner Hit " << nInnerHit << " is found. Validity: " << (*innerHit)->isValid() << std::endl;
      if( (*innerHit)->isValid() ) {
        
        DetId innerHitId = (*innerHit)->geographicalId();
        
        if ( innerHitId.det() == DetId::Tracker ) {
          /*
          std::cout << "  The inner hit is in Tracker" << std::endl;
          std::cout << "    The hit dimension: " << (*innerHit)->dimension() << std::endl;
          std::cout << "    The hit position (in local coordinates)"
          	                  << " x: " << (*innerHit)->localPosition().x()
          	                  << " y: " << (*innerHit)->localPosition().y() 
          	                  << " z: " << (*innerHit)->localPosition().z()
          	                  << std::endl;
          
          GlobalPoint innerHit_GP = globalGeometry->idToDet(innerHitId)->toGlobal( (*innerHit)->localPosition() );
          std::cout << "    The hit position (in global coordinates)"
          	                  << " x: " << innerHit_GP.x()
          	                  << " y: " << innerHit_GP.y() 
          	                  << " z: " << innerHit_GP.z()
          	                  << std::endl;
          
          // m_innerTrack_hit_X[nInnerHitValid] = innerHit_GP.x();
          // m_innerTrack_hit_Y[nInnerHitValid] = innerHit_GP.y();
          // m_innerTrack_hit_Z[nInnerHitValid] = innerHit_GP.z();
          
          float innerHit_R = sqrt( innerHit_GP.x() * innerHit_GP.x() + innerHit_GP.y() * innerHit_GP.y());
          std::cout << "    The hit global R: " << innerHit_R << std::endl;
          */
        }
        
        nInnerHitValid++;
      }
      
      nInnerHit++;
    }
    
    h_validTrackerHits = nInnerHitValid;
    //std::cout << "Number of valid tracker hits: " << nInnerHitValid << std::endl;
    
    // m_innerTrack_hit_n = nInnerHitValid;
    
    /*std::cout << "Inner track outer hit position (global coordinates)"
                    << " x:" << innerTrack->outerPosition().x()
                    << " y:" << innerTrack->outerPosition().y()
                    << " z:" << innerTrack->outerPosition().z()
                    << std::endl;
    */

    float globxtrack  = innerTrack->outerPosition().x();
    float globytrack  = innerTrack->outerPosition().y();

    float R = sqrt( globxtrack * globxtrack + globytrack * globytrack);
          
    //std::cout << "Inner track outer hit position global R: " << R << std::endl;
    
    int          innerTrack_Charge           = innerTrack->charge();
    DetId        innerTrack_OuterHit_DetId    = DetId(innerTrack->outerDetId());
    GlobalPoint  innerTrack_OuterHit_Position = GlobalPoint(innerTrack->outerPosition().x(), innerTrack->outerPosition().y(), innerTrack->outerPosition().z());
    GlobalVector innerTrack_OuterHit_Momentum = GlobalVector(innerTrack->outerMomentum().x(), innerTrack->outerMomentum().y(), innerTrack->outerMomentum().z());
    
    
    GlobalTrajectoryParameters innerTrack_OuterHit_GlobalTrajParams = GlobalTrajectoryParameters(innerTrack_OuterHit_Position, innerTrack_OuterHit_Momentum, innerTrack_Charge, &(*magneticField));
    
    const reco::Track::CovarianceMatrix innerTrack_OuterHit_CovMatr = innerTrack->outerStateCovariance();
    CurvilinearTrajectoryError innerTrack_OuterHit_CurviTrajErr = CurvilinearTrajectoryError(innerTrack_OuterHit_CovMatr);
    
    TrajectoryStateOnSurface innerTrack_OuterHit_TSOS = TrajectoryStateOnSurface(innerTrack_OuterHit_GlobalTrajParams, innerTrack_OuterHit_CurviTrajErr, trackerGeometry->idToDet(innerTrack_OuterHit_DetId)->surface());
    
    const reco::Track* outerTrack = recoMuonsSelected[i]->outerTrack().get();
    int nOuterHit = 0;
    for ( trackingRecHit_iterator outerHit = outerTrack->recHitsBegin(); outerHit != outerTrack->recHitsEnd();  ++outerHit ) {
      nOuterHit++;
      // std::cout << "Outer Hit " << nOuterHit << " is found. Validity: " << (*outerHit)->isValid() << std::endl;
      if( (*outerHit)->isValid() ) {
        DetId outerHitId = (*outerHit)->geographicalId();
        
        if ( outerHitId.det() == DetId::Muon ) {
          // std::cout << "  The outer hit is in Muon system" << std::endl;
          if ( outerHitId.subdetId() == MuonSubdetId::CSC ) {
            // std::cout << "    The hit is in CSC" << std::endl;
            const CSCDetId cscChamberId(outerHitId.rawId());
            std::cout << "    The hit is in CSC endcap " << cscChamberId.endcap() << " station " << cscChamberId.station() << " ring " << cscChamberId.ring() << std::endl;
            
            if (    cscChamberId.endcap()  != m_CSCEndcap
                 || cscChamberId.station() != m_CSCStation
                 || cscChamberId.ring()    != m_CSCRing
                 || ( m_CSCChamber != -1 && cscChamberId.chamber() != m_CSCChamber ) ) {
                 continue;
                 }
            
            //std::cout << "    The hit dimension: " << (*outerHit)->dimension() << std::endl;
            
            
            
            if ( (*outerHit)->dimension() == 4 ) { // > 1 ?
              std::vector<const TrackingRecHit*> vCSCHits2D = (*outerHit)->recHits();
              //std::cout << "      vCSCHits2D size: " << vCSCHits2D.size() << std::endl;
              
              nLayersHit = 0;
                    
              for(int i = 0; i < 6; i++) {
                  v_hitx[i] = -999.;
                  v_hity[i] = -999.;
                  v_hitz[i] = -999.;
                  v_resx[i] = -999.;
                  v_resy[i] = -999.;
              }
              
              for ( std::vector<const TrackingRecHit*>::const_iterator itCSCHits2D  = vCSCHits2D.begin();
                                                                       itCSCHits2D != vCSCHits2D.end();
                                                                     ++itCSCHits2D ) {
                                                                    
              
                  const TrackingRecHit* cscHit2D = *itCSCHits2D;
                  DetId hitId  = cscHit2D->geographicalId();
                  const CSCDetId cscDetId(hitId.rawId());
                  
                  TrajectoryStateOnSurface extrapolation = propagator->propagate(innerTrack_OuterHit_TSOS, cscGeometry->idToDet(hitId)->surface());

                  nLayersHit++;
                  
                  std::cout << "extrapolated x: " << extrapolation.localPosition().x() << std::endl;
                  std::cout << "             y: " << extrapolation.localPosition().y() << std::endl;
                  std::cout << "             z: " << extrapolation.localPosition().z() << std::endl;
                  std::cout << "hit location x: " << cscHit2D->localPosition().x() << std::endl;
                  std::cout << "             y: " << cscHit2D->localPosition().y() << std::endl;
                  std::cout << "             z: " << cscHit2D->localPosition().z() << std::endl;
                  
                  float resx = extrapolation.localPosition().x() - cscHit2D->localPosition().x();
                  float resy = extrapolation.localPosition().y() - cscHit2D->localPosition().y();
                  
                  // std::cout << "         res x: " << resx << std::endl;
                  // std::cout << "         res y: " << resy << std::endl;
                  
                  v_hitx[cscDetId.layer()-1] = cscHit2D->localPosition().x();
                  v_hity[cscDetId.layer()-1] = cscHit2D->localPosition().y();
                  v_hitz[cscDetId.layer()-1] = cscHit2D->localPosition().z();
                  
                  //residual = propagated track - actual muon hit
                  v_resx[cscDetId.layer()-1] = resx;
                  v_resy[cscDetId.layer()-1] = resy;
              
              }  
              
              h_localx = (*outerHit)->localPosition().x();
              h_localy = (*outerHit)->localPosition().y();
              h_localz = (*outerHit)->localPosition().z();
              h_nlayers = nLayersHit;
              h_charge = outerTrack->charge();
              //std::cout << h_localx << " " << h_localy << " " << h_localz << " layers: " << nLayersHit << std::endl;
         
              t_endcap = cscChamberId.endcap();
              t_station = cscChamberId.station();
              t_ring = cscChamberId.ring();
              t_chamber = cscChamberId.chamber();
              mual_ttree->Fill();
              
            }
            
            
          } else if ( outerHitId.subdetId() == MuonSubdetId::DT ) {
            // std::cout << "    The hit is in CSC" << std::endl;
          } else if ( outerHitId.subdetId() == MuonSubdetId::RPC ) {
            // std::cout << "    The hit is in RPC" << std::endl;
          } else {
            std::cout << "  Error! The Muon hit is NOT in DT, CSC or RPC" << std::endl;
          }
        } else {
          std::cout << "  Error! The outer hit is NOT in muon system" << std::endl;
        }
      }
    }
  }
}
 


// ------------ method called once each job just before starting event loop  ------------
void 
MuonTracksAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonTracksAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MuonTracksAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MuonTracksAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MuonTracksAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MuonTracksAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonTracksAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTracksAnalyzer);


//define this as a plug-in
DEFINE_FWK_MODULE(MuonTracksAnalyzer);
