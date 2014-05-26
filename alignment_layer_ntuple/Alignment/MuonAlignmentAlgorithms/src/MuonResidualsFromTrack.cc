/*
 * $Id: MuonResidualsFromTrack.cc,v 1.5 2011/10/12 23:40:24 khotilov Exp $ 
 */

#include "Alignment/MuonAlignmentAlgorithms/interface/MuonResidualsFromTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Alignment/MuonAlignmentAlgorithms/interface/MuonDT13ChamberResidual.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonDT2ChamberResidual.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonCSCChamberResidual.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonTrackDT13ChamberResidual.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonTrackDT2ChamberResidual.h"
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonTrackCSCChamberResidual.h"

// nja
#include "Alignment/MuonAlignmentAlgorithms/interface/CSCTTree.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"

#include "TDecompChol.h"
#include <math.h>

MuonResidualsFromTrack::MuonResidualsFromTrack(const edm::EventSetup& iSetup,
                                               edm::ESHandle<MagneticField> magneticField,
                                               edm::ESHandle<GlobalTrackingGeometry> globalGeometry,
                                               edm::ESHandle<Propagator> prop,
                                               const Trajectory *traj,
                                               const reco::Track* recoTrack,
                                               AlignableNavigator *navigator,
                                               double maxResidual)
  : m_recoTrack(recoTrack)
{
//nja
    //CSCLayerData dummyData;
     //MuonResidualsFromTrack(*iSetup, magneticField, globalGeometry, prop, &traj, &recoTrack, &navigator, 1000., &dummyData);

}

MuonResidualsFromTrack::MuonResidualsFromTrack(const edm::EventSetup& iSetup,
                                               edm::ESHandle<MagneticField> magneticField,
                                               edm::ESHandle<GlobalTrackingGeometry> globalGeometry,
                                               edm::ESHandle<Propagator> prop,
                                               const Trajectory *traj,
                                               const reco::Track* recoTrack,
                                               AlignableNavigator *navigator,
                                               double maxResidual,
                                               struct CSCLayerData * layerData,
                                               TTree * layerTree)
  : m_recoTrack(recoTrack)
{
  std::cout << "BEGIN MuonResidualsFromTrack" << std::endl;
  const std::string metname = " *** MuonResidualsFromTrack *** ";
  LogTrace(metname) << "Tracking Component changed!";
  
  clear();

  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  iSetup.get<TransientRecHitRecord>().get("WithTrackAngle",theTrackerRecHitBuilder);
  reco::TransientTrack track( *m_recoTrack, &*magneticField, globalGeometry );
  TransientTrackingRecHit::ConstRecHitContainer recHitsForRefit;

  layerData->eta = m_recoTrack->eta();
  layerData->pz = m_recoTrack->pz();
  layerData->pt = m_recoTrack->pt();

  int iT = 0, iM = 0;
  int iCSC = 0, iDT = 0;
  for (trackingRecHit_iterator hit = m_recoTrack->recHitsBegin(); hit != m_recoTrack->recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      DetId hitId  = (*hit)->geographicalId();
      if ( hitId.det() == DetId::Tracker ) {
        iT++;
        std::cout << "Tracker Hit " << iT << " is found. Add to refit. Dimension: " << (*hit)->dimension() << std::endl;
        
        recHitsForRefit.push_back( theTrackerRecHitBuilder->build(&**hit) );
      } else if ( hitId.det() == DetId::Muon ){
//        if ( (*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit ) {
//          LogTrace("Reco|TrackingTools|TrackTransformer") << "RPC Rec Hit discarged"; 
//          continue;
//        }
          iM++;
          std::cout << "Muon Hit " << iM << " is found. We do not add muon hits to refit. Dimension: " << (*hit)->dimension() << std::endl;
          if ( hitId.subdetId() == MuonSubdetId::DT ) {
            const DTChamberId chamberId(hitId.rawId());
            iDT++;
            std::cout << "Muon Hit in DT wheel " << chamberId.wheel() << " station " << chamberId.station() << " sector " << chamberId.sector() << "." << std::endl;
          } else if ( hitId.subdetId() == MuonSubdetId::CSC ) {
            const CSCDetId cscDetId(hitId.rawId());
            iCSC++;
            std::cout << "Muon hit in CSC endcap " << cscDetId.endcap() << " station " << cscDetId.station() << " ring " << cscDetId.ring() << " chamber " << cscDetId.chamber() << "." << std::endl;
          } else if ( hitId.subdetId() == MuonSubdetId::RPC ) {
            std::cout << "Muon Hit in RPC" << std::endl;
          } else {
            std::cout << "Warning! Muon Hit not in DT or CSC or RPC" << std::endl;
          }
//        recHitsForRefit.push_back(theMuonRecHitBuilder->build(&**hit));
      }
    }
  }
  layerData->nTracker = iT;
  layerData->nCSC = iCSC;
  layerData->nDT = iDT;
  
//  TrackTransformer trackTransformer();
//  std::vector<Trajectory> vTrackerTrajectory = trackTransformer.transform(track, recHitsForReFit);
//  std::cout << "Tracker trajectories size " << vTrackerTrajectory.size() << std::endl;

  
  TrajectoryStateOnSurface lastTrackerTsos;
  double lastTrackerTsosGlobalPositionR = 0.0;
  
  std::vector<TrajectoryMeasurement> vTrajMeasurement = traj->measurements();
  std::cout << "  Size of vector of TrajectoryMeasurements: " << vTrajMeasurement.size() << std::endl;
  int nTrajMeasurement = 0;
  for ( std::vector<TrajectoryMeasurement>::const_iterator iTrajMeasurement =  vTrajMeasurement.begin();
                                                           iTrajMeasurement != vTrajMeasurement.end();
                                                         ++iTrajMeasurement ) {
    nTrajMeasurement++;
    std::cout << "    TrajectoryMeasurement #" << nTrajMeasurement << std::endl;
    
    TrajectoryMeasurement trajMeasurement = *iTrajMeasurement;
    
    TrajectoryStateOnSurface tsos = m_tsoscomb(trajMeasurement.forwardPredictedState(), trajMeasurement.backwardPredictedState());
    TrajectoryStateOnSurface tsosF = trajMeasurement.forwardPredictedState();
    TrajectoryStateOnSurface tsosB = trajMeasurement.backwardPredictedState();
    TrajectoryStateOnSurface tsosU = trajMeasurement.updatedState();
    std::cout << "      TrajectoryMeasurement TSOS validity: " << tsos.isValid() << std::endl;
    if ( tsos.isValid() ) {
      double tsosGlobalPositionR = sqrt( tsos.globalPosition().x()*tsos.globalPosition().x() + tsos.globalPosition().y()*tsos.globalPosition().y() );
      std::cout << "         TrajectoryMeasurement TSOS localPosition"
                << " x: " << tsos.localPosition().x()
                << " y: " << tsos.localPosition().y()
                << " z: " << tsos.localPosition().z()
                << std::endl;
      std::cout << "         TrajectoryMeasurement TSOS globalPosition"
                << " x: " << tsos.globalPosition().x()
                << " y: " << tsos.globalPosition().y()
                << " R: " << tsosGlobalPositionR
                << " z: " << tsos.globalPosition().z()
                << std::endl;
      if ( tsosGlobalPositionR > lastTrackerTsosGlobalPositionR ) {
        lastTrackerTsos = tsos;
        lastTrackerTsosGlobalPositionR = tsosGlobalPositionR;
      }
    }
    
    const TransientTrackingRecHit *trajMeasurementHit = &(*trajMeasurement.recHit());
    std::cout << "      TrajectoryMeasurement hit validity: " << trajMeasurementHit->isValid() << std::endl;
    if ( trajMeasurementHit->isValid() ) {
      DetId trajMeasurementHitId  = trajMeasurementHit->geographicalId();
      int   trajMeasurementHitDim = trajMeasurementHit->dimension();
      if ( trajMeasurementHitId.det() == DetId::Tracker ) {
        std::cout << "      TrajectoryMeasurement hit Det: Tracker" << std::endl;
        std::cout << "      TrajectoryMeasurement hit dimension: " << trajMeasurementHitDim << std::endl;
        m_tracker_numHits++;
        double xresid     = tsos.localPosition().x()               - trajMeasurementHit->localPosition().x();
        double xresiderr2 = tsos.localError().positionError().xx() + trajMeasurementHit->localPositionError().xx();
        m_tracker_chi2 += xresid * xresid / xresiderr2;
        
        
        
	      if (    trajMeasurementHitId.subdetId() == StripSubdetector::TID
	           || trajMeasurementHitId.subdetId() == StripSubdetector::TEC) {
          m_contains_TIDTEC = true;
        }
        
      } else {
        std::cout << "      TrajectoryMeasurement hit det: UNKNOWN" << std::endl;
      }
    }
  }    
 
  // Should we select this track for layer debugging?
  // Make sure it meets the tracker requirements
  layerData->select = false;
  if( m_tracker_numHits >= 15 ) { //15 is minTrackerHits
      if( normalizedChi2() < 10 ) { //10 is maxTrackerRedChi2
         layerData->select = true;
      }
  }

  int iT2 = 0, iM2 = 0;
  //int trackerHits = 0, cscHits = 0, dtHits = 0; // nja
  for (trackingRecHit_iterator hit2 = m_recoTrack->recHitsBegin(); hit2 != m_recoTrack->recHitsEnd(); ++hit2) {
    if((*hit2)->isValid()) {
      DetId hitId2  = (*hit2)->geographicalId();
      if ( hitId2.det() == DetId::Tracker ) {
        iT2++;
        //trackerHits++; //nja
        //layerData->nTracker = trackerHits;
        std::cout << "Tracker Hit " << iT2 << " is found. We don't calcualte Tsos for it" << std::endl;
      } else if ( hitId2.det() == DetId::Muon ){
          iM2++;
          std::cout << "Muon Hit " << iM2 << " is found. Dimension: " << (*hit2)->dimension() << std::endl;
          if ( hitId2.subdetId() == MuonSubdetId::DT ) {
              //dtHits++; // nja
            const DTChamberId chamberId(hitId2.rawId());
            std::cout << "Muon Hit in DT wheel " << chamberId.wheel() << " station " << chamberId.station() << " sector " << chamberId.sector() << std::endl;
            
            
            
            if ( (*hit2)->dimension() > 1 ) {
            std::vector<const TrackingRecHit*> vDTSeg2D = (*hit2)->recHits();
            std::cout << "          vDTSeg2D size: " << vDTSeg2D.size() << std::endl;
            for ( std::vector<const TrackingRecHit*>::const_iterator itDTSeg2D =  vDTSeg2D.begin();
                                                                     itDTSeg2D != vDTSeg2D.end();
                                                                   ++itDTSeg2D ) {
              std::vector<const TrackingRecHit*> vDTHits1D =  (*itDTSeg2D)->recHits();
              std::cout << "            vDTHits1D size: " << vDTHits1D.size() << std::endl;
              for ( std::vector<const TrackingRecHit*>::const_iterator itDTHits1D =  vDTHits1D.begin();
                                                                       itDTHits1D != vDTHits1D.end();
                                                                     ++itDTHits1D ) {
                const TrackingRecHit* hit = *itDTHits1D;
                std::cout << "              hit dimension: " << hit->dimension() << std::endl;
                
                DetId hitId  = hit->geographicalId();
                const DTSuperLayerId superLayerId(hitId.rawId());
	              const DTLayerId layerId(hitId.rawId());
	              std::cout << "              hit superLayerId: " << superLayerId.superLayer() << std::endl;
	              std::cout << "              hit layerId: " << layerId.layer() << std::endl;
                
                if ( superLayerId.superlayer() == 2  && vDTHits1D.size() >= 3 ) {
                  if ( m_dt2.find(chamberId) == m_dt2.end() ) {
                    AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
                    m_dt2[chamberId] = new MuonDT2ChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
	                  std::cout << "              This is first appearance of the DT with hits in superlayer 2" << std::endl;
	                  
	                  // have we seen this chamber before? check if it was in dt13
	                  if ( m_dt13.find(chamberId) == m_dt13.end() ) {
	                    m_chamberIds.push_back(chamberId);
	                  }
	                  
                  }
                  
                  TrajectoryStateOnSurface extrapolation;
                  extrapolation = prop->propagate( lastTrackerTsos, globalGeometry->idToDet(hitId)->surface() );
                  
                  if ( extrapolation.isValid() ) { 
            	    std::cout << " extrapolation localPosition()"
            	              << " x: " << extrapolation.localPosition().x()
            	              << " y: " << extrapolation.localPosition().y() 
            	              << " z: " << extrapolation.localPosition().z() << std::endl;
                    m_dt2[chamberId]->addResidual(prop, &extrapolation, hit);
                  }
//            	    residualDT2IsAdded = true;
            	    
                } else if ( (superLayerId.superlayer() == 1 || superLayerId.superlayer() == 3) && vDTHits1D.size() >= 6 ) {
                  if ( m_dt13.find(chamberId) == m_dt13.end() ) {
                    AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
                    m_dt13[chamberId] = new MuonDT13ChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
                    std::cout << "              This is first appearance of the DT with hits in superlayers 1 and 3" << std::endl;
                    
                    // have we seen this chamber before? check if it was in dt2
	                  if ( m_dt2.find(chamberId) == m_dt2.end() ) {
	                    m_chamberIds.push_back(chamberId);
	                  }
            	    }
            	    
            	    TrajectoryStateOnSurface extrapolation;
                  extrapolation = prop->propagate( lastTrackerTsos, globalGeometry->idToDet(hitId)->surface() );
            	    
            	    if ( extrapolation.isValid() ) {
            	    std::cout << " extrapolation localPosition()"
            	              << " x: " << extrapolation.localPosition().x()
            	              << " y: " << extrapolation.localPosition().y() 
            	              << " z: " << extrapolation.localPosition().z() << std::endl;
                    m_dt13[chamberId]->addResidual(prop, &extrapolation, hit);
                  }
//            	    residualDT13IsAdded = true;
            	    
            	    
            	  }
              }
            }
          }
            
            


          } else if ( hitId2.subdetId() == MuonSubdetId::CSC ) {
              //cscHits++; // nja
              //layerData->nCSC = cscHits;
            const CSCDetId cscDetId2(hitId2.rawId());
            const CSCDetId chamberId(cscDetId2.endcap(), cscDetId2.station(), cscDetId2.ring(), cscDetId2.chamber());
            std::cout << "Muon hit in CSC endcap " << cscDetId2.endcap() << " station " << cscDetId2.station() << " ring " << cscDetId2.ring() << " chamber " << cscDetId2.chamber() << "." << std::endl;
         
            layerData->endcap = cscDetId2.endcap();
            layerData->station = cscDetId2.station();
            layerData->ring = cscDetId2.ring();
            layerData->chamber = cscDetId2.chamber();
            
            if ( (*hit2)->dimension() == 4 ) {
            std::vector<const TrackingRecHit*> vCSCHits2D = (*hit2)->recHits();
            std::cout << "          vCSCHits2D size: " << vCSCHits2D.size() << std::endl;
            if ( vCSCHits2D.size() >= 5 ) {
               
               
                int nLayers = 0;
            

                    for (int i = 0; i < 6; i++) { 
                        layerData->v_hitx[i] = -999.;
                        layerData->v_hity[i] = -999.;

                        layerData->v_resx[i] = -999.;
                        layerData->v_resy[i] = -999.;
                    }


              for ( std::vector<const TrackingRecHit*>::const_iterator itCSCHits2D =  vCSCHits2D.begin();
                                                                       itCSCHits2D != vCSCHits2D.end();
                                                                     ++itCSCHits2D ) {
                const TrackingRecHit* cscHit2D = *itCSCHits2D;
                std::cout << "            cscHit2D dimension: " << cscHit2D->dimension() << std::endl;
                const TrackingRecHit* hit = cscHit2D;
                std::cout << "              hit dimension: " << hit->dimension() << std::endl;
                
                DetId hitId  = hit->geographicalId();
                const CSCDetId cscDetId(hitId.rawId());
                std::cout << "              hit layer: " << cscDetId.layer() << std::endl;
                
                std::cout << " hit localPosition"
            	            << " x: " << hit->localPosition().x()
            	            << " y: " << hit->localPosition().y() 
            	            << " z: " << hit->localPosition().z()
            	            << std::endl;
            	  std::cout << " hit globalPosition"
            	            << " x: " << globalGeometry->idToDet(hitId)->toGlobal(hit->localPosition()).x()
            	            << " y: " << globalGeometry->idToDet(hitId)->toGlobal(hit->localPosition()).y()
            	            << " z: " << globalGeometry->idToDet(hitId)->toGlobal(hit->localPosition()).z()
            	            << std::endl;
                
                // not sure why we sometimes get layer == 0
                if (cscDetId.layer() == 0) continue;
  
                // have we seen this chamber before?
                std::cout << "Have we seen this chamber before?";
                if ( m_csc.find(chamberId) == m_csc.end() ) {
                  std::cout << " NO. m_csc.count() = " << m_csc.count(chamberId) << std::endl;
                  AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
                  m_csc[chamberId] = new MuonCSCChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
	                std::cout << "              This is first appearance of the CSC with hits m_csc.count() = " << m_csc.count(chamberId) << std::endl;
	                m_chamberIds.push_back(chamberId);
                  //addTrkCovMatrix(chamberId, tsos); // only for the 1st hit
            	  } else {
            	    std::cout << " YES. m_csc.count() = " << m_csc.count(chamberId) << std::endl;
            	  }
            	  
            	  std::cout << " lastTrackerTsos localPosition"
            	            << " x: " << lastTrackerTsos.localPosition().x()
            	            << " y: " << lastTrackerTsos.localPosition().y() 
            	            << " z: " << lastTrackerTsos.localPosition().z()
            	            << std::endl;
            	  std::cout << " lastTrackerTsos globalPosition"
            	            << " x: " << lastTrackerTsos.globalPosition().x()
            	            << " y: " << lastTrackerTsos.globalPosition().y() 
            	            << " z: " << lastTrackerTsos.globalPosition().z()
            	            << std::endl;
            	  std::cout << " Do extrapolation from lastTrackerTsos to hit surface" << std::endl;
            	  TrajectoryStateOnSurface extrapolation;
                extrapolation = prop->propagate( lastTrackerTsos, globalGeometry->idToDet(hitId)->surface() );
            	  std::cout << " extrapolation.isValid() = " << extrapolation.isValid() << std::endl;
            	  
            	  if ( extrapolation.isValid() ) {
            	    std::cout << " extrapolation localPosition()"
            	              << " x: " << extrapolation.localPosition().x()
            	              << " y: " << extrapolation.localPosition().y() 
            	              << " z: " << extrapolation.localPosition().z() << std::endl;
                  m_csc[chamberId]->addResidual(prop, &extrapolation, hit);
                }


                  //nja
                      layerData->charge = m_recoTrack->charge();
                      layerData->v_hitx[cscDetId.layer()-1] = hit->localPosition().x();
                      layerData->v_hity[cscDetId.layer()-1] = hit->localPosition().y();

                      if(extrapolation.isValid() ) {
                          layerData->v_resx[cscDetId.layer()-1] = extrapolation.localPosition().x() - hit->localPosition().x();
                          layerData->v_resy[cscDetId.layer()-1] = extrapolation.localPosition().y() - hit->localPosition().y();
                      }

                  nLayers++;


              } // end of loop over CSC layers
                
              layerData->nlayers = nLayers;
              //layerData->nDT = dtHits;
              //dtHits = 0;

            } else { /*(*layerData).select = false;*/ }
          }
          
          layerTree->Fill();  
            
          } else if ( hitId2.subdetId() == MuonSubdetId::RPC ) {
            std::cout << "Muon Hit in RPC" << std::endl;
          } else {
            std::cout << "Warning! Muon Hit not in DT or CSC or RPC" << std::endl;
          }
//        recHitsForRefit.push_back(theMuonRecHitBuilder->build(&**hit));
          if ( hitId2.subdetId() == MuonSubdetId::DT || hitId2.subdetId() == MuonSubdetId::CSC ) {
            
          }
      }
    }
  }
  
  
  //layerTree->Fill();
  std::cout << "END MuonResidualsFromTrack" << std::endl << std::endl;
}


MuonResidualsFromTrack::MuonResidualsFromTrack( edm::ESHandle<GlobalTrackingGeometry> globalGeometry,
                                                const reco::Muon *recoMuon,
                                                AlignableNavigator *navigator,
                                                double maxResidual )
  : m_recoMuon(recoMuon)
{
  clear();
  assert( m_recoMuon->isTrackerMuon() && m_recoMuon->innerTrack().isNonnull());
  m_recoTrack = m_recoMuon->innerTrack().get();
  
  m_tracker_chi2 = m_recoMuon->innerTrack()->chi2();
  m_tracker_numHits = m_recoMuon->innerTrack()->ndof() + 5;
  m_tracker_numHits = m_tracker_numHits > 0 ? m_tracker_numHits : 0 ;
  
  /*
  for (trackingRecHit_iterator hit = m_recoMuon->innerTrack()->recHitsBegin();  hit != m_recoMuon->innerTrack()->recHitsEnd();  ++hit)
  {
    DetId id = (*hit)->geographicalId();
    if (id.det() == DetId::Tracker)
    {
      m_tracker_numHits++;
      if (id.subdetId() == StripSubdetector::TID  ||  id.subdetId() == StripSubdetector::TEC) m_contains_TIDTEC = true;
    }
  }
  */
  
  for (std::vector<reco::MuonChamberMatch>::const_iterator chamberMatch = m_recoMuon->matches().begin();  
       chamberMatch != m_recoMuon->matches().end();  chamberMatch++)
  {
    if (chamberMatch->id.det() != DetId::Muon ) continue;
    
    for (std::vector<reco::MuonSegmentMatch>::const_iterator segMatch = chamberMatch->segmentMatches.begin();
         segMatch != chamberMatch->segmentMatches.end();  ++segMatch)
    {
      // select the only segment that belongs to track and is the best in station by dR
      if (! (segMatch->isMask(reco::MuonSegmentMatch::BestInStationByDR) &&
             segMatch->isMask(reco::MuonSegmentMatch::BelongsToTrackByDR)) ) continue;
      
      if (chamberMatch->id.subdetId() == MuonSubdetId::DT)
      {
        const DTChamberId chamberId(chamberMatch->id.rawId());

        DTRecSegment4DRef segmentDT = segMatch->dtSegmentRef;
        const DTRecSegment4D* segment = segmentDT.get();
        if (segment == 0)  continue;
        
        if ( segment->hasPhi()  &&  fabs(chamberMatch->x - segMatch->x) > maxResidual ) continue;
        if ( segment->hasZed()  &&  fabs(chamberMatch->y - segMatch->y) > maxResidual ) continue;
          
        // have we seen this chamber before?
        if (m_dt13.find(chamberId) == m_dt13.end()  &&  m_dt2.find(chamberId) == m_dt2.end()) {
          m_chamberIds.push_back(chamberId);
        }

        if (segment->hasZed())
        {
          if (m_dt2.find(chamberId) == m_dt2.end())
          {
            AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
// YP
//            m_dt2[chamberId] = new MuonTrackDT2ChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
          }
          else std::cout<<"multi segment match to tmuon: dt2  -- should not happen!"<<std::endl;
          m_dt2[chamberId]->setSegmentResidual(&(*chamberMatch), &(*segMatch));
        }
        if (segment->hasPhi())
        {
          if (m_dt13.find(chamberId) == m_dt13.end())
          {
            AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
// YP
//            m_dt13[chamberId] = new MuonTrackDT13ChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
          }
          else std::cout<<"multi segment match to tmuon: dt13  -- should not happen!"<<std::endl;
          m_dt13[chamberId]->setSegmentResidual(&(*chamberMatch), &(*segMatch));
        }
      }

      else if (chamberMatch->id.subdetId() == MuonSubdetId::CSC) 
      {
        const CSCDetId cscDetId(chamberMatch->id.rawId());
        const CSCDetId chamberId(cscDetId.chamberId());

        if ( fabs(chamberMatch->x - segMatch->x) > maxResidual ) continue;

        // have we seen this chamber before?
        if (m_csc.find(chamberId) == m_csc.end())
        {
          m_chamberIds.push_back(chamberId);
          AlignableDetOrUnitPtr chamberAlignable = navigator->alignableFromDetId(chamberId);
// YP
//          m_csc[chamberId] = new MuonTrackCSCChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable);
        }
        else std::cout<<"multi segment match to tmuon: csc  -- should not happen!"<<std::endl;
        m_csc[chamberId]->setSegmentResidual(&(*chamberMatch), &(*segMatch));
      }

    }
  }
}

// This is destructor
// It deletes all chambers residulas
MuonResidualsFromTrack::~MuonResidualsFromTrack() {
  for (std::map<DetId,MuonChamberResidual*>::const_iterator residual = m_dt13.begin(); residual != m_dt13.end(); ++residual) {
    delete residual->second;
  }
  for (std::map<DetId,MuonChamberResidual*>::const_iterator residual = m_dt2.begin();  residual != m_dt2.end();  ++residual) {
    delete residual->second;
  }
  for (std::map<DetId,MuonChamberResidual*>::const_iterator residual = m_csc.begin();  residual != m_csc.end();  ++residual) {
    delete residual->second;
  }
}


void MuonResidualsFromTrack::clear() {
  m_tracker_numHits = 0;
  m_tracker_chi2 = 0.;
  m_contains_TIDTEC = false;
  m_chamberIds.clear();
  m_dt13.clear();
  m_dt2.clear();
  m_csc.clear();
  m_trkCovMatrix.clear();
}


double MuonResidualsFromTrack::trackerRedChi2() const 
{
  if (m_tracker_numHits > 5) return m_tracker_chi2 / double(m_tracker_numHits - 5);
  else return -1.;
}


double MuonResidualsFromTrack::normalizedChi2() const 
{
  if (m_recoMuon) return m_recoTrack->normalizedChi2();
  return trackerRedChi2();
}


MuonChamberResidual * MuonResidualsFromTrack::chamberResidual(DetId chamberId, int type)
{
  if (type == MuonChamberResidual::kDT13) {
    if (m_dt13.find(chamberId) == m_dt13.end()) return NULL;
    return m_dt13[chamberId];
  }
  else if (type == MuonChamberResidual::kDT2) {
    if (m_dt2.find(chamberId) == m_dt2.end()) return NULL;
    return m_dt2[chamberId];
  }
  else if (type == MuonChamberResidual::kCSC) {
    if (m_csc.find(chamberId) == m_csc.end()) return NULL;
    return m_csc[chamberId];
  }
  else return NULL;
}

void MuonResidualsFromTrack::addTrkCovMatrix(DetId chamberId, TrajectoryStateOnSurface &tsos)
{
  const AlgebraicSymMatrix55 cov55 = tsos.localError().matrix();
  TMatrixDSym cov44(4);
  // change indices from q/p,dxdz,dydz,x,y   to   x,y,dxdz,dydz
  int subs[4] = { 3, 4, 1, 2 };
  for (int i=0;i<4;i++) for (int j=0;j<4;j++)  cov44(i,j) = cov55( subs[i], subs[j] );
  m_trkCovMatrix[chamberId] = cov44;
}

TMatrixDSym MuonResidualsFromTrack::covMatrix(DetId chamberId)
{
  TMatrixDSym result(4);
  std::cout<<"MuonResidualsFromTrack:: cov initial:"<<std::endl;
  result.Print();
  if (m_trkCovMatrix.find(chamberId) == m_trkCovMatrix.end())
  {
    std::cout<<"MuonResidualsFromTrack:: cov does not exist!"<<std::endl;
    return result;
  }
  result = m_trkCovMatrix[chamberId];

  std::cout<<"MuonResidualsFromTrack:: cov before:"<<std::endl;
  result.Print();

  // add segment's errors in quadratures to track's covariance matrix
  double r_err;
  if (m_csc.find(chamberId) == m_csc.end())
  {
    r_err = m_csc[chamberId]->residual_error();
    result(0,0) += r_err*r_err;
    r_err = m_csc[chamberId]->resslope_error();
    result(2,2) += r_err*r_err;
  }
  if (m_dt13.find(chamberId) == m_dt13.end())
  {
    r_err = m_dt13[chamberId]->residual_error();
    result(0,0) += r_err*r_err;
    r_err = m_dt13[chamberId]->resslope_error();
    result(2,2) += r_err*r_err;
  }
  if (m_dt2.find(chamberId) == m_dt2.end())
  {
    r_err = m_dt2[chamberId]->residual_error();
    result(1,1) += r_err*r_err;
    r_err = m_dt2[chamberId]->resslope_error();
    result(3,3) += r_err*r_err;
  }
  std::cout<<"MuonResidualsFromTrack:: cov after:"<<std::endl;
  result.Print();

  return result;
}



TMatrixDSym MuonResidualsFromTrack::corrMatrix(DetId chamberId)
{
  TMatrixDSym result(4);
  TMatrixDSym cov44 = covMatrix(chamberId);

  // invert it using cholesky decomposition
  TDecompChol decomp(cov44);
  bool ok = decomp.Invert(result);
  std::cout<<"MuonResidualsFromTrack:: corr after:"<<std::endl;
  result.Print();

  if (!ok){std::cout<<"MuonResidualsFromTrack:: cov inversion failed!"<<std::endl;}
  return result;
}

TMatrixD MuonResidualsFromTrack::choleskyCorrMatrix(DetId chamberId)
{
  TMatrixD result(4,4);
  TMatrixDSym corr44 = corrMatrix(chamberId);

  // get an upper triangular matrix U such that corr = U^T * U
  TDecompChol decomp(corr44);
  bool ok = decomp.Decompose();
  result = decomp.GetU();

  std::cout<<"MuonResidualsFromTrack:: corr cholesky after:"<<std::endl;
  result.Print();

  if (!ok){std::cout<<"MuonResidualsFromTrack:: corr decomposition failed!"<<std::endl;}
  return result;
}
