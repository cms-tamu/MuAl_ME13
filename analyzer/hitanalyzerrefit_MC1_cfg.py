import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


part=1
def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

import formattedfilespyMC
myfilelist = cms.untracked.vstring()
myfilelist.extend( chunks(formattedfilespyMC.fileListMC, 25)[part - 1] )
print len(myfilelist), part
process.source = cms.Source("PoolSource",
        fileNames = myfilelist
)
outFileName = "histoMC" + str(part) + ".root"




#process.options = cms.untracked.PSet(
#     SkipEvent = cms.untracked.vstring('ProductNotFound')    
#)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring("cout"),
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MC_53_V14::All"
#process.GlobalTag.globaltag = "START53_V14::All"

# this line works after 53x
process.load("Configuration.Geometry.GeometryIdeal_cff")
# this line works before 53x
#process.load("Configuration.StandardSequences.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoTracker.TrackProducer.TrackRefitter_cfi")

process.muAlGeneralTracks = process.TrackRefitter.clone()

process.muAlAncientMuonSeed = process.ancientMuonSeed.clone()

process.muAlStandAloneMuons = process.standAloneMuons.clone()
# this line switch on "old" hit based muon reconstruction
#process.muAlStandAloneMuons.STATrajBuilderParameters.BWFilterParameters.MuonTrajectoryUpdatorParameters.Granularity = 2
process.muAlStandAloneMuons.InputObjects = cms.InputTag("muAlAncientMuonSeed")

process.muAlGlobalMuons = process.globalMuons.clone()
process.muAlGlobalMuons.TrackerCollectionLabel = cms.InputTag("muAlGeneralTracks")
process.muAlGlobalMuons.MuonCollectionLabel = cms.InputTag("muAlStandAloneMuons","UpdatedAtVtx")

process.muAlTevMuons = process.tevMuons.clone()
process.muAlTevMuons.MuonCollectionLabel = cms.InputTag("muAlGlobalMuons")

process.muAlGlbTrackQual = process.glbTrackQual.clone()
process.muAlGlbTrackQual.InputCollection = cms.InputTag("muAlGlobalMuons")
process.muAlGlbTrackQual.InputLinksCollection = cms.InputTag("muAlGlobalMuons")

process.muAlMuons = process.muons1stStep.clone()
process.muAlMuons.inputCollectionTypes = cms.vstring('inner tracks','links','outer tracks','tev firstHit','tev picky','tev dyt')
#process.muAlMuons.inputCollectionTypes = cms.vstring('links','outer tracks','tev firstHit','tev picky','tev dyt') # nja
process.muAlMuons.inputCollectionLabels = cms.VInputTag( cms.InputTag("muAlGeneralTracks"),
                                                         cms.InputTag("muAlGlobalMuons"),
                                                         cms.InputTag("muAlStandAloneMuons","UpdatedAtVtx"),
                                                         cms.InputTag("muAlTevMuons","firstHit"),
                                                         cms.InputTag("muAlTevMuons","picky"),
                                                         cms.InputTag("muAlTevMuons","dyt")
                                                       )

process.muAlMuons.fillGlobalTrackQuality = cms.bool(True)
process.muAlMuons.globalTrackQualityInputTag = cms.InputTag('muAlGlbTrackQual')
#process.muAlMuons.fillGlobalTrackRefits = cms.bool(False)

# This is to load new CondDB
from CondCore.DBCommon.CondDBSetup_cfi import *

# Tracker alignment record
#process.TrackerAlignmentInputDB = cms.ESSource("PoolDBESSource",
#                                                   CondDBSetup,
#                                                   connect = cms.string('sqlite_file:alignments_MP1382.db'),
#                                                   toGet = cms.VPSet(cms.PSet(
#                                                        record = cms.string("TrackerAlignmentRcd"),
#                                                        tag = cms.string('Alignments')))
#                                                   )
#process.es_prefer_TrackerAlignmentInputDB = cms.ESPrefer("PoolDBESSource", "TrackerAlignmentInputDB")

#process.TrackerSurfaceDeformationInputDB = cms.ESSource("PoolDBESSource",
#                                                   CondDBSetup,
#                                                   connect = cms.string('sqlite_file:alignments_MP1382.db'),
#                                                   toGet = cms.VPSet(cms.PSet(cms.PSet(
#                                                        record = cms.string("TrackerSurfaceDeformationRcd"),
#                                                        tag = cms.string('Deformations'))))
#                                                   )
#process.es_prefer_TrackerSurfaceDeformationInputDB = cms.ESPrefer("PoolDBESSource", "TrackerSurfaceDeformationInputDB")

#--------------%<----------------------
# New pixel fo mp1375 Tracker alignment candidate - add the fragment below to gather_cfg.py
# --------------%<----------------------
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.myPixelTemplate = cms.ESSource("PoolDBESSource",
#                              CondDBSetup,
#                              connect = cms.string("sqlite_file:siPixelTemplates38T_CMSSW_424p1_2010_2011.db"),
#                              toGet = cms.VPSet(cms.PSet(record = cms.string("SiPixelTemplateDBObjectRcd"),
#                                #tag = cms.string("SiPixelTemplateDBObject38Tv3")
#                                tag = cms.string("SiPixelTemplateDBObject38T17")
#                              ))
#                           )
#process.es_prefer_myPixelTemplate = cms.ESPrefer("PoolDBESSource","myPixelTemplate")

process.MuonTracksAnalyzer = cms.EDAnalyzer("MuonTracksAnalyzer",
  tracksTag      = cms.InputTag("globalMuons"),
#  muonsTag       = cms.InputTag("muons"),  # nja
  muonsTag       = cms.InputTag("muAlMuons"),
  minTrackPt     = cms.double(30.),
  minTrackEta    = cms.double(-2.4),
  maxTrackEta    = cms.double(2.4),
  minTrackerHits = cms.int32(10),
  minDTHits      = cms.int32(1),
  minCSCHits     = cms.int32(1),
  dtWheel        = cms.int32(-3), # = -2,-1, 0, 1, 2; use "-3" to exclude DT
  dtStation      = cms.int32(1),
  dtSector       = cms.int32(1),
  cscEndcap      = cms.int32(1), # = 1, 2; use "0" to exclude CSC
  cscStation     = cms.int32(1),
  cscRing        = cms.int32(3),
  cscChamber     = cms.int32(-1) # = -1 to include ALL chambers
)

process.p = cms.Path(process.muAlGeneralTracks * process.muAlAncientMuonSeed * process.muAlStandAloneMuons * process.muAlGlobalMuons * process.muAlTevMuons * process.muAlGlbTrackQual * process.muAlMuons * process.MuonTracksAnalyzer)
#process.p = cms.Path(process.muAlAncientMuonSeed * process.muAlStandAloneMuons * process.muAlGlobalMuons * process.muAlTevMuons * process.muAlGlbTrackQual * process.muAlMuons * process.MuonTracksAnalyzer)
#process.p = cms.Path(process.muAlMuons * process.MuonTracksAnalyzer)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outFileName)
)

