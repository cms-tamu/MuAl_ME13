import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
        '/store/caf/user/pakhotin/singleMuonGun_RECO_v4/singleMuonGun_RECO_v4_Consolidated_v3/8f04993aadaf700f87a5ec0188006dcc/singleMuonGun_RECO_Consolidated_100_1_kAB.root'
    )
)


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))

process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring("cout"),
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string("INFO")))

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_53_V16A::All"

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.MuonTracksAnalyzer = cms.EDAnalyzer("MuonTracksAnalyzer",
  tracksTag      = cms.InputTag("globalMuons"),
  muonsTag       = cms.InputTag("muons"),
  minTrackPt     = cms.double(20.),
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
  cscRing        = cms.int32(3)
)


# process.output = cms.OutputModule("PoolOutputModule",
  # fileName = cms.untracked.string("output.root"), #SingleMu_Run2012A_MuAlCalIsolatedMu-13Jul2012-v1_MEp_1_3_1.root"),
  # outputCommands = cms.untracked.vstring( "drop *" ) 

  # #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("Path"))
# )


process.p = cms.Path(process.MuonTracksAnalyzer)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("histo.root")
  )

# process.EndPath = cms.EndPath(process.output)
