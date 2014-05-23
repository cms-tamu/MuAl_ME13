# Auto generated configuration file
# using: 
# Revision: 1.381.2.13 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/MinBias_7TeV_cfi.py -s SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO --conditions DESIGN53_V18::All --eventcontent RECO --datatier RECO -n 10 --no_exec --pileup=NoPileUp
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.AlCaRecoStreamsMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOEventContent.outputCommands,
#    outputCommands = cms.untracked.vstring('drop *', 
#        "keep *_genParticles_*_*",
#        "keep *_generator_*_*",
#        "drop edmHepMCProduct_generator_*_*",
#        'keep *_globalMuons_*_*',
#        "keep SiPixelClusteredmNewDetSetVector_siPixelClusters_*_*",
#        "keep SiStripClusteredmNewDetSetVector_siStripClusters_*_*",
#        'keep *_muonCSCDigis_*_*', 
#        'keep *_muonDTDigis_*_*', 
#        'keep *_muonRPCDigis_*_*', 
#        'keep *_dt1DRecHits_*_*', 
#        'keep *_dt2DSegments_*_*', 
#        'keep *_dt4DSegments_*_*', 
#        'keep *_csc2DRecHits_*_*', 
#        'keep *_cscSegments_*_*', 
#        'keep *_rpcRecHits_*_*', 
#        'keep L1AcceptBunchCrossings_*_*_*', 
#        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
#        'keep *_TriggerResults_*_*', 
#        'keep DcsStatuss_scalersRawToDigi_*_*'),
    fileName = cms.untracked.string('singleMuonGun_RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V14::All', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc_GRun', '')

# Generate events
process.generator = cms.EDProducer("SingleMuonGun",
  Verbosity = cms.untracked.int32(0),
  MinPt  = cms.double(30.0),
  MaxPt  = cms.double(200.0),
  MinEta = cms.double(0.8),
  MaxEta = cms.double(1.3),
  MinPhi = cms.double(2.00),
  MaxPhi = cms.double(3.55)
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.reconstruction_step,process.endjob_step,process.RECOoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
