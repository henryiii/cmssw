# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 --conditions auto:run2_data -s RAW2DIGI,RECO --process reRECO --mc --era Run2_2018 --eventcontent AOD --scenario pp --datatier AOD --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2018 -n 100
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('reRECO',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

files = """
206A6E9A-4DB2-1941-A60B-7174FA398D86.root
23BAC38C-A5CD-4844-8EA5-C7B5AA443861.root
3C001F42-8E40-7B41-BA3B-780714F6ABBD.root
4978A440-0E12-8541-86B8-81086CDC98A0.root
5163175A-C833-534E-9197-A275A38A78E2.root
60FF74A6-2B35-1243-BDCB-282EC3A1D961.root
6D6E7EAA-5553-C440-B6BD-E8FF165CADC2.root
6EEA5B7D-3D49-5F42-856F-AB11593D6EB5.root
715C11D3-94C5-5945-A295-DD120D8697F4.root
76E5A405-B666-444B-A086-C9D102519D96.root
834227F5-C683-AC43-86C8-96A0642376DF.root
83560C49-46DD-9C49-AB5A-6D72053E8CC3.root
86C39A48-D56F-3741-929C-6C72DCA015D7.root
8727F95B-0A35-D24B-A957-81391C5EE7D9.root
95F7275E-3A4D-CF42-9F6C-84F4EF4ECAFF.root
9DA43473-81B2-EB42-BF52-FE7C1DE0857B.root
9E4F2FEA-B27D-0E47-AE5B-74D01EE8CE51.root
9FBCC652-F2E5-1348-AB39-C74B46F2BB8D.root
A9FB61CD-4A96-E645-B56B-AA8037DB9118.root
AD86C759-277C-124F-95BA-F4DD7BA838AD.root
B2E83B4E-5136-6240-B6C5-C36B1728716C.root
BE332030-68F5-2A46-97D5-F851A69AB3BA.root
CD903778-801D-CA4A-B73A-ACAA54BA3D74.root
CE277B72-55FC-6A42-98DE-AB9CB048A517.root
D2FF2B84-7404-614C-961A-AC0BE5E98F65.root
D4B3768F-C95E-E249-ADDC-A608EDD26EB8.root
E0802B57-98A4-DA4E-96BF-F7D90DF1B786.root
E0968F76-5E44-0A4C-AC6F-1222662E57E0.root
E7AF7B47-7A09-A441-B255-DCE8D46C954E.root
E9740D7E-F61D-4D40-A1E7-6A8A6D647D94.root
EDD5D9BB-5E52-454A-AF1A-63391379C76A.root
F68115D9-E295-DA4E-8DB0-15520D79107F.root
""".strip().split()

name = 'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2018/RunIIAutumn18DR/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PUAvg50IdealConditions_IdealConditions_102X_upgrade2018_design_v9_ext1-v2/260000/'

# Input source - all
readFiles = cms.untracked.vstring()
readFiles.extend([name + fn for fn in files])

# Input source - classic
#readFiles = cms.untracked.vstring(name + files[14])

process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    secondaryFileNames = cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(8),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    fileName = cms.untracked.string('step3_RAW2DIGI_RECO.root'),
    outputCommands = process.AODEventContent.outputCommands
)

process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("RecoVertex.PrimaryVertexProducer.pvDumper_cfi")

process.quickTrackAssociatorByHits = cms.EDProducer("QuickTrackAssociatorByHitsProducer",
                                                    AbsoluteNumberOfHits = cms.bool(False),
                                                    Cut_RecoToSim = cms.double(0.75),
                                                    SimToRecoDenominator = cms.string('reco'), # either "sim" or "reco"
                                                    Quality_SimToReco = cms.double(0.5),
                                                    Purity_SimToReco = cms.double(0.75),
                                                    ThreeHitTracksAreSpecial = cms.bool(True),
                                                    PixelHitWeight = cms.double(1.0),
                                                    useClusterTPAssociation = cms.bool(True),
                                                    cluster2TPSrc = cms.InputTag("tpClusterProducer")
                                                )


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_design', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)

process.truth_match_step = cms.Path(process.tpClusterProducer+process.quickTrackAssociatorByHits+process.pvDumper)

myreco= cms.Sequence(process.globalreco_trackingTask)

process.reconstruction_step= cms.Path(process.reconstruction_trackingOnly)
#process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODoutput_step = cms.EndPath(process.AODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.truth_match_step,process.endjob_step,process.AODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2018

#call to customisation function customisePostEra_Run2_2018 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2018(process)

# End of customisation functions


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
