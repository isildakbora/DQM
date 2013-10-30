import FWCore.ParameterSet.Config as cms

process = cms.Process("DataScouting")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'GR_R_52_V9::All'
process.GlobalTag.globaltag = 'GR_R_53_V12::All'  


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.default.limit = 100
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                      #SkipEvent = cms.untracked.vstring('ProductNotFound')
                                    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/data/Run2012B/DataScouting/RAW/v1/000/194/535/16D13DD4-CBA2-E111-AE6F-001D09F24353.root',
                                      '/store/data/Run2012B/DataScouting/RAW/v1/000/194/533/FADCCE72-C5A2-E111-825D-003048D2BBF0.root',
                                      '/store/data/Run2012B/DataScouting/RAW/v1/000/194/533/EC8BE038-9DA2-E111-AEEC-00215AEDFD98.root',
                                      '/store/data/Run2012B/DataScouting/RAW/v1/000/194/533/E865F95C-B7A2-E111-9FFF-003048D2BC5C.root',
                                      '/store/data/Run2012B/DataScouting/RAW/v1/000/194/533/E63A7608-99A2-E111-BAD9-001D09F290BF.root')
                                       #'/afs/cern.ch/work/i/isildak/FADCCE72-C5A2-E111-825D-003048D2BBF0.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
)

process.DQMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.DQMEventContent.outputCommands,
    fileName = cms.untracked.string('DataScouting_DQM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('DQM')
    )
)
#process.DQMoutput.outputCommands.append('keep *_*_*_*')

process.dqmSaver.workflow = '/DataScouting/DQM/Test'
process.dqmsave_step = cms.Path(process.dqmSaver)

process.load('DQM.DataScouting.dataScouting_cff')

# Path and EndPath definitions
process.testmodule_step = cms.Path(process.dataScoutingDQMSequence)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(process.testmodule_step,process.dqmsave_step,process.DQMoutput_step)

