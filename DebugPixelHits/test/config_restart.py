import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.suppressError = cms.untracked.vstring("patTriggerFull")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )    

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_mc2017_realistic_v5', '')


import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi 

#old
process.HitCollectorForDebug = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('HitCollectorForDebug'),
    MaxChi2 = cms.double(30.0), ## was 30 ## TO BE TUNED
    nSigma  = cms.double(3.),   ## was 3  ## TO BE TUNED 
)


process.clusterInfo = cms.EDAnalyzer("DebugPixelHits_v1",
        
        tracker = cms.InputTag("MeasurementTrackerEvent"),
        vertices = cms.InputTag("offlinePrimaryVertices"),
        lumiScalers = cms.InputTag("scalersRawToDigi"),
        
        # configuraton for refitter
        DoPredictionsOnly = cms.bool(False),
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithAngleAndTemplate'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        RefitDirection = cms.string('oppositeToMomentum'),
        RefitRPCHits = cms.bool(True),
        Propagator = cms.string('SmartPropagatorAnyRKOpposite'),
        
        #Propagators
        PropagatorOpposite = cms.string("RungeKuttaTrackerPropagatorOpposite"),
        PropagatorAlong = cms.string("RungeKuttaTrackerPropagator"),
        Chi2MeasurementEstimator = cms.string("HitCollectorForDebug"),
        
        #Error rescaling
        rescaleError = cms.double(1),

        badComponentsFile = cms.string('/afs/cern.ch/work/g/gpetrucc/Tracking/CMSSW_9_2_3_patch2/src/RecoTracker/DebugTools/test/badComponents.txt'),
        
        debug = cms.untracked.int32(100)
)

process.tagAndProbe = cms.Path( 
    process.MeasurementTrackerEvent +
    process.clusterInfo 
)

#process.out = cms.OutputModule("PoolOutputModule",
    #fileName = cms.untracked.string("debug_Zmm_lostHits.root"),
    #outputCommands = cms.untracked.vstring("keep *", "drop *_*_*_TagProbe"),
    #SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("tagAndProbe")),
#)

#uncomment to skim the selected events
#process.end = cms.EndPath(process.out)


del process.clusterInfo.badComponentsFile
process.TFileService = cms.Service("TFileService", fileName = cms.string("debugHits_test.root"))
process.source = cms.Source("PoolSource",
                            
fileNames = cms.untracked.vstring(                                
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/E5A25746-D0CC-9947-82E8-9E6DF65E022A.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/D507C172-1BCF-124D-BA01-16BD2426D8D2.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/C8B1D249-B1D7-2A4F-A9AA-5DD639561D52.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/C722E935-09DC-0F46-BD30-7F9635330074.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/0A093C4B-1A77-1E44-AB10-2F65C04DCD84.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/682388BE-F139-AA41-AF93-249EE116437C.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/96D692DB-EBB1-D44A-B5A8-4203646AFF1E.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/A259CAE9-61EE-FB45-AFE2-77C37FC9F447.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/FD8347F0-D48B-A04C-98A2-38C646AC7C74.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/7BFE6C0C-3070-6C48-AF23-8CB455543B40.root",
"/store/relval/CMSSW_10_2_4/RelValTTbar_13/GEN-SIM-RECO/PU25ns_102X_upgrade2018_realistic_v12-v1/20000/1992FB62-C6C4-DC41-91E4-D9CF3AA8ECA6.root",
)

)

