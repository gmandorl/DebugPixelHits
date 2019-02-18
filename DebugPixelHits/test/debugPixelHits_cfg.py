import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe2")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.suppressError = cms.untracked.vstring("patTriggerFull")
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
        #'file:debug_Zmm_lostHits.root',
        #'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/031/00000/307D7C74-0E56-E711-85FD-02163E01A1BD.root',
        #'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/1077D486-DD56-E711-825A-02163E01A3F5.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/3ABA4153-5A56-E711-9177-02163E01A488.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/76E58FFB-7956-E711-A5D4-02163E01A279.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/7EFED543-1B5B-E711-A355-02163E0118F3.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/C69BEAD5-4856-E711-A2D3-02163E0140F3.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/046/00000/CCF370A9-D359-E711-8A56-02163E013475.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/047/00000/4C30A60D-4F56-E711-A660-02163E013502.root',
#'root://xrootd-cms.infn.it//store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v1/000/297/047/00000/EE9B2B82-3356-E711-95BD-02163E0133E4.root',

#altri esempi

'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/678/00000/6E1B76C0-A466-E711-A521-02163E0128F4.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/679/00000/9AD7D2B4-A766-E711-B44E-02163E011E6F.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/681/00000/B2370662-AE66-E711-9BC6-02163E01A4E3.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/009F39A0-8A68-E711-806D-02163E01A1BE.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/04870657-8968-E711-96B0-02163E019CB3.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1C9DD9FE-9E68-E711-BA29-02163E01A3B2.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1E2469ED-9768-E711-AE47-02163E01A5AC.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/1E44504F-8968-E711-9E26-02163E01A737.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/2EC43DA1-8C68-E711-BDC0-02163E01A6F7.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/30A9DAE4-9268-E711-A52C-02163E0124B2.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/30BF9297-8568-E711-B137-02163E011E01.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/44734DF8-8D68-E711-B338-02163E019E77.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/6A63CEA7-9768-E711-A772-02163E012A9F.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/6CEC663B-9068-E711-B702-02163E019B22.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/70C05ABA-8B68-E711-B4F2-02163E01A5B7.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/7C89813E-9868-E711-AD47-02163E0119C3.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/8476A489-8F68-E711-8C8C-02163E019BBE.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/9423360C-C868-E711-A5FC-02163E019C43.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/94379C77-8368-E711-9E81-02163E01A3D6.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/B0198019-8C68-E711-8DA3-02163E01340A.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/B2E5D0E3-8E68-E711-B6A7-02163E01A1CE.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/B8B20868-8768-E711-A540-02163E011F3F.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/B8FC5BEB-8968-E711-A92F-02163E01A583.root',
'/store/data/Run2017B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/298/853/00000/BA56E307-8068-E711-9C21-02163E019C2D.root'


#https://cmsweb.cern.ch/das/request?view=plain&instance=prod%2Fglobal&input=file+dataset%3D%2FSingleMuon%2FRun2017B-ZMu-PromptReco-v1%2FRAW-RECO
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/1CDCFB44-DD55-E711-8D01-02163E01450A.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/242FEC1C-D855-E711-A559-02163E01373C.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/3410C6CA-DC55-E711-A2F9-02163E012B04.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/3CEFF073-DB55-E711-9DE5-02163E011CB0.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/6686B87F-D355-E711-80A0-02163E013806.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/6EA95BA9-E055-E711-AAC9-02163E012920.root',
	#'/store/express/Run2017B/ExpressPhysics/FEVT/Express-v1/000/297/218/00000/96C34217-DB55-E711-A0F9-02163E013937.root',
    )
)
#JSON = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = JSON).getVLuminosityBlockRange()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Express_v2', '')

#process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter = cms.EDFilter("TriggerResultsFilter",
    daqPartitions = cms.uint32(1),
    hltResults = cms.InputTag("TriggerResults","","HLT"),
    l1tIgnoreMask = cms.bool(False),
    l1tResults = cms.InputTag(""),
    l1techIgnorePrescales = cms.bool(False),
    throw = cms.bool(True),
    triggerConditions = cms.vstring('HLT_IsoMu20_v*','HLT_IsoMu24_v*')
)

process.tagMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 20 && numberOfMatchedStations >= 2"+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)
process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 10 && numberOfMatchedStations >= 1 && innerTrack.isNonnull"), 
)
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('60 < mass < 140'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi 
process.HitCollectorForDebug = TrackingTools.KalmanUpdators.Chi2MeasurementEstimator_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('HitCollectorForDebug'),
    MaxChi2 = cms.double(30.0), ## was 30 ## TO BE TUNED
    nSigma  = cms.double(3.),   ## was 3  ## TO BE TUNED 
)
process.clusterInfo = cms.EDAnalyzer("DebugPixelHits_CluRef_AddHit",
        pairs = cms.InputTag("tpPairs"),
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
        #SiPixelQuality = cms.string(''),
        badComponentsFile = cms.string('/afs/cern.ch/work/g/gpetrucc/Tracking/CMSSW_9_2_3_patch2/src/RecoTracker/DebugTools/test/badComponents.txt'),
        ## https://github.com/cms-sw/cmssw/blob/9b7f92a91b55fe1bf3e38435a6afd5b97dea4c9f/RecoLocalTracker/SubCollectionProducers/src/JetCoreClusterSplitter.cc#L139-L153
        debug = cms.untracked.int32(100)
)

process.tagAndProbe = cms.Path( 
    #process.triggerResultsFilter +
    process.tagMuons +
    process.oneTag     +
    process.probeMuons +
    process.tpPairs    +
    process.onePair    +
    process.MeasurementTrackerEvent +
    process.clusterInfo 
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("debug_Zmm_lostHits.root"),
    outputCommands = cms.untracked.vstring("keep *", "drop *_*_*_TagProbe"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("tagAndProbe")),
)
#uncomment to skim the selected events
#process.end = cms.EndPath(process.out)

# this below probably not needed
process.TFileService = cms.Service("TFileService", fileName = cms.string("bon4.root"))

if False:
    del process.clusterInfo.badComponentsFile
    process.TFileService = cms.Service("TFileService", fileName = cms.string("debugHits-withBrokenCenter.root"))
    process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v2', '')
    process.source.fileNames = [
        #'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/14163DB4-AD7F-E711-A261-0025905B85FE.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/284B6A7D-AF7F-E711-9AD1-0CC47A78A3D8.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/4EA407D8-B57F-E711-853D-0CC47A4C8E66.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/6A7CA73A-B57F-E711-B8D8-0025905B85F6.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D0CE71F2-AC7F-E711-B1FD-0025905A6088.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D45679B8-B77F-E711-8929-0CC47A4D7636.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D8CDD648-B37F-E711-9850-0025905B858A.root',
#'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/F0BF8D8B-AF7F-E711-8C9A-0025905B8606.root'
        
     
     '/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/0223A19A-A885-E711-A687-0CC47A7C3610.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/069BDFC2-A785-E711-840E-0025905A6076.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/5A0E9898-AC85-E711-AF80-0025905B85D0.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/6AD9B15B-B385-E711-9D26-0CC47A74527A.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/801D68CE-A885-E711-BA7E-0CC47A4D7628.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/80A81E5A-B385-E711-B3C5-0CC47A4D762A.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/B8B890F3-AA85-E711-A8E9-0CC47A4D764A.root',
'/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/F2684D46-AC85-E711-91D7-0025905B859E.root'


#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/02AF0C46-BB7F-E711-90AE-0CC47A4D7650.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/26CB79FC-B17F-E711-B0F3-0CC47A4C8F18.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/2AC5CD5B-B57F-E711-9AF5-0025905B858A.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/409026AA-FA7F-E711-8DA2-0025905B85FC.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/62AF1FF4-B27F-E711-B93F-0CC47A4C8E56.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/8E0B8C8E-B57F-E711-A664-0CC47A78A33E.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/90354ABB-B77F-E711-95B4-0025905A60CA.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/98416BFC-E17F-E711-A85F-0CC47A7C345E.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D0AA918B-B87F-E711-BF10-0025905B8612.root

#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/24E29B4F-AE85-E711-BD92-0025905B8572.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/483E74A5-AA85-E711-9352-0CC47A7C3610.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/5016043F-AB85-E711-83C0-0025905A6094.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/5E8D5F99-AC85-E711-B6CD-0CC47A7C35D2.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/923CE28D-AC85-E711-9294-0CC47A4D75F2.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/A0E6ED96-AD85-E711-936E-0025905A60DA.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/B24F4A04-AA85-E711-9EE3-0025905B859E.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/E25D9E4F-AE85-E711-B7F7-0025905B8572.root
#/store/relval/CMSSW_9_3_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PUpmx25ns_92X_upgrade2017_realistic_v10_resub2-v1/00000/F631B9B6-AD85-E711-9447-0CC47A4D75F2.root

#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/14163DB4-AD7F-E711-A261-0025905B85FE.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/284B6A7D-AF7F-E711-9AD1-0CC47A78A3D8.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/4EA407D8-B57F-E711-853D-0CC47A4C8E66.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/6A7CA73A-B57F-E711-B8D8-0025905B85F6.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D0CE71F2-AC7F-E711-B1FD-0025905A6088.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D45679B8-B77F-E711-8929-0CC47A4D7636.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/D8CDD648-B37F-E711-9850-0025905B858A.root
#/store/relval/CMSSW_9_3_0_pre3/RelValZMM_13/GEN-SIM-RECO/PU25ns_92X_upgrade2017_realistic_v10_resub-v1/00000/F0BF8D8B-AF7F-E711-8C9A-0025905B8606.root
     
]
