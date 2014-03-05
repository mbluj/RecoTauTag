#import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.coreTools import *

#runOnMC = False
runOnMC = True
if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"


process = cms.Process("tauNTuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['TauValidationNTupleProd'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#-- Calibration tag -----------------------------------------------------------
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runOnMC:
    process.GlobalTag.globaltag = cms.string('START71_V1::All')
    #process.GlobalTag.globaltag = cms.string('START53_V7F::All')
    #process.GlobalTag.globaltag = cms.string('START53_V19D::All')
else:
    process.GlobalTag.globaltag = cms.string('GR_P_V43D::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

#-- PAT standard config -------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("RecoVertex.Configuration.RecoVertex_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Dummy output for PAT. Not used in the analysis ##
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('keep *')
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         '/store/relval/CMSSW_7_1_0_pre3/RelValQQH1352T_Tauola_13/GEN-SIM-RECO/PU50ns_POSTLS171_V2-v1/00000/04638DB7-EBA0-E311-8B24-02163E008CDD.root'
         #'/store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S8_START52_V9-v1/0000/B27ECBD4-06C2-E111-BEA5-001A9281171E.root' 
         #'/store/data/Run2012D/TauPlusX/AOD/22Jan2013-v1/30000/68EB3FB8-6187-E211-8CE1-0025904B5FBA.root'
    ),

    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop recoPFTaus_*_*_*'                      
    )
)


#comment out as already definied somwhere in PAT
#from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
#process.goodOfflinePrimaryVertices = cms.EDFilter(
#    "PrimaryVertexObjectFilter",
#    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
#    src=cms.InputTag('offlinePrimaryVertices')
#    )

process.load("PhysicsTools/PatAlgos/patSequences_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process) #create HPS Taus from the pat default sequence

if not runOnMC:
    _jetCorrections=('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])
else:
    _jetCorrections=('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])


# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process) #create pat trigger objects
process.patTriggerSequence = cms.Sequence(
    process.patTrigger*
    process.patTriggerEvent)
process.patDefaultSequence += process.patTriggerSequence

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.jetTools import *

#do not switch
#switchJetCollection(
#        process,
 #       cms.InputTag('ak5PFJets'),
 #       doJTA = True,
 #       doBTagging = True,
 #       jetCorrLabel = _jetCorrections,
 #       doType1MET = False,
 #       doJetID = True,
 #       jetIdLabel = "ak5"
#)

process.TFileService = cms.Service("TFileService",
                                               fileName = cms.string('TauValidationNTuple.root') #output file
                                                                                  )
###############################################
#############    User input  ##################
###############################################

# Enter a list with all trigger filter names you want to investigate.
# A bool with the same name will be created for each filter which denotes if a filter object is matched to the tag tau
filterName = [
               "hltPFTau35",
               "hltPFTau35Track",
               "hltPFTau35TrackPt20",
               ]

#Enter a list of HPS discriminators you want to store in the output tree for the tag tau
IDName = [  
    "decayModeFindingNewDMs",
    "decayModeFindingOldDMs",
    "decayModeFinding",
    "byLooseIsolation",
    "byVLooseCombinedIsolationDeltaBetaCorr",
    "byLooseCombinedIsolationDeltaBetaCorr",
    "byMediumCombinedIsolationDeltaBetaCorr",
    "byTightCombinedIsolationDeltaBetaCorr",
    "byCombinedIsolationDeltaBetaCorrRaw",
    "byLooseCombinedIsolationDeltaBetaCorr3Hits",
    "byMediumCombinedIsolationDeltaBetaCorr3Hits",
    "byTightCombinedIsolationDeltaBetaCorr3Hits",
    "byCombinedIsolationDeltaBetaCorrRaw3Hits",
    "chargedIsoPtSum",
    "neutralIsoPtSum",
    "puCorrPtSum",
    "byIsolationMVA3oldDMwoLTraw",
    "byVLooseIsolationMVA3oldDMwoLT",
    "byLooseIsolationMVA3oldDMwoLT",
    "byMediumIsolationMVA3oldDMwoLT",
    "byTightIsolationMVA3oldDMwoLT",
    "byVTightIsolationMVA3oldDMwoLT",
    "byVVTightIsolationMVA3oldDMwoLT",
    "byIsolationMVA3oldDMwLTraw",
    "byVLooseIsolationMVA3oldDMwLT",
    "byLooseIsolationMVA3oldDMwLT",
    "byMediumIsolationMVA3oldDMwLT",
    "byTightIsolationMVA3oldDMwLT",
    "byVTightIsolationMVA3oldDMwLT",
    "byVVTightIsolationMVA3oldDMwLT",
    "byIsolationMVA3newDMwoLTraw",
    "byVLooseIsolationMVA3newDMwoLT",
    "byLooseIsolationMVA3newDMwoLT",
    "byMediumIsolationMVA3newDMwoLT",
    "byTightIsolationMVA3newDMwoLT",
    "byVTightIsolationMVA3newDMwoLT",
    "byVVTightIsolationMVA3newDMwoLT",
    "byIsolationMVA3newDMwLTraw",
    "byVLooseIsolationMVA3newDMwLT",
    "byLooseIsolationMVA3newDMwLT",
    "byMediumIsolationMVA3newDMwLT",
    "byTightIsolationMVA3newDMwLT",
    "byVTightIsolationMVA3newDMwLT",
    "byVVTightIsolationMVA3newDMwLT",
    "againstElectronLoose",
    "againstElectronMedium",
    "againstElectronTight",
    "againstElectronMVA5raw",
    "againstElectronMVA5category",
    "againstElectronVLooseMVA5",
    "againstElectronLooseMVA5",
    "againstElectronMediumMVA5",
    "againstElectronTightMVA5",
    "againstElectronVTightMVA5",
    "againstElectronDeadECAL",
    "againstMuonLoose",
    "againstMuonMedium",
    "againstMuonTight",
    "againstMuonLoose2",
    "againstMuonMedium2",
    "againstMuonTight2",
    "againstMuonLoose3",
    "againstMuonTight3",
    "againstMuonMVAraw",
    "againstMuonLooseMVA",
    "againstMuonMediumMVA",
    "againstMuonTightMVA",
    ]

common_ntuple_branches = cms.PSet(
    index = cms.string("index"), # Index of reco object in the event
    nRecoObjects = cms.string("nTotalObjects"), # Number of reco objects in the event

    nvtx = cms.string("Nvtx"),

    PV_z = cms.string("getPV.z"),



    RunNr = cms.string("RunNr"),
    EvtNr = cms.string("EvtNr"),
    LumiSec = cms.string("LumiSec"),

    recoTauCand_Pt = cms.string("recoTauCand.pt"),
    recoTauCand_Eta = cms.string("recoTauCand.eta"),
    recoTauCand_Phi = cms.string("recoTauCand.phi"),
    recoTauCand_DecayMode = cms.string("recoTauCand.decayMode"),
    recoTauCand_Vz = cms.string("recoTauCand.vz"),

    isAltTauObjMatched = cms.string("isAltTauObjMatched"),
    altTauObj_Pt = cms.string("altTauObj.pt"),
    altTauObj_Eta = cms.string("altTauObj.eta"),
    altTauObj_Phi = cms.string("altTauObj.phi"),

    tauSeedPfJet_Pt = cms.string("recoTauCand.pfJetRef().pt"),
    tauSeedPfJet_Eta = cms.string("recoTauCand.pfJetRef().eta"),
    tauSeedPfJet_Phi = cms.string("recoTauCand.pfJetRef().phi"),
    tauSeedPfJet_Vz = cms.string("recoTauCand.pfJetRef().vz"),

    pfJet_Pt = cms.string("PfJet.pt"),
    pfJet_Eta = cms.string("PfJet.eta"),
    pfJet_Phi = cms.string("PfJet.phi"),
    pfJetVz = cms.string("PfJet.vz"),
    isPfJetMatched = cms.string("isPfJetMatched"),

)
if runOnMC:
    #common_ntuple_branches.genWeight = cms.string("? genInfo!=0 ? genInfo.weight : 1 ")

    common_ntuple_branches.isGenParticleMatched = cms.string("isGenParticelMatched")
    # Careful! Only use GenTauMatch (returns the values of the generator particle matched to the tagTau) values if "bool TauTrigMatch::GenTauMatchTest()" returns "true". Otherwise it contains (unrealsitic) default values
    common_ntuple_branches.GenParticle_Pt = cms.string("GenParticle.pt")
    common_ntuple_branches.GenParticle_Eta = cms.string("GenParticle.eta")
    common_ntuple_branches.GenParticle_Phi = cms.string("GenParticle.phi")
    common_ntuple_branches.GenParticel_pdgId = cms.string("GenParticle.pdgId")
    common_ntuple_branches.isTauGenJetMatched = cms.string("isTauGenJetMatched")
    common_ntuple_branches.isTauGenJetMatched = cms.string("isTauGenJetMatched")
    # Careful! Only use GenTauJet (returns the values of the generated tau Jet) values if "bool TauTrigMatch::GenHadTauMatch()" returns "true". Otherwise it contains (unrealsitic) default values
    common_ntuple_branches.GenTauJet_Pt = cms.string("GenTauJet.pt")
    common_ntuple_branches.GenTauJet_Eta = cms.string("GenTauJet.eta")
    common_ntuple_branches.GenTauJet_Phi = cms.string("GenTauJet.phi")
    common_ntuple_branches.GenTauJet_DecayMode = cms.string("genDecayMode") # Decay Modes encoded as: 5 * (Prong - 1) + numberOfPi0 (example: oneProng2Pi0 = 2, threeProng1Pi0 = 11) cf. recoDecayModes. Additional decays: oneProngOther = -1, threeProngOther = -10, rare = -100 and NoGenJet = -999 


process.tauNTuple = cms.EDAnalyzer('TauValidationNTupleProd',
      tauTag         = cms.InputTag("patTaus"),
      altTauTag      = cms.InputTag("patTaus"),
      trigTag        = cms.InputTag("patTriggerEvent"),
      runOnMC        = cms.bool(runOnMC),
      ntuple         = common_ntuple_branches,
      maxDR          = cms.double(0.5), #The DeltaR parameter used for the trigger matching
      filterNames    = cms.vstring(),
)

###############################################

for j in range(len(filterName)):
    setattr(common_ntuple_branches, filterName[j], cms.string( "isTrigObjMatched(%i)"%j) )

for j in range(len(IDName)):
    setattr(common_ntuple_branches, IDName[j], cms.string( "recoTauCandID(\"%s\")"%IDName[j]) )

process.tauNTuple.filterNames = cms.vstring( filterName )


process.p = cms.Path(
        process.goodOfflinePrimaryVertices*
        process.PFTau*
        process.patDefaultSequence*
        process.tauNTuple
        )
