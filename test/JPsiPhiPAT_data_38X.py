import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
			    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
			    skipEvents = cms.untracked.uint32(0),
			    fileNames = cms.untracked.vstring(
 #'file:/tmp/terhi/TestData.root'
 #       '/store/data/Run2010A/MuOnia/RECO/v4/000/144/001/0A000271-8BB1-DF11-83E4-0030487CD840.root'
	)
)

#from myAnalyzers.JPsiKsPAT.RecoInput2_cfi import *


process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR_P_V21::All')

################GLOBAL TAGS ##############

#process.GlobalTag.globaltag = cms.string('FT_53_V10_AN3::All') 
# FT_R_53_V10 is wrong
#Analysis Global tag for dataset muonia24Aug2012-v1 reco'ed in CMSSW_5_3_2_pathc4


process.GlobalTag.globaltag = cms.string('FT_53_V18_AN3::All') 
#Analysis Global tag for dataset /MuOnia/Run2012ABCD-22Jan2013-v1/AOD recoed with   CMSSW_5_3_7_pathc5!! 

#process.GlobalTag.globaltag = cms.string('FT_53_V21_AN6::All') 
#Analysis Global tag for dataset /MuOnia/Run2012ABCD-22Jan2013-v1/AOD recoed with   CMSSW_5_3_7_patch6


#process.GlobalTag.globaltag = cms.string('GR_P_V41_AN3::All')
#Analysis Global tag for dataset Run2012C-PromptReco-v2  

#process.GlobalTag.globaltag = cms.string('FT_53_V6_AN3::All') 
# FT_R_53_V6 may be wrong
#Analysis Global tag for dataset muonia13Jul2012-v1 reco'ed in CMSSW_5_3_2_pathc4

################GLOBAL TAGS #################

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.coreTools import *

process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from PhysicsTools.PatAlgos.tools.trackTools import *
#makeTrackCandidates(process, 
#        label='TrackCandsK',                            # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
#        tracks=cms.InputTag('generalTracks'),           # input track collection
#        particleType='K+',                              # particle type (for assigning a mass)
#        preselection='pt>0.3 & abs(eta)<2.5',   # preselection cut on candidates. Only methods of 'reco::Candidate' are available
#        selection='',                           # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
#        isolation={},                                   # Isolations to use ('source':deltaR; set to {} for None)
#        isoDeposits=[],
#        mcAs= None                                     # Replicate MC match as the one used for Muons
#        );                                              #  you can specify more than one collection for this

#makeTrackCandidates(process,
#        label='TrackCandsPi',                             # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
#        tracks=cms.InputTag('generalTracks'),           # input track collection
#        particleType='pi+',                              # particle type (for assigning a mass)
#        preselection='pt>0.3 & abs(eta)<2.5',   # preselection cut on candidates. Only methods of 'reco::Candidate' are available
#        selection='',                           # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
#        isolation={},                                   # Isolations to use ('source':deltaR; set to {} for None)
#        isoDeposits=[],
#        mcAs= None                                     # Replicate MC match as the one used for Muons
#        );
removeMCMatching(process,['Photons','Muons','Taus','Electrons','Jets','METs','PFAll'], outputModules=[])
#removeMCMatching(process, ['All'],"",False)
#removeMCMatching(process, ['All'])

process.allKTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('K+')
                                    )

process.allPiTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('pi+')
                                    )

process.kTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allKTracks"),
                               cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               #cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )

process.piTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allPiTracks"),
                               cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               #cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )

########## PAT ELECTRONS ##################

process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 +  process.mvaTrigNoIPV0 + process.mvaNonTrigV0 )

#Electron ID
process.patElectrons.electronIDSources = cms.PSet(
    #MVA
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
    mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
   )


process.patConversions = cms.EDProducer("PATConversionProducer",
                                        # input collection
 #                                       electronSource = cms.InputTag("gsfElectrons"),
 
                                      electronSource = cms.InputTag("cleanPatElectrons")  
                                        # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,
                                        )

################## PAT ELECTRONS ####################



###### PAT JETS ########
import PhysicsTools.PatAlgos.tools.jetTools as jetTools
jetTools.switchJetCollection(process, 
                    cms.InputTag('ak5PFJets'),   
                    doJTA            = True,            
                    doBTagging       = True,            
# MC:               jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),  
# data:                   
							jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 
                    doType1MET       = False,            
 #                   genJetCollection = cms.InputTag("ak5GenJets"),
                    doJetID      = False, #True in MC
                    jetIdLabel   = "ak5",
                    outputModules = [],
)
###### PAT JETS ########


process.bsVertexAnalysis = cms.EDAnalyzer("BsToJpsiPhiAnalysis",
                                          isMCstudy = cms.bool(False),
                                          genParticlesLabel  = cms.InputTag(""),  
                                          TrackLabel_K = cms.InputTag("kTracks"),
                                          TrackLabel_pi = cms.InputTag("piTracks"),
                                          TriggerTag = cms.InputTag("TriggerResults::HLT"),
                                          MuonTag = cms.InputTag("patMuons"),
								   					JetCollection = cms.InputTag("patJets"),
                                          StoreDeDxInfo = cms.bool( True ),
                                          JpsiMassWindowBeforeFit = cms.double(0.31), #leave this selection looser than the trigger one for the efficiency calculation
                                          JpsiMassWindowAfterFit = cms.double(0.150),
                                          JpsiPtCut      = cms.double(6),  
                                          KaonTrackPtCut = cms.double(0.6),
                                          BdKaonTrackPtCut = cms.double(0.6),
                                          PhiMassWindowBeforeFit  = cms.double(0.03), 
                                          PhiMassWindowAfterFit  = cms.double(0.02),
                                          BsLowerMassCutBeforeFit = cms.double(4.5),
                                          BsUpperMassCutBeforeFit = cms.double(6),
                                          BsLowerMassCutAfterFit  = cms.double(5),
                                          BsUpperMassCutAfterFit  = cms.double(6),
                                          KstarMassWindowBeforeFit =cms.double(0.2),
                                          KstarMassWindowAfterFit =cms.double(0.15),
                                          BdLowerMassCutBeforeFit = cms.double(4.5),
                                          BdUpperMassCutBeforeFit = cms.double(6),
                                          BdLowerMassCutAfterFit = cms.double(4.9),
                                          BdUpperMassCutAfterFit = cms.double(5.7),
                                          verbose                = cms.bool( False ), 
                                          outputFile = cms.untracked.string("DataJan13B.root"),
														BsPDGMass = cms.double(5.36677),
                                          BdPDGMass = cms.double(5.27958),
                                          BpPDGMass = cms.double(5.27926 )
                                         )

###################################################################
###################################################################
# New (easier) Onia2mumu trigger matching
#
#    # Make PAT Muons

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()


from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

#
#################################################################
#################################################################


### ==== Apply some final selection (none by default) ====
process.patMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("p>2 && abs(eta)<2.4"), 
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4),
                                           maxAbsZ = cms.double(15),
                                           maxd0 = cms.double(2)
                                           )

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('breco_JPsiPhiPAT_data_ntpl.root')
#)

#process.patDefaultSequence.remove(process.patJetCorrFactors)
#process.patDefaultSequence.remove(process.patJetCharge)
#process.patDefaultSequence.remove(process.patJetPartonMatch)
#process.patDefaultSequence.remove(process.patJetGenJetMatch)
#process.patDefaultSequence.remove(process.patJetPartons)
#process.patDefaultSequence.remove(process.patJetPartonAssociation)
#process.patDefaultSequence.remove(process.patJetFlavourAssociation)
#process.patDefaultSequence.remove(process.patJets)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJet)
#process.patDefaultSequence.remove(process.metJESCorAK5CaloJetMuons)
#process.patDefaultSequence.remove(process.patMETs)
#process.patDefaultSequence.remove(process.selectedPatJets)
#process.patDefaultSequence.remove(process.cleanPatJets)
#process.patDefaultSequence.remove(process.countPatJets)


# can I do a replace of patMuons with the sequence that includes the trigger matching?
process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)
#process.vertex = cms.Path(process.inclusiveVertexing * process.inclusiveMergedVertices * process.selectedVertices * process.bcandidates)
#process.pat = cms.Path( process.patDefaultSequence )

process.pat = cms.Path(process.mvaID +  process.patDefaultSequence +  process.patConversions)

#print(process.pat)

#patAODTrackCandsUnfiltered*patAODTrackCands*electronMatch*(patTrackCandsMCMatch*patTrackCands+patElectrons)+
#muonMatch*patMuons+pfPileUp+pfNoPileUp*(pfAllNeutralHadrons+pfAllChargedHadrons+pfAllPhotons)*tauIsoDepositPFCandidates*
#tauIsoDepositPFChargedHadrons*tauIsoDepositPFNeutralHadrons*tauIsoDepositPFGammas*tauMatch*tauGenJets*
#tauGenJetsSelectorAllHadrons*tauGenJetMatch*patTaus+photonMatch*patPhotons+patCandidateSummary*selectedPatElectrons+
#selectedPatTrackCands+selectedPatMuons+selectedPatTaus+selectedPatPhotons+selectedPatCandidateSummary*
#cleanPatMuons*(cleanPatElectrons+cleanPatTrackCands)*cleanPatPhotons*cleanPatTaus*cleanPatCandidateSummary*
#countPatElectrons+countPatMuons+countPatTaus+countPatLeptons+countPatPhotons

process.ntup = cms.Path( process.allPiTracks * process.allKTracks * process.kTracks * process.piTracks * process.bsVertexAnalysis )
#process.filter = cms.Path( process.primaryVertexFilter * process.noScraping)
process.filter = cms.Path( process.noScraping)
process.schedule = cms.Schedule( process.filter, process.pat, process.ntup )
