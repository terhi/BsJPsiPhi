import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

#-- NUMBER OF EVENTS --#
process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(20) )

#-- SOURCE FILES --#
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(

									 #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/BdToKstarJPsi_EtaPtFilter_8TeV-pythia6-evtgen/AODSIM/PU_RD2_START53_V19F-v1/00000/0002ECCA-AD4A-E311-99EE-1CC1DE050110.root'
                     		#'file:/afs/cern.ch/work/t/terhi/private/BsMCtest.root'
									'file:/afs/cern.ch/work/t/terhi/private/BuMCtest.root'
                           #'file:/lustre/cmsdata/pazzini/BsJpsiPhi/BsJpsiPhi/0618D67F-EDD8-E111-85A7-001A645C984A.root'

                            )
)

#-- LOGGER --#
process.MessageLogger = cms.Service(
                                    "MessageLogger",
                                    #destinations = cms.untracked.vstring('BuToJPsiK_log.txt'),#qui si decide il file di output
                                    #destinations = cms.untracked.vstring('BsMCtest_log.txt'),#qui si decide il file di output
                                    #destinations = cms.untracked.vstring('BdToJPsiKstar_log.txt'),#qui si decide il file di output
                                    default = cms.untracked.PSet( reportEvery = cms.untracked.int32(1000) ) #qui riporta i messaggi ogni ()
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#-- GEOMETRY + B-FIELD --#
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

#-- GLOBAL TAG --#
#process.GlobalTag.globaltag = cms.string('START53_V19F::All')
process.GlobalTag.globaltag = cms.string('START53_V7A::All')

#-- PAT LAYER 0+1 --#
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.pfTools import *
usePFIso( process )

#-- PAT OVERLAP mu/ele --#
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

#-- PAT MC MATCHING mu --#
#MUON MC-MATCHING VALUES FROM BsMuMu MUON-ID STUDIES
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi")
process.muonMatch.matched = cms.InputTag("genParticles")
process.muonMatch.maxDeltaR = cms.double(0.12)
process.muonMatch.maxDPtRel = cms.double(0.3)
process.muonMatch.checkCharge = cms.bool(True)
process.muonMatch.resolveAmbiguities = cms.bool(True)
process.muonMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT MC MATCHING ele --#
#ELECTRON MC-MATCHING VALUES FROM: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
# object        electron        photon  muon    tau to jet      jet
# maxDPtRel     0.5             1.0     0.5     3.0             3.0
# maxDeltaR     0.5             0.2     0.5     0.1             0.4
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi")
process.electronMatch.matched = cms.InputTag("genParticles")
process.electronMatch.maxDeltaR = cms.double(0.5)
process.electronMatch.maxDPtRel = cms.double(0.5)
process.electronMatch.checkCharge = cms.bool(True)
process.electronMatch.resolveAmbiguities = cms.bool(True)
process.electronMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT MC MATCHING jet --#
#JET MC-MATCHING VALUES FROM: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatching
# object        electron        photon  muon    tau to jet      jet
# maxDPtRel     0.5             1.0     0.5     3.0             3.0
# maxDeltaR     0.5             0.2     0.5     0.1             0.4
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
process.patJetPartonMatch.src     = cms.InputTag("patJets")
process.patJetPartonMatch.matched = cms.InputTag("genParticles")
process.patJetPartonMatch.maxDeltaR = cms.double(0.4)
process.patJetPartonMatch.maxDPtRel = cms.double(3.0)
process.patJetPartonMatch.checkCharge = cms.bool(False)
process.patJetPartonMatch.resolveAmbiguities = cms.bool(True)
process.patJetPartonMatch.resolveByMatchQuality = cms.bool(True)

#-- PAT TRACKS --#
from PhysicsTools.PatAlgos.tools.trackTools import *

#kaons
process.allKTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('K+')
                                    )
process.kTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allKTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )

#pions
process.allPiTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
                                    src = cms.InputTag("generalTracks"),
                                    particleType = cms.string('pi+')
                                    )
process.piTracks = cms.EDFilter("CandViewRefSelector",
                               src = cms.InputTag("allPiTracks"),
                               #cut = cms.string("pt > 0.3 & abs(eta) < 2.5")
                               cut = cms.string("pt > 0.59 & abs(eta) < 2.5")
                               )

#-- SIM HITS --#
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
                                            #src = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
                                            #setStatus  = cms.int32(8),
                                            ##   particleTypes = cms.vstring(""),
                                            #filter = cms.vstring("pt > 0.0"),  # just for testing
                                            #genParticles = cms.InputTag("genParticles")
                                            #)

#-- PAT ELECTRONS --#
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.mvaID = cms.Sequence(  process.mvaTrigV0 +  process.mvaTrigNoIPV0 + process.mvaNonTrigV0 + process.simpleEleIdSequence )

process.patElectrons.useParticleFlow  = cms.bool(True)
process.patElectrons.pfElectronSource = cms.InputTag("particleFlow")
process.patElectrons.embedTrack       = cms.bool(True)

process.patElectrons.electronIDSources = cms.PSet(
    #MVA
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
    mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
    #non-MVA
    #eidVeto = cms.InputTag("eidVeto"),
    eidTight = cms.InputTag("eidTight"),
    eidLoose = cms.InputTag("eidLoose"),
    eidRobustTight = cms.InputTag("eidRobustTight"),
    eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
    eidRobustLoose = cms.InputTag("eidRobustLoose"),
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
    )

#-- PAT CONVERSION --#
process.patConversions = cms.EDProducer("PATConversionProducer",
                                        electronSource = cms.InputTag("cleanPatElectrons")  ,
                                        #muonSource = cms.InputTag("cleanPatMuons")
                                        # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,
)

#-- PAT JETS --#
import PhysicsTools.PatAlgos.tools.jetTools as jetTools

process.patJets.addTagInfos = cms.bool(True)

jetTools.switchJetCollection(process,
                             cms.InputTag('ak5PFJets'),
                             doJTA              = True,
                             doBTagging         = True,
                             btagInfo           = cms.vstring('impactParameterTagInfos','secondaryVertexTagInfos','softPFMuonsTagInfos','softPFElectronsTagInfos'),
                             btagdiscriminators = cms.vstring('jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags'),
                             jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute']),
                             #### data:                   jetCorrLabel     = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']),
                             doType1MET       = False,
                             genJetCollection = cms.InputTag("ak5GenJets"),
                             doJetID          = True,
                             jetIdLabel       = "ak5",
                             outputModules    = [],
)

#-- ANALYZER TAGS AND PARAMETERS --#
process.bsVertexAnalysis = cms.EDAnalyzer("BsToJpsiPhiAnalysis",
                                          isMCstudy                     = cms.bool( True ),
                                          genParticlesLabel             = cms.InputTag("genParticles"),
                                          TrackLabel_K                  = cms.InputTag("kTracks"),
                                          TrackLabel_pi                 = cms.InputTag("piTracks"),
                                          TriggerTag                    = cms.InputTag("TriggerResults::HLT"),
                                          MuonTag                       = cms.InputTag("patMuons"),
                                          JetCollection                 = cms.InputTag("patJets"),
                                          #ElectronTag                  = cms.InputTag("patElectrons"),
                                          StoreDeDxInfo                 = cms.bool( False ),
                                          JpsiMassWindowBeforeFit       = cms.double(0.31), #leave this selection looser than the trigger one for the efficiency calculation
                                          JpsiMassWindowAfterFit        = cms.double(0.150),
                                          JpsiPtCut                     = cms.double(6),
                                          KaonTrackPtCut                = cms.double(0.6),
                                          BdKaonTrackPtCut              = cms.double(0.6),
                                          PhiMassWindowBeforeFit        = cms.double(0.03),
                                          PhiMassWindowAfterFit         = cms.double(0.02),
                                          BsLowerMassCutBeforeFit       = cms.double(4.5),
                                          BsUpperMassCutBeforeFit       = cms.double(6),
                                          BsLowerMassCutAfterFit        = cms.double(5),
                                          BsUpperMassCutAfterFit        = cms.double(6),
                                          KstarMassWindowBeforeFit      = cms.double(0.2),
                                          KstarMassWindowAfterFit       = cms.double(0.15),
                                          BdLowerMassCutBeforeFit       = cms.double(4.5),
                                          BdUpperMassCutBeforeFit       = cms.double(6),
                                          BdLowerMassCutAfterFit        = cms.double(4.9),
                                          BdUpperMassCutAfterFit        = cms.double(5.7),
                                          verbose                       = cms.bool( False ),
                                          TestVerbose                   = cms.bool( False ),
                                          BsPDGMass = cms.double(5.3699),
                                          BdPDGMass = cms.double(5.2794),
                                          BpPDGMass = cms.double(5.2790),
                                          outputFile                   = cms.untracked.string("BuToJPsiK.root"),
                                          #outputFile                   = cms.untracked.string("BsToJPsiPhi.root"),
                                          #outputFile                   = cms.untracked.string("BdToJPsiKstarMu.root"),
)

#-- PAT MUONS --#
#process.patMuons.useParticleFlow   = cms.bool(True)
process.patMuons.pfMuonSource      = cms.InputTag("particleFlow")
process.patMuons.embedPFCandidate  = cms.bool(True)

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.patMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone()

from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import addMCinfo, changeRecoMuonInput, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution
useL1MatchingWindowForSinglets(process)
addMCinfo(process)
#changeTriggerProcessName(process, "REDIGI36X")
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
process.muonMatchHLTL3.maxDeltaR        = 0.1
process.muonMatchHLTL3.maxDPtRel        = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR  = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel  = 10.0
process.muonMatchHLTTrackMu.maxDeltaR   = 0.1
process.muonMatchHLTTrackMu.maxDPtRel   = 10.0

process.patMuons = cms.EDFilter("PATMuonSelector",
                                src = cms.InputTag("patMuonsWithTrigger"),
                                # cut = cms.string("p>2 && abs(eta)<2.4"),
                                cut = cms.string("p>0 && abs(eta)<1000"),
)

#-- WRAPPING UP --#
process.patDefaultSequence.replace(process.patMuons,process.patMuonsWithoutTrigger * process.patTriggerMatching * process.patMuons)

process.pat = cms.Path(process.mvaID  + process.patDefaultSequence +  process.patConversions)

process.ntup = cms.Path( process.allPiTracks * process.allKTracks * process.kTracks * process.piTracks * process.bsVertexAnalysis )

process.schedule = cms.Schedule(process.pat, process.ntup )

####################################################################
#######               DUMP Completo                         ########
####################################################################
#temp = process.dumpPython()
#outputfile = file("Complete.py",'w')
#outputfile.write(temp)
#outputfile.close()
