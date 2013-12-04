// description on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BsJpsiPhi_AWG
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h" 
#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/BsToJpsiPhiAnalysis.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerReport.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector.h"
#include "TLorentzRotation.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HeavyFlavorAnalysis/BsToJpsiPhi/interface/KinematicFitInterface.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/JetReco/interface/GenJet.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
//#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include <iostream>
#include <TMath.h>

using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;

BsToJpsiPhiAnalysis::BsToJpsiPhiAnalysis(const edm::ParameterSet& iConfig):theConfig_(iConfig),
nominalJpsiMass( 3.096916 ),
nominalPhiMass(1.019 ),
nominalMuonMass(0.1056583),
nominalKaonMass(0.493677),
nominalPionMass(0.139570),
nominalKstarMass(0.892),
nominalBplusMass(5.2792) 
{
  isMCstudy_ = iConfig.getParameter<bool>("isMCstudy");
  thegenParticlesLabel_ = iConfig.getParameter<InputTag>("genParticlesLabel");
  trackLabelK_ = iConfig.getParameter<edm::InputTag>("TrackLabel_K");
  trackLabelPi_ = iConfig.getParameter<edm::InputTag>("TrackLabel_pi");
  triggerTag_ = iConfig.getParameter<edm::InputTag>("TriggerTag");
  jetCollection_ = iConfig.getParameter<edm::InputTag>("JetCollection"); 
  muonTag_ = iConfig.getParameter<edm::InputTag>("MuonTag");
 // electronTag_ = iConfig.getParameter<edm::InputTag>("ElectronTag"); 

  StoreDeDxInfo_ = iConfig.getParameter<bool>("StoreDeDxInfo");
  JpsiMassWindowBeforeFit_ = iConfig.getParameter<double>("JpsiMassWindowBeforeFit");

  BsLowerMassCutBeforeFit_  = iConfig.getParameter<double>("BsLowerMassCutBeforeFit");
  BsUpperMassCutBeforeFit_  = iConfig.getParameter<double>("BsUpperMassCutBeforeFit");
  BsLowerMassCutAfterFit_  = iConfig.getParameter<double>("BsLowerMassCutAfterFit");
  BsUpperMassCutAfterFit_  = iConfig.getParameter<double>("BsUpperMassCutAfterFit");

  JpsiMassWindowAfterFit_ = iConfig.getParameter<double>("JpsiMassWindowAfterFit");
  JpsiPtCut_ =  iConfig.getParameter<double>("JpsiPtCut");
  KaonTrackPtCut_ = iConfig.getParameter<double>("KaonTrackPtCut");
  BdKaonTrackPtCut_ = iConfig.getParameter<double>("BdKaonTrackPtCut");
  PhiMassWindowBeforeFit_ = iConfig.getParameter<double>("PhiMassWindowBeforeFit");
  PhiMassWindowAfterFit_ = iConfig.getParameter<double>("PhiMassWindowAfterFit");

  KstarMassWindowBeforeFit_ = iConfig.getParameter<double>("KstarMassWindowBeforeFit");
  KstarMassWindowAfterFit_ = iConfig.getParameter<double>("KstarMassWindowAfterFit");
  BdLowerMassCutBeforeFit_ = iConfig.getParameter<double>("BdLowerMassCutBeforeFit");
  BdUpperMassCutBeforeFit_ = iConfig.getParameter<double>("BdUpperMassCutBeforeFit");

  BdLowerMassCutAfterFit_ = iConfig.getParameter<double>("BdLowerMassCutAfterFit");
  BdUpperMassCutAfterFit_ = iConfig.getParameter<double>("BdUpperMassCutAfterFit");

  BdPDGMass_ = iConfig.getParameter<double>("BdPDGMass");
  BpPDGMass_ = iConfig.getParameter<double>("BpPDGMass");
  BsPDGMass_ = iConfig.getParameter<double>("BsPDGMass");

  verbose_                = iConfig.getParameter<bool>("verbose");
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile");
  event_counter_ = 0;
  elecounter=0;

  edm::LogInfo("RecoVertex/BsToJpsiPhiAnalysis")<< "Initializing Bs to Jpsi Phi analyser  - Output file: " << outputFile_ <<"\n";

}

BsToJpsiPhiAnalysis::~BsToJpsiPhiAnalysis() {}

void BsToJpsiPhiAnalysis::beginJob()
{

  bsRootTree_ = new BsToJpsiPhiRootTree();
  bsRootTree_->createTree(outputFile_);

}

void BsToJpsiPhiAnalysis::endJob() 
{
  bsRootTree_->writeFile();
  delete bsRootTree_;
  cout << "Total number of Events: " << event_counter_ << endl;
  cout << "Total number of tagged electrons: " <<  elecounter<< endl;	
}

void
BsToJpsiPhiAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  event_counter_++;

//  double nominalJpsiMass=3.096916;
//  double nominalPhiMass=1.019;
//  double nominalMuonMass=0.1056583;
//  double nominalKaonMass=0.493677;
//  double nominalPionMass=0.139570;
//  double nominalKstarMass=0.892;

double bestVtxProbBplus = -1;

  bsRootTree_->resetEntries();

  bsRootTree_->BpPVTrkPt_->clear();
  bsRootTree_->BpPVTrkCharge_->clear();
  bsRootTree_->BpPVTrkEta_->clear();
  bsRootTree_->BpPVTrkPhi_->clear();
  bsRootTree_->BpBJetTrkCharge_->clear();
  bsRootTree_->BpBJetTrkPt_->clear();

   bsRootTree_->PVTrkPt_->clear();
	bsRootTree_->PVTrkCharge_->clear();
	bsRootTree_->PVTrkEta_->clear();
	bsRootTree_->PVTrkPhi_->clear();
	bsRootTree_->BJetTrkCharge_->clear();
	bsRootTree_->BJetTrkPt_->clear();

  pat::CompositeCandidate BCand_best;
  TrackRef trk1Ref_best;
  TrackRef trk2Ref_best;
  TrackRef trkMu1Ref_best;
  TrackRef trkMu2Ref_best;
  RefCountedKinematicParticle bs_best;
  pat::Muon mu1_best;
  pat::Muon mu2_best;

  Vertex PVvtxCosTheta;
  Vertex BpPVvtxCosTheta;
  Vertex BdPVvtxCosTheta; 

  int BsPVVtxInd=0;
  vector<TrackRef> BpTrkRefs;
  vector<TrackRef> BsTrkRefs; 

  //double BsPDGMass=5.36677;
  //double BpPDGMass=5.27926;
  //double BdPDGMass=5.27958;

  //Save information about the pile up
  if(isMCstudy_){
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    int numInteraction = 0;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
      //sum the in-time and out-of-time interactions
      if (PVI->getBunchCrossing()==0)
        numInteraction=numInteraction+PVI->getPU_NumInteractions();
    }
    bsRootTree_->PUinteraction_=numInteraction;
  }


   edm::Handle<edm::View<pat::Jet>> myjets;
   iEvent.getByLabel(jetCollection_,myjets);
   const edm::View<pat::Jet> & jets = *myjets;
/*
	for(size_t a =0; a< jets.size(); a++){
		if( jets[a].bDiscriminator("jetProbabilityBJetTags") > 0.5  ){
		cout << a << ": " << jets[a].bDiscriminator("jetProbabilityBJetTags")<< endl;
		cout << "jet pt "<<jets[a].pt() << endl;
		}
   }	
*/
 // look https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP

  // Get primary vertices  
  int VtxIndex=-99;
  double minVtxProb=-999.;
  double MinBVtxHyp1=-999.;
  
  double BSx;
  double BSy;
  double BSz;
  double BSdx;
  double BSdy;
  double BSdz;
  double BSdxdz;
  double BSdydz;
  double BSsigmaZ;
  double BSdsigmaZ;
  double PVx;
  double PVy;
  double PVz;
  double PVerrx;
  double PVerry;
  double PVerrz;
  double BsLxyz=0;

  BSx=BSy=BSz=PVx=PVy=PVz=PVerrx=PVerry=PVerrz=-9999999;
  BSdx=BSdy=BSdz=BSsigmaZ=BSdsigmaZ=-9999999;

  reco::BeamSpot vertexBeamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
  vertexBeamSpot = *recoBeamSpotHandle;
  
  BSx = vertexBeamSpot.x0();
  BSy = vertexBeamSpot.y0();
  BSz = vertexBeamSpot.z0();
  BSdx = vertexBeamSpot.x0Error();
  BSdy = vertexBeamSpot.y0Error();
  BSdz = vertexBeamSpot.z0Error();
  BSdxdz = vertexBeamSpot.dxdz();
  BSdydz = vertexBeamSpot.dydz();
  BSsigmaZ = vertexBeamSpot.sigmaZ();
  BSdsigmaZ = vertexBeamSpot.sigmaZ0Error();
  
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
 
  bsRootTree_->NVertices_=recVtxs->size(); 
  double minPtVertex=0;
  for(size_t i = 0; i < recVtxs->size(); ++ i) {
    const Vertex &vtx = (*recVtxs)[i];
    double RecVtxProb=TMath::Prob(vtx.chi2(),(int)vtx.ndof());
    double summaPtVertex=0;
    for (reco::Vertex::trackRef_iterator trackvertex=vtx.tracks_begin(); trackvertex!=vtx.tracks_end(); trackvertex++){
      const reco::Track & RTv=*(trackvertex->get());
      summaPtVertex=summaPtVertex+abs(RTv.pt());
    }
    //cout<<"---------"<<endl;
    //cout<<summaPtVertex<<endl;
    //cout<<RecVtxProb<<endl;
    //if(RecVtxProb>minVtxProb){
    if(summaPtVertex>minPtVertex){
      VtxIndex=i;
      //minVtxProb=RecVtxProb;
      minPtVertex=summaPtVertex;
    }
  }
  
  const Vertex &RecVtx = (*recVtxs)[VtxIndex];

  if(VtxIndex!=-99)
    {
      bsRootTree_->isPV_ = 1;
      PVx = RecVtx.x(); 
      PVy= RecVtx.y();
      PVz= RecVtx.z();
      PVerrx=RecVtx.xError();
      PVerry=RecVtx.yError();
      PVerrz=RecVtx.zError();
    }
  else {
    bsRootTree_->isBS_ = 1;
    PVx=BSx; 
    PVy=BSy;
    PVz=BSz;
    PVerrx=BSdx;
    PVerry=BSdy;
    PVerrz=BSdz;
  }

  bsRootTree_->getVtx(BSx,BSy,BSz,PVx,PVy,PVz,PVerrx,PVerry,PVerrz);

  bsRootTree_->BSdx_ = BSdx;
  bsRootTree_->BSdy_ = BSdy;
  bsRootTree_->BSdz_ = BSdz;
  bsRootTree_->BSsigmaZ_ = BSsigmaZ;
  bsRootTree_->BSdsigmaZ_ = BSdsigmaZ;
  
  if(verbose_ == true){
    std::cout<<"Beam Spot (x,y,z) = ("<< BSx << ", "<< BSy << ", "<< BSz 
	     << ")  Primary Vtx = (" << PVx <<" ," << PVy << " ," <<PVz<< ")"<< std::endl;  
  }

  // run/event number 
  bsRootTree_->runNumber_ = iEvent.id().run();
  bsRootTree_->eventNumber_ = (unsigned int)iEvent.id().event();
  bsRootTree_->lumiSection_ = iEvent.luminosityBlock();

  //RunSelector
  //if (bsRootTree_->eventNumber_!=2999926893) return;


  // MC info 
  edm::Handle<GenParticleCollection> genParticles;
  if(isMCstudy_){
    iEvent.getByLabel(thegenParticlesLabel_ , genParticles );
    fillMCInfo(genParticles);
  }

  // dE/dx info
  Handle<DeDxDataValueMap> energyLossHandle;
  if(StoreDeDxInfo_)  iEvent.getByLabel("dedxHarmonic2", energyLossHandle);
  
  // HLT code for trigger bits
  edm::Handle<edm::TriggerResults>  hltresults;
  iEvent.getByLabel(triggerTag_, hltresults);
  
  const  edm::TriggerNames & triggerNames_ = iEvent.triggerNames(*hltresults);  
  int ntrigs = hltresults->size();
  for (int itrig = 0; itrig != ntrigs; ++itrig){
    TString trigName=triggerNames_.triggerName(itrig);
    if (trigName=="HLT_Mu3_TkMu0_Jpsi") bsRootTree_->triggerbit_HLTmu3Tk_ = hltresults->accept(itrig);
    if (trigName=="HLT_Mu0_Track0_Jpsi") bsRootTree_->triggerbit_HLTmu5_ = hltresults->accept(itrig); 
    if (trigName=="HLT_DoubleMu0_Quarkonium_v1") bsRootTree_->triggerbit_HLTmu7_ = hltresults->accept(itrig); 
    if (trigName=="HLT_DoubleMu3")     bsRootTree_->triggerbit_HLTdoubleMu3_      = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu0")     bsRootTree_->triggerbit_HLTdoubleMu0_      = hltresults->accept(itrig);
    if (trigName=="HLT_L1DoubleMuOpen")bsRootTree_->triggerbit_HLTL1DoubleMuOpen_ = hltresults->accept(itrig);
    if (trigName=="HLT_Mu0_TkMu0_Jpsi")bsRootTree_->triggerbit_HLTMu0Track0Jpsi_ = hltresults->accept(itrig);
    if (trigName=="HLT_L1DoubleMuOpen_Tight")bsRootTree_->triggerbit_HLTL1DoubleMuOpenTight_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3_Jpsi_v2")bsRootTree_->triggerbit_HLT_DoubleMu3_Jpsi_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3_Quarkonium_v2")bsRootTree_->triggerbit_HLT_DoubleMu3_Quarkonium_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3_Quarkonium_v1")bsRootTree_->triggerbit_HLT_DoubleMu3_Quarkonium_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon6p5_Jpsi_Displaced_v1")bsRootTree_->triggerbit_Jpsi_Displaced_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon7_Jpsi_Displaced_v1")bsRootTree_->triggerbit_7Jpsi_Displaced_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon7_Jpsi_Displaced_v2")bsRootTree_->triggerbit_7Jpsi_Displaced_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon7_Jpsi_Displaced_v3")bsRootTree_->triggerbit_7Jpsi_Displaced_v3_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3p5_Jpsi_Displaced_v2")bsRootTree_->triggerbit_3p5Jpsi_Displaced_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v1")bsRootTree_->triggerbit_4Jpsi_Displaced_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v4")bsRootTree_->triggerbit_4Jpsi_Displaced_v4_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v5")bsRootTree_->triggerbit_4Jpsi_Displaced_v5_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v1")bsRootTree_->triggerbit_5Jpsi_Displaced_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v2")bsRootTree_->triggerbit_5Jpsi_Displaced_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v4")bsRootTree_->triggerbit_5Jpsi_Displaced_v4_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v5")bsRootTree_->triggerbit_5Jpsi_Displaced_v5_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v1")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v1_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v3")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v3_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v5")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v5_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v6")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v6_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v9")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v9_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_v10")        bsRootTree_->triggerbit_Dimuon0_Jpsi_v10_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon6p5_Jpsi_Displaced_v1_5E32v8p1V5")bsRootTree_->triggerbit_Jpsi_Displaced_v1MC_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3_Quarkonium_v2_5E32v6p1V1") bsRootTree_->triggerbit_HLT_DoubleMu3_Quarkonium_v2MC_            = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu3_Jpsi_v2_5E32v6p1V1")     bsRootTree_->triggerbit_HLT_DoubleMu3_Jpsi_v2MC_                = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v1")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v2")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v3")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v5")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v6")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v9")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon10_Jpsi_Barrel_v10")        bsRootTree_->triggerbit_Dimuon10_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon13_Jpsi_Barrel_v1")        bsRootTree_->triggerbit_Dimuon13_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon13_Jpsi_Barrel_v4")        bsRootTree_->triggerbit_Dimuon13_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon13_Jpsi_Barrel_v5")        bsRootTree_->triggerbit_Dimuon13_Barrel_           = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v9") bsRootTree_->triggerbit_4Jpsi_Displaced_v9_ = hltresults->accept(itrig);
    //New for data 2012
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v10") bsRootTree_->triggerbit_4Jpsi_Displaced_v10_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v11") bsRootTree_->triggerbit_4Jpsi_Displaced_v11_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v12") bsRootTree_->triggerbit_4Jpsi_Displaced_v12_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_Muon_v15") bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_Muon_v16" && bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_==0) bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_Muon_v17" && bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_==0) bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_ = hltresults->accept(itrig);
    if (trigName=="HLT_Dimuon0_Jpsi_Muon_v18" && bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_==0) bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_ = hltresults->accept(itrig);

  }
  
  edm::Handle<pat::MuonCollection> allmuons;
  iEvent.getByLabel(muonTag_,allmuons);
 /* 
  InputTag gsfEleLabel(string("gsfElectrons"));
  Handle<GsfElectronCollection> allelectrons;
  iEvent.getByLabel(gsfEleLabel,allelectrons);
*/

  InputTag EleLabel(string("cleanPatElectrons"));
  edm::Handle<edm::View<pat::Electron> > allelectrons;
  iEvent.getByLabel(EleLabel,allelectrons);

 //  cout << "ele size " << allelectrons->size() << endl;

  if(verbose_==true){
    if(allmuons->size()>0) std::cout<<"******found number of muons= "<< allmuons->size() << std::endl;
  }

  // variables to determine minima of fit probability
  double minVtxP = -99.;   //KK hypothesis
  double minJpsiP = -99;   // jpsi alone

  if (allmuons->size() > 25)
   cout << "WARNING! muon list oversize:" << allmuons->size() << endl;      
  // loop on muons (2 muons at least)
  for(size_t i=0; i < allmuons->size(); ++i){
    const pat::Muon & mu1 = (*allmuons)[i];
    if( !mu1.isPFMuon() ) {continue;}
    if(verbose_ == true){   std::cout<<"Got one muon "<<mu1.pt()<<std::endl;}
    for (size_t j=i+1; j < allmuons->size(); ++j){
      const pat::Muon & mu2 = (*allmuons)[j];
      if( !mu2.isPFMuon() ) {continue;}
      if(verbose_ == true){   std::cout<<"Got the second muon "<<mu2.pt()<<std::endl;}
      

      if (!isMCstudy_) { if ((mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty() || mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty()) && (mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty() || mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty())) continue;}

      if(!mu1.isGlobalMuon() && !mu1.isTrackerMuon()) continue;
      if(!mu2.isGlobalMuon() && !mu2.isTrackerMuon()) continue;      
      if(verbose_==true) std::cout<<"******mu1.isGlobalMuon() == "<<mu1.isGlobalMuon()<<" mu1.isTrackerMuon()==" <<mu1.isTrackerMuon()<<" mu2.isGlobalMuon()== " <<mu2.isGlobalMuon()<<" mu2.isTrackerMuon()== "<<mu2.isTrackerMuon()<<std::endl;
      //Traker muons quality check 
      if (mu1.isTrackerMuon())
        if(!muon::isGoodMuon(mu1, muon::TrackerMuonArbitrated)) continue;
      if (mu2.isTrackerMuon())
        if(!muon::isGoodMuon(mu2, muon::TrackerMuonArbitrated)) continue;
      // *****
      // Before this point the code takes into account the T&P selections (just the "continue", the other selections must be applied offline)
      // *****
      bsRootTree_->ihaveajpsi_=1;
  
      if(mu1.charge()==mu2.charge()) continue;
      if(verbose_==true) std::cout<<"******MUONS HAVE OPPOSITE CHARGE: mu1.charge()= " <<mu1.charge()<<" mu2.charge()= "<<mu2.charge()<<std::endl;

      if(bsRootTree_->iPassedCutIdent_   < 1 ) bsRootTree_->iPassedCutIdent_ = 1 ;
      if(bsRootTree_->iPassedCutIdentBd_   < 1 ) bsRootTree_->iPassedCutIdentBd_ = 1 ;
      
      pat::CompositeCandidate Jpsi;
      Jpsi.addDaughter(mu1);
      Jpsi.addDaughter(mu2);
      AddFourMomenta addP4;
      addP4.set(Jpsi);
      
      if(verbose_==true) std::cout<<"******Di-Muon Mass= " <<Jpsi.mass()<<std::endl;

      if ( abs(Jpsi.mass()- nominalJpsiMass ) > JpsiMassWindowBeforeFit_ || Jpsi.pt() < JpsiPtCut_) continue;
      // jpsi mass windows preliminary cut
      if(bsRootTree_->iPassedCutIdent_   < 2 ) bsRootTree_->iPassedCutIdent_ = 2 ;
      if(bsRootTree_->iPassedCutIdentBd_   < 2 ) bsRootTree_->iPassedCutIdentBd_ = 2 ;
      
      edm::ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);	
      TrackRef trkMu1Ref = mu1.get<TrackRef>();
      TrackRef trkMu2Ref = mu2.get<TrackRef>();
      TrackRef muonTrkP = mu1.track();
      TrackRef muonTrkM = mu2.track();
      
      vector<TransientTrack> trk_all;
      TransientTrack mu1TT=(*theB).build(&trkMu1Ref);
      TransientTrack mu2TT=(*theB).build(&trkMu2Ref);      
      trk_all.push_back(mu1TT);
      trk_all.push_back(mu2TT);
      
      KalmanVertexFitter kvf(true);
      TransientVertex tv = kvf.vertex(trk_all);
      
      if (!tv.isValid()) continue; 
      if(verbose_==true) std::cout<<"****** MUONS HAVE VALID VERTEX FIT"<< std::endl;

      // vertex validity
      if(bsRootTree_->iPassedCutIdent_   < 3 ) bsRootTree_->iPassedCutIdent_ = 3 ;
      if(bsRootTree_->iPassedCutIdentBd_   < 3 ) bsRootTree_->iPassedCutIdentBd_ = 3 ;
       
      if (bsRootTree_->matchFilterJpsi1_!=1)
        bsRootTree_->matchFilterJpsi1_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
      if (bsRootTree_->matchFilterJpsi2_!=1)
        bsRootTree_->matchFilterJpsi2_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();

      int isCowboy=0;
      if (mu1.charge()==1) {
         float mupPhi = atan(mu1.py()/mu1.px());
         if ( mu1.px() < 0 && mu1.py() < 0 ) mupPhi -= 3.1415;
         if ( mu1.px() < 0 && mu1.py() > 0 ) mupPhi += 3.1415;
         float mumPhi = atan(mu2.py()/mu2.px());
         if ( mu2.px() < 0 && mu2.py() < 0 ) mumPhi -= 3.1415;
         if ( mu2.px() < 0 && mu2.py() > 0 ) mumPhi += 3.1415;
         if ( (mupPhi - mumPhi)>0 ) isCowboy=1;
      } else {
         float mupPhi = atan(mu2.py()/mu2.px());
         if ( mu2.px() < 0 && mu2.py() < 0 ) mupPhi -= 3.1415;
         if ( mu2.px() < 0 && mu2.py() > 0 ) mupPhi += 3.1415;
         float mumPhi = atan(mu1.py()/mu1.px());
         if ( mu1.px() < 0 && mu1.py() < 0 ) mumPhi -= 3.1415;
         if ( mu1.px() < 0 && mu1.py() > 0 ) mumPhi += 3.1415;
         if ( (mupPhi - mumPhi)>0 ) isCowboy=1;
      } 
           

      Vertex vertex = tv;
      // ***
      // Calculating variables in the closest way to the trigger
      // ***
      double vtxProb_Jpsi = TMath::Prob(vertex.chi2(),(int)vertex.ndof());
      math::XYZVector pperp(mu1.px() + mu2.px(), mu1.py() + mu2.py(), 0.); 
      reco::Vertex::Point vpoint=vertex.position();
      //translate to global point, should be improved
      GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());
      GlobalPoint displacementFromBeamspot( -1*((BSx -  secondaryVertex.x()) +  (secondaryVertex.z() - BSz) * BSdxdz),-1*((BSy - secondaryVertex.y())+  (secondaryVertex.z() - BSz) * BSdydz), 0);
      reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
      double CosAlpha = vperp.Dot(pperp)/(vperp.R()*pperp.R());
      double MuonsDCA=999;
      TrajectoryStateClosestToPoint mu1TS = mu1TT.impactPointTSCP();
      TrajectoryStateClosestToPoint mu2TS = mu2TT.impactPointTSCP();
      if (mu1TS.isValid() && mu2TS.isValid()) {
        ClosestApproachInRPhi cApp;
        cApp.calculate(mu1TS.theState(), mu2TS.theState());
        MuonsDCA=cApp.distance();
      }
      double max_Dr1=fabs( (- (mu1.vx()-BSx) * mu1.py() + (mu1.vy()-BSy) * mu1.px() ) / mu1.pt() );
      double max_Dr2=fabs( (- (mu2.vx()-BSx) * mu2.py() + (mu2.vy()-BSy) * mu2.px() ) / mu2.pt() );


      reco::Vertex::Error verr = vertex.error();
      // translate to global error, should be improved
      GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2) );
      float lxy = displacementFromBeamspot.perp();
      float lxyerr = sqrt(err.rerr(displacementFromBeamspot));

      bsRootTree_->JpsiNumberOfCandidates_++;
      
      // fill the Jpsi with highest vertex probability in the tree
      //      if(vtxProb_Jpsi > minJpsiP){
      minJpsiP = vtxProb_Jpsi;
      
      bsRootTree_->JpsiM_alone_ = Jpsi.mass();	  
      bsRootTree_->JpsiPhi_alone_ = Jpsi.phi();	  
      bsRootTree_->JpsiEta_alone_ = Jpsi.eta();	  
      bsRootTree_->JpsiPt_alone_ = Jpsi.pt();	  
      bsRootTree_->Mu1Pt_beffit_   = mu1.pt();
      bsRootTree_->Mu1Pz_beffit_   = mu1.pz();
      bsRootTree_->Mu1Eta_beffit_  = mu1.eta();
      bsRootTree_->Mu1Phi_beffit_  = mu1.phi();
      bsRootTree_->Mu2Pt_beffit_   = mu2.pt();
      bsRootTree_->Mu2Pz_beffit_   = mu2.pz();
      bsRootTree_->Mu2Eta_beffit_  = mu2.eta();
      bsRootTree_->Mu2Phi_beffit_  = mu2.phi();
      
      int pixhits1 = 0;
      const reco::HitPattern& pp1 = trkMu1Ref.get()->hitPattern();
      for (int iter=0; iter<pp1.numberOfHits(); iter++) {
	uint32_t hit = pp1.getHitPattern(iter);
	if (pp1.validHitFilter(hit) && pp1.pixelBarrelHitFilter(hit)) pixhits1++;
	if (pp1.validHitFilter(hit) && pp1.pixelEndcapHitFilter(hit)) pixhits1++;
      }
      bsRootTree_->JpsiMu1nPixHits_alone_   = pixhits1;
      
      int pixhits2 = 0;
      const reco::HitPattern& pp2 = trkMu2Ref.get()->hitPattern();
      for (int iter=0; iter<pp2.numberOfHits(); iter++) {
	uint32_t hit = pp2.getHitPattern(iter);
	if (pp2.validHitFilter(hit) && pp2.pixelBarrelHitFilter(hit)) pixhits2++;
	if (pp2.validHitFilter(hit) && pp2.pixelEndcapHitFilter(hit)) pixhits2++;
      }
      bsRootTree_->JpsiMu2nPixHits_alone_   = pixhits2;

      // Muons overlapping remover
      if ( muon::overlap(mu1,mu2,1,1,true) ) continue;

      ////           Phi
      Handle<CandidateView> allTracks;
      iEvent.getByLabel(trackLabelK_, allTracks);
      
      for (size_t k=0; k< allTracks->size(); ++k){
	for (size_t l=k+1; l< allTracks->size(); ++l){
	  
	  const Candidate & track1 = (*allTracks)[k];
	  const Candidate & track2 = (*allTracks)[l];
	  TrackRef trk1Ref = track1.get<TrackRef>();
	  TrackRef trk2Ref = track2.get<TrackRef>();
	  
	  if (track1.charge()==track2.charge()) continue;
	  if (track1.pt() < KaonTrackPtCut_) continue;
	  if (track2.pt() < KaonTrackPtCut_) continue;
          if (trk1Ref->numberOfValidHits() < 5 || trk2Ref->numberOfValidHits()<5) continue;
 
	  // passed kaon opposite sign and pt cut
	  if(bsRootTree_->iPassedCutIdent_   < 4 ) bsRootTree_->iPassedCutIdent_ = 4 ;
	  
	  // phi candidate
	  pat::CompositeCandidate PhiCand;
	  PhiCand.addDaughter(track1);
	  PhiCand.addDaughter(track2);
	  AddFourMomenta ad;
	  ad.set(PhiCand);
	  
	  if (abs(PhiCand.mass()- nominalPhiMass) > PhiMassWindowBeforeFit_) continue;

	  // passed phi mass window before fit
	  if(bsRootTree_->iPassedCutIdent_   < 5 ) bsRootTree_->iPassedCutIdent_ = 5 ;
	  bsRootTree_->PhiNumberOfCandidatesBeforeFit_++;	    
	    
	  // Muon-Track overlap check
	  if(((muonTrkP->charge() == track1.charge()) && (muonTrkP->momentum() == track1.momentum())) ||
	     ((muonTrkP->charge() == track2.charge()) && (muonTrkP->momentum() == track2.momentum()))) continue;
	  if(((muonTrkM->charge() == track1.charge()) && (muonTrkM->momentum() == track1.momentum())) ||
	     ((muonTrkM->charge() == track2.charge()) && (muonTrkM->momentum() == track2.momentum()))) continue;
	  
	  // passed muon - track overlap check
	  if(bsRootTree_->iPassedCutIdent_   < 7 ) bsRootTree_->iPassedCutIdent_ = 7 ;
	  
	  // B candidate
	  pat::CompositeCandidate BCand;
	  BCand.addDaughter(mu1);
	  BCand.addDaughter(mu2);
	  BCand.addDaughter(track1);
	  BCand.addDaughter(track2);
	  AddFourMomenta add4mom;
	  add4mom.set(BCand);
	  
	  if (BCand.mass() < BsLowerMassCutBeforeFit_ || BCand.mass() > BsUpperMassCutBeforeFit_) continue;
	  
	  // passed Bs mass cut before fit
	  if(bsRootTree_->iPassedCutIdent_   < 8 ) bsRootTree_->iPassedCutIdent_ = 8 ;
	  
	  // start fit on the B candidates
          // save the best track refs
          trk1Ref_best=trk1Ref;
          trk2Ref_best=trk2Ref;
          trkMu1Ref_best=trkMu1Ref;
          trkMu2Ref_best=trkMu2Ref;
          mu1_best=mu1;
          mu2_best=mu2;
 
	  vector<TransientTrack> t_tracks;
	  t_tracks.push_back((*theB).build(&trkMu1Ref));
	  t_tracks.push_back((*theB).build(&trkMu2Ref));
	  t_tracks.push_back((*theB).build(&trk1Ref));
	  t_tracks.push_back((*theB).build(&trk2Ref));
	  
	  if (!trkMu1Ref.isNonnull() || !trkMu2Ref.isNonnull() || !trk1Ref.isNonnull() || !trk2Ref.isNonnull()) continue;
	  // checked track references
	  if(bsRootTree_->iPassedCutIdent_   < 9 ) bsRootTree_->iPassedCutIdent_ = 9 ;
	  bsRootTree_->BsNumberOfCandidatesBeforeFit_++;	    

          vector<TransientTrack> phi_tracks;
          phi_tracks.push_back((*theB).build(&trk1Ref));
          phi_tracks.push_back((*theB).build(&trk2Ref));
          KalmanVertexFitter kvfphi;
          TransientVertex tvphi = kvfphi.vertex(phi_tracks);
          if (!tvphi.isValid()) continue;
          if(verbose_==true) std::cout<<"****** KAONS HAVE VALID VERTEX FIT"<< std::endl;
          Vertex vertexphi = tvphi;
          double vtxProb_Phi = TMath::Prob(vertexphi.chi2(),(int)vertexphi.ndof());
          if(verbose_==true) std::cout<<"Phi vertex probability:"<<vtxProb_Phi<< std::endl;

          KalmanVertexFitter kvfbs;
          TransientVertex kvfbsvertex = kvfbs.vertex(t_tracks);
          Vertex vertexbskalman = kvfbsvertex;
          if (!kvfbsvertex.isValid()) continue;
          GlobalError gigi=kvfbsvertex.positionError();
	  
	  // track info before the fit
	  bsRootTree_->K1Pt_beffit_   = track1.pt();
	  bsRootTree_->K1Pz_beffit_   = track1.pz();
	  bsRootTree_->K1Eta_beffit_  = track1.eta();
	  bsRootTree_->K1Phi_beffit_  = track1.phi();
	  bsRootTree_->K2Pt_beffit_   = track2.pt();
	  bsRootTree_->K2Pz_beffit_   = track2.pz();
	  bsRootTree_->K2Eta_beffit_  = track2.eta();
	  bsRootTree_->K2Phi_beffit_  = track2.phi();
	  
	  //call fit interface
	  KinematicFitInterface Kfitter;
	  bool fitSuccess = Kfitter.doFit(t_tracks, nominalMuonMass,  nominalKaonMass, nominalKaonMass);
	    
	  if(fitSuccess != 1) continue;
	  // Kinematic fit success
	  if(bsRootTree_->iPassedCutIdent_   < 10 ) bsRootTree_->iPassedCutIdent_ = 10 ;
	  
	  //double vtxprob_Bs = Kfitter.getProb();
          double vtxprob_Bs = TMath::Prob(vertexbskalman.chi2(),(int)vertexbskalman.ndof());
	  RefCountedKinematicParticle bs = Kfitter.getParticle();
	  RefCountedKinematicVertex bVertex = Kfitter.getVertex();
	  AlgebraicVector7 b_par = bs->currentState().kinematicParameters().vector();
	  AlgebraicSymMatrix77 bs_er = bs->currentState().kinematicParametersError().matrix();
          AlgebraicMatrix33 BVError(bVertex->error().matrix_new());
	  double fittedBsMass = b_par[6];
	 
// NEW STUFF      
/*          RefCountedKinematicTree BsTreeNo = Kfitter.getTree();
          //Apply constraint
          float BsMsigma = 0.00024;
          KinematicConstraint * bs_const = new MassKinematicConstraint( BsPDGMass, BsMsigma);
          KinematicParticleFitter constFitter;
          RefCountedKinematicTree BsTreeN1 = constFitter.fit(bs_const,BsTreeNo);
          BsTreeN1->movePointerToTheTop();
          RefCountedKinematicParticle bs = BsTreeN1->currentParticle();
          AlgebraicVector7 b_par = bs->currentState().kinematicParameters().vector();
          AlgebraicSymMatrix77 bs_er = bs->currentState().kinematicParametersError().matrix();
          RefCountedKinematicVertex bVertex = BsTreeN1->currentDecayVertex();
          AlgebraicMatrix33 BVErrorNi(bVertex->error().matrix_new());
*/
          if (!bVertex->vertexIsValid()) continue;

          if(verbose_ == true){   std::cout<<"Good kinematic vertices"<<std::endl;}
          TMatrix cova(2,2);
          cova.IsSymmetric();
           //cova(0,0)=bs_er(3,3);
           //cova(1,1)=bs_er(4,4);
          //cova(2,2)=bs_er(5,5);
           //cova(1,0)=bs_er(3,4);
          //cova(2,0)=bs_er(3,5);
          //cova(2,1)=bs_er(4,5);
           //cova(0,1)=bs_er(3,4);
          //cova(0,2)=bs_er(3,5);
          //cova(1,2)=bs_er(4,5);
          cova(0,0)=gigi.cxx();
          cova(1,1)=gigi.cyy();
          cova(0,1)=gigi.cyx();
          cova(1,0)=gigi.cyx();
//END NEW STUFF

 
	  // check if it is a valid candidate to be counted (verify passed AfterFit cuts)
	  if (abs(Jpsi.mass() - nominalJpsiMass) < JpsiMassWindowAfterFit_ && Jpsi.pt() > JpsiPtCut_ &&
	      abs(PhiCand.mass() - nominalPhiMass) < PhiMassWindowAfterFit_ &&
	      fittedBsMass > BsLowerMassCutAfterFit_ && fittedBsMass < BsUpperMassCutAfterFit_  ) bsRootTree_->BsNumberOfCandidatesAfterFit_++;
	  
	  // store values in root tree if vtx probability is higher than already stored one
	  if (vtxprob_Bs > minVtxP){

	// if several Bs candidates always clear the track vectors before saving new track quantities

	  TrackRef mu1trkref = mu1.get<TrackRef>();
	  TrackRef mu2trkref = mu2.get<TrackRef>();
	  TrackRef K1trkRef = track1.get<TrackRef>();
	  TrackRef K2trkRef = track2.get<TrackRef>();

		BsTrkRefs.clear();
		BsTrkRefs.push_back(mu1trkref);
		BsTrkRefs.push_back(mu2trkref);			
		BsTrkRefs.push_back(K1trkRef);
		BsTrkRefs.push_back(K2trkRef);
	  
	    
	    if (abs(Jpsi.mass() - nominalJpsiMass) > JpsiMassWindowAfterFit_ || Jpsi.pt() < JpsiPtCut_) continue;
	    // passed jpsi mass window after fit
	    if(bsRootTree_->iPassedCutIdent_   < 11 ) bsRootTree_->iPassedCutIdent_ = 11 ;
	    
	    if (abs(PhiCand.mass() - nominalPhiMass) > PhiMassWindowAfterFit_) continue;
	    // passed phi mass window after fit
	    if(bsRootTree_->iPassedCutIdent_   < 12 ) bsRootTree_->iPassedCutIdent_ = 12 ;
	    
	    if (fittedBsMass < BsLowerMassCutAfterFit_ || fittedBsMass > BsUpperMassCutAfterFit_) continue;
	    // passed Bs mass window after fit
	    if(bsRootTree_->iPassedCutIdent_   < 13 ) bsRootTree_->iPassedCutIdent_ = 13 ;
	    
	    // interesting only if removing best vertex P choice (test for selection) 
	    bsRootTree_->BsNumberOfCandidatesAfterBestFit_++;
	    
	    minVtxP = vtxprob_Bs;

            //Assign the best Bs candidate to store variable
            BCand_best = BCand; 

	    
            // L1/HLT match check on mu1/mu2
            bsRootTree_->matchL11_=!mu1.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
            bsRootTree_->matchL12_=!mu2.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
	    bsRootTree_->match2mu01_=!mu1.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
	    bsRootTree_->match2mu02_=!mu2.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
	    bsRootTree_->match1mu01_=!mu1.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
	    bsRootTree_->match1mu02_=!mu2.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
            bsRootTree_->matchDoubleMu31J_ = !mu1.triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
            bsRootTree_->matchDoubleMu32J_ = !mu2.triggerObjectMatchesByFilter("hltDoubleMu3JpsiL3Filtered").empty();
            bsRootTree_->matchDoubleMu31Q_ = !mu1.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
            bsRootTree_->matchDoubleMu32Q_ = !mu2.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
            bsRootTree_->matchDoubleMu71_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();
            bsRootTree_->matchDoubleMu72_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();
            bsRootTree_->matchDoubleMu41_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
            bsRootTree_->matchDoubleMu42_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
            bsRootTree_->matchDoubleMu51_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
            bsRootTree_->matchDoubleMu52_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
            bsRootTree_->matchDoubleMu01_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
            bsRootTree_->matchDoubleMu02_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
            bsRootTree_->matchDoubleMu101_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
            bsRootTree_->matchDoubleMu102_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
            bsRootTree_->matchDoubleMu131_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();
            bsRootTree_->matchDoubleMu132_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();
           
             
            //pat::TriggerObjectStandAloneCollection mu0tkMatch = patMuonP->triggerObjectMatchesByFilter("hltMu0TrackJpsiTrackMassFiltered");
            //pat::TriggerObjectStandAloneCollection mu3tkMatch = patMuonP->triggerObjectMatchesByFilter("hltMu3TrackJpsiTrackMassFiltered");
	    //-----------------------------------------------------------	       
            bool matchedMu = false, matchedTrack = false;
            pat::TriggerObjectStandAloneCollection  mu0tkmuMatch = mu1.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
            for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
              if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
              //if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
              if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
            } 
            if (matchedMu) bsRootTree_->match2mu31_=1; 
	    else if (matchedTrack) bsRootTree_->match2mu31_=1; 
	    else bsRootTree_->match2mu31_=0;	

	    matchedMu = false; matchedTrack = false;
            mu0tkmuMatch = mu2.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
            for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
             if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
              if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
            } 
            if (matchedMu) bsRootTree_->match2mu32_=1; 
	    else if (matchedTrack) bsRootTree_->match2mu32_=1; 
	    else bsRootTree_->match2mu32_=0;	
	    //-----------------------------------------------------------	       
            matchedMu = false, matchedTrack = false;
            mu0tkmuMatch = mu1.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFilteredTight");
            for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
              if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
              if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
            } 
            if (matchedMu) bsRootTree_->matchmu0tk01_=1; 
	    else if (matchedTrack) bsRootTree_->matchmu0tk01_=1; 
	    else bsRootTree_->matchmu0tk01_=0;	

            matchedMu = false, matchedTrack = false;
            mu0tkmuMatch = mu2.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFilteredTight");
            for (unsigned kt = 0; kt < mu0tkmuMatch.size(); ++kt) {
              if (mu0tkmuMatch[kt].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
              if (mu0tkmuMatch[kt].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
            } 
            if (matchedMu) bsRootTree_->matchmu0tk02_=1; 
	    else if (matchedTrack) bsRootTree_->matchmu0tk02_=1; 
	    else bsRootTree_->matchmu0tk02_=0;	
	    //-----------------------------------------------------------	       
            //pat::TriggerObjectStandAloneCollection  mu3tkmuMatch = patMuonP->triggerObjectMatchesByFilter("hltMu3TkMuJpsiTkMuMassFiltered");
            // end L1/HLT-reco matching

	    reco::Vertex reFitVertex;

	    //recalculate primary vertex without tracks from B
	    //reco::Vertex tmpFitVertex = reVertex(recVtxs, iEvent,iSetup, trkMu1Ref, trkMu2Ref, trk1Ref, trk2Ref);
	    //reco::Vertex tmpFitVertex = reVertex(recVtxs, iEvent,iSetup, mu1, mu2, trk1Ref, trk2Ref);
	    //if(tmpFitVertex.isValid()) reFitVertex = tmpFitVertex;
//	    else reFitVertex = reco::Vertex(RecVtx);   // use the original vertex if the refit vertex is invalid
//	    reFitVertex = RecVtx;   // use the original vertex if the refit vertex is invalid
	    bsRootTree_->PVx_refit_ = reFitVertex.x();
	    bsRootTree_->PVy_refit_ = reFitVertex.y();
	    bsRootTree_->PVz_refit_ = reFitVertex.z();
	    
	    bsRootTree_->PVerrx_refit_ = reFitVertex.xError();
	    bsRootTree_->PVerry_refit_ = reFitVertex.yError();
	    bsRootTree_->PVerrz_refit_ = reFitVertex.zError();
	    
	    // fill kinematic info to tree
	    bsRootTree_->BsFitChi2_  = bs->chiSquared();
	    bsRootTree_->BsFitNdof_   =(int)bs->degreesOfFreedom();
	    bsRootTree_->BsFitVtxProb_ = vtxprob_Bs;
            bsRootTree_->BsPhiVtxProb_ = vtxProb_Phi;
            // Save Jpsi Vertex probability
            bsRootTree_->JpsiVtxProb_ = vtxProb_Jpsi;
	    bsRootTree_->CosDeltaAlpha_ = CosAlpha;
            bsRootTree_->MuMuDCA_ = MuonsDCA;
            bsRootTree_->MuMuDistance_ = lxy;
            bsRootTree_->MuMuDistanceSigma_ = lxyerr;
            bsRootTree_->MuDr1_ = max_Dr1;
            bsRootTree_->MuDr2_ = max_Dr2;
  
	    GlobalVector Bsvec(b_par[3], b_par[4], b_par[5]); // the fitted momentum vector of the Bs
	    bsRootTree_->BsFitM_ = fittedBsMass;		
	    bsRootTree_->BsFitEta_ = Bsvec.eta();
	    bsRootTree_->BsFitPt_  = Bsvec.perp();
	    bsRootTree_->BsFitPz_  = Bsvec.z();
	    bsRootTree_->BsFitPhi_ = Bsvec.phi();
	    
	    RefCountedKinematicTree reftree = Kfitter.getTree();
	    setFitParKK(reftree);
            RefCountedKinematicTree jpsitree = Kfitter.getJpsiTree();
            
	    // fitted kaons
	    vector< RefCountedKinematicParticle > bs_children = reftree->finalStateParticles();
	    AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
	    AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();

            // fitted muons
	    vector< RefCountedKinematicParticle > jpsi_children = jpsitree->finalStateParticles();
	    AlgebraicVector7 bs_par3 = jpsi_children[0]->currentState().kinematicParameters().vector();
	    AlgebraicVector7 bs_par4 = jpsi_children[1]->currentState().kinematicParameters().vector();
	    double pt1 = sqrt(bs_par1[3]*bs_par1[3]+bs_par1[4]*bs_par1[4]);
	    double pt2 = sqrt(bs_par2[3]*bs_par2[3]+bs_par2[4]*bs_par2[4]);
	    bsRootTree_->K1Pt_fit_ = pt1;
	    bsRootTree_->K2Pt_fit_ = pt2;
	    
	    // fitted phi
	    TLorentzVector pK1;
	    double en1 = sqrt(bs_par1[3]*bs_par1[3]+bs_par1[4]*bs_par1[4]+bs_par1[5]*bs_par1[5]+bs_par1[6]*bs_par1[6]);
	    pK1.SetPxPyPzE(bs_par1[3],bs_par1[4],bs_par1[5],en1);
	    TLorentzVector pK2;
	    double en2 = sqrt(bs_par2[3]*bs_par2[3]+bs_par2[4]*bs_par2[4]+bs_par2[5]*bs_par2[5]+bs_par2[6]*bs_par2[6]);
	    pK2.SetPxPyPzE(bs_par2[3],bs_par2[4],bs_par2[5],en2);
	    TLorentzVector PhiFit = pK1 + pK2;
            bsRootTree_->PhiM_fit_ = PhiFit.M();

            if (mu1.isGlobalMuon()) bsRootTree_->BsMu1QualityG_=selGlobalMuon(mu1,RecVtx.position()); 
            if (mu2.isGlobalMuon()) bsRootTree_->BsMu2QualityG_=selGlobalMuon(mu2,RecVtx.position()); 
            if (mu1.isTrackerMuon()) bsRootTree_->BsMu1QualityT_=selTrackerMuon(mu1,RecVtx.position()); 
            if (mu2.isTrackerMuon()) bsRootTree_->BsMu2QualityT_=selTrackerMuon(mu2,RecVtx.position());
	    
	    bsRootTree_->BsFitVtx_x_ = bVertex->position().x();
	    bsRootTree_->BsFitVtx_y_ = bVertex->position().y();
	    bsRootTree_->BsFitVtx_z_ = bVertex->position().z();
	    
	    bsRootTree_->BsM_nofit_ = BCand.mass();
	    bsRootTree_->BsPt_nofit_ = BCand.pt();
	    bsRootTree_->BsPz_nofit_ = BCand.pz();
	    bsRootTree_->BsPhi_nofit_ = BCand.phi();
	    bsRootTree_->BsEta_nofit_ = BCand.eta();
	    
	    bsRootTree_->JpsiM_nofit_ = Jpsi.mass();
	    bsRootTree_->JpsiPhi_nofit_ = Jpsi.phi();	  
	    bsRootTree_->JpsiEta_nofit_ = Jpsi.eta();	  
	    bsRootTree_->JpsiPt_nofit_ = Jpsi.pt();	  
	    bsRootTree_->JpsiPz_nofit_ = Jpsi.pz();	  
	    
	    bsRootTree_->PhiM_nofit_ = PhiCand.mass();
	    bsRootTree_->PhiPhi_nofit_ = PhiCand.phi();	  
	    bsRootTree_->PhiEta_nofit_ = PhiCand.eta();	  
	    bsRootTree_->PhiPt_nofit_ = PhiCand.pt();	  
	    bsRootTree_->PhiPz_nofit_ = PhiCand.pz();	  
	    
	    bsRootTree_->K1Pt_nofit_   = track1.pt();
	    bsRootTree_->K1Pz_nofit_   = track1.pz();
	    bsRootTree_->K1Eta_nofit_  = track1.eta();
	    bsRootTree_->K1Phi_nofit_  = track1.phi();
	    bsRootTree_->K1Key_nofit_  = trk1Ref.key();
	    bsRootTree_->K2Pt_nofit_   = track2.pt();
	    bsRootTree_->K2Pz_nofit_   = track2.pz();
	    bsRootTree_->K2Eta_nofit_  = track2.eta();
	    bsRootTree_->K2Phi_nofit_  = track2.phi();
	    bsRootTree_->K2Key_nofit_  = trk2Ref.key();
	    
	    bsRootTree_->K1Chi2_ = trk1Ref.get()->normalizedChi2();
	    bsRootTree_->K1nHits_= trk1Ref.get()->numberOfValidHits();
	    bsRootTree_->K2Chi2_ = trk2Ref.get()->normalizedChi2();
	    bsRootTree_->K2nHits_= trk2Ref.get()->numberOfValidHits();
	    bsRootTree_->Mu1Chi2_ = trkMu1Ref.get()->normalizedChi2();
	    bsRootTree_->Mu1nHits_= trkMu1Ref.get()->numberOfValidHits();
	    bsRootTree_->Mu2Chi2_ = trkMu2Ref.get()->normalizedChi2();
	    bsRootTree_->Mu2nHits_ =trkMu2Ref.get()->numberOfValidHits();
 	    
            bsRootTree_->BsCowboy_=isCowboy;
	    bsRootTree_->Mu1d0_ = trkMu1Ref->d0();
	    bsRootTree_->Mu2d0_ = trkMu1Ref->d0();
	    bsRootTree_->Mu1dz_ = trkMu1Ref->dz();
	    bsRootTree_->Mu2dz_ = trkMu1Ref->dz();

            bsRootTree_->JpsiMu1Pt_alone_ = mu1.pt();	  
            bsRootTree_->JpsiMu2Pt_alone_ = mu2.pt();	  
            bsRootTree_->JpsiMu1Phi_alone_ = mu1.phi();	  
            bsRootTree_->JpsiMu2Phi_alone_ = mu2.phi();	  
            bsRootTree_->JpsiMu1Eta_alone_ = mu1.eta();	  
            bsRootTree_->JpsiMu2Eta_alone_ = mu2.eta();	  
            bsRootTree_->JpsiMu1d0_alone_ = trkMu1Ref->d0();
            bsRootTree_->JpsiMu2d0_alone_ = trkMu1Ref->d0();
            bsRootTree_->JpsiMu1dz_alone_ = trkMu1Ref->dz();
            bsRootTree_->JpsiMu2dz_alone_ = trkMu1Ref->dz();
            bsRootTree_->JpsiMu1chi2_alone_ = trkMu1Ref->chi2();
            bsRootTree_->JpsiMu2chi2_alone_ = trkMu1Ref->chi2();
            bsRootTree_->JpsiMu1ndof_alone_ = trkMu1Ref->ndof();
            bsRootTree_->JpsiMu2ndof_alone_ = trkMu1Ref->ndof();
            bsRootTree_->JpsiMu1nHits_alone_ =  trkMu1Ref->numberOfValidHits();
            bsRootTree_->JpsiMu2nHits_alone_ =  trkMu2Ref->numberOfValidHits();

            // muon categories:
            // 1: tracker muons
            // 2: global muons
            // 3: global + tracker muon
            // 4: neither tracker nor global
      
            if (mu1.isTrackerMuon() && !mu1.isGlobalMuon())       bsRootTree_->JpsiMuon1Cat_alone_ = 1;
            else if (!mu1.isTrackerMuon() && mu1.isGlobalMuon())  bsRootTree_->JpsiMuon1Cat_alone_ = 2;
            else if (mu1.isTrackerMuon() && mu1.isGlobalMuon())   bsRootTree_->JpsiMuon1Cat_alone_ = 3;
            else if (!mu1.isTrackerMuon() && !mu1.isGlobalMuon()) bsRootTree_->JpsiMuon1Cat_alone_ = 4;
            if (mu1.isPFMuon())       bsRootTree_->MuonCat1_ = 1;
      
            if (mu2.isTrackerMuon() && !mu2.isGlobalMuon())       bsRootTree_->JpsiMuon2Cat_alone_ = 1;
            else if (!mu2.isTrackerMuon() && mu2.isGlobalMuon())  bsRootTree_->JpsiMuon2Cat_alone_ = 2;
            else if (mu2.isTrackerMuon() && mu2.isGlobalMuon())   bsRootTree_->JpsiMuon2Cat_alone_ = 3;
            else if (!mu2.isTrackerMuon() && !mu2.isGlobalMuon()) bsRootTree_->JpsiMuon2Cat_alone_ = 4;
            if (mu2.isPFMuon())       bsRootTree_->MuonCat2_ = 1;
      
            if(mu1.isGlobalMuon()) 
      	      bsRootTree_->MuonType_ = 1;
            else if (mu1.isTrackerMuon() && !mu1.isGlobalMuon())
	      bsRootTree_->MuonType_ = 2;
            else if(mu1.isStandAloneMuon() && !mu1.isGlobalMuon() && !mu1.isTrackerMuon())
	      bsRootTree_->MuonType_ = 3;
            if (mu1.isGlobalMuon()) {
              if(muon::isGoodMuon(mu1, muon::GlobalMuonPromptTight)) {bsRootTree_->Mu1GlobalMuonPromptTight_=1;} else {bsRootTree_->Mu1GlobalMuonPromptTight_=0;}
            }
            if (mu2.isGlobalMuon()) {
              if(muon::isGoodMuon(mu2, muon::GlobalMuonPromptTight)) {bsRootTree_->Mu2GlobalMuonPromptTight_=1;} else {bsRootTree_->Mu2GlobalMuonPromptTight_=0;}
            }

            if (mu1.isTrackerMuon()) {
              if(muon::isGoodMuon(mu1, muon::TrackerMuonArbitrated)) {bsRootTree_->Mu1TrackerMuonArbitrated_=1;} else {bsRootTree_->Mu1TrackerMuonArbitrated_=0;}
              if(muon::isGoodMuon(mu1, muon::TMLastStationTight)) {bsRootTree_->Mu1TMLastStationTight_=1;} else {bsRootTree_->Mu1TMLastStationTight_=0;}       // penetration depth tight selector
              if(muon::isGoodMuon(mu1, muon::TMOneStationTight)) {bsRootTree_->Mu1TMOneStationTight_=1;} else {bsRootTree_->Mu1TMOneStationTight_=0;}       // require one well matched segment
              if(muon::isGoodMuon(mu1, muon::TMLastStationOptimizedLowPtTight)) {bsRootTree_->Mu1TMLastStationOptimizedLowPtTight_=1;} else {bsRootTree_->Mu1TMLastStationOptimizedLowPtTight_=0;} // combination of TMLastStation and TMOneStation
              if(muon::isGoodMuon(mu1, muon::TMLastStationAngTight)) {bsRootTree_->Mu1TMLastStationAngTight_=1;} else {bsRootTree_->Mu1TMLastStationAngTight_=0;}   // TMLastStationTight with additional angular cuts
              if(muon::isGoodMuon(mu1, muon::TMOneStationAngTight)) {bsRootTree_->Mu1TMOneStationAngTight_=1;} else {bsRootTree_->Mu1TMOneStationAngTight_=0;}    // TMOneStationTight with additional angular cuts
              if(muon::isGoodMuon(mu1, muon::TMLastStationOptimizedBarrelLowPtTight)) {bsRootTree_->Mu1TMLastStationOptimizedBarrelLowPtTight_=1;} else {bsRootTree_->Mu1TMLastStationOptimizedBarrelLowPtTight_=0;}   
            }
            if (mu2.isTrackerMuon()) {
              if(muon::isGoodMuon(mu2, muon::TrackerMuonArbitrated)) {bsRootTree_->Mu2TrackerMuonArbitrated_=1;} else {bsRootTree_->Mu2TrackerMuonArbitrated_=0;}
              if(muon::isGoodMuon(mu2, muon::TMLastStationTight)) {bsRootTree_->Mu2TMLastStationTight_=1;} else  {bsRootTree_->Mu2TMLastStationTight_=0;}       // penetration depth tight selector
              if(muon::isGoodMuon(mu2, muon::TMOneStationTight)) {bsRootTree_->Mu2TMOneStationTight_=1;} else  {bsRootTree_->Mu2TMOneStationTight_=0;}       // require one well matched segment
              if(muon::isGoodMuon(mu2, muon::TMLastStationOptimizedLowPtTight)) {bsRootTree_->Mu2TMLastStationOptimizedLowPtTight_=1;} else {bsRootTree_->Mu2TMLastStationOptimizedLowPtTight_=0;} // combination of TMLastStation and TMOneStation
              if(muon::isGoodMuon(mu2, muon::TMLastStationAngTight)) {bsRootTree_->Mu2TMLastStationAngTight_=1;} else {bsRootTree_->Mu2TMLastStationAngTight_=0;}   // TMLastStationTight with additional angular cuts
              if(muon::isGoodMuon(mu2, muon::TMOneStationAngTight)) {bsRootTree_->Mu2TMOneStationAngTight_=1;} else  {bsRootTree_->Mu2TMOneStationAngTight_=0;}    // TMOneStationTight with additional angular cuts
              if(muon::isGoodMuon(mu2, muon::TMLastStationOptimizedBarrelLowPtTight)) {bsRootTree_->Mu2TMLastStationOptimizedBarrelLowPtTight_=1;} else  {bsRootTree_->Mu2TMLastStationOptimizedBarrelLowPtTight_=0;}
            }

	    
	    // dedx info
	    if(StoreDeDxInfo_){
	      const DeDxDataValueMap &  eloss  = *energyLossHandle;
	      double dedxTrk = eloss[trk1Ref].dEdx();
	      double errdedxTrk = eloss[trk1Ref].dEdxError();
	      int NumdedxTrk = eloss[trk1Ref].numberOfMeasurements();
	      
	      bsRootTree_->getDeDx(dedxTrk,errdedxTrk,NumdedxTrk);
	    }
	    
	    // proper decay time and proper decay length without the refitted vertex 
	    // ctau 3D
            BsLxyz=sqrt(pow(bVertex->position().x()-PVx,2)+pow(bVertex->position().y()-PVy,2)+pow(bVertex->position().z()-PVz,2));
	    bsRootTree_->BsCt3D_ = BsPDGMass_*((bVertex->position().x()-PVx)*Bsvec.x()+
					     +(bVertex->position().y()-PVy)*Bsvec.y()+
					     +(bVertex->position().z()-PVz)*Bsvec.z())/
	      (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
	    // ctau 2d
	    bsRootTree_->BsCt2D_ = BsPDGMass_*((bVertex->position().x()-PVx)*Bsvec.x()+
					     +(bVertex->position().y()-PVy)*Bsvec.y())/
	      (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            // ctau 2d BS
            bsRootTree_->BsCt2DBS_ = BsPDGMass_*((bVertex->position().x()-BSx)*Bsvec.x()+(bVertex->position().y()-BSy)*Bsvec.y())/(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            Int_t PVCosThetaIndex = -1;
            double MinDistance = 10000000;    

            if(verbose_ == true){   std::cout<<"Try to find the min Costheta PV"<<std::endl;}

            for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++){
              const Vertex &vtx = (*recVtxs)[VtxInd];

              if (vtx.tracksSize()<11) continue;

              //cout<<" Prob kinematic= "<< Kfitter.getProb() <<" Prb kalman= "<<vtxprob_Bs<<endl;

              //Double_t PVSVvecDotBsPvec=(bVertex->position().x()-vtx.x())*Bsvec.x()+(bVertex->position().y()-vtx.y())*Bsvec.y()+(bVertex->position().z()-vtx.z())*Bsvec.z();
              Double_t PVSVvecDotBsPvec=(kvfbsvertex.position().x()-vtx.x())*Bsvec.x()+(kvfbsvertex.position().y()-vtx.y())*Bsvec.y()+(kvfbsvertex.position().z()-vtx.z())*Bsvec.z();
              //Double_t PVSVlength = TMath::Sqrt( pow((bVertex->position().x()- vtx.x()), 2.0) + pow((bVertex->position().y()-vtx.y()), 2.0) + pow((bVertex->position().z()- vtx.z()), 2.0) );
              Double_t PVSVlength = TMath::Sqrt( pow((kvfbsvertex.position().x()- vtx.x()), 2.0) + pow((kvfbsvertex.position().y()-vtx.y()), 2.0) + pow((kvfbsvertex.position().z()- vtx.z()), 2.0) );
              Double_t BsPlength = TMath::Sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z());
              Double_t BsCosTheta = PVSVvecDotBsPvec / (BsPlength * PVSVlength);
/*
              TVector3 bsmomentumx(Bsvec.x(),0,Bsvec.z());                        
              TVector3 bsmomentumy(0,Bsvec.y(),Bsvec.z());                        
              TVector3 bsflightdirx(-vtx.x()+bVertex->position().x(),0,-vtx.z()+bVertex->position().z());
              TVector3 bsflightdiry(0,-vtx.y()+bVertex->position().y(),-vtx.z()+bVertex->position().z());
              Double_t BsCosThetax=bsmomentumx.Dot(bsflightdirx)/(bsflightdirx.Mag()*bsmomentumx.Mag()); 
              Double_t BsCosThetay=bsmomentumy.Dot(bsflightdiry)/(bsflightdiry.Mag()*bsmomentumy.Mag()); 
              Double_t distance = 1-sqrt(BsCosThetax*BsCosThetax+BsCosThetay*BsCosThetay); */
              Double_t distance = 1-BsCosTheta;

              if(distance < MinDistance){
                 MinDistance = distance;
                 //cout << "mindist "<< MinDistance << endl;
                 PVCosThetaIndex = VtxInd;
              }
            }
            if(verbose_ == true){   std::cout<<"Found the min Costheta PV "<<Bsvec<<std::endl;}
            if (PVCosThetaIndex==-1) { continue;}
				BsPVVtxInd = PVCosThetaIndex;
            PVvtxCosTheta = reVertex(iSetup, vertexBeamSpot,  (*recVtxs)[PVCosThetaIndex], mu1, mu2, trk1Ref, trk2Ref); 
            //bsRootTree_->BsCt3DPVCosTheta_ = BsPDGMass*( (bVertex->position().x()-PVvtxCosTheta.x())*Bsvec.x() + (bVertex->position().y()-PVvtxCosTheta.y())*Bsvec.y() + (bVertex->position().z()-PVvtxCosTheta.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );
            bsRootTree_->BsCt3DPVCosTheta_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*Bsvec.y() + (kvfbsvertex.position().z()-PVvtxCosTheta.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );
            //bsRootTree_->BsCt2DPVCosTheta_ = BsPDGMass*( (bVertex->position().x()-PVvtxCosTheta.x())*Bsvec.x() + (bVertex->position().y()-PVvtxCosTheta.y())*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() );
            bsRootTree_->BsCt2DPVCosTheta_ = BsPDGMass_*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*BCand.px() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*BCand.py() )/( BCand.px()*BCand.px() + BCand.py()*BCand.py() );
            //bsRootTree_->BsCt2DPVCosTheta_ = BsPDGMass*( (kvfbsvertex.position().x()-PVvtxCosTheta.x())*Bsvec.x() + (kvfbsvertex.position().y()-PVvtxCosTheta.y())*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() );
            VertexDistanceXY d2Costheta;
            //Measurement1D measurement2Costheta = d2Costheta.distance(PVvtxCosTheta,bVertex->vertexState());
            Measurement1D measurement2Costheta = d2Costheta.distance(PVvtxCosTheta,kvfbsvertex);
            double error2DCostheta = measurement2Costheta.error();
            //double scale2Costheta = ((bVertex->position().x() - PVvtxCosTheta.x())*Bsvec.x()+(bVertex->position().y() - PVvtxCosTheta.y())*Bsvec.y())/
            //    (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
            //     sqrt((bVertex->position().x() - PVvtxCosTheta.x())*(bVertex->position().x() - PVvtxCosTheta.x())+(bVertex->position().y() - PVvtxCosTheta.y())*(bVertex->position().y() - PVvtxCosTheta.y())));
            double scale2Costheta = ((kvfbsvertex.position().x() - PVvtxCosTheta.x())*BCand.px()+(kvfbsvertex.position().y() - PVvtxCosTheta.y())*BCand.py())/
                (sqrt(BCand.px()*BCand.px()+BCand.py()*BCand.py())*
                 sqrt((kvfbsvertex.position().x() - PVvtxCosTheta.x())*(kvfbsvertex.position().x() - PVvtxCosTheta.x())+(kvfbsvertex.position().y() - PVvtxCosTheta.y())*(kvfbsvertex.position().y() - PVvtxCosTheta.y())));
            TVector LengthVector(2);
            //LengthVector(0)=bVertex->position().x()-PVvtxCosTheta.x();
            LengthVector(0)=kvfbsvertex.position().x()-PVvtxCosTheta.x();
            //LengthVector(1)=bVertex->position().y()-PVvtxCosTheta.y();
            LengthVector(1)=kvfbsvertex.position().y()-PVvtxCosTheta.y();
            bsRootTree_->BsCtErr2DCostheta_=sqrt(pow(BsPDGMass_*(error2DCostheta*abs(scale2Costheta))/sqrt(BCand.px()*BCand.px()+BCand.py()*BCand.py()),2)+pow(BsPDGMass_*abs(scale2Costheta)/(BCand.px()*BCand.px() + BCand.py()*BCand.py()),2)*cova.Similarity(LengthVector));


////////////// Jet tagging code BEGIN ////////////////
	/*
	vector<int> jetIndexes;
	vector<int> sharedjetTrks;
	double BprobMaxBjet=-1;
	bool ismatchedtrack;
	int sharedTrkCounter,BjetIndex=	-1;

	TrackRef mu1trkref = mu1.get<TrackRef>();
	TrackRef mu2trkref = mu2.get<TrackRef>();
	TrackRef trk1ref = track1.get<TrackRef>();
	TrackRef trk2ref = track2.get<TrackRef>();
	
	const Vertex &vtx = (*recVtxs)[PVCosThetaIndex];

	for(size_t i=0; i<jets.size(); i++){
	  const reco::TrackRefVector & jetTrackVec = jets[i].associatedTracks();
 
	  if(jets[i].isPFJet() == false) continue;

	    ismatchedtrack = false;
	    sharedTrkCounter = 0;	

	  for(size_t k=0; k<jetTrackVec.size(); k++){
		reco::TrackRef jetTrk = jetTrackVec[k];

		if(mu1trkref == jetTrk || mu2trkref == jetTrk){continue;}
		if(trk1ref == jetTrk || trk2ref == jetTrk) {continue;}
 		
		for(reco::Vertex::trackRef_iterator trackvertex=vtx.tracks_begin(); trackvertex!=vtx.tracks_end(); trackvertex++){
			
		reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
		if(jetTrk == PVtrk ) {sharedTrkCounter++; ismatchedtrack = true; }	

	     }// end of PV track loop 			
	   } // end of jet track loop

	 if(ismatchedtrack == true){
	  jetIndexes.push_back(i); 
	  sharedjetTrks.push_back(sharedTrkCounter);	
	 }
	} //end of jet loop

//	cout << "jets " << jets.size() << endl;
//	cout << " PV jet cands "<< jetIndexes.size() << endl;
//	cout << "shared tracks " << sharedTrkCounter << endl;	

	for(size_t i =0; i<jetIndexes.size(); i++ ){			
	    double Bprobability = jets[jetIndexes[i]].bDiscriminator("jetProbabilityBJetTags"); 
	    if(Bprobability <= 0) continue;
//	  cout << "Bjet prob "<< Bprobability << "current Bprobmax " << BprobMaxBjet << endl;

	    if(Bprobability > 0 && Bprobability > BprobMaxBjet){
		BjetIndex = jetIndexes[i];
		BprobMaxBjet = Bprobability;		
	  }	
	}	

//	cout << "Final Bjet prob "<<BprobMaxBjet << endl;
//        cout << "Bjet indx "<<BjetIndex << endl;
        double ptmax = -1.;

	if(BjetIndex == -1){
	 //cout << "Pt looppi" << endl;
	  for(size_t i =0; i<jetIndexes.size(); i++ ){
//	   cout << "jet pt no b jet " << jets[jetIndexes[i]].pt() << endl;
	     if( jets[jetIndexes[i]].pt() > ptmax){
		ptmax = jets[jetIndexes[i]].pt();
		BjetIndex = jetIndexes[i];
	     }
	  }
	}
//	cout << "Final Bjet index " << BjetIndex << endl;
//	loop over PV tracks	

	if(BjetIndex!=-1){

	  if( jets[BjetIndex].bDiscriminator("jetProbabilityBJetTags") > 0){
	  bsRootTree_->JetBTagProb_ =  jets[BjetIndex].bDiscriminator("jetProbabilityBJetTags"); 
//	  cout <<"Saved Bjet prob " << jets[BjetIndex].bDiscriminator("jetProbabilityBJetTags")<< endl;
	  }
		
	  bsRootTree_->BJetPt_ = jets[BjetIndex].pt();
	  bsRootTree_->BJetEta_ = jets[BjetIndex].eta();
	  bsRootTree_->BJetPhi_ = jets[BjetIndex].phi();

//	  cout<< "Jet Eta " << jets[BjetIndex].eta() << endl;
//	  cout<< "Jet Phi " << jets[BjetIndex].phi() << endl;
//	  cout<< "Jet Pt " << jets[BjetIndex].pt() << endl;

	  if(isMCstudy_ && jets[BjetIndex].genParton()){
	    const reco::GenParticle *JetParton = jets[BjetIndex].genParton();
	    Int_t PartonId = JetParton->pdgId();
//	    cout << "B jet origin " << PartonId << endl;
	    bsRootTree_->BJetParton_ = PartonId;	
	  }

	  const reco::TrackRefVector & BjetTrackVec = jets[BjetIndex].associatedTracks();
		
	  for(size_t u = 0; u < BjetTrackVec.size(); u++){
	    reco::TrackRef jetTrk = BjetTrackVec[u];	
	    const reco::Track *BjetTrack = jetTrk.get();
//	    cout << "track pt " << BjetTrack->pt() << endl;
//	    cout << "track charge " << BjetTrack->charge() << endl;
	    bsRootTree_->BJetTrkCharge_->push_back(BjetTrack->charge());
	    bsRootTree_->BJetTrkPt_->push_back(BjetTrack->pt());
	  }


	  int trkcouter = 0; 
	  bool OverlappingTrack;

	  for(reco::Vertex::trackRef_iterator trackvertex=vtx.tracks_begin(); trackvertex!=vtx.tracks_end(); trackvertex++){ //loop over Bs PV tracks

	  reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
	  OverlappingTrack = false;

	   for(size_t l = 0; l < recVtxs->size(); l++){ // loop over other PVs
	  
	    if(l != size_t(PVCosThetaIndex) ){
			 
		const Vertex &vertex = (*recVtxs)[l];	    	
		    for(reco::Vertex::trackRef_iterator trkvertex=vertex.tracks_begin(); trkvertex!=vertex.tracks_end(); trkvertex++){ //loop over trks of other PVs
		       reco::TrackRef trk = trkvertex->castTo<reco::TrackRef>();

		       if(trk == PVtrk ){OverlappingTrack = true; }	 	
		     } // end of loop over trks of other PVs
		    }	
	    } // end of loop over other PVs

		//if trk of Bs PV doesn't belong to tracks of another PV, the track is saved
		  if(OverlappingTrack == false){
		   const reco::Track *RTv = trackvertex->get();
		   bsRootTree_->PVTrkPt_->push_back(RTv->pt());
		   bsRootTree_->PVTrkCharge_->push_back(RTv->charge());
		   bsRootTree_->PVTrkEta_->push_back(RTv->eta());
		   bsRootTree_->PVTrkPhi_->push_back(RTv->phi());
		   trkcouter++;
		  }
	    } //loop over Bs PV tracks
		
	} // end of if(BjetIndex!=-1 )

*/

////////////// Jet tagging code END ////////

	    // ctau 3D MPV 
	    AlgebraicMatrix31 pB;
	    pB(0,0) = bs->currentState().globalMomentum().x();
	    pB(1,0) = bs->currentState().globalMomentum().y();
	    pB(2,0) = bs->currentState().globalMomentum().z();
	    
	    AlgebraicMatrix13 pBT;
	    pBT(0,0) = bs->currentState().globalMomentum().x();
	    pBT(0,1) = bs->currentState().globalMomentum().y();
	    pBT(0,2) = bs->currentState().globalMomentum().z();
	    
	    AlgebraicMatrix31 PV;
	    PV(0,0) = PVx;
	    PV(0,1) = PVy;
	    PV(0,2) = PVz;
	    AlgebraicMatrix31 BV;
	    BV(0,0) = bVertex->position().x();
	    BV(0,1) = bVertex->position().y();
	    BV(0,2) = bVertex->position().z();
	    AlgebraicMatrix31 lxyz = BV-PV;
	    AlgebraicMatrix33 PVError(RecVtx.error());
	    AlgebraicMatrix33 BVError(bVertex->error().matrix_new());
	    AlgebraicMatrix33 lxyzError = PVError + BVError;
	    lxyzError.Invert();
	    
	    AlgebraicMatrix11 a = pBT * lxyzError * pB ;
	    AlgebraicMatrix11 b = pBT * lxyzError * lxyz;
	    double num(b(0,0));
	    double deno(a(0,0));
	    bsRootTree_->BsCtMPV_ = (num*bs->currentState().mass())/(deno);

	    //	    cout << "value no refit (3d,2d,MPV) = (" << bsRootTree_->BsCt3D_ << "," << bsRootTree_->BsCt2D_ << "," << bsRootTree_->BsCtMPV_ << ")" << endl;
	    
	    // error on ctau 3D
	    GlobalPoint SVpos( bVertex->position().x(), bVertex->position().y(), bVertex->position().z());
	    GlobalPoint PVpos( PVx, PVy, PVz);
	    GlobalError SVerr( bVertex->error() );
	    GlobalError PVerr( RecVtx.error() );
	    VertexDistance3D d1;
	    Measurement1D measurement = d1.distance(VertexState(SVpos,SVerr),VertexState(PVpos,PVerr));
	    double error3D = measurement.error();
	    double scale1 = ((bVertex->position().x() - PVx)*Bsvec.x()+
			     (bVertex->position().y() - PVy)*Bsvec.y()+
			     (bVertex->position().z() - PVz)*Bsvec.z())/
	      (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z())*
	       sqrt((bVertex->position().x() - PVx)*(bVertex->position().x() - PVx)+
		    (bVertex->position().y() - PVy)*(bVertex->position().y() - PVy)+
		    (bVertex->position().z() - PVz)*(bVertex->position().z() - PVz)));
	    bsRootTree_->BsCtErr3D_ = BsPDGMass_*(error3D*abs(scale1))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
	    
	    // error on ctau 2D
	    VertexDistanceXY d2;
	    Measurement1D measurement2 = d2.distance(RecVtx,bVertex->vertexState());
	    double error2D = measurement2.error();
	    double scale2 = ((bVertex->position().x() - PVx)*Bsvec.x()+(bVertex->position().y() - PVy)*Bsvec.y())/
	      (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
	       sqrt((bVertex->position().x() - PVx)*(bVertex->position().x() - PVx)+(bVertex->position().y() - PVy)*(bVertex->position().y() - PVy)));
	    bsRootTree_->BsCtErr2D_ = BsPDGMass_*(error2D*abs(scale2))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            // error on ctau 2D BS
            VertexDistanceXY d2BS;
            Measurement1D measurement2BS = d2BS.distance(vertexBeamSpot,bVertex->vertexState());
            double error2DBS = measurement2BS.error();
            double scale2BS = ((bVertex->position().x() - BSx)*Bsvec.x()+(bVertex->position().y() - BSy)*Bsvec.y())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
               sqrt((bVertex->position().x() - BSx)*(bVertex->position().x() - BSx)+(bVertex->position().y() - BSy)*(bVertex->position().y() - BSy)));
	    bsRootTree_->BsCtErr2DBS_=BsPDGMass_*(error2DBS*abs(scale2BS))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());
 
	    // error on ctau 2D - 2 (first approximation)
	    bsRootTree_->BsCtErr2D2_ = sqrt((1/(Bsvec.perp()*Bsvec.perp()))*
					    (bs_er(1,1)*Bsvec.x()*Bsvec.x()+
					     bs_er(2,2)*Bsvec.y()*Bsvec.y()+
					     bs_er(1,2)*Bsvec.x()*Bsvec.y()));
	    
	    // error on ctau 3D MPV
	    bsRootTree_->BsCtErrMPV_ = BsPDGMass_/sqrt(deno);

	    //	    cout << "error no refit (3d,2d,MPV) = (" << bsRootTree_->BsCtErr3D_ << "," << bsRootTree_->BsCtErr2D_ << "," << bsRootTree_->BsCtErrMPV_ << ")" << endl;

	    // proper decay time and proper decay length with the refitted vertex
            // ctau 3D
            bsRootTree_->BsCt3Drefit_ = BsPDGMass_*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
						  +(bVertex->position().y()-reFitVertex.y())*Bsvec.y()+
						  +(bVertex->position().z()-reFitVertex.z())*Bsvec.z())/
              (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
            // ctau 2d
            bsRootTree_->BsCt2Drefit_ = BsPDGMass_*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
						  +(bVertex->position().y()-reFitVertex.y())*Bsvec.y())/
              (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());
	    
            // ctau 3D MPV
            AlgebraicMatrix31 pB2;
            pB2(0,0) = bs->currentState().globalMomentum().x();
            pB2(1,0) = bs->currentState().globalMomentum().y();
            pB2(2,0) = bs->currentState().globalMomentum().z();

            AlgebraicMatrix13 pBT2;
            pBT2(0,0) = bs->currentState().globalMomentum().x();
            pBT2(0,1) = bs->currentState().globalMomentum().y();
            pBT2(0,2) = bs->currentState().globalMomentum().z();

            AlgebraicMatrix31 PV2;
            PV2(0,0) = reFitVertex.x();
            PV2(0,1) = reFitVertex.y();
            PV2(0,2) = reFitVertex.z();
            AlgebraicMatrix31 BV2;
            BV2(0,0) = bVertex->position().x();
            BV2(0,1) = bVertex->position().y();
            BV2(0,2) = bVertex->position().z();
            AlgebraicMatrix31 lxyz2 = BV2-PV2;
            AlgebraicMatrix33 PVError2(reFitVertex.error());
            AlgebraicMatrix33 BVError2(bVertex->error().matrix_new());
            AlgebraicMatrix33 lxyzError2 = PVError2 + BVError2;
            lxyzError2.Invert();

            AlgebraicMatrix11 a2 = pBT2 * lxyzError2 * pB2 ;
            AlgebraicMatrix11 b2 = pBT2 * lxyzError2 * lxyz2;
            double num2(b2(0,0));
            double deno2(a2(0,0));
            bsRootTree_->BsCtMPVrefit_ = (num2*bs->currentState().mass())/(deno2);

	    //	    cout << "value refit (3d,2d,MPV) = (" << bsRootTree_->BsCt3Drefit_ << "," << bsRootTree_->BsCt2Drefit_ << "," << bsRootTree_->BsCtMPVrefit_ << ")" << endl;

            // error on ctau 3D
            GlobalPoint SVpos2( bVertex->position().x(), bVertex->position().y(), bVertex->position().z());
            GlobalPoint PVpos2( reFitVertex.x(), reFitVertex.y(), reFitVertex.z());
            GlobalError SVerr2( bVertex->error() );
            GlobalError PVerr2( reFitVertex.error() );
            VertexDistance3D d12;
            Measurement1D measurement12 = d12.distance(VertexState(SVpos2,SVerr2),VertexState(PVpos2,PVerr2));
            double error3D2 = measurement12.error();
            double scale12 = ((bVertex->position().x() - reFitVertex.x())*Bsvec.x()+
                             (bVertex->position().y() - reFitVertex.y())*Bsvec.y()+
                             (bVertex->position().z() - reFitVertex.z())*Bsvec.z())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z())*
               sqrt((bVertex->position().x() - reFitVertex.x())*(bVertex->position().x() - reFitVertex.x())+
                    (bVertex->position().y() - reFitVertex.y())*(bVertex->position().y() - reFitVertex.y())+
                    (bVertex->position().z() - reFitVertex.z())*(bVertex->position().z() - reFitVertex.z())));
            bsRootTree_->BsCtErr3Drefit_ = BsPDGMass_*(error3D2*abs(scale12))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
	    
            // error on ctau 2D
            VertexDistanceXY d22;
            Measurement1D measurement22 = d22.distance(reFitVertex,bVertex->vertexState());
            double error2D2 = measurement22.error();
            double scale22 = ((bVertex->position().x() - reFitVertex.x())*Bsvec.x()+
			      (bVertex->position().y() - reFitVertex.y())*Bsvec.y())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
               sqrt((bVertex->position().x() - reFitVertex.x())*(bVertex->position().x() - reFitVertex.x())+
		    (bVertex->position().y() - reFitVertex.y())*(bVertex->position().y() - reFitVertex.y())));
            bsRootTree_->BsCtErr2Drefit_ = BsPDGMass_*(error2D2*abs(scale22))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            // error on ctau 3D MPV
            bsRootTree_->BsCtErrMPVrefit_ = BsPDGMass_/sqrt(deno2);

	    //	    cout << "error refit (3d,2d,MPV) = (" << bsRootTree_->BsCtErr3Drefit_ << "," << bsRootTree_->BsCtErr2Drefit_ << "," << bsRootTree_->BsCtErrMPVrefit_ << ")" << endl;


	    VertexDistanceXY vdist;	      
	    if(Bsvec.perp()!=0) {
	      bsRootTree_->BsLxy_    = vdist.distance( reFitVertex, bVertex->vertexState() ).value(); 
	      bsRootTree_->BsLxyErr_ = vdist.distance( reFitVertex, bVertex->vertexState() ).error(); 
	      if (  (bVertex->position().x()- reFitVertex.x())*Bsvec.x()+(bVertex->position().y()-reFitVertex.y())*Bsvec.y() < 0  )
		bsRootTree_->BsLxy_ = -1.0 * bsRootTree_->BsLxy_;   // in case negative sign is necessary 
	      bsRootTree_->BsCt_     = bsRootTree_->BsLxy_     *  fittedBsMass/Bsvec.perp();
	      bsRootTree_->BsCtErr_  = bsRootTree_->BsLxyErr_  *  fittedBsMass/Bsvec.perp();
	    }
	    bsRootTree_->BsErrX_  = bs_er(1,1);
	    bsRootTree_->BsErrY_  = bs_er(2,2);
	    bsRootTree_->BsErrXY_ = bs_er(1,2); 
	    
	    VertexDistance3D vdist3d;
	    bsRootTree_->BsDist3d_    = vdist3d.distance(bVertex->vertexState(),reFitVertex).value();
	    bsRootTree_->BsDist3dErr_ = vdist3d.distance(bVertex->vertexState(),reFitVertex).error();
	    bsRootTree_->BsTime3d_    = bsRootTree_->BsDist3d_    * fittedBsMass/Bsvec.perp() * 100. /3.;
	    bsRootTree_->BsTime3dErr_ = bsRootTree_->BsDist3dErr_ * BCand.mass()/Bsvec.perp() * 100. /3.;
	    
	    bsRootTree_->BsDist2d_     = vdist.distance(bVertex->vertexState(),reFitVertex).value();
	    bsRootTree_->BsDist2dErr_ = vdist.distance(bVertex->vertexState(),reFitVertex).error();
	    bsRootTree_->BsTime2d_     = bsRootTree_->BsDist2d_ * fittedBsMass/Bsvec.perp() *100. /3.;
	    bsRootTree_->BsTime2dErr_  = bsRootTree_->BsDist2dErr_ * fittedBsMass/Bsvec.perp() * 100. /3.;
	    
            if(verbose_ == true){   std::cout<<"Calculate the angles of Bs"<<std::endl;}
	    // transversity basis angles
	    TLorentzVector pbs;
	    pbs.SetPxPyPzE(BCand.px(),BCand.py(),BCand.pz(),BCand.energy());
	   
            TLorentzVector pmuplus;
	    TLorentzVector pmuminus;
            if (jpsi_children[0]->currentState().particleCharge() == 1) {
	      pmuplus.SetXYZM(bs_par3[3],bs_par3[4],bs_par3[5],bs_par3[6]);
	      pmuminus.SetXYZM(bs_par4[3],bs_par4[4],bs_par4[5],bs_par4[6]);
            } else {
	      pmuminus.SetXYZM(bs_par3[3],bs_par3[4],bs_par3[5],bs_par3[6]);
	      pmuplus.SetXYZM(bs_par4[3],bs_par4[4],bs_par4[5],bs_par4[6]);
            } 
	    
	    TLorentzVector pkplus;
	    TLorentzVector pkminus;
            if (bs_children[0]->currentState().particleCharge() == 1) {
	      pkplus.SetXYZM(bs_par1[3],bs_par1[4],bs_par1[5],bs_par1[6]);
	      pkminus.SetXYZM(bs_par2[3],bs_par2[4],bs_par2[5],bs_par2[6]);
            } else {
	      pkminus.SetXYZM(bs_par1[3],bs_par1[4],bs_par1[5],bs_par1[6]);
	      pkplus.SetXYZM(bs_par2[3],bs_par2[4],bs_par2[5],bs_par2[6]);
            } 
	    
	    // boosting in JPsi restframe
	    TLorentzVector pjpsi;                                                                                                           
	    pjpsi = pmuplus + pmuminus;
	    TLorentzVector pphi;
	    pphi = pkplus + pkminus;
	    
	    // the betas for the boost
	    TVector3 p3_JPsi;
	    p3_JPsi = pjpsi.Vect();
	    p3_JPsi *= -1./pjpsi.E();
	    
	    // the boost matrix
	    TLorentzRotation boost_jpsi(p3_JPsi);
	    TLorentzVector p_JPsi_JPsi;
	    p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);
	    
	    // the different momenta in the new frame                                                                                                       
	    TLorentzVector p_JPsi_muplus;
	    TLorentzVector p_JPsi_Kplus;
	    TLorentzVector p_JPsi_phi;                                                                       
	    p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);                                                                      
	    p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);                                                                                              
	    p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);
	    
	    // the 3-momenta
	    TVector3 p3_JPsi_muplus;
	    p3_JPsi_muplus = p_JPsi_muplus.Vect();
	    TVector3 p3_JPsi_Kplus;
	    p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	    TVector3 p3_JPsi_phi;
	    p3_JPsi_phi = p_JPsi_phi.Vect();
	    
	    // coordinate system
	    TVector3 x,y,z;
	    x = p3_JPsi_phi.Unit();
	    y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	    y = y.Unit();
	    z = x.Cross(y);
	    
	    // Transversity Basis
	    angle_costheta = p3_JPsi_muplus.Unit() * z;
	    
	    double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costheta*angle_costheta);
	    double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costheta*angle_costheta);
	    angle_phi = TMath::ACos(cos_phi);
	    if (sin_phi < 0){
	      angle_phi =  -angle_phi;
	    }
	    
	    // boosting in phi restframe                                                                                                          
	    // the betas for the boost
	    TVector3 p3_phi;
	    p3_phi = pphi.Vect();
	    p3_phi *= -1./pphi.E();
	      
	    // the boost matrix
	    TLorentzRotation boost_phi(p3_phi);
	    TLorentzVector p_phi_phi;
	    p_phi_phi = boost_phi.VectorMultiplication(pphi);
	    
	    // the different momenta in the new frame
	    TLorentzVector p_phi_Kplus;
	    TLorentzVector p_phi_JPsi;
	    TLorentzVector p_phi_Bs;
	    p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	    p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	    p_phi_Bs = boost_phi.VectorMultiplication(pbs);
	    
	    // the 3-momenta
	    TVector3 p3_phi_Kplus;
	    p3_phi_Kplus = p_phi_Kplus.Vect();
	    TVector3 p3_phi_JPsi;
	    p3_phi_JPsi = p_phi_JPsi.Vect();
	    TVector3 p3_phi_Bs;
	    p3_phi_Bs = p_phi_Bs.Vect();
	    angle_cospsi = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
	    
	    // set cos of angle between bs momentum and decay length
	    AngleBsDecayLength = ((bVertex->position().x()-PVx) * BCand.px() + (bVertex->position().y()-PVy) * BCand.py() + 
				  (bVertex->position().z()-PVz) * BCand.pz()) / sqrt(((bVertex->position().x()-PVx) * (bVertex->position().x()-PVx) + 
										      (bVertex->position().y()-PVy) * (bVertex->position().y()-PVy) + 
										      (bVertex->position().z()-PVz) * (bVertex->position().z()-PVz)) * 
										     (BCand.px()*BCand.px() + BCand.py()*BCand.py() + 
										      BCand.pz()*BCand.pz()));

	    bsRootTree_->getAngles(angle_costheta,angle_phi,angle_cospsi,AngleBsDecayLength);
	    
	    // number of pixel/tracker/muons hits kaons
	    int pixhits1 = 0;
	    // hit pattern of the track
	    const reco::HitPattern& p = trk1Ref.get()->hitPattern();
	    // loop over the hits of the track
	    for (int iter=0; iter<p.numberOfHits(); iter++) {
		uint32_t hit = p.getHitPattern(iter);
		// if the hit is valid and in pixel barrel & endcap, print out the layer
		if (p.validHitFilter(hit) && p.pixelBarrelHitFilter(hit)) pixhits1++;
		if (p.validHitFilter(hit) && p.pixelEndcapHitFilter(hit)) pixhits1++;
	    }
	    bsRootTree_->K1pixH_ = pixhits1;
	    // count the number of valid tracker *** hits ***
	    bsRootTree_->K1trkH_= p.numberOfValidTrackerHits();
	    // count the number of tracker *** layers *** with measurement
	    bsRootTree_->K1trkLay_ =p.trackerLayersWithMeasurement();
	    bsRootTree_->K1muDTh_  =p.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->K1muCSCh_ =p.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->K1muRPCh_ =p.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
	    
	    int pixhits2=0;
	    const reco::HitPattern& p2 = trk2Ref.get()->hitPattern();
	    for (int iter=0; iter<p2.numberOfHits(); iter++) {
	      uint32_t hit = p2.getHitPattern(iter);
	      if (p2.validHitFilter(hit) && p2.pixelBarrelHitFilter(hit)) pixhits2++;
	      if (p2.validHitFilter(hit) && p2.pixelEndcapHitFilter(hit)) pixhits2++;
	    }
	    bsRootTree_->K2pixH_   = pixhits2;
	    bsRootTree_->K2trkH_   = p2.numberOfValidTrackerHits();
	    bsRootTree_->K2trkLay_ = p2.trackerLayersWithMeasurement();
	    bsRootTree_->K2muDTh_  = p2.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->K2muCSCh_ = p2.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->K2muRPCh_ = p2.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
	      
	    // number of pixel/tracker/muons hits muons
	    int pixhits3 = 0;
	    const reco::HitPattern& p3 = trkMu1Ref.get()->hitPattern();
	    for (int iter=0; iter<p3.numberOfHits(); iter++) {
	      uint32_t hit = p3.getHitPattern(iter);
	      if (p3.validHitFilter(hit) && p3.pixelBarrelHitFilter(hit)) pixhits3++;
		if (p3.validHitFilter(hit) && p3.pixelEndcapHitFilter(hit)) pixhits3++;
	    }
	    bsRootTree_->Mu1pixH_   = pixhits3;
	    bsRootTree_->Mu1trkH_   = p3.numberOfValidTrackerHits();
	    bsRootTree_->Mu1trkLay_ = p3.trackerLayersWithMeasurement();
	    bsRootTree_->Mu1muDTh_  = p3.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->Mu1muCSCh_ = p3.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->Mu1muRPCh_ = p3.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
	    
	    int pixhits4=0;
	    const reco::HitPattern& p4 = trkMu2Ref.get()->hitPattern();
	    for (int iter=0; iter<p4.numberOfHits(); iter++) {
	      uint32_t hit = p4.getHitPattern(iter);
	      if (p4.validHitFilter(hit) && p4.pixelBarrelHitFilter(hit)) pixhits4++;
	      if (p4.validHitFilter(hit) && p4.pixelEndcapHitFilter(hit)) pixhits4++;
	    }
	    bsRootTree_->Mu2pixH_   = pixhits4;
	    bsRootTree_->Mu2trkH_   = p4.numberOfValidTrackerHits();
	    bsRootTree_->Mu2trkLay_ = p4.trackerLayersWithMeasurement();
	    bsRootTree_->Mu2muDTh_  = p4.numberOfValidMuonDTHits();      // not-null, valid, muon DT
	    bsRootTree_->Mu2muCSCh_ = p4.numberOfValidMuonCSCHits();    // not-null, valid, muon CSC
	    bsRootTree_->Mu2muRPCh_ = p4.numberOfValidMuonRPCHits();     // not-null, valid, muon RPC
	    
	    
	    // deltaR matching
	    bool K1Truth = MCmatching( track1, genParticles, bsRootTree_->K1mcId_, bsRootTree_->K1momId_, bsRootTree_->K1gmomId_, 333, 531);
	    bool K2Truth = MCmatching( track2, genParticles, bsRootTree_->K2mcId_, bsRootTree_->K2momId_, bsRootTree_->K2gmomId_, 333, 531);
	    bool Mu1Truth= MCmatching( mu1,    genParticles, bsRootTree_->Mu1mcId_,bsRootTree_->Mu1momId_,bsRootTree_->Mu1gmomId_, 443, 531);
	    bool Mu2Truth= MCmatching( mu2,    genParticles, bsRootTree_->Mu2mcId_,bsRootTree_->Mu2momId_,bsRootTree_->Mu2gmomId_, 443, 531);
	    if (K1Truth==1 && K2Truth==1 && Mu1Truth==1 && Mu2Truth==1)  bsRootTree_->isMatched_ = 1;
	    else bsRootTree_->isMatched_ = 0;

            //cout<<"isMatched= "<<bsRootTree_->isMatched_<<endl;
	  }
	  
	} // trk2 loop
      } // trk1 loop
      
      
      if(verbose_ == true){   std::cout<<"Start the K* reco"<<std::endl;}
      ////           Kstar
      Handle<CandidateView> allTracksPi;
      iEvent.getByLabel(trackLabelPi_, allTracksPi);
      
      for (size_t itrack=0; itrack< allTracksPi->size(); ++itrack){
	for (size_t jtrack=itrack+1; jtrack< allTracksPi->size(); ++jtrack){
	
	  const Candidate & track1 = (*allTracksPi)[itrack];
	  const Candidate & track2 = (*allTracksPi)[jtrack];
          TrackRef trk1Ref = track1.get<TrackRef>();
          TrackRef trk2Ref = track2.get<TrackRef>();

	    if (track1.charge()==track2.charge()) continue;
	    if (track1.pt() < BdKaonTrackPtCut_) continue;
	    if (track2.pt() < BdKaonTrackPtCut_) continue;
            if (trk1Ref->numberOfValidHits() < 5 || trk2Ref->numberOfValidHits()<5) continue;
	        
	    if(bsRootTree_->iPassedCutIdentBd_   < 4 ) bsRootTree_->iPassedCutIdentBd_ = 4 ;
	        
	    // kstar candidate
	    double KaonMassSq = nominalKaonMass * nominalKaonMass;
	    double KaonE1 = sqrt(KaonMassSq+track1.px()*track1.px()+track1.py()*track1.py()+track1.pz()*track1.pz());
	    double KaonE2 = sqrt(KaonMassSq+track2.px()*track2.px()+track2.py()*track2.py()+track2.pz()*track2.pz());
	    int K1flag=0;
	    int K2flag=0;
	    double Kstmass1  = sqrt((KaonE1+track2.energy())*(KaonE1+track2.energy())
				    -(track1.px()+track2.px())*(track1.px()+track2.px())
				    -(track1.py()+track2.py())*(track1.py()+track2.py())
				    -(track1.pz()+track2.pz())*(track1.pz()+track2.pz()));
	    double Kstmass2  = sqrt((KaonE2+track1.energy())*(KaonE2+track1.energy())
				    -(track1.px()+track2.px())*(track1.px()+track2.px())
				    -(track1.py()+track2.py())*(track1.py()+track2.py())
				    -(track1.pz()+track2.pz())*(track1.pz()+track2.pz()));
	        
	    if(abs(Kstmass1 - nominalKstarMass) < abs(Kstmass2 - nominalKstarMass)){
	      if(abs(Kstmass1 - nominalKstarMass) > KstarMassWindowBeforeFit_) continue;
	      K1flag=1;
	    } else{
	      if(abs(Kstmass2 - nominalKstarMass) > KstarMassWindowBeforeFit_) continue;
	      K2flag=1;
	    }
	    if(bsRootTree_->iPassedCutIdentBd_   < 5 ) bsRootTree_->iPassedCutIdentBd_ = 5 ;
	        
	  
	    if (abs(Jpsi.mass() - nominalJpsiMass) > JpsiMassWindowBeforeFit_) continue;
	    if(bsRootTree_->iPassedCutIdentBd_   < 6 ) bsRootTree_->iPassedCutIdentBd_ = 6 ;
	        
	    // check on the overlap
	    if(((muonTrkP->charge() == track1.charge()) && (muonTrkP->momentum() == track1.momentum())) ||
	      ((muonTrkP->charge() == track2.charge()) && (muonTrkP->momentum() == track2.momentum()))) continue;
	    if(((muonTrkM->charge() == track1.charge()) && (muonTrkM->momentum() == track1.momentum())) ||
	       ((muonTrkM->charge() == track2.charge()) && (muonTrkM->momentum() == track2.momentum()))) continue;  
            // Muons overlapping remover
            if ( muon::overlap(mu1,mu2,1,1,true) ) continue;

	    if(bsRootTree_->iPassedCutIdentBd_   < 7 ) bsRootTree_->iPassedCutIdentBd_ = 7 ;

	    // B candidate
	    pat::CompositeCandidate BdCand;
	    BdCand.addDaughter(mu1);
	    BdCand.addDaughter(mu2);
	    BdCand.addDaughter(track1);
	    BdCand.addDaughter(track2);
	    AddFourMomenta add4mom;
	    add4mom.set(BdCand);
	      
	    if (BdCand.mass() < BdLowerMassCutBeforeFit_ || BdCand.mass() > BdUpperMassCutBeforeFit_) continue;
	        
	    if(bsRootTree_->iPassedCutIdentBd_   < 8 ) bsRootTree_->iPassedCutIdentBd_ = 8 ;

	    edm::ESHandle<TransientTrackBuilder> theB;
	    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	
	    TrackRef trkkst1 = track1.get<TrackRef>();
	    TrackRef trkkst2 = track2.get<TrackRef>();
	        
	    vector<TransientTrack> t_tks;
	    t_tks.push_back((*theB).build(&trkMu1Ref));
	    t_tks.push_back((*theB).build(&trkMu2Ref));
	    t_tks.push_back((*theB).build(&trkkst1));
	    t_tks.push_back((*theB).build(&trkkst2));
	        
	    if(!trkMu1Ref.isNonnull() || !trkMu2Ref.isNonnull() || !trkkst1.isNonnull() || !trkkst2.isNonnull() ) continue;
	    if(bsRootTree_->iPassedCutIdentBd_   < 9 ) bsRootTree_->iPassedCutIdentBd_ = 9 ;
	        
	    KinematicFitInterface KfitterHyp1;
            bool fitSuccessHyp1 = KfitterHyp1.doFit(t_tks, nominalMuonMass,  nominalKaonMass, nominalPionMass);
	    KinematicFitInterface KfitterHyp2;
            bool fitSuccessHyp2 = KfitterHyp2.doFit(t_tks, nominalMuonMass,  nominalPionMass, nominalKaonMass);
            if(!fitSuccessHyp1 || !fitSuccessHyp2) continue; 
	      
	    if(bsRootTree_->iPassedCutIdentBd_   < 10 ) bsRootTree_->iPassedCutIdentBd_ = 10 ;
	
	    RefCountedKinematicParticle bmesHyp1 = KfitterHyp1.getParticle();
	    RefCountedKinematicVertex bVertexHyp1 = KfitterHyp1.getVertex();
	    AlgebraicVector7 b_parHyp1 = bmesHyp1->currentState().kinematicParameters().vector();
	    AlgebraicSymMatrix77 bd_erHyp1 = bmesHyp1->currentState().kinematicParametersError().matrix();
	    double vtxProbHyp1 = TMath::Prob(bmesHyp1->chiSquared(),(int)bmesHyp1->degreesOfFreedom());

	    RefCountedKinematicParticle bmesHyp2 =  KfitterHyp2.getParticle();
	    RefCountedKinematicVertex bVertexHyp2 = KfitterHyp2.getVertex();
	    AlgebraicVector7 b_parHyp2 = bmesHyp2->currentState().kinematicParameters().vector();
	    AlgebraicSymMatrix77 bd_erHyp2 = bmesHyp2->currentState().kinematicParametersError().matrix();
	    double vtxProbHyp2 = TMath::Prob(bmesHyp2->chiSquared(),(int)bmesHyp2->degreesOfFreedom());

	
	    // 	    // temporary check
	    // 	    if( fabs(vtxProbHyp1 - vtxProbHyp2) > 0.001 ) {
	    // 	      std::cout<<"vtx probs not equal" << std::endl;
	    // 	      exit(1);
	    // 	    }
	    
	    double fittedBdMassHyp1 =  b_parHyp1[6];
	    double fittedBdMassHyp2 =  b_parHyp2[6];

	    RefCountedKinematicTree reftree1 = KfitterHyp1.getTree();
	    RefCountedKinematicTree reftree2 = KfitterHyp2.getTree() ;
            RefCountedKinematicTree jpsitree = KfitterHyp1.getJpsiTree();

            // fitted kaons
            vector< RefCountedKinematicParticle > bd_children = reftree1->finalStateParticles();
            AlgebraicVector7 bd_par1 = bd_children[0]->currentState().kinematicParameters().vector();
            AlgebraicVector7 bd_par2 = bd_children[1]->currentState().kinematicParameters().vector();

            // fitted muons
            vector< RefCountedKinematicParticle > jpsi_children = jpsitree->finalStateParticles();
            AlgebraicVector7 bd_par3 = jpsi_children[0]->currentState().kinematicParameters().vector();
            AlgebraicVector7 bd_par4 = jpsi_children[1]->currentState().kinematicParameters().vector();
	    
	    if (abs(Jpsi.mass() - nominalJpsiMass) < JpsiMassWindowAfterFit_ && Jpsi.pt() > JpsiPtCut_ &&
		(abs(Kstmass1 - nominalKstarMass)< KstarMassWindowAfterFit_  ||   
		 abs(Kstmass2 - nominalKstarMass)< KstarMassWindowAfterFit_) &&
	       ( ( fittedBdMassHyp1 > BdLowerMassCutAfterFit_ && fittedBdMassHyp1 < BdUpperMassCutAfterFit_ ) ||
		( fittedBdMassHyp2 > BdLowerMassCutAfterFit_ && fittedBdMassHyp2 < BdUpperMassCutAfterFit_ ) ) ) bsRootTree_->BdNumberOfCandidates_++;
	    
	    if(vtxProbHyp1>MinBVtxHyp1){
	      
	      if (abs(Jpsi.mass() - nominalJpsiMass) > JpsiMassWindowAfterFit_ || Jpsi.pt() < JpsiPtCut_) continue;
	      // passed jpsi mass window after fit
	      if(bsRootTree_->iPassedCutIdentBd_   < 11 ) bsRootTree_->iPassedCutIdentBd_ = 11 ;
	            
	      if( abs(Kstmass1 - nominalKstarMass)> KstarMassWindowAfterFit_  &&   
		  abs(Kstmass2 - nominalKstarMass)> KstarMassWindowAfterFit_ ) continue;
	      // if(abs(Kstmass2-0.892)> KstarMassWindowAfterFit_) continue;
	      
	      // passed jpsi kstar window after fit
	      if(bsRootTree_->iPassedCutIdentBd_   < 12 ) bsRootTree_->iPassedCutIdentBd_ = 12 ;
	      if ( ( fittedBdMassHyp1 < BdLowerMassCutAfterFit_ || fittedBdMassHyp1 > BdUpperMassCutAfterFit_ ) &&
		   ( fittedBdMassHyp2 < BdLowerMassCutAfterFit_ || fittedBdMassHyp2 > BdUpperMassCutAfterFit_ )) continue;
	      // passed Bd mass window after fit
	      if(bsRootTree_->iPassedCutIdentBd_   < 13 ) bsRootTree_->iPassedCutIdentBd_ = 13 ;


              if (mu1.isPFMuon())       bsRootTree_->BdMuonCat1_ = 1;
              if (mu2.isPFMuon())       bsRootTree_->BdMuonCat2_ = 1;
 
              bsRootTree_->BdTrack1Charge_=track1.charge();
              if (mu1.isGlobalMuon()) bsRootTree_->BdMu1QualityG_=selGlobalMuon(mu1,RecVtx.position()); 
              if (mu2.isGlobalMuon()) bsRootTree_->BdMu2QualityG_=selGlobalMuon(mu2,RecVtx.position()); 
              if (mu1.isTrackerMuon()) bsRootTree_->BdMu1QualityT_=selTrackerMuon(mu1,RecVtx.position()); 
              if (mu2.isTrackerMuon()) bsRootTree_->BdMu2QualityT_=selTrackerMuon(mu2,RecVtx.position()); 

              if (mu1.isTrackerMuon() && !mu1.isGlobalMuon())       bsRootTree_->BdJpsiMuon1Cat_alone_ = 1;
              else if (!mu1.isTrackerMuon() && mu1.isGlobalMuon())  bsRootTree_->BdJpsiMuon1Cat_alone_ = 2;
              else if (mu1.isTrackerMuon() && mu1.isGlobalMuon())   bsRootTree_->BdJpsiMuon1Cat_alone_ = 3;
              else if (!mu1.isTrackerMuon() && !mu1.isGlobalMuon()) bsRootTree_->BdJpsiMuon1Cat_alone_ = 4;

              if (mu2.isTrackerMuon() && !mu2.isGlobalMuon())       bsRootTree_->BdJpsiMuon2Cat_alone_ = 1;
              else if (!mu2.isTrackerMuon() && mu2.isGlobalMuon())  bsRootTree_->BdJpsiMuon2Cat_alone_ = 2;
              else if (mu2.isTrackerMuon() && mu2.isGlobalMuon())   bsRootTree_->BdJpsiMuon2Cat_alone_ = 3;
              else if (!mu2.isTrackerMuon() && !mu2.isGlobalMuon()) bsRootTree_->BdJpsiMuon2Cat_alone_ = 4;
 

	      MinBVtxHyp1 = vtxProbHyp1;

              // L1/HLT match check on mu1/mu2
              bsRootTree_->BdmatchL11_=!mu1.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
              bsRootTree_->BdmatchL12_=!mu2.triggerObjectMatchesByFilter("hltL1DoubleMuOpenTightL1Filtered").empty();
              bsRootTree_->Bdmatch2mu01_=!mu1.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
              bsRootTree_->Bdmatch2mu02_=!mu2.triggerObjectMatchesByFilter("hltDiMuonL3PreFiltered0").empty();
              bsRootTree_->Bdmatch1mu01_=!mu1.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
              bsRootTree_->Bdmatch1mu02_=!mu2.triggerObjectMatchesByFilter("hltDoubleMu0QuarkoniumL3PreFiltered").empty();
              bsRootTree_->BdmatchDoubleMu31Q_ = !mu1.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
              bsRootTree_->BdmatchDoubleMu32Q_ = !mu2.triggerObjectMatchesByFilter("hltDoubleMu3QuarkoniumL3Filtered").empty();
              bsRootTree_->BdmatchDoubleMu71_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();
              bsRootTree_->BdmatchDoubleMu72_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterJpsi").empty();
              bsRootTree_->BdmatchDoubleMu41_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
              bsRootTree_->BdmatchDoubleMu42_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
              bsRootTree_->BdmatchDoubleMu51_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
              bsRootTree_->BdmatchDoubleMu52_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu5Jpsi").empty();
              bsRootTree_->BdmatchDoubleMu01_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
              bsRootTree_->BdmatchDoubleMu02_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
              bsRootTree_->BdmatchDoubleMu101_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
              bsRootTree_->BdmatchDoubleMu102_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiBarrel").empty() || !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon10JpsiBarrel").empty();
              bsRootTree_->BdmatchDoubleMu131_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();
              bsRootTree_->BdmatchDoubleMu132_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterDimuon13JpsiBarrel").empty();


              bool matchedMu = false, matchedTrack = false;
              pat::TriggerObjectStandAloneCollection  mu0tkmuMatch = mu1.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
              for (unsigned k = 0; k < mu0tkmuMatch.size(); ++k) {
                if (mu0tkmuMatch[k].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
                if (mu0tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
              }
              if (matchedMu) bsRootTree_->Bdmatch2mu31_=1;
              else if (matchedTrack) bsRootTree_->Bdmatch2mu31_=1;
              else bsRootTree_->Bdmatch2mu31_=0;

              matchedMu = false; matchedTrack = false;
              mu0tkmuMatch = mu2.triggerObjectMatchesByFilter("hltMu0TkMuJpsiTkMuMassFiltered");
              for (unsigned k = 0; k < mu0tkmuMatch.size(); ++k) {
                if (mu0tkmuMatch[k].collection() == string("hltL3MuonCandidates::HLT")) matchedMu = true;
                if (mu0tkmuMatch[k].collection() == string("hltMuTkMuJpsiTrackerMuonCands::HLT")) matchedTrack = true;
              }
              if (matchedMu) bsRootTree_->Bdmatch2mu32_=1;
              else if (matchedTrack) bsRootTree_->Bdmatch2mu32_=1;
              else bsRootTree_->Bdmatch2mu32_=0;


              // end L1/HLT-reco matching
	      reco::Vertex reFitVertex;
	      //recalculate primary vertex without tracks from B
	      //reco::Vertex tmpFitVertex = reVertex(recVtxs, iEvent,iSetup, mu1, mu2, trkkst1, trkkst2);
	      //if(tmpFitVertex.isValid()) reFitVertex = tmpFitVertex;
//	      else reFitVertex = reco::Vertex(RecVtx);   // use the original vertex if the refit vertex is invalid
//	      reFitVertex = RecVtx;   // use the original vertex if the refit vertex is invalid

	      bsRootTree_->BdPVx_refit_    = reFitVertex.x();
	      bsRootTree_->BdPVy_refit_    = reFitVertex.y();
	      bsRootTree_->BdPVz_refit_    = reFitVertex.z();

	      bsRootTree_->BdPVerrx_refit_ = reFitVertex.xError();
	      bsRootTree_->BdPVerry_refit_ = reFitVertex.yError();
	      bsRootTree_->BdPVerrz_refit_ = reFitVertex.zError();
	      
              bsRootTree_->BdJpsiM_nofit_ = Jpsi.mass();
              bsRootTree_->BdJpsiPhi_nofit_ = Jpsi.phi();
              bsRootTree_->BdJpsiEta_nofit_ = Jpsi.eta();
              bsRootTree_->BdJpsiPt_nofit_ = Jpsi.pt();
              bsRootTree_->BdJpsiPz_nofit_ = Jpsi.pz();

              bsRootTree_->BdMu1Pt_beffit_   = mu1.pt();
              bsRootTree_->BdMu1Pz_beffit_   = mu1.pz();
              bsRootTree_->BdMu1Eta_beffit_  = mu1.eta();
              bsRootTree_->BdMu1Phi_beffit_  = mu1.phi();
              bsRootTree_->BdMu2Pt_beffit_   = mu2.pt();
              bsRootTree_->BdMu2Pz_beffit_   = mu2.pz();
              bsRootTree_->BdMu2Eta_beffit_  = mu2.eta();
              bsRootTree_->BdMu2Phi_beffit_  = mu2.phi();

	      bsRootTree_->BdFitChi2_Hyp1_  = bmesHyp1->chiSquared();
	      bsRootTree_->BdFitNdof_Hyp1_   =(int)bmesHyp1->degreesOfFreedom();
	      bsRootTree_->BdFitChi2_Hyp2_  = bmesHyp2->chiSquared();
	      bsRootTree_->BdFitNdof_Hyp2_   =(int)bmesHyp2->degreesOfFreedom();

	      bsRootTree_->BdFitVtxProb_Hyp1_ = vtxProbHyp1;
	      bsRootTree_->BdFitVtxProb_Hyp2_ = vtxProbHyp2;
	      bsRootTree_->BdFitM_Hyp1_ = b_parHyp1[6];		
	      bsRootTree_->BdFitM_Hyp2_ = b_parHyp2[6];		

	      GlobalVector BdvecHyp1(b_parHyp1[3], b_parHyp1[4], b_parHyp1[5]); // the fitted momentum vector 	
	      bsRootTree_->BdFitEta_Hyp1_ = BdvecHyp1.eta();
	      bsRootTree_->BdFitPt_Hyp1_  = BdvecHyp1.perp();
	      bsRootTree_->BdFitPz_Hyp1_  = BdvecHyp1.z();
	      bsRootTree_->BdFitPhi_Hyp1_ = BdvecHyp1.phi();

	      GlobalVector BdvecHyp2(b_parHyp2[3], b_parHyp2[4], b_parHyp2[5]); // the fitted momentum vector 	
	      bsRootTree_->BdFitEta_Hyp2_ = BdvecHyp2.eta();
	      bsRootTree_->BdFitPt_Hyp2_  = BdvecHyp2.perp();
	      bsRootTree_->BdFitPz_Hyp2_  = BdvecHyp2.z();
	      bsRootTree_->BdFitPhi_Hyp2_ = BdvecHyp2.phi();

	      setFitParHyp1( reftree1 );
	      setFitParHyp2( reftree2 );

	      bsRootTree_->BdFitVtx_x_Hyp1_ = bVertexHyp1->position().x();
	      bsRootTree_->BdFitVtx_y_Hyp1_ = bVertexHyp1->position().y();
	      bsRootTree_->BdFitVtx_z_Hyp1_ = bVertexHyp1->position().z();
	 
	      bsRootTree_->BdFitVtx_x_Hyp2_ = bVertexHyp2->position().x();
	      bsRootTree_->BdFitVtx_y_Hyp2_ = bVertexHyp2->position().y();
	      bsRootTree_->BdFitVtx_z_Hyp2_ = bVertexHyp2->position().z();

              bsRootTree_->BdCowboy_=isCowboy;

	      bsRootTree_->BdM_nofit_ = BdCand.mass();
	      bsRootTree_->BdPt_nofit_ = BdCand.pt();
	      bsRootTree_->BdPz_nofit_ = BdCand.pz();
	      bsRootTree_->BdPhi_nofit_ = BdCand.phi();
	      bsRootTree_->BdEta_nofit_ = BdCand.eta();
	      
	      bsRootTree_->KstarMass_nofit_Hyp1_ = Kstmass1;
	      bsRootTree_->KstarMass_nofit_Hyp2_ = Kstmass2;	      

	      bsRootTree_->BdK1Pt_nofit_   = track1.pt();
	      bsRootTree_->BdK1Pz_nofit_   = track1.pz();
	      bsRootTree_->BdK1Eta_nofit_  = track1.eta();
	      bsRootTree_->BdK1Phi_nofit_  = track1.phi();
	      bsRootTree_->BdK1Key_nofit_  = trkkst1.key();
	      bsRootTree_->BdK2Pt_nofit_   = track2.pt();
	      bsRootTree_->BdK2Pz_nofit_   = track2.pz();
	      bsRootTree_->BdK2Eta_nofit_  = track2.eta();
	      bsRootTree_->BdK2Phi_nofit_  = track2.phi();
	      bsRootTree_->BdK2Key_nofit_  = trkkst2.key();	      

              // Save Jpsi Vertex probability
              bsRootTree_->BdJpsiVtxProb_ = vtxProb_Jpsi;
              bsRootTree_->BdCosDeltaAlpha_ = CosAlpha;
              bsRootTree_->BdMuMuDCA_ = MuonsDCA;
              bsRootTree_->BdMuMuDistance_ = lxy;
              bsRootTree_->BdMuMuDistanceSigma_ = lxyerr;
              bsRootTree_->BdMuDr1_ = max_Dr1;
              bsRootTree_->BdMuDr2_ = max_Dr2;

	      RefCountedKinematicVertex bdVertex;
	      AlgebraicSymMatrix77 bd_er;
	      GlobalVector Bdvec;
	      double Bdmass;
	      if(K1flag==1)       {bdVertex =  bVertexHyp1; bd_er = bd_erHyp1; Bdvec = BdvecHyp1; Bdmass = fittedBdMassHyp1; }
	      else if (K2flag==1) {bdVertex =  bVertexHyp2; bd_er = bd_erHyp2; Bdvec = BdvecHyp2; Bdmass = fittedBdMassHyp2; }
	      else {std::cout<<"error flag" << std::endl;  exit(1);}


	      Int_t PVCosThetaIndex = -1;
	      Double_t MinDistance = 10000000; 	
	      Double_t distance = 0;	

	      for(unsigned int BdVtxInd=0; BdVtxInd<recVtxs->size(); BdVtxInd++){
	        const Vertex &vtx = (*recVtxs)[BdVtxInd];

	        Double_t PVSVvecDotBdPvec = ( bdVertex->position().x()- vtx.x() )*Bdvec.x() + (bdVertex-> position().y() - vtx.y())*Bdvec.y() + (bdVertex-> position().z() - vtx.z() )*Bdvec.z();
	        Double_t PVSVlength = TMath::Sqrt( pow((bdVertex->position().x()- vtx.x()), 2.0) + pow((bdVertex->position().y()- vtx.y()), 2.0) + pow((bdVertex->position().z()- vtx.z()), 2.0) );

    	        Double_t BplusPlength = TMath::Sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y() + Bdvec.z()*Bdvec.z());
	        Double_t BplusCosTheta = PVSVvecDotBdPvec / (BplusPlength * PVSVlength);
	        distance = 1-BplusCosTheta;

	        if(distance < MinDistance){
		  MinDistance = distance;
		  PVCosThetaIndex = BdVtxInd;
	        }			
	      }

	      if(PVCosThetaIndex == -1){ continue;}
           
              BdPVvtxCosTheta = reVertex(iSetup, vertexBeamSpot,  (*recVtxs)[PVCosThetaIndex], mu1, mu2, trkkst1 , trkkst2);



	      // proper decay time and proper decay length with the refitted vertex
	      VertexDistanceXY vdist;	      
	      if(Bdvec.perp()!=0) {
		bsRootTree_->BdLxy_    = vdist.distance( reFitVertex, bdVertex->vertexState() ).value(); 
		bsRootTree_->BdLxyErr_ = vdist.distance( reFitVertex, bdVertex->vertexState() ).error(); 
		if (  (bdVertex->position().x()- reFitVertex.x())*Bdvec.x()+(bdVertex->position().y()-reFitVertex.y())*Bdvec.y() < 0  )
		  bsRootTree_->BdLxy_ = -1.0 * bsRootTree_->BdLxy_;   // in case negative sign is necessary 
		bsRootTree_->BdCt_     = bsRootTree_->BdLxy_     *  Bdmass/Bdvec.perp();
		bsRootTree_->BdCtErr_  = bsRootTree_->BdLxyErr_  *  Bdmass/Bdvec.perp();
	      }
  
              // ctau 2d BS
              bsRootTree_->BdCt2DBS_ = BdPDGMass_*((bdVertex->position().x()-BSx)*Bdvec.x()+(bdVertex->position().y()-BSy)*Bdvec.y())/(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y());
              // error on ctau 2D BS
              VertexDistanceXY d2BS;
              Measurement1D measurement2BS = d2BS.distance(vertexBeamSpot,bdVertex->vertexState());
              double error2DBS = measurement2BS.error();
              double scale2BS = ((bdVertex->position().x() - BSx)*Bdvec.x()+(bdVertex->position().y() - BSy)*Bdvec.y())/(sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y())*sqrt((bdVertex->position().x() - BSx)*(bdVertex->position().x() - BSx)+(bdVertex->position().y() - BSy)*(bdVertex->position().y() - BSy)));
              bsRootTree_->BdCtErr2DBS_=BdPDGMass_*(error2DBS*abs(scale2BS))/sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y());


// 	      if(BdCand.pt()!=0) {
//                 bsRootTree_->BdLxy_ = ((bdVertex->position().x()-PVx)*Bdvec.x()+(bdVertex->position().y()-PVy)*Bdvec.y())/Bdvec.perp();
//                 bsRootTree_->BdCt_  = bsRootTree_->BdLxy_*Bdmass/Bdvec.perp();
//               }

              bsRootTree_->BdErrX_  = bd_er(1,1);
              bsRootTree_->BdErrY_  = bd_er(2,2);
              bsRootTree_->BdErrXY_ = bd_er(1,2); 
              
              VertexDistance3D vdist3d;
              bsRootTree_->BdDist3d_    = vdist3d.distance(bdVertex->vertexState(),RecVtx).value();
              bsRootTree_->BdDist3dErr_ = vdist3d.distance(bdVertex->vertexState(),RecVtx).error();
              bsRootTree_->BdTime3d_    = bsRootTree_->BdDist3d_    * Bdmass/Bdvec.perp() * 100. /3.;
              bsRootTree_->BdTime3dErr_ = bsRootTree_->BdDist3dErr_ * Bdmass/Bdvec.perp() * 100. /3.;
              
                       
              bsRootTree_->BdDist2d_     = vdist.distance(bdVertex->vertexState(),RecVtx).value();
              bsRootTree_->BdDist2dErr_ = vdist.distance(bdVertex->vertexState(),RecVtx).error();
              bsRootTree_->BdTime2d_     = bsRootTree_->BdDist2d_ * Bdmass/Bdvec.perp() *100. /3.;
              bsRootTree_->BdTime2dErr_  = bsRootTree_->BdDist2dErr_ * Bdmass/Bdvec.perp() * 100. /3.;
        
	    // transversity basis angles
            TLorentzVector pmuplus;
	    TLorentzVector pmuminus;
            if (jpsi_children[0]->currentState().particleCharge() == 1) {
	      pmuplus.SetXYZM(bd_par3[3],bd_par3[4],bd_par3[5],bd_par3[6]);
	      pmuminus.SetXYZM(bd_par4[3],bd_par4[4],bd_par4[5],bd_par4[6]);
            } else {
	      pmuminus.SetXYZM(bd_par3[3],bd_par3[4],bd_par3[5],bd_par3[6]);
	      pmuplus.SetXYZM(bd_par4[3],bd_par4[4],bd_par4[5],bd_par4[6]);
            } 
	    
	    TLorentzVector pkplus;
	    TLorentzVector pkminus;
            if (bd_children[0]->currentState().particleCharge() == 1) {
	      pkplus.SetXYZM(bd_par1[3],bd_par1[4],bd_par1[5],bd_par1[6]);
	      pkminus.SetXYZM(bd_par2[3],bd_par2[4],bd_par2[5],bd_par2[6]);
            } else {
	      pkminus.SetXYZM(bd_par1[3],bd_par1[4],bd_par1[5],bd_par1[6]);
	      pkplus.SetXYZM(bd_par2[3],bd_par2[4],bd_par2[5],bd_par2[6]);
            } 
	    
	    // boosting in JPsi restframe
	    TLorentzVector pjpsi;                                                                                                           
	    pjpsi = pmuplus + pmuminus;
	    TLorentzVector pphi;
	    pphi = pkplus + pkminus;
	    
	    // the betas for the boost
	    TVector3 p3_JPsi;
	    p3_JPsi = pjpsi.Vect();
	    p3_JPsi *= -1./pjpsi.E();
	    
	    // the boost matrix
	    TLorentzRotation boost_jpsi(p3_JPsi);
	    TLorentzVector p_JPsi_JPsi;
	    p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);
	    
	    // the different momenta in the new frame                                                                                                       
	    TLorentzVector p_JPsi_muplus;
	    TLorentzVector p_JPsi_Kplus;
	    TLorentzVector p_JPsi_phi;                                                                       
	    p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);                                                                      
	    p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);                                                                                              
	    p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);
	    
	    // the 3-momenta
	    TVector3 p3_JPsi_muplus;
	    p3_JPsi_muplus = p_JPsi_muplus.Vect();
	    TVector3 p3_JPsi_Kplus;
	    p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	    TVector3 p3_JPsi_phi;
	    p3_JPsi_phi = p_JPsi_phi.Vect();
	    
	    // coordinate system
	    TVector3 x,y,z;
	    x = p3_JPsi_phi.Unit();
	    y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	    y = y.Unit();
	    z = x.Cross(y);
	    
	    // Transversity Basis
	    double Bdangle_costheta = p3_JPsi_muplus.Unit() * z;
	    
	    double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	    double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	    double Bdangle_phi = TMath::ACos(cos_phi);
	    if (sin_phi < 0){
	      Bdangle_phi =  -Bdangle_phi;
	    }
	    
	    // boosting in phi restframe                                                                                                          
	    // the betas for the boost
	    TVector3 p3_phi;
	    p3_phi = pphi.Vect();
	    p3_phi *= -1./pphi.E();
	      
	    // the boost matrix
	    TLorentzRotation boost_phi(p3_phi);
	    TLorentzVector p_phi_phi;
	    p_phi_phi = boost_phi.VectorMultiplication(pphi);
	    
	    // the different momenta in the new frame
	    TLorentzVector p_phi_Kplus;
	    TLorentzVector p_phi_JPsi;
	    p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	    p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	    
	    // the 3-momenta
	    TVector3 p3_phi_Kplus;
	    p3_phi_Kplus = p_phi_Kplus.Vect();
	    TVector3 p3_phi_JPsi;
	    p3_phi_JPsi = p_phi_JPsi.Vect();

            bsRootTree_->Bdcostheta_=Bdangle_costheta;
            bsRootTree_->Bdcospsi_= -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
            bsRootTree_->Bdphi_=Bdangle_phi;

              if(verbose_ == true){   std::cout<<"Bd MC matching"<<std::endl;}
	      // deltaR matching
	      bool K1Truth = MCmatching( track1, genParticles, bsRootTree_->BdK1mcId_, bsRootTree_->BdK1momId_, bsRootTree_->BdK1gmomId_, 313, 511);
	      bool K2Truth = MCmatching( track2, genParticles, bsRootTree_->BdK2mcId_, bsRootTree_->BdK2momId_, bsRootTree_->BdK2gmomId_, 313, 511);
	      bool Mu1Truth= MCmatching( mu1,    genParticles, bsRootTree_->BdMu1mcId_,bsRootTree_->BdMu1momId_,bsRootTree_->BdMu1gmomId_, 443, 511);
	      bool Mu2Truth= MCmatching( mu2,    genParticles, bsRootTree_->BdMu2mcId_,bsRootTree_->BdMu2momId_,bsRootTree_->BdMu2gmomId_, 443, 511);
	      
	      if (K1Truth==1 && K2Truth==1 && Mu1Truth==1 && Mu2Truth==1)  bsRootTree_->isMatchedBd_ = 1;
	      else bsRootTree_->isMatchedBd_ = 0; 

	    }	      

	    
	  } // trk2 end loop
	} // trk1 end loop
  
      //////////////// 
		// B+ meson
	  /////////////// 


	Handle<CandidateView> allTracksK;
        iEvent.getByLabel(trackLabelK_, allTracksK);		    

        for(size_t itracks = 0; itracks < allTracksK->size(); itracks++){

	  const Candidate & KplusTrack = (*allTracksK)[itracks];
	  TrackRef trkKplusRef = KplusTrack.get<TrackRef>();
	  if (KplusTrack.pt() < KaonTrackPtCut_) continue; 						
	  //if (KplusTrack.charge()!= 1) continue;

	  TrackRef muonTrack1 = mu1.track();
     TrackRef muonTrack2 = mu2.track();

	  if( muonTrack1->momentum() == KplusTrack.momentum() && muonTrack1->charge() == KplusTrack.charge() ) {continue; } 
	  if( muonTrack2->momentum() == KplusTrack.momentum() && muonTrack2->charge() == KplusTrack.charge() ) {continue; }


          // loose cuts to reduce tracks 
	  if(mu1.pt()<2.5 || mu2.pt()<2.5){ continue; }
	  if(KplusTrack.pt() < 1.5 ){ continue;	} 

          pat::CompositeCandidate Bplus;		

	  Bplus.addDaughter(mu1);
	  Bplus.addDaughter(mu2); 	
	  Bplus.addDaughter(KplusTrack);			
	  AddFourMomenta add4mom;		
	  add4mom.set(Bplus);					

	  if( Bplus.mass() > 6 || Bplus.mass() < 4 ) continue; 

	  edm::ESHandle<TransientTrackBuilder> theBplusBuilder; 
	  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theBplusBuilder);
     TrackRef trkBpMu1Ref = mu1.get<TrackRef>();
	  TrackRef trkBpMu2Ref = mu2.get<TrackRef>();

	  reco::TransientTrack JpsiBpMu1 = (*theBplusBuilder).build(&trkBpMu1Ref);
	  reco::TransientTrack JpsiBpMu2 = (*theBplusBuilder).build(&trkBpMu2Ref);

	  KinematicParticleFactoryFromTransientTrack pFactory;

	  //JPSI FIT AND MASS CONSTRAINT 

	  //if a jpsi formed of two muons, passes the jpsi mass constraint fit, the muon tracks really form a real Jpsi and that jpsi can be used as a Jpsi candidate for B+    	
          //The mass of a muon and the insignificant mass sigma 
          //to avoid singularities in the covariance matrix.
          ParticleMass muon_mass = 0.10565837; //pdg mass
          //ParticleMass pion_mass = 0.13957018;
          ParticleMass kaon_mass = 0.493677;
          //ParticleMass psi_mass = 3.096916;
          //ParticleMass bc_mass = 6.277;
          //float bc_sigma = 0.006;
          float muon_sigma = muon_mass*1.e-6;
          //float pion_sigma = pion_mass*1.e-6;
          float kaon_sigma = 0.000016;

          //initial chi2 and ndf before kinematic fits.
          float chi = 0.;
          float ndf = 0.;

	  std::vector<RefCountedKinematicParticle> JpsiMuons;

  	  JpsiMuons.push_back(pFactory.particle (JpsiBpMu1, nominalMuonMass, chi, ndf, muon_sigma));
	  JpsiMuons.push_back(pFactory.particle (JpsiBpMu2, nominalMuonMass, chi, ndf, muon_sigma));

	  KinematicParticleVertexFitter Fitter;
	  RefCountedKinematicTree JpsiTree = Fitter.fit(JpsiMuons);
  
	  // if the fit fails, do not consider this as candidate
	  if(JpsiTree->isEmpty()) continue;
  
	  KinematicParticleFitter constFitter;
  
  	  float jpsiMsigma = 0.00004;
  	  KinematicConstraint * jpsi_const = new MassKinematicConstraint(nominalJpsiMass,jpsiMsigma);
  
	  JpsiTree = constFitter.fit(jpsi_const,JpsiTree);
 	  if(JpsiTree->isEmpty()) continue;

	  JpsiTree->movePointerToTheTop();
	  RefCountedKinematicParticle JpsiBp = JpsiTree->currentParticle();

	  //when to Jpsi is succesfully reconstructed one can combine kaon track with reconstructed Jpsi candidate		
	  reco::TransientTrack TransTrackKplus = (*theBplusBuilder).build(&trkKplusRef);

	  vector<RefCountedKinematicParticle> allParticlesTrk;
	  allParticlesTrk.push_back(JpsiBp);	
	  allParticlesTrk.push_back(pFactory.particle (TransTrackKplus, nominalKaonMass, chi, ndf, kaon_sigma));

	  KinematicParticleVertexFitter Bplusfitter; 
	  RefCountedKinematicTree Bplus_Tree = Bplusfitter.fit(allParticlesTrk); 

	  if (Bplus_Tree->isEmpty()) continue;
     					
	  Bplus_Tree->movePointerToTheTop();
	  RefCountedKinematicParticle bplusmes = Bplus_Tree->currentParticle();
	  RefCountedKinematicVertex bplusVertex = Bplus_Tree->currentDecayVertex();

	  VertexDistanceXY vtxdist;

          vector<TransientTrack> t_tracksbp;
          t_tracksbp.push_back(JpsiBpMu1);
          t_tracksbp.push_back(JpsiBpMu2);
          t_tracksbp.push_back(TransTrackKplus);

          KalmanVertexFitter kvfbp;
          TransientVertex kvfbpvertex = kvfbp.vertex(t_tracksbp); 
          Vertex vertexbpkalman = kvfbpvertex;
          if (!kvfbpvertex.isValid()) continue;
          GlobalError gigibp=kvfbpvertex.positionError();

          TMatrix covabp(2,2);
          covabp.IsSymmetric();
          covabp(0,0)=gigibp.cxx();
          covabp(1,1)=gigibp.cyy();
          covabp(0,1)=gigibp.cyx();
          covabp(1,0)=gigibp.cyx();

	  //Double_t vtxProbBplus = TMath::Prob(bplusmes->chiSquared(),(int)bplusmes->degreesOfFreedom());
          double vtxProbBplus = TMath::Prob(vertexbpkalman.chi2(),(int)vertexbpkalman.ndof());
          

	  if(bestVtxProbBplus < vtxProbBplus){ 
        bestVtxProbBplus = vtxProbBplus;
	//		cout << "B+ candidate found, vtx prob "<< bestVtxProbBplus << endl;

			BpTrkRefs.clear(); 		
			BpTrkRefs.push_back(trkBpMu1Ref);
			BpTrkRefs.push_back(trkBpMu2Ref);	
		   BpTrkRefs.push_back(trkKplusRef);

	     AlgebraicVector7 bplus_par = bplusmes->currentState().kinematicParameters().vector();

		 
       bsRootTree_->BplusCharge_ = KplusTrack.charge();
	    bsRootTree_->BplusM_fit_ = bplus_par[6];
	    bsRootTree_->BplusVtxProb_ = vtxProbBplus;
	    bsRootTree_->BplusChi2_ = bplusmes->chiSquared();
	    bsRootTree_->BplusPt_ =  sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) );
	    bsRootTree_->BplusPtot_ = sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) + pow(bplus_par[5],2.0) );
	    bsRootTree_->KplusPt_ = KplusTrack.pt();
	    bsRootTree_->KplusPtot_ = KplusTrack.p();
	    bsRootTree_->BplusMu1Pt_ = mu1.pt(); 
	    bsRootTree_->BplusMu2Pt_ = mu2.pt();
	    bsRootTree_->BplusMu1Eta_ = mu1.eta(); 
	    bsRootTree_->BplusMu2Eta_ = mu2.eta(); 

            // Save Jpsi Vertex probability
            bsRootTree_->BpJpsiVtxProb_ = vtxProb_Jpsi;
            bsRootTree_->BpCosDeltaAlpha_ = CosAlpha;
            bsRootTree_->BpMuMuDCA_ = MuonsDCA;
            bsRootTree_->BpMuMuDistance_ = lxy;
            bsRootTree_->BpMuMuDistanceSigma_ = lxyerr;
            bsRootTree_->BpMuDr1_ = max_Dr1;
            bsRootTree_->BpMuDr2_ = max_Dr2;

            if (mu1.isPFMuon())       bsRootTree_->BpMuonCat1_ = 1;
            if (mu2.isPFMuon())       bsRootTree_->BpMuonCat2_ = 1;

	    bsRootTree_->BplusMu1Ptot_ = mu1.p(); 
	    bsRootTree_->BplusMu2Ptot_ = mu2.p();

	    GlobalVector BplusVec(bplus_par[3], bplus_par[4], bplus_par[5]);
     	    bsRootTree_->BplusEta_ = BplusVec.eta();
	    bsRootTree_->BplusPhi_ = BplusVec.phi();

	    pat::CompositeCandidate Jpsi_bplus;
	    Jpsi_bplus.addDaughter(mu1);
	    Jpsi_bplus.addDaughter(mu2);
      	    AddFourMomenta addP4;
      	    addP4.set(Jpsi_bplus);

	    bsRootTree_->JpsiMass_bplus_ = Jpsi_bplus.mass();
	    //Jpsi pt from the JpsiTree		
	    AlgebraicVector7 jpsi_par = JpsiBp->currentState().kinematicParameters().vector();
	    bsRootTree_->JpsiPtbplus_fit_ = sqrt( pow(jpsi_par[3],2.0) + pow(jpsi_par[4],2.0) );	
	    bsRootTree_->JpsiPt_bplus_ = Jpsi_bplus.pt();	

            GlobalPoint BplusVtxPos( bplusVertex->position().x(), bplusVertex->position().y(), bplusVertex->position().z());
	    GlobalError BplusVtxErr( bplusVertex->error() );	

	    Int_t PVCosThetaIndex = -1;
	    Double_t MinDistance = 10000000; 	
	    Double_t distance = 0;	

	    for(unsigned int BpVtxInd=0; BpVtxInd<recVtxs->size(); BpVtxInd++){
	      const Vertex &vtx = (*recVtxs)[BpVtxInd];

	      Double_t PVSVvecDotBplusPvec = ( bplusVertex-> position().x()- vtx.x() )*BplusVec.x() + (bplusVertex-> position().y() - vtx.y())*BplusVec.y() + (bplusVertex-> position().z() - vtx.z() )*BplusVec.z();
	      Double_t PVSVlength = TMath::Sqrt( pow((bplusVertex->position().x()- vtx.x()), 2.0) + pow((bplusVertex->position().y()- vtx.y()), 2.0) + pow((bplusVertex->position().z()- vtx.z()), 2.0) );

    	      Double_t BplusPlength = TMath::Sqrt(BplusVec.x()*BplusVec.x()+BplusVec.y()*BplusVec.y() + BplusVec.z()*BplusVec.z());
	      Double_t BplusCosTheta = PVSVvecDotBplusPvec / (BplusPlength * PVSVlength);
	      distance = 1-BplusCosTheta;

	      if(distance < MinDistance){
		MinDistance = distance;
		PVCosThetaIndex = BpVtxInd;
	      }			
	    }

       bsRootTree_->BplusPVindex_ = PVCosThetaIndex;
	    if(PVCosThetaIndex == -1){ continue;}
           
	    reco::Vertex BplusPVVtx = (*recVtxs)[PVCosThetaIndex]; 
	    double BplusPVx = BplusPVVtx.x(); 
	    double BplusPVy = BplusPVVtx.y();
	    double BplusPVz= BplusPVVtx.z();
            double BplusPVerrx=BplusPVVtx.xError();
	    double BplusPVerry=BplusPVVtx.yError();
      	    double BplusPVerrz=BplusPVVtx.zError();
            reco::TrackRef nulltrack;

            BpPVvtxCosTheta = reVertex(iSetup, vertexBeamSpot,  (*recVtxs)[PVCosThetaIndex], mu1, mu2, trkKplusRef , nulltrack);
            bsRootTree_->BpCt2DPVCosTheta_ = BpPDGMass_*( (kvfbpvertex.position().x()-BpPVvtxCosTheta.x())*Bplus.px() + (kvfbpvertex.position().y()-BpPVvtxCosTheta.y())*Bplus.py() )/( Bplus.px()*Bplus.px() + Bplus.py()*Bplus.py() );
   
            VertexDistanceXY d2Costheta;
            Measurement1D measurement2Costheta = d2Costheta.distance(BpPVvtxCosTheta,kvfbpvertex);
            double error2DCostheta = measurement2Costheta.error();
            double scale2Costheta = ((kvfbpvertex.position().x()-BpPVvtxCosTheta.x())*Bplus.px()+(kvfbpvertex.position().y()-BpPVvtxCosTheta.y())*Bplus.py())/
                (sqrt(Bplus.px()*Bplus.px()+Bplus.py()*Bplus.py())*
                 sqrt((kvfbpvertex.position().x() - BplusPVx)*(kvfbpvertex.position().x() - BplusPVx)+(kvfbpvertex.position().y() -BpPVvtxCosTheta.y())*(kvfbpvertex.position().y() - BplusPVy)));
            TVector LengthVector(2);
            LengthVector(0)=kvfbpvertex.position().x()-BpPVvtxCosTheta.x();
            LengthVector(1)=kvfbpvertex.position().y()-BpPVvtxCosTheta.y();
            bsRootTree_->BpCtErr2DCostheta_=sqrt(pow(BpPDGMass_*(error2DCostheta*abs(scale2Costheta))/sqrt(Bplus.px()*Bplus.px()+Bplus.py()*Bplus.py()),2)+pow(BpPDGMass_*abs(scale2Costheta)/(Bplus.px()*Bplus.px() + Bplus.py()*Bplus.py()),2)*covabp.Similarity(LengthVector));


            vector<TransientTrack> trk_BpJpsi;
            trk_BpJpsi.push_back( (*theBplusBuilder).build(&trkBpMu1Ref) );
            trk_BpJpsi.push_back( (*theBplusBuilder).build(&trkBpMu2Ref) );      
     
            KalmanVertexFitter kvfBpJpsi;
            TransientVertex tvBpJpsi = kvfBpJpsi.vertex(trk_BpJpsi);

	    std::pair<bool,Measurement1D> ImpactPar3DKandJpsiVtx = IPTools::absoluteImpactParameter3D(TransTrackKplus, tvBpJpsi); 

	    if(ImpactPar3DKandJpsiVtx.first){
              bsRootTree_->IP3DKandJpsiVtx_ = ImpactPar3DKandJpsiVtx.second.value();
  	      bsRootTree_->IP3DKandJpsiVtxErr_ = ImpactPar3DKandJpsiVtx.second.error();
    	    }

	    bsRootTree_->BpmatchDoubleMu01_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
            bsRootTree_->BpmatchDoubleMu02_ = !mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsi").empty();
            bsRootTree_->BpmatchDoubleMu41_ = !mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
            bsRootTree_->BpmatchDoubleMu42_ = !mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty();
	    bsRootTree_->BpmatchDoubleMu01DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
	    bsRootTree_->BpmatchDoubleMu02DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();


	    // deltaR matching for Bplus 

	    bool BplusKTruth = MCmatchingBplusK(KplusTrack, genParticles, bsRootTree_->BplusKmcId_, bsRootTree_->BplusKmomId_,521);
	    bool BplusMu1Truth = MCmatching(mu1, genParticles, bsRootTree_->BplusMu1mcId_, bsRootTree_->BplusMu1momId_,bsRootTree_->BplusMu1gmomId_,443,521);
	    bool BplusMu2Truth= MCmatching( mu2,genParticles, bsRootTree_->BplusMu2mcId_,bsRootTree_->BplusMu2momId_,bsRootTree_->BplusMu2gmomId_, 443, 521);
	    if (BplusKTruth && BplusMu1Truth && BplusMu2Truth) bsRootTree_->isMatchedBplus_ = 1;
	    else bsRootTree_->isMatchedBplus_ = 0; 
          } //end of if sentence
	} //end of B+ loop 
      
    } // loop on muons2
  } // loop on muons1

  // Loop again on tracks to search for the Bc
  Handle<CandidateView> allTracksPi;
  iEvent.getByLabel(trackLabelPi_, allTracksPi);
  //iEvent.getByLabel(trackLabelK_, allTracksPi);

  double BcAngStore=99999;
  double minBcProb=-999.;
  double maxBcP=-999999.;

  if (minVtxP>0 && bsRootTree_->BsFitM_>5.00 && bsRootTree_->BsFitM_<5.43) {
  //if (minVtxP>0) {
   for (size_t i=0; i< allTracksPi->size(); ++i){
      const Candidate & track1 = (*allTracksPi)[i];
      //if (track1.pt() < 0.4) continue;
      pat::CompositeCandidate BcCand;
      BcCand.addDaughter(BCand_best);
      BcCand.addDaughter(track1);
      AddFourMomenta add4mom;
      add4mom.set(BcCand);
      //CHECK
      //cout<<"Pion: "<<track1.energy()<<" - "<<track1.px()<<" - "<<track1.py()<<" - "<<track1.pz()<<endl;
      //cout<<"Bs  : "<<BCand_best.energy()<<" - "<<BCand_best.px()<<" - "<<BCand_best.py()<<" - "<<BCand_best.pz()<<endl;
      //cout<<"Bc  : "<<BcCand.energy()<<" - "<<BcCand.px()<<" - "<<BcCand.py()<<" - "<<BcCand.pz()<<endl;
      // PDG Bc mass 6.277
      if (BcCand.mass() < 5.50 || BcCand.mass() > 7) continue;
    //  cout<<"------------------- Bc mass = "<<BcCand.mass()<<endl;
       
      ESHandle<MagneticField> bFieldHandle;
      iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

      // TransientTrack kaon1TT(trk1Ref_best, &(*bFieldHandle) );
      // TransientTrack kaon2TT(trk2Ref_best, &(*bFieldHandle) );
      // TransientTrack muon1TT(trkMu1Ref_best, &(*bFieldHandle) );
      // TransientTrack muon2TT(trkMu2Ref_best, &(*bFieldHandle) );
      TrackRef trk3Ref = track1.get<TrackRef>();
      TransientTrack pion1TT(trk3Ref, &(*bFieldHandle) );

	  

      if (trk3Ref->numberOfValidHits() < 5) continue;


      if (trk3Ref==trk1Ref_best) { continue;}
      if (trk3Ref==trk2Ref_best) { continue;}
      if (trkMu1Ref_best==trk3Ref) { continue;}
      if (trkMu2Ref_best==trk3Ref) { continue;}

      edm::ESHandle<TransientTrackBuilder> theB;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

      //TransientTrack pion1TT;
      //pion1TT=(*theB).build(&trk3Ref);
      //trk_all.push_back((*theB).build(&trkMu2Ref));

      vector<TransientTrack> t_tracks;
      t_tracks.push_back((*theB).build(&trkMu1Ref_best));
      t_tracks.push_back((*theB).build(&trkMu2Ref_best));
      t_tracks.push_back((*theB).build(&trk1Ref_best));
      t_tracks.push_back((*theB).build(&trk2Ref_best));

      //call fit interface
      KinematicFitInterface Kfitter2;
      bool fitSuccess = Kfitter2.doFit(t_tracks, nominalMuonMass,  nominalKaonMass, nominalKaonMass);

      if(fitSuccess != 1) {
	//cout<<"Bs fit fail"<<endl;	 
	continue;}
      RefCountedKinematicParticle bs_best = Kfitter2.getParticle();
      RefCountedKinematicVertex BsVertex = Kfitter2.getVertex();
      RefCountedKinematicTree Bs4BcTree = Kfitter2.getTree();

      // if the fit fails, do not consider this as candidate
      if(Bs4BcTree->isEmpty()) continue;

      KinematicParticleFitter constBsFitter;

      float BsMsigma = 0.0002;
      KinematicConstraint * bs_const = new MassKinematicConstraint( BsPDGMass_, BsMsigma);
      KinematicParticleFitter constFitter;
      RefCountedKinematicTree BsTreeN1 = constFitter.fit(bs_const,Bs4BcTree);

      // if the fit fails, do not consider this as candidate
      if(BsTreeN1->isEmpty()) continue;
      BsTreeN1->movePointerToTheTop();
      RefCountedKinematicParticle bsConstr = BsTreeN1->currentParticle();

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack p_Factory;
	       
      //The mass of a muon and the insignificant mass sigma 
      //to avoid singularities in the covariance matrix.
      //ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass pion_mass = 0.13957018;
      //ParticleMass kaon_mass = 0.493677;
      //ParticleMass psi_mass = 3.096916;
      //ParticleMass bc_mass = 6.277;
      //float bc_sigma = 0.006;
      //float muon_sigma = muon_mass*1.e-6;
      float pion_sigma = pion_mass*1.e-6;
      //float kaon_sigma = 0.000016;
	       
      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      // vector<RefCountedKinematicParticle> kaonParticles;
      // vector<RefCountedKinematicParticle> muonParticles;
      vector<RefCountedKinematicParticle> bcParticles;
      // kaonParticles.push_back(pFactory.particle(kaon1TT,kaon_mass,chi,ndf,kaon_sigma));
      // kaonParticles.push_back(pFactory.particle(kaon2TT,kaon_mass,chi,ndf,kaon_sigma));
      // muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
      // muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
      bcParticles.push_back(p_Factory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
      bcParticles.push_back(bsConstr);
      //bcParticles.push_back(bs_best);
 
      //KinematicParticleVertexFitter Bcfitter;   
      KinematicConstrainedVertexFitter Bcfitter;   
      RefCountedKinematicTree BcVertexFitTree = Bcfitter.fit(bcParticles);
      //KinematicParticleFitter fitter;
      //KinematicConstraint * Bc_c = new MassKinematicConstraint(bc_mass,bc_sigma);
      // add mass constraint to the ks fit to do a constrained fit:  
      //BcVertexFitTree = fitter.fit(Bc_c,BcVertexFitTree);

      if (BcVertexFitTree->isEmpty()) {
        //std::cout << "Empty vertex from the Bc vertex fit" << std::endl;
        continue;
       }
      if (!BcVertexFitTree->isValid()) {
        //std::cout << "invalid vertex from the Bc vertex fit" << std::endl;
        continue; 
       }

      BcVertexFitTree->movePointerToTheTop();
      RefCountedKinematicParticle BcFromFit = BcVertexFitTree->currentParticle();
      RefCountedKinematicVertex BcVertex = BcVertexFitTree->currentDecayVertex();

                int ClosestPVindex = -1;
	        Double_t MinDistance = 10000000; 	
                double BcPVx;	
                double BcPVy;	
                double BcPVz;
		//PV selection with minimum z distance between PV and SV
		for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++){
			const Vertex &vtx = (*recVtxs)[VtxInd];
			Double_t PVSVdistance = TMath::Abs(BcVertex->position().z()-vtx.z());
			if(PVSVdistance < MinDistance){
				MinDistance = PVSVdistance;
				ClosestPVindex = VtxInd;	
			}			
		}
                if (ClosestPVindex==-1){ continue;}
	        const reco::Vertex &BcRecVtx = (*recVtxs)[ClosestPVindex]; 
		BcPVx = BcRecVtx.x(); 
	        BcPVy = BcRecVtx.y();
	        BcPVz= BcRecVtx.z();
                std::pair<bool,Measurement1D> ImpactPar3D = IPTools::absoluteImpactParameter3D(pion1TT, BcRecVtx);
              //  cout<<"Impact parameter 3D: "<<ImpactPar3D.second.value()<<endl;

      AlgebraicVector7 b_par = BcFromFit->currentState().kinematicParameters().vector();
      GlobalVector Bcvec(b_par[3], b_par[4], b_par[5]);



      math::XYZVector Bcpperp(b_par[3], b_par[4], b_par[5]);
      reco::Vertex::Point Bcvperp(BcVertex->position().x(),BcVertex->position().y(),BcVertex->position().z());


      double vtxprob_Bc = TMath::Prob(BcVertex->chiSquared(), (int)BcVertex->degreesOfFreedom());
      double AngBcBs = reco::deltaR(BCand_best.eta(),BCand_best.phi(),Bcvec.eta(),Bcvec.phi());
      double AngBcPi = reco::deltaR(Bcvec.eta(),Bcvec.phi(),track1.eta(),track1.phi());
      double BcCt;

      VertexDistanceXY vdist;
      VertexDistance3D dist3;
      double BcLxy    = vdist.distance( BcRecVtx, BcVertex->vertexState() ).value();
      double BsLxy    = vdist.distance( BcRecVtx, BsVertex->vertexState() ).value();
      double BcLxyz    = dist3.distance( BcRecVtx, BcVertex->vertexState() ).value();
      double BsLxyz    = dist3.distance( BcRecVtx, BsVertex->vertexState() ).value();
      if(sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z())!=0) {
        double BcL_ =(BcVertex->position().x()-BcPVx)*Bcvec.x()+(BcVertex->position().y()-BcPVy)*Bcvec.y();
        //BcLxy =sqrt(pow(BcVertex->position().x()-PVx,2)+pow(BcVertex->position().y()-PVy,2));
        //BsLxy =sqrt(pow(BsVertex->position().x()-PVx,2)+pow(BsVertex->position().y()-PVy,2));
        //BcLxyz =sqrt(pow(BcVertex->position().x()-PVx,2)+pow(BcVertex->position().y()-PVy,2));
        //BsLxyz =sqrt(pow(BsVertex->position().x()-PVx,2)+pow(BsVertex->position().y()-PVy,2));
        //double BcLxyErr_ = vdist.distance( reFitVertex, bVertex->vertexState() ).error();
        if ( BcL_ < 0  )
          BcL_ = -1.0 * BcL_;   // in case negative sign is necessary 
        BcCt=BcL_*6.277/(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y());
      }
      else {continue;}
 //     cout<<"Distanza Bs "<<BsLxy<<" Distanza Bc "<<BcLxy<<endl;
 //     cout<<"3D Distanza Bs "<<BsLxyz<<" Distanza Bc "<<BcLxyz<<endl;
      if (BsLxy<BcLxy) continue;
      //if (vtxprob_Bc<0.05) continue;
      if (AngBcPi>1.) continue;
      if (TMath::Abs(Bcvperp.Dot(Bcpperp)/(Bcvperp.R()*Bcpperp.R()))<0.9) continue;
      if (ImpactPar3D.second.value()>0.5) continue;
      //if (BcCt<0.01) continue;
      if (track1.pt()<0.3) continue;
      //if (AngBcBs<BcAngStore) {
      if (energyLossHandle.isValid() && track1.p()<1) {
        const DeDxDataValueMap &  eloss  = *energyLossHandle;
        double dedxTrk = eloss[trk3Ref].dEdx();
        //double errdedxTrk = eloss[trk3Ref].dEdxError();
        //int NumdedxTrk = eloss[trk3Ref].numberOfMeasurements();
        double trackmass=TMath::Sqrt((dedxTrk-2.557)/2.579)*track1.p();
//        cout<<"Track mass: "<<trackmass<<endl;
//        if (trackmass<0.35 || trackmass!=trackmass) cout<<"Probable pion presence!!!!!!!!!!!!!!!!!"<<endl;
      }
 //     cout<<"AngBcPi: "<<AngBcPi<<endl;
      

      if (vtxprob_Bc > minBcProb){
      //if (sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z()) > maxBcP){
//        cout<<"Bc mass (fitted) = "<<BcFromFit->currentState().mass()<<endl;
//        cout<<"Pion momentum= "<<track1.p()<<endl;
        minBcProb=vtxprob_Bc;
        maxBcP=sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z());
        BcAngStore=AngBcBs;
        bsRootTree_->BcCt_=BcCt;
        bsRootTree_->BcCosAlpha_ = Bcvperp.Dot(Bcpperp)/(Bcvperp.R()*Bcpperp.R());
        //bsRootTree_->BcP_=sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y()+Bcvec.z()*Bcvec.z());
        bsRootTree_->BcIP3D_=ImpactPar3D.second.value();
        bsRootTree_->BcP_=sqrt(Bcvec.x()*Bcvec.x()+Bcvec.y()*Bcvec.y());
        bsRootTree_->BctrackPt_ = track1.pt(); 
        bsRootTree_->BcM_ = BcFromFit->currentState().mass(); 
        bsRootTree_->BcAng_ = AngBcBs;
        bsRootTree_->BcAngPi_ = AngBcPi;
        bsRootTree_->BcProb_ = vtxprob_Bc;
      } 
    }
  } 

	////////////////////// JET TAGGING Bs BEGIN //////////////////

	vector<int> jetIndexesBs;
	vector<int> sharedjetTrksBs;
	double BprobMaxBjetBs=-1;
	bool ismatchedtrackBs;
	int sharedTrkCounterBs,BjetIndexBs=-1;
		
   if( PVvtxCosTheta.isValid() &&  BsTrkRefs.size() == 4 ){
/*
			cout << "Bs track ref vec size "<< BsTrkRefs.size() << endl;
		   cout << "Bs track ref 0 " << BsTrkRefs[0]->pt() << endl;
			cout << "Bs track ref 1 " << BsTrkRefs[1]->pt() << endl;
			cout << "Bs track ref 2 " << BsTrkRefs[2]->pt() << endl;
			cout << "Bs track ref 2 " << BsTrkRefs[3]->pt() << endl;
*/
	for(size_t i=0; i<jets.size(); i++){ //loop over jets
	  const reco::TrackRefVector & jetTrackVec = jets[i].associatedTracks();
 
	  if(jets[i].isPFJet() == false) continue;

	    ismatchedtrackBs = false;
	    sharedTrkCounterBs = 0;	

	  for(size_t k=0; k<jetTrackVec.size(); k++){ //loop ovet jet tracks
		reco::TrackRef jetTrk = jetTrackVec[k];

		if(BsTrkRefs[0] == jetTrk || BsTrkRefs[1] == jetTrk){continue; cout << "jet trk overlap" << endl;}

		if(BsTrkRefs[2] == jetTrk || BsTrkRefs[3] == jetTrk){continue; cout << "jet trk overlap" << endl;} 	
	
		for(reco::Vertex::trackRef_iterator trackvertex=PVvtxCosTheta.tracks_begin(); trackvertex!=PVvtxCosTheta.tracks_end(); trackvertex++){ //loop over PV tracks
			
		 reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
		 if(jetTrk == PVtrk ) {sharedTrkCounterBs++; ismatchedtrackBs = true; }	

	     }// end of PV track loop 			
	   } // end of jet track loop
		int counter =0;
	 if(ismatchedtrackBs == true){
	  jetIndexesBs.push_back(i); 
	  sharedjetTrksBs.push_back(sharedTrkCounterBs);	
	 }
	} //end of jet loop

//	cout << "B jets Bs " << jets.size() << endl;
//	cout << " PV jet cands Bs "<< jetIndexesBs.size() << endl;


	for(size_t i =0; i<jetIndexesBs.size(); i++ ){			
	    double Bprobability = jets[jetIndexesBs[i]].bDiscriminator("jetProbabilityBJetTags"); 
		
	    if(Bprobability <= 0) continue;
//		   cout << "Bjet prob Bs "<< Bprobability << " current Bprobmax " << BprobMaxBjetBs << endl;

	    if(Bprobability > 0 && Bprobability > BprobMaxBjetBs){
		BjetIndexBs = jetIndexesBs[i];
		BprobMaxBjetBs = Bprobability;				
	  }	
	}	

//	cout << "Final Bjet prob "<<BprobMaxBjetsBs << endl;
 //       cout << "Bjet indx "<<BjetIndexBs << endl;
        double ptmax = -1.;

	if(BjetIndexBs == -1){
	 //cout << "Pt looppi" << endl;
	  for(size_t i =0; i<jetIndexesBs.size(); i++ ){
//	   cout << "jet pt, no b jet Bp " << jets[jetIndexesBp[i]].pt() << endl;
	     if( jets[jetIndexesBs[i]].pt() > ptmax){
		ptmax = jets[jetIndexesBs[i]].pt();
		BjetIndexBs = jetIndexesBs[i];
	     }
	  }
	}
//	cout << "Final Bjet index " << BjetIndexBs << endl;
//	loop over PV tracks	

	if(BjetIndexBs!=-1){

//	  if( jets[BjetIndexBp].bDiscriminator("jetProbabilityBJetTags") > 0){
	  bsRootTree_->JetBTagProb_ =  jets[BjetIndexBs].bDiscriminator("jetProbabilityBJetTags"); 
//	  cout <<"Saved Bjet prob Bs " << jets[BjetIndexBs].bDiscriminator("jetProbabilityBJetTags")<< endl;
//	cout << "B jet index" << BjetIndexBp << endl;
//   cout << "jet coll. size " <<jets.size() << endl;
//	  }
		
	  bsRootTree_->BJetPt_ = jets[BjetIndexBs].pt();
	  bsRootTree_->BJetEta_ = jets[BjetIndexBs].eta();
	  bsRootTree_->BJetPhi_ = jets[BjetIndexBs].phi();

//	  cout<< "Jet Eta " << jets[BjetIndexBs].eta() << endl;
//	  cout<< "Jet Phi " << jets[BjetIndexBs].phi() << endl;
//	  cout<< "Jet Pt " << jets[BjetIndexBs].pt() << endl;

	  if(isMCstudy_ && jets[BjetIndexBs].genParton()){
	    const reco::GenParticle *JetParton = jets[BjetIndexBs].genParton();
	    Int_t PartonId = JetParton->pdgId();
//	    cout << "B jet origin " << PartonId << endl;
	    bsRootTree_->BJetParton_ = PartonId;	
	  }

	  const reco::TrackRefVector & BjetTrackVec = jets[BjetIndexBs].associatedTracks();
	
	  for(size_t u = 0; u < BjetTrackVec.size(); u++){
	    reco::TrackRef jetTrk = BjetTrackVec[u];	
	    const reco::Track *BjetTrack = jetTrk.get();
//    cout << "track pt " << BjetTrack->pt() << endl;
//	    cout << "track charge " << BjetTrack->charge() << endl;
	    bsRootTree_->BJetTrkCharge_->push_back(BjetTrack->charge());
	    bsRootTree_->BJetTrkPt_->push_back(BjetTrack->pt());
	  }


	  int trkcouter = 0; 
	  bool OverlappingTrack;
	
	  for(reco::Vertex::trackRef_iterator trackvertex=PVvtxCosTheta.tracks_begin(); trackvertex!=PVvtxCosTheta.tracks_end(); trackvertex++){ //loop over Bs PV tracks

	  reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
	  OverlappingTrack = false;

	   for(size_t l = 0; l < recVtxs->size(); l++){ // loop over other PVs
	  
	    if(l != size_t(BsPVVtxInd) ){	 
		 const Vertex &vertex = (*recVtxs)[l];	    	

		    for(reco::Vertex::trackRef_iterator trkvertex=vertex.tracks_begin(); trkvertex!=vertex.tracks_end(); trkvertex++){ //loop over trks of other PVs
		       reco::TrackRef trk = trkvertex->castTo<reco::TrackRef>();

		       if(trk == PVtrk ){OverlappingTrack = true; }	 	
		     } // end of loop over trks of other PVs
		    }	
	    } // end of loop over other PVs

	//	 if(PVtrk == BsTrkRefs[0] || PVtrk == BsTrkRefs[1] || PVtrk == BsTrkRefs[2] || PVtrk ==BsTrkRefs[3]){cout << "PV track overlap with Bs tracks" << endl;}	

		//if trk of Bs PV doesn't belong to tracks of another PV, the track is saved
		  if(OverlappingTrack == false){
			if(PVtrk != BsTrkRefs[0] && PVtrk != BsTrkRefs[1] && PVtrk != BsTrkRefs[2] && PVtrk != BsTrkRefs[3]){
		   const reco::Track *RTv = trackvertex->get();
		   bsRootTree_->PVTrkPt_->push_back(RTv->pt());
		   bsRootTree_->PVTrkCharge_->push_back(RTv->charge());
		   bsRootTree_->PVTrkEta_->push_back(RTv->eta());
		   bsRootTree_->PVTrkPhi_->push_back(RTv->phi());
			}
		  }
	    } //end of loop over Bs PV tracks
		
	} // end of if(BjetIndex!=-1 )

  } //end of if BsPVindex and BsTrkRef vec size

	////////////////////////// JET TAGGING Bs END /////////////////




	////////////////////// JET TAGGING B+ BEGIN //////////////////

	vector<int> jetIndexesBp;
	vector<int> sharedjetTrksBp;
	double BprobMaxBjetBp=-1;
	bool ismatchedtrack;
	int sharedTrkCounterBp,BjetIndexBp=	-1;
		
//	const Vertex &vtx = (*recVtxs)[PVCosThetaIndex];
   if( bsRootTree_->BplusPVindex_ != -999 &&  BpTrkRefs.size() == 3 ){

 /*  cout << "Bp track ref vec size, jee "<< BpTrkRefs.size() << endl;
		   cout << "Bp track ref 0 " << BpTrkRefs[0]->pt() << endl;
			cout << "Bp track ref 1 " << BpTrkRefs[1]->pt() << endl;
			cout << "Bp track ref 2 " << BpTrkRefs[2]->pt() << endl;
*/
	for(size_t i=0; i<jets.size(); i++){
	  const reco::TrackRefVector & jetTrackVec = jets[i].associatedTracks();
 
	  if(jets[i].isPFJet() == false) continue;

	    ismatchedtrack = false;
	    sharedTrkCounterBp = 0;	

	  for(size_t k=0; k<jetTrackVec.size(); k++){
		reco::TrackRef jetTrk = jetTrackVec[k];

		if(BpTrkRefs[0] == jetTrk || BpTrkRefs[1] == jetTrk){continue; cout << "jet trk overlap" << endl;}
		if(BpTrkRefs[2] == jetTrk ) {continue;}
 		
		for(reco::Vertex::trackRef_iterator trackvertex=BpPVvtxCosTheta.tracks_begin(); trackvertex!=BpPVvtxCosTheta.tracks_end(); trackvertex++){
			
		reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
		if(jetTrk == PVtrk ) {sharedTrkCounterBp++; ismatchedtrack = true; }	

	     }// end of PV track loop 			
	   } // end of jet track loop

	 if(ismatchedtrack == true){
	  jetIndexesBp.push_back(i); 
	  sharedjetTrksBp.push_back(sharedTrkCounterBp);	
	 }
	} //end of jet loop

//	cout << "B jets Bp " << jets.size() << endl;
//	cout << " PV jet cands Bp "<< jetIndexesBp.size() << endl;
//	cout << "shared tracks " << sharedTrkCounterBp << endl;	

	for(size_t i =0; i<jetIndexesBp.size(); i++ ){			
	    double Bprobability = jets[jetIndexesBp[i]].bDiscriminator("jetProbabilityBJetTags"); 
		
	    if(Bprobability <= 0) continue;
//		   cout << "Bjet prob Bp "<< Bprobability << " current Bprobmax " << BprobMaxBjetBp << endl;

	    if(Bprobability > 0 && Bprobability > BprobMaxBjetBp){
		BjetIndexBp = jetIndexesBp[i];
		BprobMaxBjetBp = Bprobability;		
	  }	
	}	

//	cout << "Final Bjet prob "<<BprobMaxBjetBp << endl;
//        cout << "Bjet indx "<<BjetIndexBp << endl;
        double ptmax = -1.;

	if(BjetIndexBp == -1){
	 //cout << "Pt looppi" << endl;
	  for(size_t i =0; i<jetIndexesBp.size(); i++ ){
//	   cout << "jet pt, no b jet Bp " << jets[jetIndexesBp[i]].pt() << endl;
	     if( jets[jetIndexesBp[i]].pt() > ptmax){
		ptmax = jets[jetIndexesBp[i]].pt();
		BjetIndexBp = jetIndexesBp[i];
	     }
	  }
	}
//	cout << "Final Bjet index " << BjetIndexBp << endl;
//	loop over PV tracks	

	if(BjetIndexBp!=-1){

//	  if( jets[BjetIndexBp].bDiscriminator("jetProbabilityBJetTags") > 0){
	  bsRootTree_->BpJetBTagProb_ =  jets[BjetIndexBp].bDiscriminator("jetProbabilityBJetTags"); 
	//  cout <<"Saved Bjet prob Bp " << jets[BjetIndexBp].bDiscriminator("jetProbabilityBJetTags")<< endl;
//	cout << "B jet index" << BjetIndexBp << endl;
//   cout << "jet coll. size " <<jets.size() << endl;
//	  }
		
	  bsRootTree_->BpBJetPt_ = jets[BjetIndexBp].pt();
	  bsRootTree_->BpBJetEta_ = jets[BjetIndexBp].eta();
	  bsRootTree_->BpBJetPhi_ = jets[BjetIndexBp].phi();

//	  cout<< "Jet Eta " << jets[BjetIndexBp].eta() << endl;
//	  cout<< "Jet Phi " << jets[BjetIndexBp].phi() << endl;
//	  cout<< "Jet Pt " << jets[BjetIndexBp].pt() << endl;

	  if(isMCstudy_ && jets[BjetIndexBp].genParton()){
	    const reco::GenParticle *JetParton = jets[BjetIndexBp].genParton();
	    Int_t PartonId = JetParton->pdgId();
//	    cout << "B jet origin " << PartonId << endl;
	    bsRootTree_->BpBJetParton_ = PartonId;	
	  }

	  const reco::TrackRefVector & BjetTrackVec = jets[BjetIndexBp].associatedTracks();
	
	
	  for(size_t u = 0; u < BjetTrackVec.size(); u++){
	    reco::TrackRef jetTrk = BjetTrackVec[u];	
	    const reco::Track *BjetTrack = jetTrk.get();
//    cout << "track pt " << BjetTrack->pt() << endl;
//	    cout << "track charge " << BjetTrack->charge() << endl;
	    bsRootTree_->BpBJetTrkCharge_->push_back(BjetTrack->charge());
	    bsRootTree_->BpBJetTrkPt_->push_back(BjetTrack->pt());
	  }


	  int trkcouter = 0; 
	  bool OverlappingTrack;
	 
       size_t BpPVind = bsRootTree_->BplusPVindex_;
//		 cout << "Bp vtx index "<< BpPVind << endl;

// 		const Vertex &Pvvtx = (*recVtxs)[BpPVind];
//     size_t PVTracks = recVtxs[BpPVind]->nTracks();
//	    cout << "PV tracks " << Pvvtx.nTracks()  << endl;

	  for(reco::Vertex::trackRef_iterator trackvertex=BpPVvtxCosTheta.tracks_begin(); trackvertex!=BpPVvtxCosTheta.tracks_end(); trackvertex++){ //loop over Bp PV tracks

	  reco::TrackRef PVtrk = trackvertex->castTo<reco::TrackRef>();
	  OverlappingTrack = false;

	   for(size_t l = 0; l < recVtxs->size(); l++){ // loop over other PVs
	  
	    if(l != size_t(bsRootTree_->BplusPVindex_) ){
			 
		const Vertex &vertex = (*recVtxs)[l];	    	
		    for(reco::Vertex::trackRef_iterator trkvertex=vertex.tracks_begin(); trkvertex!=vertex.tracks_end(); trkvertex++){ //loop over trks of other PVs
		       reco::TrackRef trk = trkvertex->castTo<reco::TrackRef>();

		       if(trk == PVtrk ){OverlappingTrack = true; }	 	
		     } // end of loop over trks of other PVs
		    }	
	    } // end of loop over other PVs

	//	 if(PVtrk == BpTrkRefs[0] || PVtrk == BpTrkRefs[1] || PVtrk == BpTrkRefs[2]){cout << "PV track overlap with Bp tracks" << endl;}	
		//if trk of Bp PV doesn't belong to tracks of another PV, the track is saved
		  if(OverlappingTrack == false){
			if(PVtrk != BpTrkRefs[0] && PVtrk != BpTrkRefs[1] && PVtrk != BpTrkRefs[2]){
		   const reco::Track *RTv = trackvertex->get();
		   bsRootTree_->BpPVTrkPt_->push_back(RTv->pt());
		   bsRootTree_->BpPVTrkCharge_->push_back(RTv->charge());
		   bsRootTree_->BpPVTrkEta_->push_back(RTv->eta());
		   bsRootTree_->BpPVTrkPhi_->push_back(RTv->phi());
		   trkcouter++;
			}
		  }
	    } //end of loop over Bp PV tracks
		
	} // end of if(BjetIndex!=-1 )

  } //end of if BpPVindex and BpTrkRef vec size

	////////////////////////// JET TAGGING B+ END /////////////////


  ////////////////////// Muon Tagging begin//////////////////////////

  if(verbose_ == true){   std::cout<<"Begin muon tagging"<<std::endl;}
 // cout<<"Begin tagging"<<std::endl;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);

  bool bptagged=false;
  bool bstagged=false;
  bool bdtagged=false;

//  cout << "Number of muons " << allmuons->size() << endl;	
  for(size_t iter=0; iter < allmuons->size(); ++iter){ 
    const pat::Muon & tagMuon = (*allmuons)[iter];
    if( !tagMuon.isPFMuon() ) { continue;}
    if (!tagMuon.innerTrack().isNonnull()) continue;

  bool bpoverlap=false;
  bool bsoverlap=false;
  bool bdoverlap=false;


//	cout << "Tag muon pt " << tagMuon.pt() << endl;
    if (tagMuon.pt()== bsRootTree_->BplusMu1Pt_  || tagMuon.pt()== bsRootTree_->BplusMu2Pt_ || tagMuon.pt() == bsRootTree_->KplusPt_ ){ bpoverlap=true;}
/*
	cout << "Bplus mass " <<bsRootTree_->BplusM_fit_  << endl;
	cout << "Bplus K pt " << bsRootTree_->KplusPt_ <<endl;
	cout << "Bplus Mu1 pt " << bsRootTree_->BplusMu1Pt_ <<endl;
 	cout << "Bplus Mu2 pt " << bsRootTree_->BplusMu2Pt_ <<endl;	
*/  
    if (tagMuon.pt() == bsRootTree_->BdMu1Pt_beffit_ || tagMuon.pt() == bsRootTree_->BdMu2Pt_beffit_ || tagMuon.pt() == bsRootTree_->BdK1Pt_nofit_ || tagMuon.pt() == bsRootTree_->BdK2Pt_nofit_ ){ bdoverlap=true;}


    if(tagMuon.pt()== bsRootTree_->JpsiMu1Pt_alone_  || tagMuon.pt()== bsRootTree_->JpsiMu2Pt_alone_ || tagMuon.pt()==bsRootTree_->K1Pt_nofit_ || tagMuon.pt()==bsRootTree_->K2Pt_nofit_){bsoverlap=true;}

    TrackRef trkTagMuRef = tagMuon.get<TrackRef>();
    TransientTrack trkTagMuTT = (*ttrackBuilder).build(&trkTagMuRef);
    std::pair<bool,Measurement1D> ImpactPar3D;
 //   if(verbose_ == true){   std::cout<<"Analyzing Tag muon, pt:"<< tagMuon.pt()<<"Overlap Bs,Bp,Bd:"<<bsoverlap<<bpoverlap<<bdoverlap<<std::endl;} 

  //  std::cout<<"Analyzing Tag muon, pt:"<< tagMuon.pt()<<" Overlap Bs,Bp,Bd:"<<bsoverlap<<bpoverlap<<bdoverlap<<std::endl;

    if (bpoverlap==false && bptagged==false && BpPVvtxCosTheta.isValid()) {
      bptagged=true;

      bsRootTree_->BpTagIP_ = (IPTools::absoluteImpactParameter3D(trkTagMuTT, BpPVvtxCosTheta)).second.value();  
      bsRootTree_->BpTagPt_ = tagMuon.pt();
      bsRootTree_->BpTagCharge_ = tagMuon.charge();
      bsRootTree_->BpTagEta_ = tagMuon.eta();
      bsRootTree_->BpTagPhi_ = tagMuon.phi();
      bsRootTree_->BpTagP_ = tagMuon.p();
      bsRootTree_->BpTagPF_ = tagMuon.isPFMuon();
      if (isMCstudy_ && tagMuon.genParticlesSize()>0){ bsRootTree_->BpTagGENID_=tagMuon.genParticle()->pdgId();}
    }

    if (bsoverlap==false && bstagged==false && PVvtxCosTheta.isValid()  ) {

      bstagged=true;

      bsRootTree_->BsTagIP_ = (IPTools::absoluteImpactParameter3D(trkTagMuTT, PVvtxCosTheta)).second.value();
      bsRootTree_->BsTagPt_ = tagMuon.pt();
      bsRootTree_->BsTagCharge_ = tagMuon.charge();
      bsRootTree_->BsTagEta_ = tagMuon.eta();
      bsRootTree_->BsTagPhi_ = tagMuon.phi();
      bsRootTree_->BsTagP_ = tagMuon.p();
      bsRootTree_->BsTagPF_ = tagMuon.isPFMuon();
      if (isMCstudy_ && tagMuon.genParticlesSize()>0){ bsRootTree_->BsTagGENID_=tagMuon.genParticle()->pdgId(); }
    }
    if (bdoverlap==false && bdtagged==false && BdPVvtxCosTheta.isValid()) {
      bdtagged=true;

      bsRootTree_->BdTagIP_ = (IPTools::absoluteImpactParameter3D(trkTagMuTT, BdPVvtxCosTheta)).second.value();
      bsRootTree_->BdTagPt_ = tagMuon.pt();
      bsRootTree_->BdTagCharge_ = tagMuon.charge();
      bsRootTree_->BdTagEta_ = tagMuon.eta();
      bsRootTree_->BdTagPhi_ = tagMuon.phi();
      bsRootTree_->BdTagP_ = tagMuon.p();
      bsRootTree_->BdTagPF_ = tagMuon.isPFMuon();
      if (isMCstudy_ && tagMuon.genParticlesSize()>0){ bsRootTree_->BdTagGENID_=tagMuon.genParticle()->pdgId(); }
    }
  }
  
   ////////////////////// Muon Tagging end //////////////////////////

   ////////////////////// Electron Tagging begin ////////////////////


if(verbose_ == true){   std::cout<<"Begin electron tagging"<<std::endl;}
 // cout<<"Begin tagging"<<std::endl;
 // edm::ESHandle<TransientTrackBuilder> ttrackBuilder;
 // iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);

  bool bptaggedEle=false;
  bool bstaggedEle=false;
  bool bdtaggedEle=false;

//  cout << "Number of muons " << allmuons->size() << endl;	
  for(size_t iter=0; iter < allelectrons->size(); ++iter){ 
    const pat::Electron& tagEle = (*allelectrons)[iter];
  //  if( !tagEle.isPFElectron() ) { continue;}
  //  if (!tagEle.innerTrack().isNonnull()) continue;

 // cout << "inside e-loop!" << endl;

  bool bpoverlapEle=false;
  bool bsoverlapEle=false;
  bool bdoverlapEle=false;

  reco::GsfTrackRef trkTagEleRef = tagEle.gsfTrack(); 
  TransientTrack trkTagEleTT = (*ttrackBuilder).build(&trkTagEleRef);
  std::pair<bool,Measurement1D> ImpactPar3D;


    if (bpoverlapEle==false && bptaggedEle==false && BpPVvtxCosTheta.isValid()) {
      bptaggedEle=true;
//	   cout << "Ele found ! " << endl;
		 elecounter++;
      bsRootTree_->BpTagEIP_ = (IPTools::absoluteImpactParameter3D(trkTagEleTT, BpPVvtxCosTheta)).second.value();  

//	cout << "tagEle id " << tagEle.electronID("mvaNonTrigV0") << endl;
//   cout << "tag Ele pt " <<tagEle.pt() << endl;

		bsRootTree_->BpTagEIDNonTrig_ = tagEle.electronID("mvaNonTrigV0");
		bsRootTree_->BpTagEIDTrig_ = tagEle.electronID("mvaTrigV0");
		bsRootTree_->BpTagEIDTrigNoIP_ = tagEle.electronID("mvaTrigNoIPV0");
      bsRootTree_->BpTagEPt_ = tagEle.pt();
      bsRootTree_->BpTagECharge_ = tagEle.charge();
      bsRootTree_->BpTagEEta_ = tagEle.eta();
      bsRootTree_->BpTagEPhi_ = tagEle.phi();
      bsRootTree_->BpTagEP_ = tagEle.p();

      if (isMCstudy_ && tagEle.genParticlesSize()>0){ bsRootTree_->BpTagEGENID_=tagEle.genParticle()->pdgId();}

    }


    if (bsoverlapEle==false && bstaggedEle==false && PVvtxCosTheta.isValid()  ) {

      bstaggedEle=true;
//	cout << "Tagged e- found!" << endl;
      bsRootTree_->BsTagEIP_ = (IPTools::absoluteImpactParameter3D(trkTagEleTT, PVvtxCosTheta)).second.value();

		bsRootTree_->BsTagEIDNonTrig_ = tagEle.electronID("mvaNonTrigV0");
		bsRootTree_->BsTagEIDTrig_ = tagEle.electronID("mvaTrigV0");
		bsRootTree_->BsTagEIDTrigNoIP_ = tagEle.electronID("mvaTrigNoIPV0");
      bsRootTree_->BsTagEPt_ = tagEle.pt();
      bsRootTree_->BsTagECharge_ = tagEle.charge();
      bsRootTree_->BsTagEEta_ = tagEle.eta();
      bsRootTree_->BsTagEPhi_ = tagEle.phi();
      bsRootTree_->BsTagEP_ = tagEle.p();
 //     bsRootTree_->BsTagPF_ = tagMuon.isPFMuon();
      if (isMCstudy_ && tagEle.genParticlesSize()>0){ bsRootTree_->BsTagEGENID_=tagEle.genParticle()->pdgId(); 
	//	 cout << "tag e id " << tagEle.genParticle()->pdgId() << endl; 
}

    }


    if (bdoverlapEle==false && bdtaggedEle==false && BdPVvtxCosTheta.isValid()) {
      bdtaggedEle=true;

      bsRootTree_->BdTagEIP_ = (IPTools::absoluteImpactParameter3D(trkTagEleTT, BdPVvtxCosTheta)).second.value();
		bsRootTree_->BdTagEIDNonTrig_ = tagEle.electronID("mvaNonTrigV0");
		bsRootTree_->BdTagEIDTrig_ = tagEle.electronID("mvaTrigV0");
		bsRootTree_->BdTagEIDTrigNoIP_ = tagEle.electronID("mvaTrigNoIPV0");
      bsRootTree_->BdTagEPt_ = tagEle.pt();
      bsRootTree_->BdTagECharge_ = tagEle.charge();
      bsRootTree_->BdTagEEta_ = tagEle.eta();
      bsRootTree_->BdTagEPhi_ = tagEle.phi();
      bsRootTree_->BdTagEP_ = tagEle.p();
      //bsRootTree_->BdTagPF_ = tagEle.isPFMuon();
      if (isMCstudy_ && tagEle.genParticlesSize()>0){ bsRootTree_->BdTagEGENID_=tagEle.genParticle()->pdgId(); }

    }
  }

  ////////////////////// Electron Tagging end ///////////////////////

//cout << "Bplus charge " << bsRootTree_->BplusCharge_ << endl; 
//cout << "Bp matched " <<  bsRootTree_->isMatchedBplus_ << endl;	
  //if ( MinBVtxHyp1>0 || minVtxP>0 )
  //if ((isMCstudy_ && (bsRootTree_->ChannelID_==1 || bsRootTree_->ChannelID_==2 || bsRootTree_->ChannelID_==3 ||  bsRootTree_->ChannelID_==4)) || !isMCstudy_)
    if(verbose_ == true){   std::cout<<"Fill the TTree"<<std::endl;}
    bsRootTree_->fill();
}








GlobalVector BsToJpsiPhiAnalysis::flightDirection(const reco::Vertex &pv, reco::Vertex &sv){
  GlobalVector res(sv.position().X() - pv.position().X(),
                    sv.position().Y() - pv.position().Y(),
                    sv.position().Z() - pv.position().Z());
  return res;
}

void BsToJpsiPhiAnalysis::fillMCInfo( edm::Handle<GenParticleCollection> & genParticles){

  int iNumberOfBdecays = 0;
  Int_t SVindex =0;
  // this is a list of all the PDG ids of B mesons
  std::set<int> listOfBmesonIds;
  listOfBmesonIds.insert(511 );   // Bd
  listOfBmesonIds.insert(521 );   // B+
  listOfBmesonIds.insert(10511 );    // B_0*0
  listOfBmesonIds.insert(10521 );    // B_0*+
  listOfBmesonIds.insert(513 );   // B*d
  listOfBmesonIds.insert(523 );   // B*d+
  listOfBmesonIds.insert(10513 );   // B1(L)0
  listOfBmesonIds.insert(10523 );   // B1(L)+
  listOfBmesonIds.insert(20513 );   // B1(H)0
  listOfBmesonIds.insert(20523 );   // B1(H)+
  listOfBmesonIds.insert(515 );    // B2*_0
  listOfBmesonIds.insert(525 );    // B2*_+
  listOfBmesonIds.insert(531 );   // Bs
  listOfBmesonIds.insert(10531 );    // B_s0*_0
  listOfBmesonIds.insert(533 );   // B*s
  listOfBmesonIds.insert(10533 );   // Bs1(L)0
  listOfBmesonIds.insert(20533 );   // Bs1(H)0
  listOfBmesonIds.insert(535 );    // Bs2*_0
  listOfBmesonIds.insert(541 );   // Bc+
  listOfBmesonIds.insert(10541 );   // B*c0+
  listOfBmesonIds.insert(543 );   // B*c+
  listOfBmesonIds.insert(10543 );   // Bc1(L)+
  listOfBmesonIds.insert(20543 );   // Bc1(H)+
  listOfBmesonIds.insert(545 );    // Bc2*_0
  
  listOfBmesonIds.insert(551 );   // etab(1S)
  listOfBmesonIds.insert(10551 );   // chib(1P)
  listOfBmesonIds.insert(100551 );   // etab(2S)
  listOfBmesonIds.insert(110551 );   // chib(2P)
  listOfBmesonIds.insert(200551 );   // etab(3S)
  listOfBmesonIds.insert(210551 );   // chib(3P)
  listOfBmesonIds.insert(553 );   // upsilon(1S)
  listOfBmesonIds.insert(10553 );   // hb(1P)
  listOfBmesonIds.insert(20553 );   // chib1(1P)
  listOfBmesonIds.insert(30553 );   // upsilon1(1D)
  listOfBmesonIds.insert(100553 );   // upsilon(2S)
  listOfBmesonIds.insert(110553 );   // hb(2P)
  listOfBmesonIds.insert(120553 );   // chib1(2P)
  listOfBmesonIds.insert(130553 );   // upsilon1(2D)
  listOfBmesonIds.insert(200553 );   // upsilon(3S)
  listOfBmesonIds.insert(210553 );   // hb(3P)
  listOfBmesonIds.insert(220553 );   // chib1(3P)
  listOfBmesonIds.insert(300553 );   // upsilon(4S)
  listOfBmesonIds.insert(9000553 );   // upsilon(10860)
  listOfBmesonIds.insert(9010553 );   // upsilon(11020)
  listOfBmesonIds.insert(555 );   // chib2(1P)
  listOfBmesonIds.insert(10555 );   // etab2(1D)
  listOfBmesonIds.insert(20555 );   // upsilon2(1D)
  listOfBmesonIds.insert(100555 );   // chib2(2P)
  listOfBmesonIds.insert(110555 );   // etab2(2D)
  listOfBmesonIds.insert(120555 );   // upsilon2(2D)
  listOfBmesonIds.insert(200555 );   // chib2(3P)
  listOfBmesonIds.insert(557 );   // upsilon3(1D)
  listOfBmesonIds.insert(100557 );   // upsilon3(2D)
  
  listOfBmesonIds.insert(5122 );   // lambda_b0
  listOfBmesonIds.insert(5112 );   // sigma_b-
  listOfBmesonIds.insert(5212 );   // sigma_b0
  listOfBmesonIds.insert(5222 );   // sigma_b+
  listOfBmesonIds.insert(5114 );   // sigma*_b-
  listOfBmesonIds.insert(5214 );   // sigma*_b0
  listOfBmesonIds.insert(5224 );   // sigma*_b+
  listOfBmesonIds.insert(5132 );   // Xi_b-
  listOfBmesonIds.insert(5232 );   // Xi_b0
  listOfBmesonIds.insert(5312 );   // Xi'_b-
  listOfBmesonIds.insert(5322 );   // Xi'_b0
  listOfBmesonIds.insert(5314 );   // Xi*_b-
  listOfBmesonIds.insert(5324 );   // Xi*_b0
  listOfBmesonIds.insert(5332 );   // Omega_b-
  listOfBmesonIds.insert(5334 );   // Omega*_b-
  listOfBmesonIds.insert(5142 );   // Xi_bc0
  listOfBmesonIds.insert(5242 );   // Xi_bc+
  listOfBmesonIds.insert(5412 );   // Xi'_bc0
  listOfBmesonIds.insert(5422 );   // Xi'_bc+
  listOfBmesonIds.insert(5414 );   // Xi*_bc0
  listOfBmesonIds.insert(5424 );   // Xi*_bc+
  listOfBmesonIds.insert(5342 );   // Omega_bc0
  listOfBmesonIds.insert(5432 );   // Omega'_bc0
  listOfBmesonIds.insert(5434 );   // Omega*_bc0
  listOfBmesonIds.insert(5442 );   // Omega_bcc+
  listOfBmesonIds.insert(5444 );   // Omega*_bcc+
  listOfBmesonIds.insert(5512 );   // Xi_bb-
  listOfBmesonIds.insert(5522 );   // Xi_bb0
  listOfBmesonIds.insert(5514 );   // Xi*_bb-
  listOfBmesonIds.insert(5524 );   // Xi*_bb0
  listOfBmesonIds.insert(5532 );   // Omega_bb-
  listOfBmesonIds.insert(5524 );   // Omega*_bb-
  listOfBmesonIds.insert(5542 );   // Omega_bbc0
  listOfBmesonIds.insert(5544 );   // Omega*_bbc0
  listOfBmesonIds.insert(554 );   // Omega_bbb-

  const Candidate * Jpsi = 0;
  const Candidate * Phi = 0;
  const Candidate * mup = 0;
  const Candidate * mum = 0;
  const Candidate * Kp = 0;
  const Candidate * Km = 0;
  
  // loop over all particles
  cout<<"NewEvent-------"<<endl;
  for( size_t i = 0; i < genParticles->size(); ++ i ) {
    const GenParticle & genBsCand = (*genParticles)[ i ];
    int MC_particleID=genBsCand.pdgId();
    int absMC_particleID = abs(MC_particleID);


/////////////////////// Bp Matching //////////////////
        bool isPosMu=0, isNegMu=0, isJpsi=0, isKplus=0, isNegPi=0, isPosPi=0, isKstar = 0, isNeutralPi=0; 

      if( abs(MC_particleID) == 521 ){

		
	int numBplusDaughters = genBsCand.numberOfDaughters();
	int jpsiIndex = -5;
	int kstarIndex = -5;	
   int kpmIndex = -5;
	//checking of decay channels B+ -> Jpsi(mu+,mu-) K+ and B+ -> Jpsi K*(K+, pi0)  
	if(numBplusDaughters == 2){ //check if Bplus has 2 daughters
          for(int k = 0; k < numBplusDaughters; k++ ){ //check if the daughter's are J/psi and (K+ or K*)
	    if( genBsCand.daughter(k)->pdgId() == 443){ 
	      isJpsi = 1;
	      jpsiIndex = k;
     	    }  
	    else if( abs(genBsCand.daughter(k)->pdgId()) == 321 ){ isKplus = 1; kpmIndex = k; }
	    else if( genBsCand.daughter(k)->pdgId() == 323 ){ 
	      isKstar = 1; 
	      kstarIndex = k;	
	    }
	  }

	  if( isJpsi == 1 && ( isKplus == 1 || isKstar == 1 ) ){	
	    //check if J/psi has two muon daughters 
	    const Candidate *JpsiCand = genBsCand.daughter(jpsiIndex);
      	    for(unsigned int j = 0; j < JpsiCand->numberOfDaughters(); j++ ){ //check if Jpsi daughters are mu+ and mu-
	      const Candidate * Jpsidau = JpsiCand->daughter(j);
	      if(Jpsidau->pdgId() == 13 ){ isNegMu = 1; }
              else if (Jpsidau->pdgId() == -13 ){ isPosMu = 1; }
	    } 
				//if(isPosMu != 1 &&  isNegMu!= 1) continue;
         if( isPosMu == 1 && isNegMu == 1 && isKplus == 1 ){
				if(kpmIndex >=0) {
					const Candidate *KaonCand = genBsCand.daughter(kpmIndex);
					if( KaonCand->pdgId() == 321) {bsRootTree_->BplusDecayChannel_ = 1; cout << "B+ DecayChannel = 1" <<endl; }
					if( KaonCand->pdgId() == -321) {bsRootTree_->BplusDecayChannel_ = -1; cout << "B+ DecayChannel = -1" <<endl;}
				}
	     	// channel 1 = B+ -> Jpsi(mu+,mu-) K+
       	    }						
	   isKplus = 0; 
	   if(isKstar == 1){
	     const Candidate *Kstar = genBsCand.daughter(kstarIndex);
  	     if (Kstar->numberOfDaughters() == 2){							
  	       for(unsigned int j = 0; j < Kstar->numberOfDaughters(); j++ ){ //check if Kstar daughters are Kplus and neutral pion
	         const Candidate * Kstardau = Kstar->daughter(j);
		 if(Kstardau->pdgId() == 321 ){ isKplus = 1; }
		 else if (Kstardau->pdgId() == 111 ){ isNeutralPi = 1; }
	       } 

	       if( isNegMu == 1 && isPosMu == 1 && isNeutralPi == 1 && isKplus == 1 ){
	         bsRootTree_->BplusDecayChannel_ = 3; // channel 3 = B+ -> Jpsi K*(K+ pi0) 
	       }
	     }
	   }
         }
       }
     // end of B+ -> Jpsi(mu+,mu-) K+ channel search					

	//checking of decay channel B+ -> Jpsi K+ pi+ pi- 	
       if(numBplusDaughters == 4){	
         for(int j = 0; j < numBplusDaughters; j++ ){ //check if Bplus has 4 daughters
	   if(genBsCand.daughter(j)->pdgId() == 443 ){ isJpsi = 1; } //check if the daughter's are J/psi, pi+, pi- and K+
           else if(genBsCand.daughter(j)->pdgId() == 321 ){ isKplus = 1; }
	   else if(genBsCand.daughter(j)->pdgId() == 211 ){ isPosPi = 1; }
	   else if(genBsCand.daughter(j)->pdgId() == -211 ){ isNegPi = 1; }  			
	 }  

	 if(isJpsi == 1 && isKplus == 1 && isPosPi == 1 && isNegPi == 1 ){
	   bsRootTree_->BplusDecayChannel_ = 2;	
         } 	
       } // end of  B+ -> Jpsi K+ pi+ pi- channel search	
     } // end of Bplus search      
 /////////////////////// Bp Matching //////////////////


    // if this particle id is in the list (i.e. if it is a B meson)
    if( listOfBmesonIds.find( absMC_particleID ) != listOfBmesonIds.end()){
      
      // check if this particle has no daughter which is a B meson (cascade decays)
      // loop over daughters
      bool hasBDaughter=0;
      int numBsDaughters = genBsCand.numberOfDaughters();      
      for(int idau=0; idau < numBsDaughters; idau++) 
	if( listOfBmesonIds.find( abs(genBsCand.daughter(idau)->pdgId())) != listOfBmesonIds.end() ) hasBDaughter=1;
/*
      if( abs(MC_particleID)==511) {
        cout<<"Mother of "<<MC_particleID<<":"<<genBsCand.mother(0)->pdgId()<<endl;
      }
      if( abs(MC_particleID)==531) {
        cout<<"Mother of "<<MC_particleID<<":"<<genBsCand.mother(0)->pdgId()<<endl;
      }
 */  
      if( abs(genBsCand.mother(0)->pdgId()) != 511 && abs(MC_particleID)==511) {
        const Candidate * genBsCand2 =genBsCand.daughter(0);
        bool isBdKstarKmPip=false;
        bool isBdJpsiMC=false;
        bool isBdKstarKpPim=false;
        double bdsvx=0;
        double bdsvy=0;
        double bdsvz=0;
        double bdpvx=genBsCand.vx();
        double bdpvy=genBsCand.vy();
        double bdpvz=genBsCand.vz();
        double bdmomx=genBsCand.px();
        double bdmomy=genBsCand.py();
        double bdmomz=genBsCand.pz();

        const Candidate *bdmuplus=0;
        const Candidate *bdmuminus=0;
        const Candidate *bdkplus=0;
        const Candidate *bdkminus=0;
        int BdHasMixed=-1;

        if (MC_particleID==(-1)*genBsCand.daughter(0)->pdgId()) {
          //bsRootTree_->BdEndFlavour_=(int)(genBsCand.daughter(0)->pdgId()/511);
          BdHasMixed=1;
          if (genBsCand2->numberOfDaughters()==2) {
            if (abs(genBsCand2->daughter(0)->pdgId())==443 && abs(genBsCand2->daughter(1)->pdgId())==313) {
              if (genBsCand2->daughter(0)->numberOfDaughters()>1) 
                if (abs(genBsCand2->daughter(0)->daughter(0)->pdgId())==13 || abs(genBsCand2->daughter(0)->daughter(1)->pdgId())==13) {
                  isBdJpsiMC=true;
                  if (genBsCand2->daughter(0)->daughter(0)->charge()==1) {bdmuplus=genBsCand2->daughter(0)->daughter(0); bdmuminus=genBsCand2->daughter(0)->daughter(1); }
                  if (genBsCand2->daughter(0)->daughter(1)->charge()==1) {bdmuplus=genBsCand2->daughter(0)->daughter(1); bdmuminus=genBsCand2->daughter(0)->daughter(0); }
                }
              if (genBsCand2->daughter(1)->numberOfDaughters()==2) {
                if (genBsCand2->daughter(1)->daughter(0)->pdgId()==-321 && genBsCand2->daughter(1)->daughter(1)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand2->daughter(1)->daughter(0)->pdgId()==321 && genBsCand2->daughter(1)->daughter(1)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand2->daughter(1)->daughter(1)->pdgId()==-321 && genBsCand2->daughter(1)->daughter(0)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand2->daughter(1)->daughter(1)->pdgId()==321 && genBsCand2->daughter(1)->daughter(0)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand2->daughter(1)->daughter(0)->pdgId()<0) {bdkminus=genBsCand2->daughter(1)->daughter(0); bdkplus=genBsCand2->daughter(1)->daughter(1);}
                if (genBsCand2->daughter(1)->daughter(1)->pdgId()<0) {bdkminus=genBsCand2->daughter(1)->daughter(1); bdkplus=genBsCand2->daughter(1)->daughter(0);}
              }
            }
            if (abs(genBsCand2->daughter(1)->pdgId())==443 && abs(genBsCand2->daughter(0)->pdgId())==313) {
              if (genBsCand2->daughter(1)->numberOfDaughters()>1)
                if (abs(genBsCand2->daughter(1)->daughter(0)->pdgId())==13 || abs(genBsCand2->daughter(1)->daughter(1)->pdgId())==13) {
                  isBdJpsiMC=true;
                  if (genBsCand2->daughter(1)->daughter(0)->charge()==1) {bdmuplus=genBsCand2->daughter(1)->daughter(0); bdmuminus=genBsCand2->daughter(1)->daughter(1); }
                  if (genBsCand2->daughter(1)->daughter(1)->charge()==1) {bdmuplus=genBsCand2->daughter(1)->daughter(1); bdmuminus=genBsCand2->daughter(1)->daughter(0); }
                }
              if (genBsCand2->daughter(0)->numberOfDaughters()==2) {
                if (genBsCand2->daughter(0)->daughter(0)->pdgId()==-321 && genBsCand2->daughter(0)->daughter(1)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand2->daughter(0)->daughter(0)->pdgId()==321 && genBsCand2->daughter(0)->daughter(1)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand2->daughter(0)->daughter(1)->pdgId()==-321 && genBsCand2->daughter(0)->daughter(0)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand2->daughter(0)->daughter(1)->pdgId()==321 && genBsCand2->daughter(0)->daughter(0)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand2->daughter(0)->daughter(0)->pdgId()<0) {bdkminus=genBsCand2->daughter(0)->daughter(0); bdkplus=genBsCand2->daughter(0)->daughter(1);}
                if (genBsCand2->daughter(0)->daughter(1)->pdgId()<0) {bdkminus=genBsCand2->daughter(0)->daughter(1); bdkplus=genBsCand2->daughter(0)->daughter(0);}
              }
            }
            bdsvx=genBsCand2->daughter(0)->vx();
            bdsvy=genBsCand2->daughter(0)->vy();
            bdsvz=genBsCand2->daughter(0)->vz();
          }
        } else {
          //bsRootTree_->BdEndFlavour_=(int)(MC_particleID/511);
          if (genBsCand.numberOfDaughters()==2) {
            if (abs(genBsCand.daughter(0)->pdgId())==443 && abs(genBsCand.daughter(1)->pdgId())==313) {
              if (genBsCand.daughter(0)->numberOfDaughters()>1)
                if (abs(genBsCand.daughter(0)->daughter(0)->pdgId())==13 || abs(genBsCand.daughter(0)->daughter(1)->pdgId())==13) {
                  isBdJpsiMC=true;
                  if (genBsCand.daughter(0)->daughter(0)->charge()==1) {bdmuplus=genBsCand.daughter(0)->daughter(0); bdmuminus=genBsCand.daughter(0)->daughter(1); }
                  if (genBsCand.daughter(0)->daughter(1)->charge()==1) {bdmuplus=genBsCand.daughter(0)->daughter(1); bdmuminus=genBsCand.daughter(0)->daughter(0); }
                }
              if (genBsCand.daughter(1)->numberOfDaughters()==2) {
                if (genBsCand.daughter(1)->daughter(0)->pdgId()==-321 && genBsCand.daughter(1)->daughter(1)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand.daughter(1)->daughter(0)->pdgId()==321 && genBsCand.daughter(1)->daughter(1)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand.daughter(1)->daughter(1)->pdgId()==-321 && genBsCand.daughter(1)->daughter(0)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand.daughter(1)->daughter(1)->pdgId()==321 && genBsCand.daughter(1)->daughter(0)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand.daughter(1)->daughter(0)->pdgId()<0) {bdkminus=genBsCand.daughter(1)->daughter(0); bdkplus=genBsCand.daughter(1)->daughter(1);}
                if (genBsCand.daughter(1)->daughter(1)->pdgId()<0) {bdkminus=genBsCand.daughter(1)->daughter(1); bdkplus=genBsCand.daughter(1)->daughter(0);}
              }
            }
            if (abs(genBsCand.daughter(1)->pdgId())==443 && abs(genBsCand.daughter(0)->pdgId())==313) {
              if (genBsCand.daughter(1)->numberOfDaughters()>1)
                if (abs(genBsCand.daughter(1)->daughter(0)->pdgId())==13 || abs(genBsCand.daughter(1)->daughter(1)->pdgId())==13) {
                  isBdJpsiMC=true;
                  if (genBsCand.daughter(1)->daughter(0)->charge()==1) {bdmuplus=genBsCand.daughter(1)->daughter(0); bdmuminus=genBsCand.daughter(1)->daughter(1); }
                  if (genBsCand.daughter(1)->daughter(1)->charge()==1) {bdmuplus=genBsCand.daughter(1)->daughter(1); bdmuminus=genBsCand.daughter(1)->daughter(0); }
                }
              if (genBsCand.daughter(0)->numberOfDaughters()==2) {
                if (genBsCand.daughter(0)->daughter(0)->pdgId()==-321 && genBsCand.daughter(0)->daughter(1)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand.daughter(0)->daughter(0)->pdgId()==321 && genBsCand.daughter(0)->daughter(1)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand.daughter(0)->daughter(1)->pdgId()==-321 && genBsCand.daughter(0)->daughter(0)->pdgId()==211) isBdKstarKmPip=true;
                if (genBsCand.daughter(0)->daughter(1)->pdgId()==321 && genBsCand.daughter(0)->daughter(0)->pdgId()==-211) isBdKstarKpPim=true;
                if (genBsCand.daughter(0)->daughter(0)->pdgId()<0) {bdkminus=genBsCand.daughter(0)->daughter(0); bdkplus=genBsCand.daughter(0)->daughter(1);}
                if (genBsCand.daughter(0)->daughter(1)->pdgId()<0) {bdkminus=genBsCand.daughter(0)->daughter(1); bdkplus=genBsCand.daughter(0)->daughter(0);}
              }
            }
            bdsvx=genBsCand.daughter(0)->vx();
            bdsvy=genBsCand.daughter(0)->vy();
            bdsvz=genBsCand.daughter(0)->vz();
          } 
        }
        if (isBdJpsiMC==true && isBdKstarKmPip==true) bsRootTree_->BdChannelID_=1;
        if (isBdJpsiMC==true && isBdKstarKpPim==true) bsRootTree_->BdChannelID_=2;
        if (isBdJpsiMC==true && ( isBdKstarKmPip==true || isBdKstarKpPim==true )) {
          // calculate gen ctau 2D
          //double Lxy2D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy);
          double Lxy2D = sqrt(pow(bdsvx-bdpvx,2)+pow(bdsvy-bdpvy,2));
          bsRootTree_->BdCt2DMC_ = Lxy2D * 5.2795/sqrt(bdmomx*bdmomx+bdmomy*bdmomy);
          // calculate gen ctau 3D
          //double Lxy3D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy+(bssvz-bspvz)*bsmomz);
          double Lxy3D = sqrt(pow(bdsvx-bdpvx,2)+pow(bdsvy-bdpvy,2)+pow(bdsvz-bdpvz,2));
          bsRootTree_->BdCt3DMC_ = Lxy3D * 5.2795/sqrt(bdmomx*bdmomx+bdmomy*bdmomy+bdmomz*bdmomz);

          TLorentzVector pmuplus;
          TLorentzVector pmuminus;	      
	  TLorentzVector pkplus;
	  TLorentzVector pkminus;
          pmuplus.SetXYZM(bdmuplus->px(),bdmuplus->py(),bdmuplus->pz(),bdmuplus->mass());
          pmuminus.SetXYZM(bdmuminus->px(),bdmuminus->py(),bdmuminus->pz(),bdmuminus->mass());
          pkplus.SetXYZM(bdkplus->px(),bdkplus->py(),bdkplus->pz(),bdkplus->mass());
	  pkminus.SetXYZM(bdkminus->px(),bdkminus->py(),bdkminus->pz(),bdkminus->mass());
	     
	  // boosting in JPsi restframe
	  TLorentzVector pjpsi;
	  pjpsi = pmuplus + pmuminus;
	  TLorentzVector pphi;
	  pphi = pkplus + pkminus;
	      
	  // the betas for the boost
	  TVector3 p3_JPsi;
	  p3_JPsi = pjpsi.Vect();
	  p3_JPsi *= -1./pjpsi.E();
	    
	  // the boost matrix
	  TLorentzRotation boost_jpsi(p3_JPsi);
	  TLorentzVector p_JPsi_JPsi;
	  p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);
	     
	  // the different momenta in the new frame
	  TLorentzVector p_JPsi_muplus;
	  TLorentzVector p_JPsi_Kplus;
	  TLorentzVector p_JPsi_phi;
	      
	  p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
	  p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);                                                                                
	  p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);	      
	      
	  // the 3-momenta
	  TVector3 p3_JPsi_muplus;
	  p3_JPsi_muplus = p_JPsi_muplus.Vect();
	  TVector3 p3_JPsi_Kplus;
	  p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	  TVector3 p3_JPsi_phi;
	  p3_JPsi_phi = p_JPsi_phi.Vect();
	      
	  // coordinate system
	  TVector3 x,y,z;
	  x = p3_JPsi_phi.Unit();
	  y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	  y = y.Unit();
	  z = x.Cross(y);
	   
          double Bdangle_costheta=p3_JPsi_muplus.Unit() * z;   
	  bsRootTree_->BdcosthetaMC_ = Bdangle_costheta;
           
	  double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	  double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - Bdangle_costheta*Bdangle_costheta);
	  double Bdangle_phi = TMath::ACos(cos_phi);
	  if (sin_phi < 0){
	    Bdangle_phi =  -Bdangle_phi;
          }
	  bsRootTree_->BdphiMC_ = Bdangle_phi;
	      
	  // boosting in phi restframe
	  TVector3 p3_phi;
	  p3_phi = pphi.Vect();
	  p3_phi *= -1./pphi.E();
	      
	  // the boost matrix
	  TLorentzRotation boost_phi(p3_phi);
	  TLorentzVector p_phi_phi;
	  p_phi_phi = boost_phi.VectorMultiplication(pphi);
	      
	  // the different momenta in the new frame
	  TLorentzVector p_phi_Kplus;
	  TLorentzVector p_phi_JPsi;
	  TLorentzVector p_phi_Bs;
	      
	  p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	  p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	      
	  // the 3-momenta
	  TVector3 p3_phi_Kplus;
	  p3_phi_Kplus = p_phi_Kplus.Vect();
	  TVector3 p3_phi_JPsi;
	  p3_phi_JPsi = p_phi_JPsi.Vect();
	  bsRootTree_->BdcospsiMC_ = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
        }	
        bsRootTree_->BdIniFlavour_=(int)(MC_particleID/511);
      }
/*
      bool isPosMu=0, isNegMu=0, isJpsi=0, isKplus=0, isNegPi=0, isPosPi=0, isKstar = 0, isNeutralPi=0; 

      if( abs(MC_particleID) == 521 ){

	cout << "Bplus MC matching! " << endl; 	
	int numBplusDaughters = genBsCand.numberOfDaughters();
	int jpsiIndex = -5;
	int kstarIndex = -5;	
   int kpmIndex = -5;
	//checking of decay channels B+ -> Jpsi(mu+,mu-) K+ and B+ -> Jpsi K*(K+, pi0)  
	if(numBplusDaughters == 2){ //check if Bplus has 2 daughters
          for(int k = 0; k < numBplusDaughters; k++ ){ //check if the daughter's are J/psi and (K+ or K*)
	    if( genBsCand.daughter(k)->pdgId() == 443){ 
	      isJpsi = 1;
	      jpsiIndex = k;
     	    }  
	    else if( abs(genBsCand.daughter(k)->pdgId()) == 321 ){ isKplus = 1; kpmIndex = k; }
	    else if( genBsCand.daughter(k)->pdgId() == 323 ){ 
	      isKstar = 1; 
	      kstarIndex = k;	
	    }
	  }

	  if( isJpsi == 1 && ( isKplus == 1 || isKstar == 1 ) ){	
	    //check if J/psi has two muon daughters 
	    const Candidate *JpsiCand = genBsCand.daughter(jpsiIndex);
      	    for(unsigned int j = 0; j < JpsiCand->numberOfDaughters(); j++ ){ //check if Jpsi daughters are mu+ and mu-
	      const Candidate * Jpsidau = JpsiCand->daughter(j);
	      if(Jpsidau->pdgId() == 13 ){ isNegMu = 1; }
              else if (Jpsidau->pdgId() == -13 ){ isPosMu = 1; }
	    } 
				//if(isPosMu != 1 &&  isNegMu!= 1) continue;
            if( isPosMu == 1 && isNegMu == 1 && isKplus == 1 ){
	      bsRootTree_->BplusDecayChannel_ = 1;	// channel 1 = B+ -> Jpsi(mu+,mu-) K+
			cout << "Jpsi Cand id " << genBsCand.daughter(jpsiIndex) << endl;
			cout << "Kaon Cand id " << genBsCand.daughter(kpmIndex) << endl;
       	    }						
	   isKplus = 0; 
	   if(isKstar == 1){
	     const Candidate *Kstar = genBsCand.daughter(kstarIndex);
  	     if (Kstar->numberOfDaughters() == 2){							
  	       for(unsigned int j = 0; j < Kstar->numberOfDaughters(); j++ ){ //check if Kstar daughters are Kplus and neutral pion
	         const Candidate * Kstardau = Kstar->daughter(j);
		 if(Kstardau->pdgId() == 321 ){ isKplus = 1; }
		 else if (Kstardau->pdgId() == 111 ){ isNeutralPi = 1; }
	       } 

	       if( isNegMu == 1 && isPosMu == 1 && isNeutralPi == 1 && isKplus == 1 ){
	         bsRootTree_->BplusDecayChannel_ = 3; // channel 3 = B+ -> Jpsi K*(K+ pi0) 
	       }
	     }
	   }
         }
       }
     // end of B+ -> Jpsi(mu+,mu-) K+ channel search					

	//checking of decay channel B+ -> Jpsi K+ pi+ pi- 	
       if(numBplusDaughters == 4){	
         for(int j = 0; j < numBplusDaughters; j++ ){ //check if Bplus has 4 daughters
	   if(genBsCand.daughter(j)->pdgId() == 443 ){ isJpsi = 1; } //check if the daughter's are J/psi, pi+, pi- and K+
           else if(genBsCand.daughter(j)->pdgId() == 321 ){ isKplus = 1; }
	   else if(genBsCand.daughter(j)->pdgId() == 211 ){ isPosPi = 1; }
	   else if(genBsCand.daughter(j)->pdgId() == -211 ){ isNegPi = 1; }  			
	 }  

	 if(isJpsi == 1 && isKplus == 1 && isPosPi == 1 && isNegPi == 1 ){
	   bsRootTree_->BplusDecayChannel_ = 2;	
         } 	
       } // end of  B+ -> Jpsi K+ pi+ pi- channel search	
     } // end of Bplus search      
*/
 
      if( abs(genBsCand.mother(0)->pdgId()) != 531 && abs(MC_particleID)==531) {
        vector<int> figlie;
        vector<int> request;
        vector<int> requestDauPhi;
        vector<int> requestDauPsi;
        vector<int> requestDauPsi2;
        vector<int> requestDauPsi3;
        vector<int> requestDauPsi4;
        request.push_back(333);
        request.push_back(443);
        requestDauPhi.push_back(321);
        requestDauPhi.push_back(321);
        requestDauPsi.push_back(13);
        requestDauPsi.push_back(13);
        requestDauPsi2.push_back(13);
        requestDauPsi2.push_back(13);
        requestDauPsi2.push_back(22);
        requestDauPsi3.push_back(13);
        requestDauPsi3.push_back(13);
        requestDauPsi3.push_back(22);
        requestDauPsi3.push_back(22);
        requestDauPsi4.push_back(13);
        requestDauPsi4.push_back(13);
        requestDauPsi4.push_back(22);
        requestDauPsi4.push_back(22);
        requestDauPsi4.push_back(22);
        bool MCpsi=false;
        bool MCphi=false;
        bool MCpsi2=false;
        bool MCpsi3=false;
        bool MCpsi4=false;
        bool MCbs=false;

        double bssvx=0;
        double bssvy=0;
        double bssvz=0;
        double bspvx=genBsCand.vx();
        double bspvy=genBsCand.vy();
        double bspvz=genBsCand.vz();
        double bsmomx=genBsCand.px();
        double bsmomy=genBsCand.py();
        double bsmomz=genBsCand.pz();
	TLorentzVector pmuplus;
	TLorentzVector pmuminus;	      
	TLorentzVector pkplus;
	TLorentzVector pkminus;
        int BsHasMixed=-1; 
	

        if (MC_particleID==(-1)*genBsCand.daughter(0)->pdgId()) { 
          BsHasMixed=1;
          //bsRootTree_->BsEndFlavour_=(int)(genBsCand.daughter(0)->pdgId()/531);
          //cout<<"Mixing: "<<bsRootTree_->BsEndFlavour_<<endl;
          const Candidate * genBsCand2 =genBsCand.daughter(0);
          cout<<"Mixed - First vertex:"<<genBsCand.vx()<<" "<<genBsCand.vy()<<" "<<genBsCand.vz()<<" Second Vertex:"<<genBsCand2->vx()<<" "<<genBsCand2->vy()<<" "<<genBsCand2->vz()<<endl;		
          int numBsDau = genBsCand2->numberOfDaughters();
          for (int ghepensimi=0; ghepensimi<numBsDau; ghepensimi++) {
            int numBsDauDau = genBsCand2->daughter(ghepensimi)->numberOfDaughters();
            figlie.push_back(genBsCand2->daughter(ghepensimi)->pdgId());
            //cout<<"Figlie di Bs:"<<genBsCand2->daughter(ghepensimi)->pdgId()<<endl;
            vector<int> figlieFiglie;
            for (int ghepensimi2=0; ghepensimi2<numBsDauDau; ghepensimi2++) {
              cout<<"             >"<<genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()<<endl;
              figlieFiglie.push_back(abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()));
              if (abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->charge()==-1)  pmuminus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (abs(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->charge()==1)  pmuplus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
	      if (genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==321)  pkplus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
	      if (genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==-321)  pkminus.SetXYZM(genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand2->daughter(ghepensimi)->daughter(ghepensimi2)->mass());
            }
            sort(figlieFiglie.begin(),figlieFiglie.end());
            if (figlieFiglie.size()==requestDauPhi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPhi.begin())) MCphi=true;
            if (figlieFiglie.size()==requestDauPsi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi.begin())) MCpsi=true;
            if (figlieFiglie.size()==requestDauPsi2.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi2.begin())) MCpsi2=true;
            if (figlieFiglie.size()==requestDauPsi3.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi3.begin())) MCpsi3=true;
            if (figlieFiglie.size()==requestDauPsi4.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi4.begin())) MCpsi4=true;
            if (abs(genBsCand2->daughter(ghepensimi)->pdgId())==443) {
              cout<<"Last Vertex:"<<genBsCand2->daughter(ghepensimi)->vx()<<" "<<genBsCand2->daughter(ghepensimi)->vy()<<" "<<genBsCand2->daughter(ghepensimi)->vz()<<endl;
              bssvx=genBsCand2->daughter(ghepensimi)->vx();
              bssvy=genBsCand2->daughter(ghepensimi)->vy();
              bssvz=genBsCand2->daughter(ghepensimi)->vz();
	
            }
          }
          sort(figlie.begin(),figlie.end());
          if (figlie.size()==request.size() && equal(figlie.begin(),figlie.end(),request.begin())) MCbs=true;
          if ((MCpsi || MCpsi2 || MCpsi3) && MCbs) {cout<<"Channel!"<<endl; bsRootTree_->ChannelID_=0;}
          if (MCphi && MCpsi && MCbs) {cout<<"Signal!"<<endl; bsRootTree_->ChannelID_=1;}
          if (MCphi && MCpsi2 && MCbs) {cout<<"Strange!"<<endl; bsRootTree_->ChannelID_=2;}
          if (MCphi && MCpsi3 && MCbs) {cout<<"Strange2!"<<endl; bsRootTree_->ChannelID_=3;}
          if (MCphi && MCpsi4 && MCbs) {cout<<"Strange3!"<<endl; bsRootTree_->ChannelID_=4;}
        }
        else {
          //bsRootTree_->BsEndFlavour_=(int)(MC_particleID/531);
          //cout<<"No mixing: "<<bsRootTree_->BsEndFlavour_<<endl;
          int numBsDau = genBsCand.numberOfDaughters();
          for (int ghepensimi=0; ghepensimi<numBsDau; ghepensimi++) {
            int numBsDauDau = genBsCand.daughter(ghepensimi)->numberOfDaughters();
            figlie.push_back(genBsCand.daughter(ghepensimi)->pdgId());
            cout<<"Figlie di Bs:"<<genBsCand.daughter(ghepensimi)->pdgId()<<endl;
            vector<int> figlieFiglie;
            for (int ghepensimi2=0; ghepensimi2<numBsDauDau; ghepensimi2++) {
              cout<<"             >"<<genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()<<endl;
              figlieFiglie.push_back(abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()));
              if (abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->charge()==-1)  pmuminus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
              if (abs(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId())==13 && genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->charge()==1)  pmuplus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
	      if (genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==321)  pkplus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
	      if (genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pdgId()==-321)  pkminus.SetXYZM(genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->px(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->py(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->pz(),genBsCand.daughter(ghepensimi)->daughter(ghepensimi2)->mass());
            }
            sort(figlieFiglie.begin(),figlieFiglie.end());
            if (figlieFiglie.size()==requestDauPhi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPhi.begin())) MCphi=true;
            if (figlieFiglie.size()==requestDauPsi.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi.begin())) MCpsi=true;
            if (figlieFiglie.size()==requestDauPsi2.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi2.begin())) MCpsi2=true;
            if (figlieFiglie.size()==requestDauPsi3.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi3.begin())) MCpsi3=true;
            if (figlieFiglie.size()==requestDauPsi4.size() && equal(figlieFiglie.begin(),figlieFiglie.end(),requestDauPsi4.begin())) MCpsi4=true;
            if (abs(genBsCand.daughter(ghepensimi)->pdgId())==443) {
              cout<<"Last Vertex:"<<genBsCand.daughter(ghepensimi)->vx()<<" "<<genBsCand.daughter(ghepensimi)->vy()<<" "<<genBsCand.daughter(ghepensimi)->vz()<<endl;
              cout<<"Bs Vertex:"<<genBsCand.vx()<<" "<<genBsCand.vy()<<" "<<genBsCand.vz()<<endl;
              cout<<"Vertex:"<<genBsCand.mother(0)->vx()<<" "<<genBsCand.mother(0)->vy()<<" "<<genBsCand.mother(0)->vz()<<endl;
              bssvx=genBsCand.daughter(ghepensimi)->vx();
              bssvy=genBsCand.daughter(ghepensimi)->vy();
              bssvz=genBsCand.daughter(ghepensimi)->vz();
				
            }
          }
          sort(figlie.begin(),figlie.end());
          if (figlie.size()==request.size() && equal(figlie.begin(),figlie.end(),request.begin())) MCbs=true;
          if (MCpsi && MCbs) {cout<<"Channel!"<<endl; bsRootTree_->ChannelID_=0;}
          if (MCphi && MCpsi && MCbs) {cout<<"Signal!"<<endl; bsRootTree_->ChannelID_=1;}
          if (MCphi && MCpsi2 && MCbs) {cout<<"Strange!"<<endl; bsRootTree_->ChannelID_=2;}
          if (MCphi && MCpsi3 && MCbs) {cout<<"Strange2!"<<endl; bsRootTree_->ChannelID_=3;}
          if (MCphi && MCpsi4 && MCbs) {cout<<"Strange3!"<<endl; bsRootTree_->ChannelID_=4;}
        }
        if (MC_particleID == 531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) { bsRootTree_->BsIniFlavour_=1; bsRootTree_->BsEndFlavour_=(-1)*BsHasMixed;}
        if (MC_particleID == -531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) { bsRootTree_->BsIniFlavour_=-1; bsRootTree_->BsEndFlavour_=BsHasMixed;}
        //cout<<"Ini flavour: "<< bsRootTree_->BsIniFlavour_<<endl;
        //if (MC_particleID == 531) bsRootTree_->BsIniFlavour_=1;
        //if (MC_particleID == -531) bsRootTree_->BsIniFlavour_=-1;
        if (MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) {
		bsRootTree_->SVZpos_[SVindex] = Double_t(bssvz);
		bsRootTree_->FirstBsMCZpos_ = genBsCand.vz(); 
	// vertex z pos where the first Bs generated: genBsCand.vz() i.e primary vertex z
		SVindex++;
          // calculate gen ctau 2D
          //double Lxy2D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy);
          double Lxy2D = sqrt(pow(bssvx-bspvx,2)+pow(bssvy-bspvy,2));
	  bsRootTree_->BsLxy2DMC_ = Lxy2D;
	
          bsRootTree_->BsCt2DMC_ = Lxy2D * BsPDGMass_/sqrt(bsmomx*bsmomx+bsmomy*bsmomy);
	  bsRootTree_->BsPtMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy);	
          // calculate gen ctau 3D
          //double Lxy3D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy+(bssvz-bspvz)*bsmomz);
          double Lxy3D = sqrt(pow(bssvx-bspvx,2)+pow(bssvy-bspvy,2)+pow(bssvz-bspvz,2));
          bsRootTree_->BsCt3DMC_ = Lxy3D * BsPDGMass_/sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);
	  bsRootTree_->BsLxy3DMC_ = Lxy3D;
	  bsRootTree_->BsPMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);

	      // boosting in JPsi restframe
	      TLorentzVector pjpsi;
	      pjpsi = pmuplus + pmuminus;
	      TLorentzVector pphi;
	      pphi = pkplus + pkminus;
	      
	      // the betas for the boost
	      TVector3 p3_JPsi;
	      p3_JPsi = pjpsi.Vect();
	      p3_JPsi *= -1./pjpsi.E();
	      
	      // the boost matrix
	      TLorentzRotation boost_jpsi(p3_JPsi);
	      TLorentzVector p_JPsi_JPsi;
	      p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);
	      
	      // the different momenta in the new frame
	      TLorentzVector p_JPsi_muplus;
	      TLorentzVector p_JPsi_Kplus;
	      TLorentzVector p_JPsi_phi;
	      
	      p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
	      p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);                                                                                
	      p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);	      
	      
	      // the 3-momenta
	      TVector3 p3_JPsi_muplus;
	      p3_JPsi_muplus = p_JPsi_muplus.Vect();
	      TVector3 p3_JPsi_Kplus;
	      p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	      TVector3 p3_JPsi_phi;
	      p3_JPsi_phi = p_JPsi_phi.Vect();
	      
	      // coordinate system
	      TVector3 x,y,z;
	      x = p3_JPsi_phi.Unit();
	      y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	      y = y.Unit();
	      z = x.Cross(y);
	   
              double angle_costhetaMC=p3_JPsi_muplus.Unit() * z;   
	      bsRootTree_->BscosthetaMC_= angle_costhetaMC;
         //     cout<<"angle_costhetaMC"<<angle_costhetaMC<<endl; 
	      double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
	      double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
	      double angle_phiMC = TMath::ACos(cos_phi);
	      if (sin_phi < 0){
		angle_phiMC =  -angle_phiMC;
	      }
	      bsRootTree_->BsphiMC_ = angle_phiMC;
	      
	      // boosting in phi restframe
	      TVector3 p3_phi;
	      p3_phi = pphi.Vect();
	      p3_phi *= -1./pphi.E();
	      
	      // the boost matrix
	      TLorentzRotation boost_phi(p3_phi);
	      TLorentzVector p_phi_phi;
	      p_phi_phi = boost_phi.VectorMultiplication(pphi);
	      
	      // the different momenta in the new frame
	      TLorentzVector p_phi_Kplus;
	      TLorentzVector p_phi_JPsi;
	      TLorentzVector p_phi_Bs;
	      
	      p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	      p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	      
	      // the 3-momenta
	      TVector3 p3_phi_Kplus;
	      p3_phi_Kplus = p_phi_Kplus.Vect();
	      TVector3 p3_phi_JPsi;
	      p3_phi_JPsi = p_phi_JPsi.Vect();
	      bsRootTree_->BscospsiMC_ = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();
        }
      }      

      // if this is a real B decay (no B mesons as daughters
      if(hasBDaughter == 0){
	//count the number of B decays, should be equal two, for bbbar events
	iNumberOfBdecays++;
	int arrayIndex = iNumberOfBdecays - 1; // array index starts at zero

	// protect array bounds
	if(arrayIndex>=9) break;

		
	if(MC_particleID==(-1)*genBsCand.mother(0)->pdgId()){
//	  cout << "Mixing! " << endl; 		
	  bsRootTree_->BmesonsId_[arrayIndex] = genBsCand.mother(0)->pdgId();
//	  cout << "Original flavour " << genBsCand.mother(0)->pdgId() << endl;		
	} 		
	else bsRootTree_->BmesonsId_[arrayIndex] = MC_particleID;
   
   bsRootTree_->GenNumberOfDaughters_[arrayIndex] = numBsDaughters;
	// generator variables
	bsRootTree_->BMMC_[arrayIndex] = genBsCand.mass();
	bsRootTree_->BPtMC_[arrayIndex] = genBsCand.pt();
	bsRootTree_->BPxMC_[arrayIndex] = genBsCand.px();
	bsRootTree_->BPyMC_[arrayIndex] = genBsCand.py();
	bsRootTree_->BPzMC_[arrayIndex] = genBsCand.pz();
	bsRootTree_->BEtaMC_[arrayIndex] = genBsCand.eta();
	bsRootTree_->BPhiMC_[arrayIndex] = genBsCand.phi();
	bsRootTree_->BVtxMC_x_[arrayIndex] = genBsCand.mother(0)->vx();
	bsRootTree_->BVtxMC_y_[arrayIndex] = genBsCand.mother(0)->vy();
	bsRootTree_->BVtxMC_z_[arrayIndex] = genBsCand.mother(0)->vz();
	//generated primary vertex
	if (abs(MC_particleID)== 531){
	  bsRootTree_->genBsVtx_x_= genBsCand.mother(0)->vx();
	  bsRootTree_->genBsVtx_y_= genBsCand.mother(0)->vy();
	  bsRootTree_->genBsVtx_z_= genBsCand.mother(0)->vz();

	  /*
	  for(int i=0; i <10 ; i++)
	    cout << "gen PV x: " << genBsCand.mother(0)->vx() << 
	      " , gen PV y: " << genBsCand.mother(0)->vy() << 
	      " , gen PV z: " << genBsCand.mother(0)->vz() << endl;  
	  */
	}
	
        for(int j = 0; j < numBsDaughters; ++ j) {
	  if(j>=14) break; // protect array bounds
	  const Candidate * Bsdau = genBsCand.daughter( j );
	  
	  if (abs(Bsdau->pdgId()) == 443) Jpsi = Bsdau;
	  if (abs(Bsdau->pdgId()) == 333) Phi = Bsdau;
	  
	  // generator variables
	  bsRootTree_->BDauIdMC_[arrayIndex][j] = Bsdau->pdgId();
	  bsRootTree_->BDauMMC_[arrayIndex][j] = Bsdau->mass();
	  bsRootTree_->BDauPtMC_[arrayIndex][j] = Bsdau->pt();
	  bsRootTree_->BDauPzMC_[arrayIndex][j] = Bsdau->pz();
	  bsRootTree_->BDauEtaMC_[arrayIndex][j] = Bsdau->eta();
	  bsRootTree_->BDauPhiMC_[arrayIndex][j] = Bsdau->phi();
	  bsRootTree_->BSVtxMC_x_[arrayIndex]   = Bsdau->vx(); 
	  bsRootTree_->BSVtxMC_y_[arrayIndex]   = Bsdau->vy(); 
	  bsRootTree_->BSVtxMC_z_[arrayIndex]   = Bsdau->vz(); 
	  //Generated secondary vertex.
	  if ( abs(Bsdau->pdgId())== 443){
	    bsRootTree_->genBsSVtx_x_= Bsdau->vx();
	    bsRootTree_->genBsSVtx_y_= Bsdau->vy();
	    bsRootTree_->genBsSVtx_z_= Bsdau->vz();
	    /*	    
	    for(int i=0; i <10 ; i++)
	      cout << "gen SV x: " << Bsdau->vx() <<
		" , gen SV y: " << Bsdau->vy() <<
		" , gen SV z: " << Bsdau->vz() << endl;
	    */
	  }

	  // daughter of daughter (muons, kaons in case of jpsi phi)
	  int numBsDaughtersDaughters = Bsdau->numberOfDaughters();
	  bsRootTree_->GenNumberOfDaughtersDaughters_[arrayIndex][j] = numBsDaughtersDaughters;
	  for(int k=0; k< numBsDaughtersDaughters; k++){
	    if(k>=9) break; //protect array bounds
	    const Candidate * Bsdaudau = Bsdau->daughter(k);

	    if ( Bsdaudau->pdgId() == -13) mup = Bsdaudau;
	    if ( Bsdaudau->pdgId() == 13) mum = Bsdaudau;
	    if ( Bsdaudau->pdgId() == 321) Kp = Bsdaudau;
	    if ( Bsdaudau->pdgId() == -321) Km = Bsdaudau;
	    
	    // generator variables
	    bsRootTree_->BDauDauIdMC_[arrayIndex][j][k] = Bsdaudau->pdgId();
	    bsRootTree_->BDauDauMMC_[arrayIndex][j][k] = Bsdaudau->mass();
	    bsRootTree_->BDauDauPtMC_[arrayIndex][j][k] = Bsdaudau->pt();
	    bsRootTree_->BDauDauPzMC_[arrayIndex][j][k] = Bsdaudau->pz();
	    bsRootTree_->BDauDauEtaMC_[arrayIndex][j][k] = Bsdaudau->eta();
	    bsRootTree_->BDauDauPhiMC_[arrayIndex][j][k] = Bsdaudau->phi();

	    if (abs(MC_particleID)== 531 && numBsDaughters == 2 && Jpsi && Phi && mup && mum && Kp && Km){

	      TLorentzVector pmuplus;
	      TLorentzVector pmuminus;	      
	      TLorentzVector pkplus;
	      TLorentzVector pkminus;
	      // extra check
	      if (mup)
		pmuplus.SetXYZM(mup->px(),mup->py(),mup->pz(),mup->mass());
	      if (mum)
		pmuminus.SetXYZM(mum->px(),mum->py(),mum->pz(),mum->mass());
	      if (Kp)
		pkplus.SetXYZM(Kp->px(),Kp->py(),Kp->pz(),Kp->mass());
	      if (Km)
		pkminus.SetXYZM(Km->px(),Km->py(),Km->pz(),Km->mass());
	     

 
	      /*
	      cout << "Mu+ " << mup->pdgId() << " (px,py,pz,m): (" << mup->px() << "," << mup->py() << "," << mup->pz() << "," << mup->mass() << ")" << endl;    	    
	      cout << "Mu- " << mum->pdgId() << " (px,py,pz,m): (" << mum->px() << "," << mum->py() << "," << mum->pz() << "," << mum->mass() << ")" << endl;    	    
	      cout << "K+ " << Kp->pdgId() << " (px,py,pz,m): (" << Kp->px() << "," << Kp->py() << "," << Kp->pz() << "," << Kp->mass() << ")" << endl;    	    
	      cout << "K- " << Km->pdgId() << " (px,py,pz,m): (" << Km->px() << "," << Km->py() << "," << Km->pz() << "," << Km->mass() << ")" << endl;    	    
	      */


	      // boosting in JPsi restframe
	      TLorentzVector pjpsi;
	      pjpsi = pmuplus + pmuminus;
	      TLorentzVector pphi;
	      pphi = pkplus + pkminus;
	      
	      // the betas for the boost
	      TVector3 p3_JPsi;
	      p3_JPsi = pjpsi.Vect();
	      p3_JPsi *= -1./pjpsi.E();
	      
	      // the boost matrix
	      TLorentzRotation boost_jpsi(p3_JPsi);
	      TLorentzVector p_JPsi_JPsi;
	      p_JPsi_JPsi = boost_jpsi.VectorMultiplication(pjpsi);
	      
	      // the different momenta in the new frame
	      TLorentzVector p_JPsi_muplus;
	      TLorentzVector p_JPsi_Kplus;
	      TLorentzVector p_JPsi_phi;
	      
	      p_JPsi_muplus = boost_jpsi.VectorMultiplication(pmuplus);
	      p_JPsi_Kplus = boost_jpsi.VectorMultiplication(pkplus);                                                                                
	      p_JPsi_phi = boost_jpsi.VectorMultiplication(pphi);	      
	      
	      // the 3-momenta
	      TVector3 p3_JPsi_muplus;
	      p3_JPsi_muplus = p_JPsi_muplus.Vect();
	      TVector3 p3_JPsi_Kplus;
	      p3_JPsi_Kplus = p_JPsi_Kplus.Vect();
	      TVector3 p3_JPsi_phi;
	      p3_JPsi_phi = p_JPsi_phi.Vect();
	      
	      // coordinate system
	      TVector3 x,y,z;
	      x = p3_JPsi_phi.Unit();
	      y = p3_JPsi_Kplus.Unit() - (x * (x * p3_JPsi_Kplus.Unit()));
	      y = y.Unit();
	      z = x.Cross(y);
	   
              double angle_costhetaMC=p3_JPsi_muplus.Unit() * z;   
	      bsRootTree_->costhetaMC_[arrayIndex] = angle_costhetaMC;
           
	      double cos_phi = p3_JPsi_muplus.Unit() * x / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
	      double sin_phi = p3_JPsi_muplus.Unit() * y / TMath::Sqrt(1 - angle_costhetaMC*angle_costhetaMC);
	      double angle_phiMC = TMath::ACos(cos_phi);
	      if (sin_phi < 0){
		angle_phiMC =  -angle_phiMC;
	      }
	      bsRootTree_->phiMC_[arrayIndex] = angle_phiMC;
	      
	      // boosting in phi restframe
	      TVector3 p3_phi;
	      p3_phi = pphi.Vect();
	      p3_phi *= -1./pphi.E();
	      
	      // the boost matrix
	      TLorentzRotation boost_phi(p3_phi);
	      TLorentzVector p_phi_phi;
	      p_phi_phi = boost_phi.VectorMultiplication(pphi);
	      
	      // the different momenta in the new frame
	      TLorentzVector p_phi_Kplus;
	      TLorentzVector p_phi_JPsi;
	      TLorentzVector p_phi_Bs;
	      
	      p_phi_Kplus = boost_phi.VectorMultiplication(pkplus);
	      p_phi_JPsi = boost_phi.VectorMultiplication(pjpsi);
	      
	      // the 3-momenta
	      TVector3 p3_phi_Kplus;
	      p3_phi_Kplus = p_phi_Kplus.Vect();
	      TVector3 p3_phi_JPsi;
	      p3_phi_JPsi = p_phi_JPsi.Vect();
	      bsRootTree_->cospsiMC_[arrayIndex] = -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit();

	      //	      cout << "(costheta,phi,cospsi): "<< p3_JPsi_muplus.Unit() * z << "," << angle_phi << "," << -1 * p3_phi_Kplus.Unit() * p3_phi_JPsi.Unit() << ")" << endl;
	      
	    }
	    
	  }// loop Bs daughters daughters
	} // loop Bs daughters
	
	
	// calculate gen ctau 2D
	double Lxy2D = ((bsRootTree_->BSVtxMC_x_[arrayIndex] - bsRootTree_->BVtxMC_x_[arrayIndex])*bsRootTree_->BPxMC_[arrayIndex]+
			(bsRootTree_->BSVtxMC_y_[arrayIndex] - bsRootTree_->BVtxMC_y_[arrayIndex])*bsRootTree_->BPyMC_[arrayIndex]);
	bsRootTree_->BCt_MC2D_[arrayIndex] = Lxy2D * bsRootTree_->BMMC_[arrayIndex]/(bsRootTree_->BPtMC_[arrayIndex]*bsRootTree_->BPtMC_[arrayIndex]); 	
	
	// calculate gen ctau 3D
	double Lxy3D = ((bsRootTree_->BSVtxMC_x_[arrayIndex] - bsRootTree_->BVtxMC_x_[arrayIndex])*bsRootTree_->BPxMC_[arrayIndex]+
			(bsRootTree_->BSVtxMC_y_[arrayIndex] - bsRootTree_->BVtxMC_y_[arrayIndex])*bsRootTree_->BPyMC_[arrayIndex]+
			(bsRootTree_->BSVtxMC_z_[arrayIndex] - bsRootTree_->BVtxMC_z_[arrayIndex])*bsRootTree_->BPzMC_[arrayIndex]);
	bsRootTree_->BCt_MC3D_[arrayIndex] = Lxy3D * bsRootTree_->BMMC_[arrayIndex]/
	  (bsRootTree_->BPtMC_[arrayIndex]*bsRootTree_->BPtMC_[arrayIndex]+
	   bsRootTree_->BPzMC_[arrayIndex]*bsRootTree_->BPzMC_[arrayIndex]); 	
	
	// calculate gen ctau 
	double deltaX =  bsRootTree_->BSVtxMC_x_[arrayIndex] - 	bsRootTree_->BVtxMC_x_[arrayIndex];
	double deltaY =  bsRootTree_->BSVtxMC_y_[arrayIndex] - 	bsRootTree_->BVtxMC_y_[arrayIndex];
	bsRootTree_->BLxy_MC_[arrayIndex] = sqrt( deltaX*deltaX + deltaY*deltaY);
	if(deltaX * genBsCand.px() + deltaY * genBsCand.py() < 0 )  bsRootTree_->BLxy_MC_[arrayIndex] = -1. *  bsRootTree_->BLxy_MC_[arrayIndex];
	bsRootTree_->BCt_MC_[arrayIndex] = bsRootTree_->BLxy_MC_[arrayIndex] * bsRootTree_->BMMC_[arrayIndex] / bsRootTree_->BPtMC_[arrayIndex];
      }
    }
  
    // check if there is a Jpsi (prompt or non-prompt) in the event
    if(absMC_particleID == 443 ) bsRootTree_->isGenJpsiEvent_ = 1;
    
  }

  bsRootTree_->GenNumberOfBdecays_ = iNumberOfBdecays;

}



bool BsToJpsiPhiAnalysis::selGlobalMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx) {

  TrackRef iTrack = aMuon.innerTrack();
  const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(aMuon.originalObject());
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

  TrackRef gTrack = aMuon.globalTrack();
  const reco::HitPattern& q = gTrack->hitPattern();

  bool trackOK = false;
  // cooler way of cutting on tracks
//  if (_applyExpHitcuts) {
//    float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
//    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
//  } else  
  trackOK = (iTrack->found() > 10);

  return (// isMuonInAccept(aMuon) &&
          aMuon.isPFMuon() &&
          gTrack->chi2()/gTrack->ndof() < 10.0 &&
          gTrack->hitPattern().numberOfValidMuonHits()>0 &&
          rmu1->numberOfMatchedStations()>1 &&
          p.numberOfValidMuonHits() > 0 &&
          p.trackerLayersWithMeasurement()>5);
//	  trackOK &&
//  	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
//  	  aMuon.muonID("TrackerMuonArbitrated") &&
//  	  aMuon.muonID("TMOneStationTight") &&
//          p.pixelLayersWithMeasurement() > 1 &&
//	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
//          fabs(iTrack->dz(RefVtx)) < 15.0);

}


bool BsToJpsiPhiAnalysis::selTrackerMuon(const pat::Muon aMuon, const math::XYZPoint RefVtx) {
  
  TrackRef iTrack = aMuon.innerTrack();
  const reco::HitPattern& p = iTrack->hitPattern();
  const reco::HitPattern& ei = iTrack->trackerExpectedHitsInner();
  const reco::HitPattern& eo = iTrack->trackerExpectedHitsOuter();

  bool trackOK = false;
  // cooler way of cutting on tracks
//  if (_applyExpHitcuts) {
//    float fHits = iTrack->found() / (iTrack->found() + iTrack->lost() + ei.numberOfHits() + eo.numberOfHits());
//    trackOK = (fHits >= 0.8 && (p.hasValidHitInFirstPixelBarrel() || p.hasValidHitInFirstPixelEndcap() ));
  // old way of cutting on tracks  
//  } else
 trackOK = (iTrack->found() > 10);

  return (// isMuonInAccept(aMuon) &&
  	  trackOK &&
   	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
  	  aMuon.muonID("TrackerMuonArbitrated") &&
 	  aMuon.muonID("TMOneStationTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 && 
          fabs(iTrack->dz(RefVtx)) < 15.0 );
}


//------------------------------------------
void BsToJpsiPhiAnalysis::setFitParKK(RefCountedKinematicTree& myTree)
{
  
  
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
  
  // first particle: kaon 1   
  
  AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++)  bsRootTree_->K1Fit_par_[i] = bs_par1[i];
  
  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->K1Fit_sigX_ = sqrt(bs_err1(0,0));
  bsRootTree_->K1Fit_sigY_ = sqrt(bs_err1(1,1));
  bsRootTree_->K1Fit_sigZ_ = sqrt(bs_err1(2,2));
  bsRootTree_->K1Fit_sigPX_ = sqrt(bs_err1(3,3));
  bsRootTree_->K1Fit_sigPY_ = sqrt(bs_err1(4,4));
  bsRootTree_->K1Fit_sigPZ_ = sqrt(bs_err1(5,5));
  
  // first particle: kaon 2  
  
    
  AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++) bsRootTree_->K2Fit_par_[i] = bs_par2[i];
  
  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->K2Fit_sigX_ = sqrt(bs_err2(0,0));
  bsRootTree_->K2Fit_sigY_ = sqrt(bs_err2(1,1));
  bsRootTree_->K2Fit_sigZ_ = sqrt(bs_err2(2,2));
  bsRootTree_->K2Fit_sigPX_ = sqrt(bs_err2(3,3));
  bsRootTree_->K2Fit_sigPY_ = sqrt(bs_err2(4,4));
  bsRootTree_->K2Fit_sigPZ_ = sqrt(bs_err2(5,5));
  
}


//------------------------------------------
void BsToJpsiPhiAnalysis::setFitParHyp1(RefCountedKinematicTree& myTree)
{
  
  
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
  
  // first particle: kaon 1   
  
  AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++)  bsRootTree_->BdK1_kpi_par_Hyp1_[i] = bs_par1[i];
  
  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK1_kpi_sigX_Hyp1_ = sqrt(bs_err1(0,0));
  bsRootTree_->BdK1_kpi_sigY_Hyp1_ = sqrt(bs_err1(1,1));
  bsRootTree_->BdK1_kpi_sigZ_Hyp1_ = sqrt(bs_err1(2,2));
  bsRootTree_->BdK1_kpi_sigPX_Hyp1_ = sqrt(bs_err1(3,3));
  bsRootTree_->BdK1_kpi_sigPY_Hyp1_ = sqrt(bs_err1(4,4));
  bsRootTree_->BdK1_kpi_sigPZ_Hyp1_ = sqrt(bs_err1(5,5));
  
  // first particle: kaon 2  
  
    
  AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++) bsRootTree_->BdK2_kpi_par_Hyp1_[i] = bs_par2[i];
  
  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK2_kpi_sigX_Hyp1_ = sqrt(bs_err2(0,0));
  bsRootTree_->BdK2_kpi_sigY_Hyp1_ = sqrt(bs_err2(1,1));
  bsRootTree_->BdK2_kpi_sigZ_Hyp1_ = sqrt(bs_err2(2,2));
  bsRootTree_->BdK2_kpi_sigPX_Hyp1_ = sqrt(bs_err2(3,3));
  bsRootTree_->BdK2_kpi_sigPY_Hyp1_ = sqrt(bs_err2(4,4));
  bsRootTree_->BdK2_kpi_sigPZ_Hyp1_ = sqrt(bs_err2(5,5));
  
}

//------------------------------------------
void BsToJpsiPhiAnalysis::setFitParHyp2(RefCountedKinematicTree& myTree)
{
  
  
  vector< RefCountedKinematicParticle > bs_children = myTree->finalStateParticles();
  
  // first particle: kaon 1   
  
  AlgebraicVector7 bs_par1 = bs_children[0]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++)  bsRootTree_->BdK1_kpi_par_Hyp2_[i] = bs_par1[i];
  
  AlgebraicSymMatrix77 bs_err1 = bs_children[0]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK1_kpi_sigX_Hyp2_ = sqrt(bs_err1(0,0));
  bsRootTree_->BdK1_kpi_sigY_Hyp2_ = sqrt(bs_err1(1,1));
  bsRootTree_->BdK1_kpi_sigZ_Hyp2_ = sqrt(bs_err1(2,2));
  bsRootTree_->BdK1_kpi_sigPX_Hyp2_ = sqrt(bs_err1(3,3));
  bsRootTree_->BdK1_kpi_sigPY_Hyp2_ = sqrt(bs_err1(4,4));
  bsRootTree_->BdK1_kpi_sigPZ_Hyp2_ = sqrt(bs_err1(5,5));
  
  // first particle: kaon 2  
  
    
  AlgebraicVector7 bs_par2 = bs_children[1]->currentState().kinematicParameters().vector();
  for(int i=0; i<7; i++) bsRootTree_->BdK2_kpi_par_Hyp2_[i] = bs_par2[i];
  
  AlgebraicSymMatrix77 bs_err2 = bs_children[1]->currentState().kinematicParametersError().matrix();
  bsRootTree_->BdK2_kpi_sigX_Hyp2_ = sqrt(bs_err2(0,0));
  bsRootTree_->BdK2_kpi_sigY_Hyp2_ = sqrt(bs_err2(1,1));
  bsRootTree_->BdK2_kpi_sigZ_Hyp2_ = sqrt(bs_err2(2,2));
  bsRootTree_->BdK2_kpi_sigPX_Hyp2_ = sqrt(bs_err2(3,3));
  bsRootTree_->BdK2_kpi_sigPY_Hyp2_ = sqrt(bs_err2(4,4));
  bsRootTree_->BdK2_kpi_sigPZ_Hyp2_ = sqrt(bs_err2(5,5));
  
}



bool  BsToJpsiPhiAnalysis::MCmatching(const Candidate & track1,  edm::Handle<GenParticleCollection> & genParticles,
				      int &K1mcId, int &K1momId, int &K1gmomId,
				      int condMom, int condGMom){
  if(!isMCstudy_ ) return 0;
  bool K1Truth = 0;
  double MinDRK=999.;
  
  K1mcId = -9999999;
  K1momId = -9999999;
  K1gmomId = -9999999;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    double DeltaRK1 = deltaR(p.eta(), p.phi(), track1.eta(), track1.phi() );
   
    if(DeltaRK1<MinDRK && DeltaRK1<0.05){
      MinDRK=DeltaRK1;
      K1mcId=p.pdgId();
		
      if(p.mother()!=0) K1momId=p.mother()->pdgId();
      if(p.mother()!=0 && p.mother()->mother()!=0) {K1gmomId=p.mother()->mother()->pdgId(); }
		
      if (abs(K1momId)==condMom && abs(K1gmomId)==condGMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }
//	cout << "Muon id " << K1mcId << endl;
//	cout << "Muon gmom id " << K1gmomId << endl;

  return K1Truth;
  
}

bool BsToJpsiPhiAnalysis::MCmatchingBplusK(const Candidate & track1,  edm::Handle<GenParticleCollection> & genParticles,
				      int &K1mcId, int &K1momId,
				      int condMom){
  if(!isMCstudy_ ) return 0;
  bool K1Truth = 0;
  double MinDRK=999.;
  
  K1mcId = -9999999;
  K1momId = -9999999;


  for(size_t i = 0; i < genParticles->size(); ++ i){
    const GenParticle & p = (*genParticles)[i];
    double DeltaRK1 = deltaR(p.eta(), p.phi(), track1.eta(), track1.phi() );
   
    if(DeltaRK1<MinDRK && DeltaRK1<0.05){
      MinDRK=DeltaRK1;
      K1mcId=p.pdgId();
	
      if(p.mother()!=0) {K1momId=p.mother()->pdgId(); } 
      if (abs(K1momId)==condMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }

//	cout << "KmcId " <<  K1mcId << endl;
//   cout << "Kmom Id " << K1momId << endl;
  return K1Truth;
  
}

reco::Vertex BsToJpsiPhiAnalysis::reVertex(const edm::EventSetup& iSetup, reco::BeamSpot vertexBeamSpot, reco::Vertex RecVtx, pat::Muon mu1, pat::Muon mu2, TrackRef trk3, TrackRef trk4){

        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
        std::vector<reco::TransientTrack> vrtxRefit;
        bool removedTracks = false;
        for (reco::Vertex::trackRef_iterator trackvertex=RecVtx.tracks_begin(); trackvertex!=RecVtx.tracks_end(); trackvertex++){
          //const reco::Track & RTv=*(trackvertex->get());
          TrackRef tref = trackvertex->castTo<TrackRef>();
          if (!tref.isNonnull()) continue;
          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mu1.originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mu2.originalObject());
          if (rmu1->track().key()==trackvertex->key() || rmu2->track().key()==trackvertex->key() || trk3.key()==trackvertex->key() || trk4.key()==trackvertex->key()) {
            removedTracks = true;
          } else {
            TransientTrack tTrk = (*theB).build(&tref); //  fpTTB->build(*(*itTBR));i
            if (tTrk.isValid())  vrtxRefit.push_back(tTrk);
          }
        }
        AdaptiveVertexFitter avf;
        TransientVertex newVtx = avf.vertex(vrtxRefit, vertexBeamSpot);
        if (newVtx.isValid()) {
          return reco::Vertex(newVtx);
        } else {
          return RecVtx;
        }
}
	    
