// description on https://twiki.cern.ch/twiki/bin/viewauth/CMS/BsJpsiPhi_AWG
#include "TrackingTools/IPTools/interface/IPTools.h"  
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
#include <vector>
#include <TMath.h>

using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;


BsToJpsiPhiAnalysis::BsToJpsiPhiAnalysis(const edm::ParameterSet& iConfig):theConfig_(iConfig), 
//matcher(iConfig),
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
  muonTag_ = iConfig.getParameter<edm::InputTag>("MuonTag");
  jetCollection_ = iConfig.getParameter<edm::InputTag>("JetCollection");


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


  verbose_                = iConfig.getParameter<bool>("verbose");
  outputFile_ = iConfig.getUntrackedParameter<std::string>("outputFile");
  event_counter_ = 0;
  mcNomatch_counter_ = 0;
  tagmuonNum_  =0;

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
  cout << "Total number of unmatched MC muons " << mcNomatch_counter_ << endl;
  cout << "Total number of MC muons " << tagmuonNum_ << endl;
}

void
BsToJpsiPhiAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  event_counter_++;

  bsRootTree_->resetEntries();

  pat::CompositeCandidate BCand_best;
  TrackRef trk1Ref_best;
  TrackRef trk2Ref_best;
  TrackRef trkMu1Ref_best;
  TrackRef trkMu2Ref_best;
  RefCountedKinematicParticle bs_best;
  pat::Muon mu1_best;
  pat::Muon mu2_best;

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


  Float_t chi = 0.;
  Float_t ndf = 0.;
  Float_t muon_sigma = 0.0000000001;
  Float_t kaon_sigma = 0.000016; 

  // Get primary vertices
  Double_t bestVtxProbBplus=-1; 

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

  double BplusPVx;
  double BplusPVy;
  double BplusPVz;
  double BplusPVerrx;
  double BplusPVerry;
  double BplusPVerrz;


  double BsLxyz=0;
  Int_t NoPVcounter=0;

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
  
//  edm::Handle<reco::PFJetCollection> jets;
//  iEvent.getByLabel("ak5PFJets",jets);

   edm::Handle<edm::View<pat::Jet>> myjets;
   iEvent.getByLabel(jetCollection_,myjets);
   const edm::View<pat::Jet> & jets = *myjets;

  bsRootTree_->PtJetTrkCharge_->clear();
  bsRootTree_->PtJetTrkPt_->clear();
  

  int PtJetIndex = -1;
 
  for(size_t i = 0; i< jets.size(); i++){
//   cout << "i " << i << " pt jet pt " <<  jets[i].pt() << endl;
	if(jets[i].pt() > 10){

	PtJetIndex = i;
 	break;
       } 
 }


  if(PtJetIndex!=-1){

	bsRootTree_->PtJetPt_  = jets[PtJetIndex].pt();
	const reco::TrackRefVector & jetTrackVec = jets[PtJetIndex].associatedTracks();

	for(size_t k =0; k < jetTrackVec.size(); k++){
		TrackRef jetTrk = jetTrackVec[k];
//		cout << "trackref p " << jetTrk->momentum() << endl;
//		cout << "trackref charge " << jetTrk->charge() << endl;

		const reco::Track *jetTrack = jetTrk.get();
//		cout << "track pt " << jetTrack->pt() << endl;
//		cout << "track charge " << jetTrack->charge() << endl;
		bsRootTree_->PtJetTrkCharge_->push_back(jetTrack->charge());
		bsRootTree_->PtJetTrkPt_->push_back(jetTrack->pt());
	 } 
  }
 

// jetProbabilityBJetTags
//   const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);

  double BprobMax = 0;
  int BjetIndex=-1;
  for(size_t i =0; i < jets.size(); i++){
	double Bprobability = jets[i].bDiscriminator("jetProbabilityBJetTags"); 
	
	if(TMath::Abs(Bprobability) > BprobMax && jets[i].pt() > 10){
		BjetIndex = i;
		BprobMax = Bprobability;
//		cout << Bprobability << endl; 	
	}
  }
 
	

  bsRootTree_->BJetTrkCharge_->clear();
  bsRootTree_->BJetTrkPt_->clear();

  if(BjetIndex!=-1){
	bsRootTree_->JetBTagProb_ = jets[BjetIndex].bDiscriminator("jetProbabilityBJetTags"); 
	bsRootTree_->BJetPt_ = jets[BjetIndex].pt();

//	cout << "Max Bjet prob " << jets[BjetIndex].bDiscriminator("jetProbabilityBJetTags") << endl;
//	cout << "Bjet pt " << jets[BjetIndex].pt() << endl;

	const reco::TrackRefVector & jetTrackVec = jets[BjetIndex].associatedTracks();

	 for(size_t k =0; k < jetTrackVec.size(); k++){
		TrackRef jetTrk = jetTrackVec[k];
		const reco::Track *jetTrack = jetTrk.get();
		bsRootTree_->BJetTrkPt_->push_back(jetTrack->pt());
		bsRootTree_->BJetTrkCharge_->push_back(jetTrack->charge());
	}	
 }

/*
 const reco::TrackRefVector & 	jetTrackVec = jets[closestJetIndex].associatedTracks();
//  cout << "tracks in jet " << jetTrackVec.size() << endl;

  vector<int> *trackCharge = new vector<int>();
  vector<float> *trackPt = new vector<float>();
  
 	 for(size_t k =0; k < jetTrackVec.size(); k++){
		TrackRef jetTrk = jetTrackVec[k];
//		cout << "trackref p " << jetTrk->momentum() << endl;
//		cout << "trackref charge " << jetTrk->charge() << endl;

		const reco::Track *jetTrack = jetTrk.get();
//		cout << "track pt " << jetTrack->pt() << endl;
//		cout << "track charge " << jetTrack->charge() << endl;
		trackCharge->push_back( jetTrack->charge() );
		trackPt->push_back(jetTrack->pt());
	 }

//  bsRootTree_->JetTrkCharge_->push_back(trackCharge);
//  bsRootTree_->JetTrkPt_->push_back(trackPt);

	cout << "tracks in jet: " << jetTrackVec.size() << " = " << bsRootTree_->JetTrkCharge_->at(NumMu)->size() << " = " << bsRootTree_->JetTrkPt_->at(NumMu)->size() << endl; 
*/

/*
  edm::Handle<reco::JetTagCollection> bTagHandle;
iEvent.getByLabel("jetProbabilityBJetTags", bTagHandle);
 const reco::JetTagCollection & bTags = *(bTagHandle.product());
//const JetTagCollection &bTags = *(bTagHandle);	

Double_t dRMin = 1000000;


if(jets->size()!=0  ){
  const PFJet &jet = (*jets)[0];	

  for(size_t k = 0; k< bTags.size(); k++){


  Double_t dR = deltaR(bTags[k].first->eta(), bTags[k].first->phi(), jet.eta(), jet.phi() );
	if(dR < 0.02 && dR < dRMin ){
		dRMin = dR;
		bsRootTree_->JetBTagProb_= bTags[k].second;
		//cout << " Jet B tag Prob " <<bTags[k].second << endl;
		//cout << "dRMin " << dRMin << endl; 	
	}	
  }	
}
*/

  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices", recVtxs);
 
  bsRootTree_->NVertices_=recVtxs->size(); 
  double minPtVertex=0;
 
  
  for(size_t i = 0; i < recVtxs->size(); ++i) {
    const Vertex &vtx = (*recVtxs)[i];
    double RecVtxProb=TMath::Prob(vtx.chi2(),(int)vtx.ndof());
 		
    Double_t summaPtVertex=0;
    Int_t trackNumberCounter =0;	
    for (reco::Vertex::trackRef_iterator trackvertex=vtx.tracks_begin(); trackvertex!=vtx.tracks_end(); trackvertex++){
      const reco::Track & RTv=*(trackvertex->get());
      summaPtVertex=summaPtVertex+abs(RTv.pt());
	trackNumberCounter++;
	
    }
	//cout <<" trackNumberCounter "<< trackNumberCounter  << endl;
	//cout <<" summaPtVertex "<< summaPtVertex  << endl;

	if(i <= 29){

		bsRootTree_->PVZpos_[i] = Double_t( vtx.z() );
		bsRootTree_->NTracksInPV_[i] = trackNumberCounter;
		bsRootTree_->PVAbsPt_[i] = summaPtVertex;
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
	
	//selection of highest pt PV (old)
	 //const Vertex RecVtx muutettu Vertex RecVtx alla oleva koodi kääntyy ja toimii
   Vertex RecVtx;
   Vertex BplusPVVtx;
   Vertex BdPVVtx;  
   Vertex PVvtxCosTheta; 
   Vertex PVvtxHighestPt; 
   Bool_t isHighestPtPV = false;

   if(VtxIndex != -99){	
   PVvtxHighestPt = (*recVtxs)[VtxIndex];
   isHighestPtPV = true;
/*
   cout <<"Highest pt PVx "<< PVvtxHighestPt.x() << endl;
   cout <<"Highest pt PVy "<< PVvtxHighestPt.y() << endl;
   cout <<"Highest pt PVz "<< PVvtxHighestPt.z() << endl;
*/
   }
/*

  if(VtxIndex!=-99) //if the PV with highest pt is found, it's selected to be a "right" PV
    {
      bsRootTree_->isPV_ = 1;
      PVx = RecVtx.x(); 
      PVy= RecVtx.y();
      PVz= RecVtx.z();
      PVerrx=RecVtx.xError();
      PVerry=RecVtx.yError();
      PVerrz=RecVtx.zError();
    }
  else { //if there is not reconstructed PV's the PV is chosen to be the beam spot
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
  } */

  // run/event number 
  bsRootTree_->runNumber_ = iEvent.id().run();
  bsRootTree_->eventNumber_ = (unsigned int)iEvent.id().event();
  bsRootTree_->lumiSection_ = iEvent.luminosityBlock();

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

    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v9")bsRootTree_->triggerbit_4Jpsi_Displaced_v9_ = hltresults->accept(itrig);
//New trigger for the summer12 BsKKJPsi MC

if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v10")bsRootTree_->triggerbit_4Jpsi_Displaced_v10_ = hltresults->accept(itrig);
//New for data 2012

if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v11")bsRootTree_->triggerbit_4Jpsi_Displaced_v11_ = hltresults->accept(itrig);
//New for data 2012

if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v12")bsRootTree_->triggerbit_4Jpsi_Displaced_v12_ = hltresults->accept(itrig);
//New for data 2012

    if (trigName=="HLT_DoubleMu4_Jpsi_Displaced_v5")bsRootTree_->triggerbit_4Jpsi_Displaced_v5_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v1")bsRootTree_->triggerbit_5Jpsi_Displaced_v1_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v2")bsRootTree_->triggerbit_5Jpsi_Displaced_v2_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v4")bsRootTree_->triggerbit_5Jpsi_Displaced_v4_ = hltresults->accept(itrig);
    if (trigName=="HLT_DoubleMu5_Jpsi_Displaced_v5")bsRootTree_->triggerbit_5Jpsi_Displaced_v5_ = hltresults->accept(itrig);

//new trigger for MC 2012 and data:
 if (trigName=="HLT_Dimuon0_Jpsi_Muon_v15")bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v15_ = hltresults->accept(itrig);
 if (trigName=="HLT_Dimuon0_Jpsi_Muon_v16")bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v16_ = hltresults->accept(itrig);
 if (trigName=="HLT_Dimuon0_Jpsi_Muon_v17")bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v17_ = hltresults->accept(itrig);
if (trigName=="HLT_Dimuon0_Jpsi_Muon_v18")bsRootTree_->triggerbit_Dimuon0_Jpsi_Muon_v18_ = hltresults->accept(itrig);

//associated filter: hltVertexmumuFilterJpsiMuon

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
  }
  
  edm::Handle<pat::MuonCollection> allmuons;
  iEvent.getByLabel(muonTag_,allmuons);
  
  if(verbose_==true){
    if(allmuons->size()>0) std::cout<<"******found number of muons= "<< allmuons->size() << std::endl;
  }
  
  vector<Vertex> PVvec;
  // variables to determine minima of fit probability
  double minVtxP = -99.;   //KK hypothesis
  double minJpsiP = -99;   // jpsi alone


//	cout << "LOOPPI ALKAA" <<endl;


  // loop on muons (2 muons at least)
  for(size_t i=0; i < allmuons->size(); ++i){
    const pat::Muon & mu1 = (*allmuons)[i];
    
    for (size_t j=i+1; j < allmuons->size(); ++j){
      const pat::Muon & mu2 = (*allmuons)[j];
      
   //   if (allmuons->size() > 25) cout << "WARNING! muon list oversize:" << allmuons->size() << endl;      

	//new part: only muons that fire the trigger are used to construct the B candidates 
      if ((mu1.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty() || mu2.triggerObjectMatchesByFilter("hltDisplacedmumuFilterDoubleMu4Jpsi").empty()) && (mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty() || mu2.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty())) continue;

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
/*
      bsRootTree_->RecoMuListSize_ = allmuons->size();      
      bsRootTree_->MuRecoPt1_[i] = mu1.pt();
      bsRootTree_->MuRecoEta1_[i] = mu1.eta();
      bsRootTree_->MuRecoPhi1_[i] = mu1.phi();
      bsRootTree_->MuRecoChg1_[i] = mu1.charge();
      
      bsRootTree_->MuRecoPt2_[j] = mu2.pt();
      bsRootTree_->MuRecoEta2_[j] = mu2.eta();
      bsRootTree_->MuRecoPhi2_[j] = mu2.phi();
      bsRootTree_->MuRecoChg2_[j] = mu2.charge();
*/      


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
   	 
      KalmanVertexFitter kvf;
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
          // Muons overlapping remover
          if ( muon::overlap(mu1,mu2,1,1,true) ) continue;
	  
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
	  
	  double vtxprob_Bs = Kfitter.getProb();
	  //TEST double vtxprob_Bs = BCand.pt();
	  RefCountedKinematicParticle bs = Kfitter.getParticle();
          //bs_best = bs;
	  RefCountedKinematicVertex bVertex = Kfitter.getVertex();
	  AlgebraicVector7 b_par = bs->currentState().kinematicParameters().vector();
	  AlgebraicSymMatrix77 bs_er = bs->currentState().kinematicParametersError().matrix();
	  
	  double fittedBsMass = b_par[6];
	  
	  // check if it is a valid candidate to be counted (verify passed AfterFit cuts)
	  if (abs(Jpsi.mass() - nominalJpsiMass) < JpsiMassWindowAfterFit_ && Jpsi.pt() > JpsiPtCut_ &&
	      abs(PhiCand.mass() - nominalPhiMass) < PhiMassWindowAfterFit_ &&
	      fittedBsMass > BsLowerMassCutAfterFit_ && fittedBsMass < BsUpperMassCutAfterFit_  ) bsRootTree_->BsNumberOfCandidatesAfterFit_++;

	  		  
	  // store values in root tree if vtx probability is higher than already stored one
	  if (vtxprob_Bs > minVtxP){
	    
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

	bsRootTree_->matchDoubleMu01DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
	bsRootTree_->matchDoubleMu02DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();

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
	   
/*	    cout << "SV x refit"<<reFitVertex.x() <<endl;
	    cout << "SV x fit"<<bVertex->position().x() <<endl;
 
	    cout << "SV y refit"<<reFitVertex.y() <<endl;
	    cout << "SV y fit"<<bVertex->position().y() <<endl;

	    cout << "SV z refit"<<reFitVertex.z() <<endl;
	    cout << "SV z fit"<<bVertex->position().z() <<endl;	*/	

	    bsRootTree_->PVerrx_refit_ = reFitVertex.xError();
	    bsRootTree_->PVerry_refit_ = reFitVertex.yError();
	    bsRootTree_->PVerrz_refit_ = reFitVertex.zError();
	    
	    // fill kinematic info to tree
	    bsRootTree_->BsFitChi2_  = bs->chiSquared();
	    bsRootTree_->BsFitNdof_   =(int)bs->degreesOfFreedom();
	    bsRootTree_->BsFitVtxProb_ = Kfitter.getProb();
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
	    bsRootTree_->BsFitM_ = b_par[6];		
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

		//Onko tan just tarvetta olla tassa?
/*            if (mu1.isGlobalMuon()) bsRootTree_->BsMu1QualityG_=selGlobalMuon(mu1,RecVtx.position()); 
            if (mu2.isGlobalMuon()) bsRootTree_->BsMu2QualityG_=selGlobalMuon(mu2,RecVtx.position()); 
            if (mu1.isTrackerMuon()) bsRootTree_->BsMu1QualityT_=selTrackerMuon(mu1,RecVtx.position()); 
            if (mu2.isTrackerMuon()) bsRootTree_->BsMu2QualityT_=selTrackerMuon(mu2,RecVtx.position());
*/	    
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
      
            if (mu2.isTrackerMuon() && !mu2.isGlobalMuon())       bsRootTree_->JpsiMuon2Cat_alone_ = 1;
            else if (!mu2.isTrackerMuon() && mu2.isGlobalMuon())  bsRootTree_->JpsiMuon2Cat_alone_ = 2;
            else if (mu2.isTrackerMuon() && mu2.isGlobalMuon())   bsRootTree_->JpsiMuon2Cat_alone_ = 3;
            else if (!mu2.isTrackerMuon() && !mu2.isGlobalMuon()) bsRootTree_->JpsiMuon2Cat_alone_ = 4;
      
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
	   
	     //New part where the closest PV with respect to SV is selected 

/*	    Int_t ClosestPVindex = -1;
	    Double_t MinDistance = 10000000; 		
					
		//PV selection with minimum z distance between PV and SV
		for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++){
			const Vertex &vtx = (*recVtxs)[VtxInd];
			Double_t PVSVdistance = TMath::Abs(bVertex->position().z()-vtx.z());

			if(PVSVdistance < MinDistance){
				MinDistance = PVSVdistance;
				ClosestPVindex = VtxInd;	
			}			
		}*/

		//PV selection with miniminal angle between PV-SV vec and Bs pvec


		Int_t PVCosThetaIndex = -1;
		Double_t MinDistance = 10000000; 	
		
//		Int_t PVCosThetaIndex = -1;
//		MinDistance = 10000000; 	
		Bool_t isPVCosTheta = false;


		for(unsigned int VtxInd=0; VtxInd<recVtxs->size(); VtxInd++){
			const Vertex &vtx = (*recVtxs)[VtxInd];

			Double_t PVSVvecDotBsPvec = ( bVertex-> position().x()- vtx.x() )*Bsvec.x() + (bVertex-> position().y() - vtx.y())*Bsvec.y() + (bVertex-> position().z() - vtx.z() )*Bsvec.z();

			Double_t PVSVlength = TMath::Sqrt( pow((bVertex->position().x()- vtx.x()), 2.0) + pow((bVertex->position().y()- vtx.y()), 2.0) + pow((bVertex->position().z()- vtx.z()), 2.0) );
			Double_t BsPlength = TMath::Sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z());
					
			Double_t BsCosTheta = PVSVvecDotBsPvec / (BsPlength * PVSVlength);
			Double_t distance = 1-BsCosTheta;

				if(distance < MinDistance){
					MinDistance = distance;
					//cout << "mindist "<< MinDistance << endl;
					PVCosThetaIndex = VtxInd;
				}			
		}
			
	
		bsRootTree_->PVindex_ = PVCosThetaIndex;


  // JET Stuff
  if(jets.size()!=0){	

  bsRootTree_->BsJetPx_= jets[0].px();
  bsRootTree_->BsJetPy_= jets[0].py();
  bsRootTree_->BsJetPz_= jets[0].pz();
	
  Double_t dRBsAndJet = deltaR(Bsvec.eta(), Bsvec.phi(), jets[0].eta(), jets[0].phi() );	
  bsRootTree_->BsPtJetdR_= dRBsAndJet;
  bsRootTree_->BsJetCharge_ = jets[0].charge();
  }

  if(BjetIndex!=-1){
   	
   bsRootTree_->BsBJetdR_ = deltaR(Bsvec.eta(), Bsvec.phi(), jets[BjetIndex].eta(), jets[BjetIndex].phi() );
  }


  

/*		if(dRBsAndJet >2.95 && dRBsAndJet < 3.05 ){


		cout << "Bs px "<<Bsvec.x()<< endl;
		cout << "Bs py " <<Bsvec.y() << endl;
		cout << "Bs pz "<<Bsvec.z() << endl;

		cout << "jet px "<<jet.px() << endl;
		cout << "jet py " <<jet.py() << endl;
		cout << "jet pz "<<jet.pz() << endl;
		
		Double_t BsPdotJetP = Bsvec.x()*jet.px() + Bsvec.y()*jet.py() + Bsvec.z()*jet.pz();
		Double_t denominator = jet.p()*TMath::Sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z());

			if(denominator!=0) {
			Double_t CosAngleBsPJetP = BsPdotJetP/(denominator);	
			cout << "angle between jet p and Bs p " << 180.0/TMath::Pi()* TMath::ACos(CosAngleBsPJetP) << endl;
		
			}
		}*/
//	}


/*	if(PVCosThetaIndex != -1){
		isPVCosTheta = true;
		PVvtxCosTheta = (*recVtxs)[PVCosThetaIndex];
	
		cout <<"CosTheta PVx "<< PVvtxCosTheta.x() << endl;
		cout <<"CosTheta PVy "<< PVvtxCosTheta.y() << endl;
		cout <<"CosTheta PVz "<< PVvtxCosTheta.z() << endl;	
	}*/

			

		if(PVCosThetaIndex == -1){ //if there is no reconstructed PV's the PV is chosen to be the beam spot
			bsRootTree_->isBS_ = 1;
			PVx=BSx; 
    			PVy=BSy;
    			PVz=BSz;
    			PVerrx=BSdx;
    			PVerry=BSdy;
    			PVerrz=BSdz;

			reco::Vertex::Error err;
  			err(1,1)=BSdx;
  			err(2,2)=BSdy;
  			err(3,3)=BSdz;
								
			RecVtx = Vertex(reco::Vertex::Point(BSx,BSy,BSz),err); //Beam spot converted to vertex
    			NoPVcounter++;	
		}


		else{
			bsRootTree_->isPV_ = 1;
			RecVtx = (*recVtxs)[PVCosThetaIndex]; 
			PVx = RecVtx.x(); 
		        PVy = RecVtx.y();
		        PVz= RecVtx.z();
		        PVerrx=RecVtx.xError();
	 	        PVerry=RecVtx.yError();
      			PVerrz=RecVtx.zError();

			//Int_t TrackNro=-1;
			//if(PVCosThetaIndex<30){
			//TrackNro = bsRootTree_-> NTracksInPV_[PVCosThetaIndex];
			//cout << "Number of tracks in PV Bs "<<TrackNro << endl;
			//}
			
		/*	cout <<"cos theta PVx "<< PVx << endl;
			cout <<"cos theta PVy "<< PVy << endl;

			cout <<"cos theta PVz "<< PVz << endl;
		*/
		}
	
		bsRootTree_->getVtx(BSx,BSy,BSz,PVx,PVy,PVz,PVerrx,PVerry,PVerrz);

		bsRootTree_->BSdx_ = BSdx;
  		bsRootTree_->BSdy_ = BSdy;
  		bsRootTree_->BSdz_ = BSdz;
  		bsRootTree_->BSsigmaZ_ = BSsigmaZ;
  		bsRootTree_->BSdsigmaZ_ = BSdsigmaZ;

	
	//voiko tan siirtaa tahan rivilta 858?
	if (mu1.isGlobalMuon()) bsRootTree_->BsMu1QualityG_=selGlobalMuon(mu1,RecVtx.position()); 
            if (mu2.isGlobalMuon()) bsRootTree_->BsMu2QualityG_=selGlobalMuon(mu2,RecVtx.position()); 
            if (mu1.isTrackerMuon()) bsRootTree_->BsMu1QualityT_=selTrackerMuon(mu1,RecVtx.position()); 
            if (mu2.isTrackerMuon()) bsRootTree_->BsMu2QualityT_=selTrackerMuon(mu2,RecVtx.position());


	    // proper decay time and proper decay length without the refitted vertex 
	    
            BsLxyz=sqrt(pow(bVertex->position().x()-PVx,2)+pow(bVertex->position().y()-PVy,2)+pow(bVertex->position().z()-PVz,2));

	   // ctau 3D

	bsRootTree_->BsCt3D_ = 5.3667*( (bVertex->position().x()-PVx)*Bsvec.x() + (bVertex->position().y()-PVy)*Bsvec.y() + (bVertex->position().z()-PVz)*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );

	    // ctau 2d


	bsRootTree_->BsCt2D_ = 5.3667*( (bVertex->position().x()-PVx)*Bsvec.x() + (bVertex->position().y()-PVy)*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() );

	   // BsCtau 2D and 3D calculated for PV selected with CosTheta 	

/*	if(isPVCosTheta == true){

	bsRootTree_->BsCt3DPVCosTheta_ = 5.3667*( (bVertex->position().x()-PVvtxCosTheta.x())*Bsvec.x() + (bVertex->position().y()-PVvtxCosTheta.y())*Bsvec.y() + (bVertex->position().z()-PVvtxCosTheta.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );

	bsRootTree_->BsCt2DPVCosTheta_ = 5.3667*( (bVertex->position().x()-PVvtxCosTheta.x())*Bsvec.x() + (bVertex->position().y()-PVvtxCosTheta.y())*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() );

          }*/

	//BsCtau 2D and 3D calculated for highest pt PV

       if(isHighestPtPV == true){

	bsRootTree_->BsCt3DPVHighestPt_ = 5.3667*( (bVertex->position().x()-PVvtxHighestPt.x())*Bsvec.x() + (bVertex->position().y()-PVvtxHighestPt.y())*Bsvec.y() + (bVertex->position().z()-PVvtxHighestPt.z())*Bsvec.z() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y() + Bsvec.z()*Bsvec.z() );

	bsRootTree_->BsCt2DPVHighestPt_ = 5.3667*( (bVertex->position().x()-PVvtxHighestPt.x())*Bsvec.x() + (bVertex->position().y()-PVvtxHighestPt.y())*Bsvec.y() )/( Bsvec.x()*Bsvec.x() + Bsvec.y()*Bsvec.y());
      }


            // ctau 2d BS
            bsRootTree_->BsCt2DBS_ = 5.3667*( (bVertex->position().x()-BSx)*Bsvec.x() +(bVertex->position().y()-BSy)*Bsvec.y())/ (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());


	    
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
	    bsRootTree_->BsCtErr3D_ = b_par[6]*(error3D*abs(scale1))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
	    
	    // error on ctau 2D
	    VertexDistanceXY d2;
	    Measurement1D measurement2 = d2.distance(RecVtx,bVertex->vertexState());
	    double error2D = measurement2.error();
	    double scale2 = ((bVertex->position().x() - PVx)*Bsvec.x()+(bVertex->position().y() - PVy)*Bsvec.y())/
	      (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
	       sqrt((bVertex->position().x() - PVx)*(bVertex->position().x() - PVx)+(bVertex->position().y() - PVy)*(bVertex->position().y() - PVy)));
	    bsRootTree_->BsCtErr2D_ = b_par[6]*(error2D*abs(scale2))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            // error on ctau 2D BS
            VertexDistanceXY d2BS;
            Measurement1D measurement2BS = d2BS.distance(vertexBeamSpot,bVertex->vertexState());
            double error2DBS = measurement2BS.error();
            double scale2BS = ((bVertex->position().x() - BSx)*Bsvec.x()+(bVertex->position().y() - BSy)*Bsvec.y())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
               sqrt((bVertex->position().x() - BSx)*(bVertex->position().x() - BSx)+(bVertex->position().y() - BSy)*(bVertex->position().y() - BSy)));
	    bsRootTree_->BsCtErr2DBS_=5.3663*(error2DBS*abs(scale2BS))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());
 
	    // error on ctau 2D - 2 (first approximation)
	    bsRootTree_->BsCtErr2D2_ = sqrt((1/(Bsvec.perp()*Bsvec.perp()))*
					    (bs_er(1,1)*Bsvec.x()*Bsvec.x()+
					     bs_er(2,2)*Bsvec.y()*Bsvec.y()+
					     bs_er(1,2)*Bsvec.x()*Bsvec.y()));
	    
	    // error on ctau 3D MPV
	    bsRootTree_->BsCtErrMPV_ = b_par[6]/sqrt(deno);

	    //	    cout << "error no refit (3d,2d,MPV) = (" << bsRootTree_->BsCtErr3D_ << "," << bsRootTree_->BsCtErr2D_ << "," << bsRootTree_->BsCtErrMPV_ << ")" << endl;

	    // proper decay time and proper decay length with the refitted vertex
            // ctau 3D
            bsRootTree_->BsCt3Drefit_ = b_par[6]*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
						  +(bVertex->position().y()-reFitVertex.y())*Bsvec.y()+
						  +(bVertex->position().z()-reFitVertex.z())*Bsvec.z())/
              (Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
            // ctau 2d
            bsRootTree_->BsCt2Drefit_ = b_par[6]*((bVertex->position().x()-reFitVertex.x())*Bsvec.x()+
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
            bsRootTree_->BsCtErr3Drefit_ = b_par[6]*(error3D2*abs(scale12))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y()+Bsvec.z()*Bsvec.z());
	    
            // error on ctau 2D
            VertexDistanceXY d22;
            Measurement1D measurement22 = d22.distance(reFitVertex,bVertex->vertexState());
            double error2D2 = measurement22.error();
            double scale22 = ((bVertex->position().x() - reFitVertex.x())*Bsvec.x()+
			      (bVertex->position().y() - reFitVertex.y())*Bsvec.y())/
              (sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y())*
               sqrt((bVertex->position().x() - reFitVertex.x())*(bVertex->position().x() - reFitVertex.x())+
		    (bVertex->position().y() - reFitVertex.y())*(bVertex->position().y() - reFitVertex.y())));
            bsRootTree_->BsCtErr2Drefit_ = b_par[6]*(error2D2*abs(scale22))/sqrt(Bsvec.x()*Bsvec.x()+Bsvec.y()*Bsvec.y());

            // error on ctau 3D MPV
 	bsRootTree_->BsCtErrMPVrefit_ = b_par[6]/sqrt(deno2);

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

	    bsRootTree_->BsPVDist2d_     = vdist.distance(bVertex->vertexState(),RecVtx).value();
	    bsRootTree_->BsPVDist3d_     = vdist3d.distance(bVertex->vertexState(),RecVtx).value();
	    
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
	  }
	
	  
	} // trk2 loop
      } // trk1 loop
      
   
		//////////////// 
		// B+ meson
		/////////////// 

	
	Handle<CandidateView> allTracksK;
      		    iEvent.getByLabel(trackLabelK_, allTracksK);		    

		for(size_t itracks = 0; itracks < allTracksK->size(); itracks++){
			
			const Candidate & KplusTrack = (*allTracksK)[itracks];
			TrackRef trkKplusRef = KplusTrack.get<TrackRef>();
			if (KplusTrack.pt() < KaonTrackPtCut_) continue; 						
			if (KplusTrack.charge()!= 1) continue;

			TrackRef muonTrack1 = mu1.track();
      			TrackRef muonTrack2 = mu2.track();

	if( muonTrack1->momentum() == KplusTrack.momentum() && muonTrack1->charge() == KplusTrack.charge() ) {
//	cout << "Track overlap "<< muonTrack1->momentum() << " " << KplusTrack.momentum() << endl;	 
	continue; } 

	if( muonTrack2->momentum() == KplusTrack.momentum() && muonTrack2->charge() == KplusTrack.charge() ) {
//	cout << "Track overlap " << muonTrack2->momentum() << " " << KplusTrack.momentum() << endl;	 	
	continue; 
	}


		// loose cuts to reduce tracks 
	if(mu1.pt()<2.5 || mu2.pt()<2.5){ 
		//cout << "JEE mu1.pt() " << mu1.pt() <<endl;
		//cout << "JEE mu2.pt() " << mu2.pt() <<endl;
		continue;
	}
	if(KplusTrack.pt() < 1.5 ){
		//cout << "JEE KplusTrack.pt() " << KplusTrack.pt() <<endl;
		continue;
	} 
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
	 									
	Double_t vtxProbBplus = TMath::Prob(bplusmes->chiSquared(),(int)bplusmes->degreesOfFreedom());
	VertexDistanceXY vtxdist;

	if(bestVtxProbBplus < vtxProbBplus){ 
		bestVtxProbBplus = vtxProbBplus; 		
		AlgebraicVector7 bplus_par = bplusmes->currentState().kinematicParameters().vector();

		bsRootTree_->BplusM_fit_ = bplus_par[6];
		bsRootTree_->BplusVtxProb_ = vtxProbBplus;
		bsRootTree_->BplusChi2_ = bplusmes->chiSquared();
		bsRootTree_->BplusPt_ =  sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) );
		bsRootTree_->BplusPtot_ = sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) + pow(bplus_par[5],2.0) );
		bsRootTree_->KplusPt_ = KplusTrack.pt();
		bsRootTree_->KplusPtot_ = KplusTrack.p();
		bsRootTree_->BplusMu1Pt_ = mu1.pt(); 
		bsRootTree_->BplusMu2Pt_ = mu2.pt();
		  
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
		bsRootTree_->JpsiPt_bplus_ = sqrt( pow(jpsi_par[3],2.0) + pow(jpsi_par[4],2.0) );	

		GlobalPoint BplusVtxPos( bplusVertex->position().x(), bplusVertex->position().y(), bplusVertex->position().z());
		GlobalError BplusVtxErr( bplusVertex->error() );	

		//JET stuff
		if(jets.size()!=0){
		   
		 bsRootTree_->BplusJetPx_= jets[0].px();
		 bsRootTree_->BplusJetPy_= jets[0].py();	
		 bsRootTree_->BplusJetPz_= jets[0].pz();
		 bsRootTree_->BplusPtJetdR_ = deltaR(BplusVec.eta(), BplusVec.phi(), jets[0].eta(), jets[0].phi() );	
		
		
		
//		Double_t dRBpAndJet = deltaR(BplusVec.eta(), BplusVec.phi(), jet.eta(), jet.phi() );	
		//cout << "Bplus dR " << deltaR(BplusVec.eta(), BplusVec.phi(), jet.eta(), jet.phi() )  << endl;

 /* 		bsRootTree_->BplusJetCharge_ = jet.charge();
		cout << "Carhge " << jet.charge() << endl;
	 */
		}
		if(BjetIndex!=-1){
		 bsRootTree_->BplusBJetdR_ = deltaR(BplusVec.eta(), BplusVec.phi(), jets[BjetIndex].eta(), jets[BjetIndex].phi() );
		}	

/*		
		Double_t MytestJpsiP = sqrt( pow(jpsi_par[3],2.0) + pow(jpsi_par[4],2.0) +  pow(jpsi_par[5],2.0));	
		Double_t MytestMuP1dotP2 = mu1.px()*mu2.px() + mu1.py()*mu2.py() + mu1.pz()*mu2.pz();
		Double_t Mu12Ptot = sqrt(mu1.p()*mu1.p() + mu2.p()*mu2.p() + 2*MytestMuP1dotP2);

		cout <<  "My Bp Jpsi Ptot " <<MytestJpsiP << endl;
		cout <<  "My Bp mu12 Ptot " <<  Mu12Ptot << endl;
		cout <<  "My mu1 ptot " << mu1.p() << endl;
		cout <<  "My mu2 ptot "  << mu2.p() << endl;
*/				
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
					//cout << "mindist "<< MinDistance << endl;
					PVCosThetaIndex = BpVtxInd;
				}			
		}

		bsRootTree_->BplusPVindex_ = PVCosThetaIndex;
		if(PVCosThetaIndex == -1){ //if there is no reconstructed PV's the PV is chosen to be the beam spot
			bsRootTree_->BplusIsBS_ =1; 
			BplusPVx=BSx; 
    			BplusPVy=BSy;
    			BplusPVz=BSz;
    			BplusPVerrx=BSdx;
    			BplusPVerry=BSdy;
    			BplusPVerrz=BSdz;

			reco::Vertex::Error err;
  			err(1,1)=BSdx;
  			err(2,2)=BSdy;
  			err(3,3)=BSdz;
								
			BplusPVVtx = Vertex(reco::Vertex::Point(BSx,BSy,BSz),err); //Beam spot converted to vertex   			
		}


		else{
			bsRootTree_->BplusIsPV_ = 1; 
			BplusPVVtx = (*recVtxs)[PVCosThetaIndex]; 
			BplusPVx = BplusPVVtx.x(); 
		        BplusPVy = BplusPVVtx.y();
		        BplusPVz= BplusPVVtx.z();
		        BplusPVerrx=BplusPVVtx.xError();
	 	        BplusPVerry=BplusPVVtx.yError();
      			BplusPVerrz=BplusPVVtx.zError();
			
			//Int_t TrackNroBp=-1;
			//if(PVCosThetaIndex < 30){
			//TrackNroBp = bsRootTree_-> NTracksInPV_[PVCosThetaIndex];
			//cout << "Number of tracks in PV Bp "<<TrackNroBp << endl;
			//}	
			//cout << "PV index " << PVCosThetaIndex << endl;
		}

		double BplusLxy = vtxdist.distance( BplusPVVtx, bplusVertex->vertexState() ).value(); 
		bsRootTree_->BplusLxy_ = BplusLxy;
			
		double BplusLxyErr  = vtxdist.distance( BplusPVVtx, bplusVertex->vertexState() ).error();
		bsRootTree_->BplusLxyErr_ = BplusLxyErr;
				
				//bsRootTree_->BplusCt2D_ = (nominalBplusMass*(vtxdist.distance( RecVtx, bplusVertex->vertexState() ).value()) )/sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) );



		double BplusLxyz = sqrt( pow(bplusVertex->position().x()-BplusPVx,2) + pow(bplusVertex->position().y()-BplusPVy,2) + pow(bplusVertex->position().z()-BplusPVz,2) ); 

		bsRootTree_->BplusLxyz_ = BplusLxyz; 


		GlobalPoint BplusPVpos(BplusPVx,BplusPVy,BplusPVz);	
		GlobalError BplusPVerr(BplusPVVtx.error() );			
					
		VertexDistance3D bplusD;
		Measurement1D bplusMeasurement = bplusD.distance(VertexState(BplusVtxPos,BplusVtxErr),VertexState(BplusPVpos,BplusPVerr));

		 bsRootTree_->BplusLxyzErr_ = bplusMeasurement.error(); // error in the 3d distance between vertices

		double BplusMomlength = sqrt( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) + pow(bplus_par[5],2.0) ) ;
				//double BplusCt = (nominalBplusMass*BplusLxyz) / ( BplusMomlength );
				
				//bsRootTree_-> BplusCt3D_ = BplusCt;			
				// (PVBp - BpVtx) . pvec /()
		double SdotP = (bplusVertex-> position().x()- BplusPVx)*bplus_par[3] + (bplusVertex-> position().y() - BplusPVy)*bplus_par[4] + 
(bplusVertex-> position().z() - BplusPVz)*bplus_par[5];
								  			  		
		bsRootTree_-> BplusCosTheta_ = SdotP / (BplusMomlength * sqrt( pow((bplusVertex->position().x()- BplusPVx), 2.0) + pow((bplusVertex->position().y()- BplusPVy), 2.0) + pow((bplusVertex->position().z()- BplusPVz), 2.0) ) );
	
				
		double BplusCtErr3D =bplus_par[6]*((bplusMeasurement.error()*abs(bsRootTree_->BplusCosTheta_) ) / BplusMomlength );
		
		bsRootTree_->BplusCtErr3D_ = BplusCtErr3D; //error in 3d decay length

//AlgebraicSymMatrix77 bplus_er = bplusmes->currentState().kinematicParametersError().matrix();
				
		double BplusCt2D_2 = ( nominalBplusMass*((bplusVertex->position().x()- BplusPVx)*bplus_par[3] + (bplusVertex->position().y()- BplusPVy)*bplus_par[4]) )/( pow(bplus_par[3],2.0) + pow(bplus_par[4],2.0) );
 	
		bsRootTree_->BplusCt2D_2_ = BplusCt2D_2;
	
		double BplusCt3D_2 = nominalBplusMass*( (bplusVertex->position().x()- BplusPVx)*bplus_par[3] + (bplusVertex->position().y()- BplusPVy)*bplus_par[4] + (bplusVertex->position().z()- BplusPVz)*bplus_par[5] )/( BplusMomlength*BplusMomlength );
				
		bsRootTree_->BplusCt3D_2_ = BplusCt3D_2;
			
		VertexDistanceXY d2;
	    	Measurement1D measurement2 = d2.distance(RecVtx,bplusVertex->vertexState());

	    	double error2D = measurement2.error();
	   			
double scale2D = ((bplusVertex->position().x() - BplusPVx)*bplus_par[3] + (bplusVertex->position().y()- BplusPVy)*bplus_par[4])/ (sqrt( pow(bplus_par[3], 2.0) + pow(bplus_par[4], 2.0) )*sqrt((bplusVertex->position().x() - BplusPVx)*(bplusVertex->position().x() - BplusPVx) + (bplusVertex->position().y() - BplusPVy)*(bplusVertex->position().y() - BplusPVy)));

double BplusCtErr2D = bplus_par[6]*(error2D*abs(scale2D))/sqrt( pow(bplus_par[3], 2.0) + pow(bplus_par[4], 2.0) );	
	    
	bsRootTree_->BplusCtErr2D_ = BplusCtErr2D;
		
	bsRootTree_-> BplusCtPerSigma2D_ = ( BplusCt2D_2/BplusCtErr2D );
	bsRootTree_-> BplusCtPerSigma3D_ = ( BplusCt3D_2/BplusCtErr3D ); 
	bsRootTree_-> LxyPerSigma2D_ = ( BplusLxy/BplusLxyErr ); 
				
	std::pair<bool,Measurement1D> ImpactPar3DBplusK = IPTools::absoluteImpactParameter3D(TransTrackKplus, BplusPVVtx); 
	
	if(ImpactPar3DBplusK.first){
  		bsRootTree_->BplusKIP3D_ = ImpactPar3DBplusK.second.value();
  		bsRootTree_->BplusKIP3DErr_  = ImpactPar3DBplusK.second.error();
    	}

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

	bsRootTree_->BpmatchDoubleMu01DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
	bsRootTree_->BpmatchDoubleMu02DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();


	// deltaR matching for Bplus 

	bool BplusKTruth = MCmatchingBplusK(KplusTrack, genParticles, bsRootTree_->BplusKmcId_, bsRootTree_->BplusKmomId_,521);

	bool BplusMu1Truth = MCmatching(mu1, genParticles, bsRootTree_->BplusMu1mcId_, bsRootTree_->BplusMu1momId_,bsRootTree_->BplusMu1gmomId_,443,521);

	bool BplusMu2Truth= MCmatching( mu2,genParticles, bsRootTree_->BplusMu2mcId_,bsRootTree_->BplusMu2momId_,bsRootTree_->BplusMu2gmomId_, 443, 521);
	      
	if (BplusKTruth==1 && BplusMu1Truth==1 && BplusMu2Truth==1)  bsRootTree_->isMatchedBplus_ = 1;
	else bsRootTree_->isMatchedBplus_ = 0; 		
			} //end of if sentence
			
//B+ ID = 521
//mu- ID = 13 
//K+ ID = 321 

//B0 ID = 511
//J/Psi ID = 443
//K* ID = 313

	} //end of B+ loop     


	

      ////           Kstar
      Handle<CandidateView> allTracksPi;
      iEvent.getByLabel(trackLabelPi_, allTracksPi);
      //cout<< "INSIDE THE Bd K LOOPS!" << endl;
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
	      
	//cout << "Bd mass "<<BdCand.mass() << endl;

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

	//cout<< "INSIDE THE Bd LOOPS!" << endl;
	      
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
	
/*	 bsRootTree_->BdmatchDiMuon0Mu01_=!mu1.triggerObjectMatchesByFilter("HLT_Dimuon0_Jpsi_v*").empty();	
	 bsRootTree_->BdmatchDiMuon0Mu02_=!mu2.triggerObjectMatchesByFilter("HLT_Dimuon0_Jpsi_v*").empty();
*/
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

		bsRootTree_->BdmatchDoubleMu01DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();
		bsRootTree_->BdmatchDoubleMu02DiMuon0_ = !mu1.triggerObjectMatchesByFilter("hltVertexmumuFilterJpsiMuon").empty();



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


		Int_t BdPVCosThetaIndex = -1;
		Double_t MinDistance = 10000000; 	
		Double_t distance = 0;	

		for(unsigned int BdVtxInd=0; BdVtxInd<recVtxs->size(); BdVtxInd++){
			const Vertex &vtx = (*recVtxs)[BdVtxInd];

		Double_t PVSVvecDotBdPvec = ( bdVertex-> position().x()- vtx.x() )*Bdvec.x() + (bdVertex-> position().y() - vtx.y())*Bdvec.y() + (bdVertex-> position().z() - vtx.z() )*Bdvec.z();

		Double_t PVSVlengthBd = TMath::Sqrt( pow((bdVertex->position().x()- vtx.x()), 2.0) + pow((bdVertex->position().y()- vtx.y()), 2.0) + pow((bdVertex->position().z()- vtx.z()), 2.0) );

		Double_t BdPlength = TMath::Sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y() + Bdvec.z()*Bdvec.z());
					
		Double_t BdCosTheta = PVSVvecDotBdPvec / (BdPlength * PVSVlengthBd);
		distance = 1-BdCosTheta;

				if(distance < MinDistance){
					MinDistance = distance;
					//cout << "mindist Bd "<< MinDistance << endl;
					BdPVCosThetaIndex = BdVtxInd;
				}			
		}

		
		if(BdPVCosThetaIndex == -1){ //if there is no reconstructed PV's the PV is chosen to be the beam spot
			
			reco::Vertex::Error err;
  			err(1,1)=BSdx;
  			err(2,2)=BSdy;
  			err(3,3)=BSdz;
								
			BdPVVtx = Vertex(reco::Vertex::Point(BSx,BSy,BSz),err); //Beam spot converted to vertex   			
		}


		else{
	
			BdPVVtx = (*recVtxs)[BdPVCosThetaIndex]; 
			
			//cout << "Min dist " <<MinDistance << endl;		
			//Int_t TrackNroBp=-1;
			//if(PVCosThetaIndex < 30){
			//TrackNroBp = bsRootTree_-> NTracksInPV_[PVCosThetaIndex];
			//cout << "Number of tracks in PV Bp "<<TrackNroBp << endl;
			//}	
			//cout << "PV index " << PVCosThetaIndex << endl;
		}


		// end L1/HLT-reco matching
	      reco::Vertex reFitVertex;
	      //recalculate primary vertex without tracks from B
	      //reco::Vertex tmpFitVertex = reVertex(recVtxs, iEvent,iSetup, mu1, mu2, trkkst1, trkkst2);
	      //if(tmpFitVertex.isValid()) reFitVertex = tmpFitVertex;
//	      else reFitVertex = reco::Vertex(RecVtx);   // use the original vertex if the refit vertex is invalid
//	      reFitVertex = RecVtx;   // use the original vertex if the refit vertex is invalid
		
	//	cout << "refit x "<<reFitVertex.x() << endl;
	//	cout << "refit y " << reFitVertex.y() << endl;
	//	cout << "refit z " << reFitVertex.z() << endl;

	      bsRootTree_->BdPVx_refit_    = reFitVertex.x();
	      bsRootTree_->BdPVy_refit_    = reFitVertex.y();
	      bsRootTree_->BdPVz_refit_    = reFitVertex.z();

	      bsRootTree_->BdPVerrx_refit_ = reFitVertex.xError();
	      bsRootTree_->BdPVerry_refit_ = reFitVertex.yError();
	      bsRootTree_->BdPVerrz_refit_ = reFitVertex.zError();
			      


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
              bsRootTree_->BdCt2DBS_ = 5.2795*((bdVertex->position().x()-BSx)*Bdvec.x()+(bdVertex->position().y()-BSy)*Bdvec.y())/(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y());
              // error on ctau 2D BS
              VertexDistanceXY d2BS;
              Measurement1D measurement2BS = d2BS.distance(vertexBeamSpot,bdVertex->vertexState());
              double error2DBS = measurement2BS.error();
              double scale2BS = ((bdVertex->position().x() - BSx)*Bdvec.x()+(bdVertex->position().y() - BSy)*Bdvec.y())/(sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y())*sqrt((bdVertex->position().x() - BSx)*(bdVertex->position().x() - BSx)+(bdVertex->position().y() - BSy)*(bdVertex->position().y() - BSy)));
              bsRootTree_->BdCtErr2DBS_=5.2795*(error2DBS*abs(scale2BS))/sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y());

	
	   // ctau 3D

	bsRootTree_->BdCt3D_ = 5.2795*( (bdVertex->position().x()-BdPVVtx.x())*Bdvec.x() + (bdVertex->position().y()-BdPVVtx.y())*Bdvec.y() + (bdVertex->position().z()-BdPVVtx.z())*Bdvec.z() )/( Bdvec.x()*Bdvec.x() + Bdvec.y()*Bdvec.y() + Bdvec.z()*Bdvec.z() );

	    // ctau 2d

	bsRootTree_->BdCt2D_ = 5.2795*( (bdVertex->position().x()-BdPVVtx.x())*Bdvec.x() + (bdVertex->position().y()-BdPVVtx.y())*Bdvec.y() )/( Bdvec.x()*Bdvec.x() + Bdvec.y()*Bdvec.y() );


	// error on ctau 2D BdPVVtx
              VertexDistanceXY d2BdPV;
              Measurement1D measurement2BdPV = d2BdPV.distance(BdPVVtx,bdVertex->vertexState());
              double error2DBdPV = measurement2BdPV.error();
              double scale2BdPV = ((bdVertex->position().x() - BdPVVtx.x())*Bdvec.x()+(bdVertex->position().y() - BdPVVtx.y())*Bdvec.y())/( sqrt( Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y() )* sqrt( (bdVertex->position().x() - BdPVVtx.x())*(bdVertex->position().x() - BdPVVtx.x())+(bdVertex->position().y() - BdPVVtx.y())*(bdVertex->position().y() - BdPVVtx.y()) ) );
              bsRootTree_->BdCtErr2D_=5.2795*(error2DBdPV*abs(scale2BdPV))/sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y());

	// error on ctau 3D BdPVVtx

	GlobalPoint SVpos( bdVertex->position().x(), bdVertex->position().y(), bdVertex->position().z());
	    GlobalPoint PVpos( BdPVVtx.x(), BdPVVtx.y(), BdPVVtx.z());
	    GlobalError SVerr( bdVertex->error() );
	    GlobalError PVerr( BdPVVtx.error() );
	    VertexDistance3D d1;
	    Measurement1D measurement = d1.distance(VertexState(SVpos,SVerr),VertexState(PVpos,PVerr));
	    double error3D = measurement.error();
	    double scale1 = ((bdVertex->position().x() - BdPVVtx.x())*Bdvec.x()+
			     (bdVertex->position().y() - BdPVVtx.y())*Bdvec.y()+
			     (bdVertex->position().z() - BdPVVtx.z())*Bdvec.z())/
	      (sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y()+Bdvec.z()*Bdvec.z())*
	       sqrt((bdVertex->position().x() - BdPVVtx.x())*(bdVertex->position().x() - BdPVVtx.x()) + (bdVertex->position().y() - BdPVVtx.y())*(bdVertex->position().y() - BdPVVtx.y()) + (bdVertex->position().z() - BdPVVtx.z())*(bdVertex->position().z() - BdPVVtx.z() )));

	    bsRootTree_->BdCtErr3D_ = 5.2795*(error3D*abs(scale1))/sqrt(Bdvec.x()*Bdvec.x()+Bdvec.y()*Bdvec.y()+Bdvec.z()*Bdvec.z());


// 	      if(BdCand.pt()!=0) {
//                 bsRootTree_->BdLxy_ = ((bdVertex->position().x()-PVx)*Bdvec.x()+(bdVertex->position().y()-PVy)*Bdvec.y())/Bdvec.perp();
//                 bsRootTree_->BdCt_  = bsRootTree_->BdLxy_*Bdmass/Bdvec.perp();
//               }

              bsRootTree_->BdErrX_  = bd_er(1,1);
              bsRootTree_->BdErrY_  = bd_er(2,2);
              bsRootTree_->BdErrXY_ = bd_er(1,2); 
              
              VertexDistance3D vdist3d;
              bsRootTree_->BdDist3d_    = vdist3d.distance(bdVertex->vertexState(),BdPVVtx).value();
              bsRootTree_->BdDist3dErr_ = vdist3d.distance(bdVertex->vertexState(),BdPVVtx).error();
              bsRootTree_->BdTime3d_    = bsRootTree_->BdDist3d_    * Bdmass/Bdvec.perp() * 100. /3.;
              bsRootTree_->BdTime3dErr_ = bsRootTree_->BdDist3dErr_ * Bdmass/Bdvec.perp() * 100. /3.;
              
                       
              bsRootTree_->BdDist2d_     = vdist.distance(bdVertex->vertexState(),BdPVVtx).value();
              bsRootTree_->BdDist2dErr_ = vdist.distance(bdVertex->vertexState(),BdPVVtx).error();
              bsRootTree_->BdTime2d_     = bsRootTree_->BdDist2d_ * Bdmass/Bdvec.perp() *100. /3.;
              bsRootTree_->BdTime2dErr_  = bsRootTree_->BdDist2dErr_ * Bdmass/Bdvec.perp() * 100. /3.;
 

 //JET STUFF       
/* if(jets->size()!=0){	
  const PFJet &jet = (*jets)[0];
  bsRootTree_->BdJetPx_= jet.px();
  bsRootTree_->BdJetPy_= jet.py();
  bsRootTree_->BdJetPz_= jet.pz();
	
  Double_t dRBdAndJet = deltaR(Bdvec.eta(), Bdvec.phi(), jet.eta(), jet.phi() );	
  bsRootTree_->BdJetdR_= dRBdAndJet;
  bsRootTree_->BdJetCharge_ = jet.charge();
  }

*/		



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
  
      
    } // loop on muons2
  } // loop on muons1

bool isBplusStudy = true;
bool isBdStudy = false;

Int_t NumMu = 0;
bool nomatchflag = false;

edm::ESHandle<TransientTrackBuilder> ttrackBuilder;
iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);
 
for(size_t iter=0; iter < allmuons->size(); ++iter){ 
    const pat::Muon & tagMuon = (*allmuons)[iter];
    if( !tagMuon.isPFMuon() ) {continue;}  //tag muon MUST BE PF muon!

 if(isBplusStudy == true){

	if (tagMuon.pt()== bsRootTree_->BplusMu1Pt_  || tagMuon.pt()== bsRootTree_->BplusMu2Pt_ ){ 
	//cout<<"Muon overlap"<<endl; 
	continue;
        }
    
        if (tagMuon.pt() == bsRootTree_->KplusPt_ ){ 
	//cout<<"Kaon overlap"<<endl; 
	continue; 
	}

 } 	

else if(isBdStudy == true){
	if (tagMuon.pt() == bsRootTree_->BdMu1Pt_beffit_ || tagMuon.pt() == bsRootTree_->BdMu2Pt_beffit_){continue;}

	
	if (tagMuon.pt() == bsRootTree_->BdK1Pt_nofit_ || tagMuon.pt() == bsRootTree_->BdK2Pt_nofit_){ continue;}
 }


 else{

    if(tagMuon.pt()== bsRootTree_->JpsiMu1Pt_alone_  || tagMuon.pt()== bsRootTree_->JpsiMu2Pt_alone_ ){ 
	//cout<<"Muon overlap"<<endl; 
	continue;
    }    
   
    if(tagMuon.pt()==bsRootTree_->K1Pt_nofit_ || tagMuon.pt()==bsRootTree_->K2Pt_nofit_){ 
	//cout<<"Kaon overlap"<<endl; 
	continue;
    }
 }
	
	Double_t MuonIP3D = -999;
	Double_t MuonIP3DErr = -999; 

//  if(tagMuon.isTrackerMuon() || tagMuon.isGlobalMuon()){

    

//    if(tagMuon.isTrackerMuon()) bsRootTree_->TrackerTagMuon_[NumMu] = 1;
//    if(tagMuon.isGlobalMuon()) bsRootTree_->GlobalTagMuon_[NumMu] = 1;

    if(tagMuon.isPFMuon()) bsRootTree_->PFTagMuon_[NumMu] = 1;
 	 
    if(tagMuon.innerTrack().isNonnull()){	
    	TrackRef trkTagMuRef = tagMuon.get<TrackRef>();
    	TransientTrack trkTagMuTT = (*ttrackBuilder).build(&trkTagMuRef);

	std::pair<bool,Measurement1D> ImpactPar3D;

    if(isBplusStudy==true){

    	ImpactPar3D = IPTools::absoluteImpactParameter3D(trkTagMuTT, BplusPVVtx); 	
    }	

    
    if(isBdStudy==true){

    	ImpactPar3D = IPTools::absoluteImpactParameter3D(trkTagMuTT, BdPVVtx); 	
    }	

    else{
	    		  
  	ImpactPar3D = IPTools::absoluteImpactParameter3D(trkTagMuTT, RecVtx); 	
	}

    	if(ImpactPar3D.first){
    		MuonIP3D = ImpactPar3D.second.value();
    		MuonIP3DErr = ImpactPar3D.second.error();
    	}	
    } 

    bsRootTree_->TagMuRecoIP3D_[NumMu] = MuonIP3D;
    bsRootTree_->TagMuRecoIP3DErr_[NumMu] = MuonIP3DErr;	
     
    bsRootTree_->TagMuRecoPt_[NumMu] = tagMuon.pt();
    bsRootTree_->TagMuRecoP_[NumMu] = tagMuon.p();
    bsRootTree_->TagMuRecoEta_[NumMu] = tagMuon.eta();
    bsRootTree_->TagMuRecoPhi_[NumMu] = tagMuon.phi();
    bsRootTree_->TagMuRecoChg_[NumMu] = tagMuon.charge();
   
 double dRjetMin = 9999;
 double closestJetIndex = -1;

  
 for(size_t k=0; k < jets.size(); k++){
  double dR = deltaR(jets[k].eta(), jets[k].phi(), tagMuon.eta(), tagMuon.phi() );
	
//   if(jets[k].isPFJet() == true ){ cout << "We have a PF jet! " << endl;}
 //  if(jets[k].isCaloJet() == true) {cout << "We have a Calo jet! " << endl;}
//   if(jets[k].isBasicJet() == true) {cout << "We have a Basic jet! " << endl;} 

	if(dR < dRjetMin ){
	  dRjetMin = dR;
	  closestJetIndex = k;	  		
	}

 } //end of jet loop
	
 // cout << "closest jet wrt tag muon " <<dRjetMin << endl;

   if(closestJetIndex !=-1){
//    const PFJet &closestjet = (*jets)[closestJetIndex];
 	
    Double_t PjetDotPmu = ( tagMuon.px()*jets[closestJetIndex].px() + tagMuon.py()*jets[closestJetIndex].py() + tagMuon.pz()*jets[closestJetIndex].pz() ) / ( jets[closestJetIndex].px() * jets[closestJetIndex].px() + jets[closestJetIndex].py() * jets[closestJetIndex].py() + jets[closestJetIndex].pz() * jets[closestJetIndex].pz() );

    Double_t MuonPtRel = ( tagMuon.px() - PjetDotPmu * jets[closestJetIndex].px() )* ( tagMuon.px() - PjetDotPmu * jets[closestJetIndex].px() ) + ( tagMuon.py() - PjetDotPmu * jets[closestJetIndex].py() ) * ( tagMuon.py() - PjetDotPmu * jets[closestJetIndex].py() ) + ( tagMuon.pz() - PjetDotPmu * jets[closestJetIndex].pz() ) * ( tagMuon.pz() - PjetDotPmu * jets[closestJetIndex].pz() ) ;

  bsRootTree_->TagMuRecoPtRel_[NumMu] = TMath::Sqrt(MuonPtRel);


 //   cout << "dR jet and tag muon " << deltaR(tagMuon.eta(),tagMuon.phi(), jets[closestJetIndex].eta(),jets[closestJetIndex].phi()) << endl;    
    
 //   bsRootTree_->TagMuJetdR_[NumMu] = deltaR(tagMuon.eta(),tagMuon.phi(), jet.eta(),jet.phi());	
 //  cout << "muon pt rel "<< TMath::Sqrt(MuonPtRel) << endl;
	
   }	
   	
	//pdgid (mu-) = 13, (mu+) = -13
    double MinDRK = 999;
    int MinDRKId = -1;
    double closestMudR = 0.15;
    vector<int> closestParticles;
    vector<double> closestParticlePtDiffs;
    vector<double> closestMuPtDiffs;
    vector<double> minDeltaRs;
    int closestMuIndex = -1;
    int closestParticleIndex = -1;

    if (isMCstudy_){
//	cout << "tag muon pt " << tagMuon.pt() << " NumMu=" << NumMu << "isPF muon" << bsRootTree_->PFTagMuon_[NumMu] << endl;
//	cout << "gen particle size " <<genParticles->size() << endl;
      for(size_t i = 0; i < genParticles->size(); ++ i){
        const GenParticle & p = (*genParticles)[i];
        double DeltaRK1 = deltaR(p.eta(), p.phi(), tagMuon.eta(), tagMuon.phi() );
	double ptdiff = abs(tagMuon.pt() - p.pt());
	//	cout << "dR: " << DeltaRK1 << "  pt diff: " << ptdiff << endl;
	
	if(DeltaRK1 < MinDRK){MinDRK = DeltaRK1; MinDRKId = p.pdgId(); }

	if(DeltaRK1 < 0.15){
		closestParticles.push_back(i);
		closestParticlePtDiffs.push_back(ptdiff);
		minDeltaRs.push_back(DeltaRK1);	
				
	}
      } //end of gen particle loop


      double MinPtdiff = 9999; 
      double MindR=-5;		
      for(size_t k=0; k< closestParticles.size(); k++){

//	 cout << " min dR: " << minDeltaRs[k] << " Min Pt diff: " << closestParticlePtDiffs[k] << endl;
	 if(closestParticlePtDiffs[k] < MinPtdiff){
		MinPtdiff = closestParticlePtDiffs[k];
		MindR = minDeltaRs[k];
		closestParticleIndex = closestParticles[k];
	 } 
      }	

	
      if(closestParticleIndex!=-1){

   //     cout << "closest index: " << closestParticleIndex << " min dR: " << MindR << " Min Pt diff: " << MinPtdiff << endl;

	const GenParticle & ClosestP = (*genParticles)[closestParticleIndex];

	bsRootTree_->TagMuSimu_[NumMu]=ClosestP.pdgId();

//	cout << "closest particle's ID " << ClosestP.pdgId() << endl;

	if(ClosestP.mother()!=0){
	   bsRootTree_->MuRecoMCmother_[NumMu]=ClosestP.mother()->pdgId();
	}

	if(ClosestP.mother(0)->mother(0)!=0){
	   bsRootTree_->MuRecoMCMotherMother_[NumMu]=ClosestP.mother(0)->mother(0)->pdgId(); 
	}		
		
      }	
	else { if(NumMu == 0 && tagMuon.pt() > 2.9 ){mcNomatch_counter_++ ;}
//	cout << "Event without MC matching! Closest dR " << MinDRK << " particle ID " << MinDRKId << endl; 
	 }	


/////////////////// old MC matching 
/*	
	double MinDR = 0.15;
     for(size_t i = 0; i < genParticles->size(); ++ i){
        const GenParticle & p = (*genParticles)[i];
        double DeltaRK1 = deltaR(p.eta(), p.phi(), tagMuon.eta(), tagMuon.phi() );

	if(DeltaRK1 < MinDRK){MinDRK = DeltaRK1; MinDRKId = p.pdgId(); }

	//find closest muon
        if(DeltaRK1<closestMudR && abs(p.pdgId())==13 ){ 	
		closestMuIndex = i;
		closestMudR = DeltaRK1;
	   	
	}
	//find closest particle 
        if(DeltaRK1 < MinDR){
          MinDR = DeltaRK1;
          closestParticleIndex = i;	   	          	  	
        }	
      } //end of gen particle loop

	if(closestParticleIndex!=-1){
	const GenParticle & ClosestP = (*genParticles)[closestParticleIndex];

		//if closest particle is gen muon, save gen muon info
		if( abs(ClosestP.pdgId() ) == 13 ){
		  bsRootTree_->TagMuSimu_[NumMu]=ClosestP.pdgId();

		  if(ClosestP.mother()!=0) {bsRootTree_->MuRecoMCmother_[NumMu]=ClosestP.mother()->pdgId();}

		  if(ClosestP.mother(0)->mother(0)!=0) {bsRootTree_->MuRecoMCMotherMother_[NumMu]=ClosestP.mother(0)->mother(0)->pdgId(); }

		cout << "closest muon dR " << MinDRK << endl;  		
		}

		//if closest particle is not muon, but gen muon is still found nearby reco muon, save gen muon info
		else if(closestMuIndex != -1){
		  const GenParticle & ClosestM = (*genParticles)[closestMuIndex];
		  bsRootTree_->TagMuSimu_[NumMu]=ClosestM.pdgId();

		  if(ClosestM.mother()!=0) {bsRootTree_->MuRecoMCmother_[NumMu]=ClosestM.mother()->pdgId();}

		  if(ClosestM.mother(0)->mother(0)!=0) {bsRootTree_->MuRecoMCMotherMother_[NumMu]=ClosestM.mother(0)->mother(0)->pdgId(); }
		 cout << "closest muon dR " << MinDRK << endl; 		
		}

		//gen muon misidentified
		else{

		 bsRootTree_->TagMuSimu_[NumMu]=ClosestP.pdgId();

		  if(ClosestP.mother()!=0) {bsRootTree_->MuRecoMCmother_[NumMu]=ClosestP.mother()->pdgId();}

		  if(ClosestP.mother(0)->mother(0)!=0) {bsRootTree_->MuRecoMCMotherMother_[NumMu]=ClosestP.mother(0)->mother(0)->pdgId(); }		
		cout << "closest mis ID'd particle dR " << MinDRK << endl; 
		}
	} 
	else {mcNomatch_counter_++ ; cout << "closest no matched particle dR " << MinDRK << endl; }
*/
///////////////////// old MC matching

    } // end of if(MCstudy_) 

    NumMu=NumMu+1;	   
  
} //end of tag muon loop

bsRootTree_->TagMuListSize_ = NumMu;
tagmuonNum_ = tagmuonNum_ + NumMu+1;

// cout << "Number of PF tag muons " << NumMu << endl;			
 // cout <<"Number of events with no PV " << NoPVcounter << endl;
  //if ( MinBVtxHyp1>0 || minVtxP>0 )
  //if ((isMCstudy_ && (bsRootTree_->ChannelID_==1 || bsRootTree_->ChannelID_==2 || bsRootTree_->ChannelID_==3 ||  bsRootTree_->ChannelID_==4)) || !isMCstudy_)
    
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
 // cout<<"NewEvent-------"<<endl;
  for( size_t i = 0; i < genParticles->size(); ++ i ) {
    const GenParticle & genBsCand = (*genParticles)[ i ];
    int MC_particleID=genBsCand.pdgId();
    int absMC_particleID = abs(MC_particleID);


	//B+ meson decay channel ID
    bool isPosMu=0, isNegMu=0, isJpsi=0, isKplus=0, isNegPi=0, isPosPi=0, isKstar = 0, isNeutralPi=0; 
   		
			
    if( MC_particleID == 521 ){	
	
	int numBplusDaughters = genBsCand.numberOfDaughters();
	int jpsiIndex = -5;
	int kstarIndex = -5;	
	//checking of decay channels B+ -> Jpsi(mu+,mu-) K+ and B+ -> Jpsi K*(K+, pi0)  
	
	if(numBplusDaughters == 2){ //check if Bplus has 2 daughters

			for(int k = 0; k < numBplusDaughters; k++ ){ //check if the daughter's are J/psi and (K+ or K*)
				if( genBsCand.daughter(k)->pdgId() == 443){ 
					isJpsi = 1;
					jpsiIndex = k;
						 
				}  
				else if( genBsCand.daughter(k)->pdgId() == 321 ){ isKplus = 1; }
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
	}// end of B+ -> Jpsi(mu+,mu-) K+ channel search					
	
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



 //End of B+ meson decay channel ID



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

        if (MC_particleID==(-1)*genBsCand.daughter(0)->pdgId()) {
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
        if (MC_particleID == 511) bsRootTree_->BdIniFlavour_=1;
        if (MC_particleID == -511) bsRootTree_->BdIniFlavour_=-1;
      }
      
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
  
	

        if (MC_particleID==(-1)*genBsCand.daughter(0)->pdgId()) { 
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
        if (MC_particleID == 531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) bsRootTree_->BsIniFlavour_=1;
        if (MC_particleID == -531 && MCphi && (MCpsi || MCpsi2 || MCpsi3 || MCpsi4) && MCbs) bsRootTree_->BsIniFlavour_=-1;
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
	
          bsRootTree_->BsCt2DMC_ = Lxy2D * 5.3667/sqrt(bsmomx*bsmomx+bsmomy*bsmomy);
	  bsRootTree_->BsPtMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy);	
          // calculate gen ctau 3D
          //double Lxy3D = ((bssvx-bspvx)*bsmomx+(bssvy-bspvy)*bsmomy+(bssvz-bspvz)*bsmomz);
          double Lxy3D = sqrt(pow(bssvx-bspvx,2)+pow(bssvy-bspvy,2)+pow(bssvz-bspvz,2));
          bsRootTree_->BsCt3DMC_ = Lxy3D * 5.3667/sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);
	  bsRootTree_->BsLxy3DMC_ = Lxy3D;
	  bsRootTree_->BsPMC_ = TMath::Sqrt(bsmomx*bsmomx+bsmomy*bsmomy+bsmomz*bsmomz);

	cout<<"Line 3455 "<<endl;

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
              cout<<"angle_costhetaMC"<<angle_costhetaMC<<endl; 
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

        bsRootTree_->BmesonsId_[arrayIndex] = MC_particleID;
		
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
	  trackOK &&
         gTrack->chi2()/gTrack->ndof() < 20.0 &&
          q.numberOfValidMuonHits() > 0 &&
  	  iTrack->chi2()/iTrack->ndof() < 1.8 &&
  	  aMuon.muonID("TrackerMuonArbitrated") &&
  	  aMuon.muonID("TMOneStationTight") &&
          p.pixelLayersWithMeasurement() > 1 &&
	  fabs(iTrack->dxy(RefVtx)) < 3.0 &&
          fabs(iTrack->dz(RefVtx)) < 15.0);

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



bool BsToJpsiPhiAnalysis::MCmatching(const Candidate & track1,  edm::Handle<GenParticleCollection> & genParticles,
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
      if(p.mother()!=0 && p.mother()->mother()!=0) K1gmomId=p.mother()->mother()->pdgId(); 
      if (abs(K1momId)==condMom && abs(K1gmomId)==condGMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }

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
      if(p.mother()!=0) K1momId=p.mother()->pdgId(); 
      if (abs(K1momId)==condMom) K1Truth = 1;
      else K1Truth = 0;
    }
  }

  return K1Truth;
  
}


	    
reco::Vertex BsToJpsiPhiAnalysis::reVertex(const edm::Handle<reco::VertexCollection> &handle, const edm::Event &iEvent, const edm::EventSetup& iSetup,
					   pat::Muon mu1,pat::Muon mu2,TrackRef trk3,TrackRef trk4){
//OLD					   TrackRef trk1, TrackRef trk2, TrackRef trk3, TrackRef trk4){


  //copied from Onia2MuMu/VertexReProducer
  const edm::Provenance *prov = handle.provenance();
  if (prov == 0) throw cms::Exception("CorruptData") << "Vertex handle doesn't have provenience.";
  edm::ParameterSetID psid = prov->psetID();
  
  edm::pset::Registry *psregistry = edm::pset::Registry::instance();
  edm::ParameterSet psetFromProvenance;
  if (!psregistry->getMapped(psid, psetFromProvenance)) 
    throw cms::Exception("CorruptData") << "Vertex handle parameter set ID id = " << psid;
  
  if (prov->moduleName() != "PrimaryVertexProducer") 
    throw cms::Exception("Configuration") << "Vertices to re-produce don't come from a PrimaryVertexProducer, but from a " << prov->moduleName() <<".\n";

	  const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mu1.originalObject());
	  const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mu2.originalObject());
 
  edm::InputTag tracksTag_   = psetFromProvenance.getParameter<edm::InputTag>("TrackLabel");
  edm::InputTag   beamSpotTag_ = psetFromProvenance.getParameter<edm::InputTag>("beamSpotLabel");
  
  Handle<reco::TrackCollection> pvtracks;   iEvent.getByLabel(tracksTag_,   pvtracks);
  Handle<reco::BeamSpot>        pvbeamspot; iEvent.getByLabel(beamSpotTag_, pvbeamspot);
  
  TrackCollection muonLess;
  muonLess.reserve(pvtracks->size()-4);

  for (size_t i = 0; i < pvtracks->size(); ++i) {
    if (i == rmu1->track().key()){
      double dpt = rmu1->pt() - (*pvtracks)[i].pt();   double deta =  rmu1->eta() - (*pvtracks)[i].eta(); double dphi= rmu1->phi() - (*pvtracks)[i].phi();
      double delta = dpt*dpt + deta*deta + dphi*dphi;
      if( delta > 0.0001){
	std::cout<<"BsToJpsiPhiAnalysis::reVertex: ERROR" << std::endl;
	exit(1);
      }
    } 
    if (i == rmu2->track().key()){
      double dpt = rmu2->pt() - (*pvtracks)[i].pt();   double deta =  rmu2->eta() - (*pvtracks)[i].eta(); double dphi= rmu2->phi() - (*pvtracks)[i].phi();
      double delta = dpt*dpt + deta*deta + dphi*dphi;
      if( delta > 0.0001){
	std::cout<<"BsToJpsiPhiAnalysis::reVertex: ERROR" << std::endl;
	exit(1);
      }      
    } 
    if (i == trk3.key()){
      double dpt = trk3->pt() - (*pvtracks)[i].pt();   double deta =  trk3->eta() - (*pvtracks)[i].eta(); double dphi= trk3->phi() - (*pvtracks)[i].phi();
      double delta = dpt*dpt + deta*deta + dphi*dphi;
      if( delta > 0.0001){
	std::cout<<"BsToJpsiPhiAnalysis::reVertex: ERROR" << std::endl;
	exit(1);
      }
    } 
  }


 for (size_t i = 0; i < pvtracks->size(); ++i) {
   if (i == rmu1->track().key()) continue;
   if (i == rmu2->track().key()) continue;
   if (i == trk3.key()) continue;
   if (i == trk4.key()) continue;
   muonLess.push_back((*pvtracks)[i]);
 }


 edm::ESHandle<TransientTrackBuilder> theB;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
 
 vector<reco::TransientTrack> t_tks; t_tks.reserve(muonLess.size());
 
 for (reco::TrackCollection::const_iterator it = muonLess.begin(); it != muonLess.end(); ++it) {
   t_tks.push_back((*theB).build(*it));
   t_tks.back().setBeamSpot(*pvbeamspot);
 }

 vector<TransientVertex> pvs = PrimaryVertexProducerAlgorithm(psetFromProvenance).vertices(t_tks, *pvbeamspot);
 
 if(pvs.size() > 0) return reco::Vertex(pvs.front());
 return reco::Vertex();
}
