#ifndef HeavyFlavorAnalysis_BsToJpsiPhi_BsToJpsiPhiRootTree_h
#define HeavyFlavorAnalysis_BsToJpsiPhi_BsToJpsiPhiRootTree_h

#include <string>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>

class BsToJpsiPhiRootTree {
public:
  
	BsToJpsiPhiRootTree();
	
	~BsToJpsiPhiRootTree();
	
	void resetEntries(); 
	void writeFile();
	void createTree(const std::string filename);

	// read tree from single file
	void readTree(const std::string filename); 

	// read tree as chain from multiple files
	void readTree(std::vector<std::string> filenames);

 
	void getAngles(const double aa, const double bb, const double cc, const double dd);
	void getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff, const double gg, 
		    const double hh, const double ii);


	void getDeDx(const double f1, const double f2, const int f3);
  
	void fill();  

	void setBranchAddresses();
public:
        int BplusCharge_;
        double BplusMu1Eta_;
        double BplusMu2Eta_;
        double JpsiPtbplus_fit_;
        double BpTagIP_;
        double BpTagPt_;
        int BpTagCharge_;
        double BpTagEta_;
        double BpTagPhi_;
        double BpTagP_;
        int BpTagPF_;
        int BpTagGENID_;
		  int BpTagMomGENID_ ;
		  int BpTagGmomGENID_ ;	

		  int TrackMultiplicity_ ;
		  int TrackMultiplicityBp_ ;
		  int TrackMultiplicityBd_ ;	
		  int	MuonMultiplicity_ ;
		  int ElectronMultiplicity_;

		  double BpTagEIDNonTrig_;
  		  double BpTagEIDTrig_; 
  		  double BpTagEIDTrigNoIP_; 	
		  double BpTagEIP_;
		  double BpTagEPt_;
        int BpTagECharge_;
        double BpTagEEta_;
        double BpTagEPhi_;
        double BpTagEP_;
        int BpTagEGENID_;
		  int BpTagEGmomGENID_;
		  int BpTagEMomGENID_;		

		  double BsTagEIDNonTrig_;
  		  double BsTagEIDTrig_; 
  		  double BsTagEIDTrigNoIP_;
		  double BsTagEIP_;
		  double BsTagEPt_;
        int BsTagECharge_;
        double BsTagEEta_;
        double BsTagEPhi_;
        double BsTagEP_;
        int BsTagEGENID_;
		  int BsTagEGmomGENID_;
		  int BsTagEMomGENID_;	


		  double BdTagEIDNonTrig_;
  		  double BdTagEIDTrig_; 
  		  double BdTagEIDTrigNoIP_; 	
		  double BdTagEIP_;
		  double BdTagEPt_;
        int BdTagECharge_;
        double BdTagEEta_;
        double BdTagEPhi_;
        double BdTagEP_;
        int BdTagEGENID_;
		  int BdTagEGmomGENID_;
		  int BdTagEMomGENID_;	

        double BsTagIP_;
        double BsTagPt_;
        int BsTagCharge_;
        double BsTagEta_;
        double BsTagPhi_;
        double BsTagP_;
        int BsTagPF_;
        int BsTagGENID_;
		  int BsTagMomGENID_ ;
		  int BsTagGmomGENID_ ;	

        double BdTagIP_;
        double BdTagPt_;
        int BdTagCharge_;
        double BdTagEta_;
        double BdTagPhi_;
        double BdTagP_;
        int BdTagPF_;
        int BdTagGENID_;
		  int BdTagMomGENID_ ;
		  int BdTagGmomGENID_ ;	
		

        int  BdMuonCat1_;
        int  BdMuonCat2_;
        double  BplusM_fit_;
        double  BplusVtxProb_;
        double  BplusChi2_;
        double  BplusPt_;
        double  BplusPtot_;
        double  KplusPt_;
        double  KplusPtot_;
        double  BplusMu1Pt_;
        double  BplusMu2Pt_;
        double  BpJpsiVtxProb_;
        double  BpCosDeltaAlpha_;
        double  BpMuMuDCA_;
        double  BpMuMuDistance_;
        double  BpMuMuDistanceSigma_;
        double  BpMuDr1_;
        double  BpMuDr2_;
        int  BpMuonCat1_;
        int  BpMuonCat2_;
        double  BplusMu1Ptot_;
        double  BplusMu2Ptot_;
        double  BplusEta_;
        double  BplusPhi_;
        double  JpsiMass_bplus_;
        double  JpsiPt_bplus_;
        int  BplusPVindex_;
        double  BpCt2DPVCosTheta_;
        double  BpCtErr2DCostheta_;
        double  IP3DKandJpsiVtx_;
        double  IP3DKandJpsiVtxErr_;
        int  BpmatchDoubleMu01_;
        int  BpmatchDoubleMu02_;
        int  BpmatchDoubleMu41_;
        int  BpmatchDoubleMu42_;
        int  BpmatchDoubleMu01DiMuon0_;
        int  BpmatchDoubleMu02DiMuon0_;
        int  BplusKmcId_;
        int  BplusKmomId_;
        int  BplusMu1mcId_;
        int  BplusMu1momId_;
        int  BplusMu1gmomId_;
        int  BplusMu2mcId_;
        int  BplusMu2momId_;
        int  BplusMu2gmomId_;
        int  isMatchedBplus_;
        int  BplusDecayChannel_;

        double  BcP_ ;
        double  BcCosAlpha_ ;
        double  BcIP3D_ ;
        double  BcCt_ ;
        double  BcM_;
        double  BcMC_;
        double  BcAng_;
        double  BcAngPi_;
        double  BcProb_;
        double  BctrackPt_;

	int runNumber_;
	int PUinteraction_;
	int EleRecoMCmother1_;
	int EleRecoMC1_;
	unsigned int eventNumber_;
	int lumiSection_;

        int triggerbit_HLTmu3Tk_;
        int triggerbit_HLTmu5_;
        int triggerbit_HLTdoubleIsoMu3_;
        int triggerbit_HLTdoubleMu3_;
	int triggerbit_HLTdoubleMu0_;
	int triggerbit_HLTL1DoubleMuOpen_;
        int triggerbit_HLTmu7_;
        int triggerbit_HLTMu0Track0Jpsi_;
        int triggerbit_HLTL1DoubleMuOpenTight_;
        int triggerbit_HLT_DoubleMu3_Jpsi_v2_;
        int triggerbit_HLT_DoubleMu3_Jpsi_v2MC_;
        int triggerbit_HLT_DoubleMu3_Quarkonium_v2_;
        int triggerbit_HLT_DoubleMu3_Quarkonium_v2MC_;
        int triggerbit_HLT_DoubleMu3_Quarkonium_v1_;
        int triggerbit_Jpsi_Displaced_v1_;
        int triggerbit_Jpsi_Displaced_v1MC_;
        int triggerbit_7Jpsi_Displaced_v1_;
        int triggerbit_7Jpsi_Displaced_v2_;
        int triggerbit_7Jpsi_Displaced_v3_;
        int triggerbit_3p5Jpsi_Displaced_v2_;
        int triggerbit_4Jpsi_Displaced_v1_;
        int triggerbit_4Jpsi_Displaced_v4_;
        int triggerbit_4Jpsi_Displaced_v5_;
        int triggerbit_5Jpsi_Displaced_v1_;
        int triggerbit_5Jpsi_Displaced_v2_;
        int triggerbit_5Jpsi_Displaced_v4_;
        int triggerbit_5Jpsi_Displaced_v5_;
        int triggerbit_Dimuon0_Jpsi_v1_;
        int triggerbit_Dimuon0_Jpsi_v3_;
        int triggerbit_Dimuon0_Jpsi_v5_;
        int triggerbit_Dimuon0_Jpsi_v6_;
        int triggerbit_Dimuon0_Jpsi_v9_;
        int triggerbit_Dimuon0_Jpsi_v10_;
		  int triggerbit_Dimuon0_Jpsi_v14_;
        int triggerbit_Dimuon10_Barrel_;
        int triggerbit_Dimuon13_Barrel_;
        int triggerbit_Dimuon0_Jpsi_Muon_v15_;
        int triggerbit_4Jpsi_Displaced_v9_;
        int triggerbit_4Jpsi_Displaced_v10_;
        int triggerbit_4Jpsi_Displaced_v11_;
        int triggerbit_4Jpsi_Displaced_v12_;

	int matchL11_ ;
	int matchL12_ ;
	int match2mu01_ ;
	int match2mu02_ ;
	int match1mu01_ ;
	int match1mu02_ ;
	int matchDoubleMu31J_ ;
	int matchDoubleMu32J_ ;
	int matchDoubleMu31Q_ ;
	int matchDoubleMu32Q_ ;
	int matchDoubleMu71_ ;
	int matchDoubleMu72_ ;
	int matchDoubleMu41_ ;
	int matchDoubleMu42_ ;
	int matchDoubleMu51_ ;
	int matchDoubleMu52_ ;
	int matchDoubleMu01_ ;
	int matchDoubleMu02_ ;
	int matchDoubleMu101_ ;
	int matchDoubleMu102_ ;
	int matchDoubleMu131_ ;
	int matchDoubleMu132_ ;
	int match2mu31_ ;
	int match2mu32_ ;
        int matchmu0tk01_;
        int matchmu0tk02_;
        int matchFilterJpsi1_;
        int matchFilterJpsi2_;
	int BdmatchL11_ ;
	int BdmatchL12_ ;
	int Bdmatch2mu01_ ;
	int Bdmatch2mu02_ ;
	int Bdmatch2mu31_ ;
	int Bdmatch2mu32_ ;
	int Bdmatch1mu01_ ;
	int Bdmatch1mu02_ ;
	int BdmatchDoubleMu31Q_ ;
	int BdmatchDoubleMu32Q_ ;
	int BdmatchDoubleMu71_ ;
	int BdmatchDoubleMu72_ ;
	int BdmatchDoubleMu41_ ;
	int BdmatchDoubleMu42_ ;
	int BdmatchDoubleMu51_ ;
	int BdmatchDoubleMu52_ ;
	int BdmatchDoubleMu01_ ;
	int BdmatchDoubleMu02_ ;
	int BdmatchDoubleMu101_ ;
	int BdmatchDoubleMu102_ ;
	int BdmatchDoubleMu131_ ;
	int BdmatchDoubleMu132_ ;
	int matchDoubleMu121_ ;
	int matchDoubleMu122_ ;
	

	int MuonType_;
	
	int isPV_;
	int isBS_;
        int NVertices_;

        int  ihaveajpsi_;
        int  BsCowboy_;
        int  BdCowboy_;
        double  BsPhiVtxProb_;
        int  BsMu1QualityG_;
        int  BsMu2QualityG_;
        int  BsMu1QualityT_;
        int  BsMu2QualityT_;
        int  BdMu1QualityG_;
        int  BdMu2QualityG_;
        int  BdMu1QualityT_;
        int  BdMu2QualityT_;
        int  BsSoftMuon1_;
        int  BsSoftMuon2_;
        int  BpSoftMuon1_;
        int  BpSoftMuon2_;
        int  BdSoftMuon1_;
        int  BdSoftMuon2_;

	double	BSx_ ;
	double	BSy_ ;
	double	BSz_ ;
        double  BSdx_ ;
        double  BSdy_ ;
        double  BSdz_ ;
        double  BSsigmaZ_ ;
        double  BSdsigmaZ_ ;
	double	PVx_ ;
	double	PVy_ ;
	double	PVz_ ;
	double	PVerrx_ ;
	double	PVerry_ ;
	double	PVerrz_ ;

        double JpsiVtxProb_;
	double CosDeltaAlpha_;
        double MuMuDCA_;
        double MuMuDistance_;
        double MuMuDistanceSigma_;
        double MuDr1_;
        double MuDr2_;
        double BdJpsiVtxProb_;
	double BdCosDeltaAlpha_;
        double BdMuMuDCA_;
        double BdMuMuDistance_;
        double BdMuMuDistanceSigma_;
        double BdMuDr1_;
        double BdMuDr2_;


	double PVx_refit_;
	double PVy_refit_;
	double PVz_refit_;
	double PVerrx_refit_;
	double PVerry_refit_;
	double PVerrz_refit_;

    
        double JpsiM_alone_;
        double JpsiPhi_alone_;
        double JpsiEta_alone_;
        double JpsiPt_alone_;
        double JpsiMu1Pt_alone_;
        double JpsiMu2Pt_alone_;
        double JpsiMu1Phi_alone_;
        double JpsiMu2Phi_alone_;
        double JpsiMu1Eta_alone_;
        double JpsiMu2Eta_alone_;
        int    MuonCat1_;
        int    MuonCat2_;
        int    JpsiMuon1Cat_alone_;
        int    JpsiMuon2Cat_alone_;
        int    BdJpsiMuon1Cat_alone_;
        int    BdJpsiMuon2Cat_alone_;

        int    Mu1GlobalMuonPromptTight_;
        int    Mu2GlobalMuonPromptTight_;
        int    Mu1TrackerMuonArbitrated_;
        int    Mu1TMLastStationTight_;
        int    Mu1TMOneStationTight_;
        int    Mu1TMLastStationOptimizedLowPtTight_;
        int    Mu1TMLastStationAngTight_;
        int    Mu1TMOneStationAngTight_;
        int    Mu1TMLastStationOptimizedBarrelLowPtTight_;
        int    Mu2TrackerMuonArbitrated_;
        int    Mu2TMLastStationTight_;
        int    Mu2TMOneStationTight_;
        int    Mu2TMLastStationOptimizedLowPtTight_;
        int    Mu2TMLastStationAngTight_;
        int    Mu2TMOneStationAngTight_;
        int    Mu2TMLastStationOptimizedBarrelLowPtTight_;



        double JpsiMu1d0_alone_;
        double JpsiMu2d0_alone_;
        double JpsiMu1dz_alone_;
        double JpsiMu2dz_alone_;
        int JpsiMu1ndof_alone_;
        int JpsiMu2ndof_alone_;
        double JpsiMu1chi2_alone_;
        double JpsiMu2chi2_alone_;
        int JpsiMu1nHits_alone_;
        int JpsiMu2nHits_alone_;
        int JpsiMu1nPixHits_alone_;
        int JpsiMu2nPixHits_alone_;
  
        double  K1Pt_beffit_;
        double  K1Pz_beffit_;
        double  K1Eta_beffit_;
        double  K1Phi_beffit_;
        double  K2Pt_beffit_;
        double  K2Pz_beffit_;
        double  K2Eta_beffit_;
        double  K2Phi_beffit_;

        double  Mu1Pt_beffit_;
        double  Mu1Pz_beffit_;
        double  Mu1Eta_beffit_;
        double  Mu1Phi_beffit_;
        double  Mu2Pt_beffit_;
        double  Mu2Pz_beffit_;
        double  Mu2Eta_beffit_;
        double  Mu2Phi_beffit_;
        double  BdMu1Pt_beffit_;
        double  BdMu1Pz_beffit_;
        double  BdMu1Eta_beffit_;
        double  BdMu1Phi_beffit_;
        double  BdMu2Pt_beffit_;
        double  BdMu2Pz_beffit_;
        double  BdMu2Eta_beffit_;
        double  BdMu2Phi_beffit_;


        double BsFitChi2_;
        int    BsFitNdof_;

        double BsFitVtxProb_;

        double BsFitM_;
        double K1Pt_fit_;
        double K2Pt_fit_;
        double PhiM_fit_;
    	double BsFitEta_;
	double BsFitPt_;
	double BsFitPz_;
	double BsFitPhi_;

	double BsFitVtx_x_;
	double BsFitVtx_y_;
	double BsFitVtx_z_;

        double BsM_nofit_;
        double BsPhi_nofit_;
        double BsEta_nofit_;
        double BsPt_nofit_;
        double BsPz_nofit_;

        double JpsiM_nofit_;
        double JpsiPhi_nofit_;
        double JpsiEta_nofit_;
        double JpsiPt_nofit_;
        double JpsiPz_nofit_;
        double BdJpsiM_nofit_;
        double BdJpsiPhi_nofit_;
        double BdJpsiEta_nofit_;
        double BdJpsiPt_nofit_;
        double BdJpsiPz_nofit_;

        double PhiM_nofit_;
        double PhiPhi_nofit_;
        double PhiEta_nofit_;
        double PhiPt_nofit_;
        double PhiPz_nofit_;

	double costhetaMC_[10];
	double phiMC_[10];
	double cospsiMC_[10];
	double BscosthetaMC_;
	double BsphiMC_;
	double BscospsiMC_;

        double  K1Pt_nofit_;
        double  K1Pz_nofit_;
        double  K1Eta_nofit_;
        double  K1Phi_nofit_;
	int     K1Key_nofit_;
        double  K2Eta_nofit_;
        double  K2Pt_nofit_;
        double  K2Pz_nofit_;
        double  K2Phi_nofit_;
	int     K2Key_nofit_;

        double  K1Chi2_;
        int     K1nHits_;
        double  K2Chi2_;
        int     K2nHits_;
        int     K1pixH_;
        int     K1trkH_;
        int     K2pixH_;
        int     K2trkH_;

        double  Mu1Chi2_;
        int     Mu1nHits_;
        double  Mu2Chi2_;
        int     Mu2nHits_;
        int     Mu1pixH_;
        int     Mu1trkH_;
        int     Mu2pixH_;
        int     Mu2trkH_;

        double Mu1d0_;
        double Mu2d0_;
        double Mu1dz_;
        double Mu2dz_;	

	double costheta_;
	double phi_;
	double cospsi_;
	double Bdcostheta_;
	double Bdphi_;
	double Bdcospsi_;
	double BdcosthetaMC_;
	double BdphiMC_;
	double BdcospsiMC_;
	double AngleBsDecayLength_;

	int isMatched_;
	int isMatchedBd_;

	double BsLxy_;
	double BsLxyErr_;

	double BsCt_;
	double BsCt3D_;
	double BsCt2D_;
	double BsCt2DBS_;
	double BdCt2DBS_;
	double BdCt2DMC_;
	double BdCt3DMC_;
	double BsCtMPV_;
	double BsCt3Drefit_;
	double BsCt2Drefit_;
	double BsCtMPVrefit_;
	double BsCtErr_;
	double BsCtErr3D_;
	double BsCtErr2D_;
	double BsCtErr2DBS_;
	double BdCtErr2DBS_;
	double BsCtErr2D2_;
	double BsCtErrMPV_;
	double BsCtErr3Drefit_;
	double BsCtErr2Drefit_;
	double BsCtErrMPVrefit_;

	double BsErrX_;
	double BsErrY_;
	double BsErrXY_;

	int JpsiNumberOfCandidates_;
	int PhiNumberOfCandidatesBeforeFit_;
	int BsNumberOfCandidatesBeforeFit_;
	int BsNumberOfCandidatesAfterFit_;
	int BsNumberOfCandidatesAfterBestFit_;

        int     K1trkLay_;
        int     K1muDTh_;
        int     K1muCSCh_;
        int     K1muRPCh_;
        int     K2trkLay_;
        int     K2muDTh_;
        int     K2muCSCh_;
        int     K2muRPCh_;
        int     Mu1trkLay_;
        int     Mu1muDTh_;
        int     Mu1muCSCh_;
        int     Mu1muRPCh_;
        int     Mu2trkLay_;
        int     Mu2muDTh_;
        int     Mu2muCSCh_;
        int     Mu2muRPCh_;

        int K1mcId_;
        int K1momId_;
        int K1gmomId_;
        int K2mcId_;
        int K2momId_;
        int K2gmomId_;
        int Mu1mcId_;
        int Mu1momId_;
        int Mu1gmomId_;
        int Mu2mcId_;
        int Mu2momId_;
        int Mu2gmomId_;
     
        double BsCtErr2DCostheta_;
        Double_t BsCt3DPVCosTheta_;
        Double_t BsCt2DPVCosTheta_;
	double BsDist3d_;
	double BsDist3dErr_;
	double BsTime3d_;
	double BsTime3dErr_;
	double BsDist2d_;
	double BsDist2dErr_;
	double BsTime2d_;
	double BsTime2dErr_;

	double dedxTrk_;
	double errdedxTrk_;
	int numdedxTrk_;

	int iPassedCutIdent_;
	int tagEleCharge_;
	double tagEleP_;
        double tagElePhi_;
	double tagEleEta_;
	double tagEleHcalOverEcal_;
	double tagEleSigmaIetaIeta_;
	int iPassedCutIdentBd_;
	int BdTrack1Charge_;


        double K1Fit_par_[7];
        double K2Fit_par_[7];
        double K1Fit_sigX_;
        double K1Fit_sigY_;
        double K1Fit_sigZ_;
        double K2Fit_sigX_;
        double K2Fit_sigY_;
        double K2Fit_sigZ_;
        double K1Fit_sigPX_;
        double K1Fit_sigPY_;
        double K1Fit_sigPZ_;
        double K2Fit_sigPX_;
        double K2Fit_sigPY_;
        double K2Fit_sigPZ_;

        Double_t PVZpos_[30];
        Double_t SVZpos_[30];
        Int_t NTracksInPV_[30];
        Double_t PVAbsPt_[30];
        Double_t BsPtMC_;
        Double_t FirstBsMCZpos_;

        Double_t BsLxy3DMC_;
        Double_t BsLxy2DMC_;
        Double_t BsPMC_;

        double BsCt2DMC_; 
        double BsCt3DMC_; 
        int BsIniFlavour_; 
        int BsEndFlavour_; 
        int BdEndFlavour_; 
        int BdIniFlavour_; 
		  int BdKstarKaon_;	
		  int BdKstarPion_;	
        int ChannelID_; 
        int BdChannelID_; 
	int GenNumberOfBdecays_;
	int BmesonsId_[10];
	int BDauIdMC_[10][15];
	int BDauDauIdMC_[10][15][10];
    	int GenNumberOfDaughters_[10];
	int GenNumberOfDaughtersDaughters_[10][15];

	double BDauMMC_[10][15];
	double BDauPtMC_[10][15];
	double BDauPzMC_[10][15];
	double BDauEtaMC_[10][15];
	double BDauPhiMC_[10][15];

	double BDauDauMMC_[10][15][10];
	double BDauDauPtMC_[10][15][10];
	double BDauDauPzMC_[10][15][10];
	double BDauDauEtaMC_[10][15][10];
	double BDauDauPhiMC_[10][15][10];

	double BMMC_[10];
	double BPtMC_[10];
	double BPxMC_[10];
	double BPyMC_[10];
	double BPzMC_[10];
	double BEtaMC_[10];
	double BPhiMC_[10];

	double BVtxMC_x_[10];
	double BVtxMC_y_[10];
	double BVtxMC_z_[10];
	double BSVtxMC_x_[10];
	double BSVtxMC_y_[10];
	double BSVtxMC_z_[10];
	double BLxy_MC_[10];
	double BCt_MC_[10];
	double BCt_MC2D_[10];
	double BCt_MC3D_[10];

        double genBsVtx_z_, genBsVtx_y_, genBsVtx_x_ ;
        double genBsSVtx_z_, genBsSVtx_y_, genBsSVtx_x_ ;

	int isGenJpsiEvent_;


	// for the Bd->Kstar analysis
	double BdFitChi2_Hyp1_;
   int    BdFitNdof_Hyp1_;

   double BdFitVtxProb_Hyp1_;
   double BdFitVtxProb_;
	int BdNumberOfCandidates_;

	double BdPVx_refit_   ;
	double BdPVy_refit_   ;
	double BdPVz_refit_   ;
	double BdPVerrx_refit_;
	double BdPVerry_refit_;
	double BdPVerrz_refit_;

        double BdFitM_Hyp1_;
  	double BdFitEta_Hyp1_;
	double BdFitPt_Hyp1_;
	double BdFitPz_Hyp1_;
	double BdFitPhi_Hyp1_;

	double BdFitVtx_x_Hyp1_;
	double BdFitVtx_y_Hyp1_;
	double BdFitVtx_z_Hyp1_;

        double BdM_nofit_;
        double BdPhi_nofit_;
        double BdEta_nofit_;
        double BdPt_nofit_;
        double BdPz_nofit_;

	double KstarMass_nofit_Hyp1_ ;
	double KstarMass_nofit_Hyp2_ ; 

	double BdK1_kpi_par_Hyp1_[7];
	double BdK2_kpi_par_Hyp1_[7];
	double BdK1_kpi_sigX_Hyp1_;
        double BdK1_kpi_sigY_Hyp1_;
        double BdK1_kpi_sigZ_Hyp1_;
	double BdK2_kpi_sigX_Hyp1_;
        double BdK2_kpi_sigY_Hyp1_;
        double BdK2_kpi_sigZ_Hyp1_;
	double BdK1_kpi_sigPX_Hyp1_;
	double BdK1_kpi_sigPY_Hyp1_;
        double BdK1_kpi_sigPZ_Hyp1_;
	double BdK2_kpi_sigPX_Hyp1_;
        double BdK2_kpi_sigPY_Hyp1_;
        double BdK2_kpi_sigPZ_Hyp1_;

	double BdFitChi2_Hyp2_;
        int    BdFitNdof_Hyp2_;

        double BdFitVtxProb_Hyp2_;
   

        double BdFitM_Hyp2_;
  	double BdFitEta_Hyp2_;
	double BdFitPt_Hyp2_;
	double BdFitPz_Hyp2_;
	double BdFitPhi_Hyp2_;

	double BdFitVtx_x_Hyp2_;
	double BdFitVtx_y_Hyp2_;
	double BdFitVtx_z_Hyp2_;

    
	double BdK1_kpi_par_Hyp2_[7];
	double BdK2_kpi_par_Hyp2_[7];
	double BdK1_kpi_sigX_Hyp2_;
        double BdK1_kpi_sigY_Hyp2_;
        double BdK1_kpi_sigZ_Hyp2_;
	double BdK2_kpi_sigX_Hyp2_;
        double BdK2_kpi_sigY_Hyp2_;
        double BdK2_kpi_sigZ_Hyp2_;
	double BdK1_kpi_sigPX_Hyp2_;
	double BdK1_kpi_sigPY_Hyp2_;
        double BdK1_kpi_sigPZ_Hyp2_;
	double BdK2_kpi_sigPX_Hyp2_;
        double BdK2_kpi_sigPY_Hyp2_;
        double BdK2_kpi_sigPZ_Hyp2_;


	double BdK1Pt_nofit_  ; 
	double BdK1Pz_nofit_  ; 
	double BdK1Eta_nofit_ ; 
	double BdK1Phi_nofit_ ; 
	int BdK1Key_nofit_  ; 
	double BdK2Pt_nofit_  ; 
	double BdK2Pz_nofit_  ; 
	double BdK2Eta_nofit_ ; 
	double BdK2Phi_nofit_ ; 
	int BdK2Key_nofit_  ; 

	double BdLxy_;
	double BdLxyErr_;
	double BdErrX_;
	double BdErrY_;
	double BdErrXY_;
	double BdCt_;
	double BdCtErr_;

   double BdCt2DPVCosTheta_;
   double Bdt2DPVCosTheta_;
   double BdCtErr2DCostheta_;
   double BdtErr2DCostheta_;
	double BdCt2DPVCosThetaUnfitP_ ;
	double BdDist3d_;
	double BdDist3dErr_;
	double BdTime3d_;
	double BdTime3dErr_;
	double BdDist2d_;
	double BdDist2dErr_;
	double BdTime2d_;
	double BdTime2dErr_;

	int BdK1mcId_;
        int BdK1momId_;
        int BdK1gmomId_;
        int BdK2mcId_;
        int BdK2momId_;
        int BdK2gmomId_;
        int BdMu1mcId_;
        int BdMu1momId_;
        int BdMu1gmomId_;
        int BdMu2mcId_;
        int BdMu2momId_;
        int BdMu2gmomId_;

	//b jet tagging  
	int BJetParton_;
 	double BJetEta_;
	double BJetPhi_;
	double BJetPt_;
	double JetBTagProb_; 

   int BpBJetParton_;
	double BpBJetEta_;
	double BpBJetPhi_;
	double BpBJetPt_;
   double BpJetBTagProb_;
 
   std::vector<int> *PVTrkCharge_;
	std::vector<float> *PVTrkPt_;
	std::vector<float> *PVTrkEta_;
	std::vector<float> *PVTrkPhi_;   

   std::vector<int> *BpPVTrkCharge_;
	std::vector<float> *BpPVTrkPt_;
	std::vector<float> *BpPVTrkEta_;
	std::vector<float> *BpPVTrkPhi_;   


	std::vector<int> *BJetTrkCharge_;
	std::vector<float> *BJetTrkPt_;   
	std::vector<int> *BpBJetTrkCharge_;
	std::vector<float> *BpBJetTrkPt_;   
	//b jet tagging 


	TFile* bsFile_;
	TTree* bsTree_; 
};

#endif

