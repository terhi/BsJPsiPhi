How to run in the 5_X releases:

-create a directory HeavyFlavorAnalysis in the src dir of the release

-cd HeavyFlavorAnalysis

-clone the BsJpsiPhi package here

-cd ..

-cvs co  RecoVertex/KinematicFitPrimitives 

-cvs co -r V01-12-05 FWCore/Utilities/src/isFinite.cpp  

-addpkg MuonAnalysis/MuonAssociators

to include the MVA electron ID you need further steps:
https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Recipe_for_53X

Recipe for CMSSW_4_4_2_patch8:

cvs co -r V02-00-00 MuonAnalysis/MuonAssociators
cvs co -r V00-00-16 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools
rm /src/EGamma/EGammaAnalysisTools/test/ElectronIsoAnalyzer.cc
cvs co -r V00-04-00 CondFormats/EgammaObjects  
scramv1 b -j 4
