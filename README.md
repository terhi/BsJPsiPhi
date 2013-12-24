How to run in the 5_X releases:

-create a directory HeavyFlavorAnalysis in the src dir of the release

-cd HeavyFlavorAnalysis

-clone the BsJpsiPhi package here

-cd ..

-export CVSROOT=":ext:<cern-user-account>@lxplus5.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW"
(to get the packages with CVS, git doesn't work so well to me)

-cvs co  RecoVertex/KinematicFitPrimitives 

-cvs co -r V01-12-05 FWCore/Utilities/src/isFinite.cpp  

-git-cms-addpkg MuonAnalysis/MuonAssociators

to include the MVA electron ID you need further steps:
https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#Recipe_for_53X

-cvs co -r V09-00-01    RecoEgamma/EgammaTools

-cvs co -r V00-00-09    EgammaAnalysis/ElectronTools

-cd EgammaAnalysis/ElectronTools/data/

-cat download.url | xargs wget 

-Maybe some plugins of the ElectronTools need to be removed, the default PatAlgos are ok, some xml files have problem I removed the relative ele ID in the python

--------------------------------------------------------------------------------
Recipe for CMSSW_4_4_2_patch8:

cvs co -r V02-00-00 MuonAnalysis/MuonAssociators

cvs co -r V00-00-16 -d EGamma/EGammaAnalysisTools UserCode/EGamma/EGammaAnalysisTools

rm /src/EGamma/EGammaAnalysisTools/test/ElectronIsoAnalyzer.cc

cvs co -r V00-04-00 CondFormats/EgammaObjects  

scramv1 b -j 4
