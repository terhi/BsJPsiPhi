How to run in the 5_X releases:

create a directory HeavyFlavorAnalysis in the src dir of the release
#>cd HeavyFlavorAnalysis
clone the BsJpsiPhi package here
#>cd ..
#>cvs co  RecoVertex/KinematicFitPrimitives 
#>cvs co -r V01-12-05 FWCore/Utilities/src/isFinite.cpp  
#>addpkg MuonAnalysis/MuonAssociators

to include the MVA electron ID yuo need further steps:
