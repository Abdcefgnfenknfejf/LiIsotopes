
#ifndef NucleiAnalysisTree_hh
#define NucleiAnalysisTree_hh

#include "TreeInterface.hh"

class TFile;

class NucleiAnalysisTree : public IO::TreeInterface {
public:
  NucleiAnalysisTree();
  void ReadAdditionalWeights(TFile* fileWeights);

  // Event information
  IO::TreeBranch<UInt_t>   Run                { "Run",                  0 };
  IO::TreeBranch<UInt_t>   Event              { "Event",                0 };
  IO::TreeBranch<UInt_t>   Time               { "Time",                 0 };
  IO::TreeBranch<UChar_t>  TriggerFlags	      { "TriggerFlags",         0 };
  IO::TreeBranch<Float_t>  Rigidity           { "Rigidity" ,            0 };

  // RTI
  IO::TreeBranch<Float_t>  IGRFCutOff25PN        { "IGRFCutOff25PN",       -99999.0 };
  IO::TreeBranch<Float_t>  IGRFCutOff30PN        { "IGRFCutOff30PN",       -99999.0 };
  IO::TreeBranch<Float_t>  IGRFCutOff35PN        { "IGRFCutOff35PN",       -99999.0 };
  IO::TreeBranch<Float_t>  IGRFCutOff40PN        { "IGRFCutOff40PN",       -99999.0 };

  // RTI
  IO::TreeBranch<Float_t>  GeomagneticCutOff25PN{ "GeomagneticCutOff25PN",       -99999.0 };
  IO::TreeBranch<Float_t>  GeomagneticCutOff30PN{ "GeomagneticCutOff30PN",       -99999.0 };
  IO::TreeBranch<Float_t>  GeomagneticCutOff35PN{ "GeomagneticCutOff35PN",       -99999.0 };
  IO::TreeBranch<Float_t>  GeomagneticCutOff40PN{ "GeomagneticCutOff40PN",       -99999.0 };

 // Monte-Carlo generator information
  //  IO::TreeBranch<Float_t>  McSimulateRigidityCutOff			        { "McSimulateRigidityCutOff",                           -99999.0 };
  IO::TreeBranch<UChar_t>  McParticleId						{ "McParticleId",               IO::ValueLimitMode::HighestValue };
  IO::TreeBranch<Float_t>  McGeneratedMomentum					{ "McGeneratedMomentum",                                     0.0 };
  IO::TreeBranch<Float_t>  McGeneratedMomentumAtUpperTof			{ "McGeneratedMomentumAtUpperTof",                           0.0 };
  IO::TreeBranch<Float_t>  McGeneratedMomentumAtLowerTof			{ "McGeneratedMomentumAtLowerTof",                           0.0 };
  IO::TreeBranch<Float_t>  McGeneratedMomentumAtRich				{ "McGeneratedMomentumAtRich",                               0.0 };
  IO::TreeBranch<Float_t>  McGeneratedBeta					{ "McGeneratedBeta",                                      -999.0 };
  IO::TreeBranch<Float_t>  McGeneratedBetaAtRich				{ "McGeneratedBetaAtRich",                                -999.0 };
  IO::TreeBranch<Double_t> McEventWeight					{ "McEventWeight",                                           0.0 };
  IO::TreeBranch<Float_t>  McGeneratedStartPosX					{ "McGeneratedStartPosX",                                 -999.0 };
  IO::TreeBranch<Float_t>  McGeneratedStartPosY					{ "McGeneratedStartPosY",                                 -999.0 };
  IO::TreeBranch<Float_t>  McGeneratedStartPosZ					{ "McGeneratedStartPosZ",                                 -999.0 };
  IO::TreeBranch<Float_t>  McGeneratedDirectionTheta				{ "McGeneratedDirectionTheta",                            -999.0 };
  IO::TreeBranch<Float_t>  McGeneratedDirectionPhi				{ "McGeneratedDirectionPhi",                              -999.0 };
  IO::TreeBranch<Float_t>  McPrimaryParticleEndPosX				{ "McPrimaryParticleEndPosX",                             100000 };
  IO::TreeBranch<Float_t>  McPrimaryParticleEndPosY				{ "McPrimaryParticleEndPosY",                             100000 };
  IO::TreeBranch<Float_t>  McPrimaryParticleEndPosZ				{ "McPrimaryParticleEndPosZ",                             100000 };
  IO::TreeBranch<UChar_t>  McPrimaryParticleDecayProcess			{ "McPrimaryParticleDecayProcess",                             0 };

  
  // Monte-Carlo additional reweighting information: total weight = event weight x additional weight.
  IO::InMemoryTreeBranch<Double_t> McTotalWeight				{ "McTotalWeight",                                           1.0 };
							
  // TOF information
  IO::TreeBranch<UChar_t> TofNumberOfLayers                 { "TofNumberOfLayers",                           0 };
  IO::TreeBranch<Bool_t>  HasTofBeta                        { "HasTofBeta",                              false };
  IO::TreeBranch<Float_t> TofBeta                           { "TofBeta",                               - 999.0 };
  IO::TreeBranch<Float_t> TofInverseBeta                    { "TofInverseBeta",                         -999.0 };
  IO::TreeBranch<Bool_t>  HasTofCharge                      { "HasTofCharge",                            false };
  IO::TreeBranch<Float_t> TofCharge                         { "TofCharge",                                 0.0 };
  IO::TreeBranch<Float_t> TofLowerCharge                    { "TofLowerCharge",                            0.0 };
  IO::TreeBranch<Float_t> TofUpperCharge                    { "TofUpperCharge",                            0.0 };
  IO::TreeBranch<Char_t>  TofBetaNumberOfLayers             { "TofBetaNumberOfLayers",                      -1 };
  IO::TreeBranch<Float_t> TofBetaNormalizedChi2Coo          { "TofBetaNormalizedChi2Coo",                   -1 };
  IO::TreeBranch<Float_t> TofBetaNormalizedChi2Time         { "TofBetaNormalizedChi2Time",                  -1 };
  IO::TreeBranch<Float_t> TofInverseBetaError               { "TofInverseBetaError",                         0 };
  IO::TreeBranch<Float_t> TofTrdMatchingNorm                { "TofTrdMatchingNorm",                     -999.0 };
  IO::TreeBranch<Float_t> TofTimeDifference                 { "TofTimeDifference",                      -999.0 };
  IO::TreeBranch<std::vector<UChar_t>> TofGoodQPath          { "TofGoodQPath",            IO::TreeVectorSize(4) };   //do i need this and Tof raw charge?
  IO::TreeBranch<std::vector<Float_t>> TofRawChargesFromDeDx { "TofRawChargesFromDeDx",   IO::TreeVectorSize(4) }; 
  IO::TreeBranch<std::vector<Float_t>> TofBetaFromDeDx       { "TofBetaFromDeDx",         IO::TreeVectorSize(4) };
  //  IO::TreeBranch<Float_t> TofRawCharge                       { "TofRawCharge",                              0.0 };
  //  IO::TreeBranch<Float_t> TofRawLowerCharge                  { "TofRawLowerCharge",                         0.0 };
  //  IO::TreeBranch<Float_t> TofRawUpperCharge                  { "TofRawUpperCharge",                         0.0 };

  // ACC infromation
  IO::TreeBranch<UChar_t>              AccNumberOfHits                  { "AccNumberOfHits",                     0 };
  IO::TreeBranch<UChar_t>              AccNumberOfStoredClusters        { "AccNumberOfStoredClusters",           0 };

  // Tracker information
  IO::TreeBranch<UChar_t>              TrackerNumberOfTracks                                 { "TrackerNumberOfTracks",                                   0 };
  IO::TreeBranch<UShort_t>             TrackerNumberOfHitsX                                  { "TrackerNumberOfHitsX",                                    0 };
  IO::TreeBranch<UShort_t>             TrackerNumberOfHitsY                                  { "TrackerNumberOfHitsY",                                    0 };
  IO::TreeBranch<UShort_t>             TrackerNumberOfHitsYInner                             { "TrackerNumberOfHitsYInner",                               0 };
  IO::TreeBranch<UShort_t>             TrackerNumberOfHitsXYInner                            { "TrackerNumberOfHitsXYInner",                              0 };

  IO::TreeBranch<Bool_t>               HasTrackerCharge                                      { "HasTrackerCharge",                                    false };
  IO::TreeBranch<Float_t>              InnerTrackerCharge                                    { "InnerTrackerCharge",                                    0.0 };
  IO::TreeBranch<Float_t>              TrackerChargeHL                                       { "TrackerChargeHL",                                       0.0 }; // other charge estimator (by HuLiu)
  IO::TreeBranch<std::vector<Float_t>> TrackerCharges                                        { "TrackerCharges",                      IO::TreeVectorSize(9) };
  IO::TreeBranch<std::vector<Float_t>> TrackerChargesX                                       { "TrackerChargesX",                     IO::TreeVectorSize(9) };
  IO::TreeBranch<std::vector<Float_t>> TrackerChargesY                                       { "TrackerChargesY",                     IO::TreeVectorSize(9) };
  
  IO::TreeBranch<Bool_t>               HasTrackerTrack                                       { "HasTrackerTrack",                                     false };
  IO::TreeBranch<Bool_t>               HasTrackerTrackChoutkoInnerFit                        { "HasTrackerTrackChoutkoInnerFit",                      false };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerRigidity                      { "TrackerTrackChoutkoInnerRigidity",                      0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerChiSquareX                    { "TrackerTrackChoutkoInnerChiSquareX",                    0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerChiSquareY                    { "TrackerTrackChoutkoInnerChiSquareY",                    0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerRigidityRelError              { "TrackerTrackChoutkoInnerRigidityRelError",              0.0 };

  IO::TreeBranch<Bool_t>               HasTrackerTrackChoutkoInnerPlusL1Fit                  { "HasTrackerTrackChoutkoInnerPlusL1Fit",                false };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerPlusL1Rigidity                { "TrackerTrackChoutkoInnerPlusL1Rigidity",                0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerPlusL1ChiSquareX              { "TrackerTrackChoutkoInnerPlusL1ChiSquareX",              0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerPlusL1ChiSquareY              { "TrackerTrackChoutkoInnerPlusL1ChiSquareY",              0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoInnerPlusL1RigidityRelError        { "TrackerTrackChoutkoInnerPlusL1RigidityRelError",        0.0 };

  IO::TreeBranch<Bool_t>               HasTrackerTrackChoutkoAllFit                          { "HasTrackerTrackChoutkoAllFit",                        false };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllRigidity                        { "TrackerTrackChoutkoAllRigidity",                        0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllChiSquareX                      { "TrackerTrackChoutkoAllChiSquareX",                      0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllChiSquareY                      { "TrackerTrackChoutkoAllChiSquareY",                      0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllRigidityRelError                { "TrackerTrackChoutkoAllRigidityRelError",                0.0 };

  IO::TreeBranch<Bool_t>               HasTrackerTrackKalmanAllFit                           { "HasTrackerTrackKalmanAllFit",                         false };
  IO::TreeBranch<Float_t>              TrackerTrackKalmanAllRigidity                         { "TrackerTrackKalmanAllRigidity",                         0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackKalmanAllChiSquareX                       { "TrackerTrackKalmanAllChiSquareX",                       0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackKalmanAllChiSquareY                       { "TrackerTrackKalmanAllChiSquareY",                       0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackKalmanAllRigidityRelError                 { "TrackerTrackKalmanAllRigidityRelError",                 0.0 };

  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaTofRigidity                 { "TrackerTrackChoutkoAllBetaTofRigidity",                 0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaTofChiSquareX               { "TrackerTrackChoutkoAllBetaTofChiSquareX",               0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaTofChiSquareY               { "TrackerTrackChoutkoAllBetaTofChiSquareY",               0.0 };

  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaRichRigidity                { "TrackerTrackChoutkoAllBetaRichRigidity",                0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaRichChiSquareX              { "TrackerTrackChoutkoAllBetaRichChiSquareX",              0.0 };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoAllBetaRichChiSquareY              { "TrackerTrackChoutkoAllBetaRichChiSquareY",              0.0 };

  IO::TreeBranch<std::vector<Bool_t>>  TrackerLayerYPattern                                  { "TrackerLayerYPattern",                IO::TreeVectorSize(9) };
  IO::TreeBranch<std::vector<Bool_t>>  TrackerLayerXYPattern                                 { "TrackerLayerXYPattern",               IO::TreeVectorSize(9) };
  IO::TreeBranch<Char_t>               TrackerPattern                                        { "TrackerPattern",                                         -1 };
  IO::TreeBranch<Bool_t>               TrackerTrackIsNotInSolarArrayShadow                   { "TrackerTrackIsNotInSolarArrayShadow",                  true };
  IO::TreeBranch<std::vector<Float_t>> TrackerFitResidualsVector                             { "TrackerFitResidualsVector",           IO::TreeVectorSize(9) };
  IO::TreeBranch<std::vector<Float_t>> TrackerFitResidualsVectorInner                        { "TrackerFitResidualsVectorInner",      IO::TreeVectorSize(9) };

  IO::TreeBranch<Float_t>              TrackerRawChargeFromDeDx                              { "TrackerRawChargeFromDeDx",                              0.0 };
  IO::TreeBranch<Float_t>              TrackerBetaFromDeDx                                   { "TrackerBetaFromDeDx",                                   0.0 };

  // TrackerTrack information at selected Z = middle of AMS
  IO::TreeBranch<Float_t>              TrackerTrackX                 { "TrackerTrackX",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackY                 { "TrackerTrackY",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackZ                 { "TrackerTrackZ",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackTheta             { "TrackerTrackTheta",            -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackPhi               { "TrackerTrackPhi",              -999.0 };

  IO::TreeBranch<std::vector<Bool_t>>  TrackerTrackIsInFiducialLayer {"TrackerTrackIsInFiducialLayer", IO::TreeVectorSize(9)};
  IO::InMemoryTreeBranch<Float_t>      DeltaYOverDeltaZ              { "DeltaYOverDeltaZ",            -999.0 };    //What does this mean? do I need it? 

  // TrackerTrack information at selected Z = UpperToF (z = 63.65 cm)
  IO::TreeBranch<Float_t>              TrackerTrackUpTofX                 { "TrackerTrackUpTofX",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackUpTofY                 { "TrackerTrackUpTofY",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackUpTofTheta             { "TrackerTrackUpTofTheta",            -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackUpTofPhi               { "TrackerTrackUpTofPhi",              -999.0 };

  // TrackerTrack information at selected Z = LowerToF (z = -63.65 cm)
  IO::TreeBranch<Float_t>              TrackerTrackLowTofX                 { "TrackerTrackLowTofX",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackLowTofY                 { "TrackerTrackLowTofY",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackLowTofTheta             { "TrackerTrackLowTofTheta",            -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackLowTofPhi               { "TrackerTrackLowTofPhi",              -999.0 };

  // TrackerTrack information at selected Z = AglTop   (z = -72.8 cm) !TEMPTATIVE! maybe pmts z also useful?
  IO::TreeBranch<Float_t>              TrackerTrackAglTopX                 { "TrackerTrackAglTopX",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackAglTopY                 { "TrackerTrackAglTopY",                -999.0 };

  // TrackerTrack information at selected Z = NaFTop   (z = -75.3 cm)
  IO::TreeBranch<Float_t>              TrackerTrackNaFTopX                 { "TrackerTrackNaFTopX",                -999.0 };
  IO::TreeBranch<Float_t>              TrackerTrackNaFTopY                 { "TrackerTrackNaFTopY",                -999.0 };

  // TRD information
  IO::TreeBranch<UShort_t> TrdNumberOfSegmentsXZ                { "TrdNumberOfSegmentsXZ",                    0 };
  IO::TreeBranch<UShort_t> TrdNumberOfSegmentsYZ                { "TrdNumberOfSegmentsYZ",                    0 };
  IO::TreeBranch<UShort_t> TrdNumberOfVertices                  { "TrdNumberOfVertices",                      0 };
  IO::TreeBranch<UShort_t> TrdNumberOfVerticesXZ                { "TrdNumberOfVerticesXZ",                    0 };
  IO::TreeBranch<UShort_t> TrdNumberOfVerticesYZ                { "TrdNumberOfVerticesYZ",                    0 };
  IO::TreeBranch<UShort_t> TrdNumberOfHitsOnTrack               { "TrdNumberOfHitsOnTrack",                   0 };
  IO::TreeBranch<UShort_t> TrdNumberOfActiveLayers              { "TrdNumberOfActiveLayers",                  0 };
  IO::TreeBranch<UShort_t> TrdNumberOfRawHits                   { "TrdNumberOfRawHits",                       0 };

  // RICH information
  IO::TreeBranch<Bool_t>  HasRichRing                           { "HasRichRing",                          false };
  IO::TreeBranch<Short_t> RichNumberOfHits                      { "RichNumberOfHits",                        -3 };
  IO::TreeBranch<UChar_t> RichNumberOfRings                     { "RichNumberOfRings",                        0 };
  IO::TreeBranch<Float_t> RichBeta                              { "RichBeta",                             999.0 };
  IO::TreeBranch<Float_t> TofRichBetaMatch                      { "TofRichBetaMatch",                     999.0 };
  IO::TreeBranch<Float_t> RichBetaError                         { "RichBetaError",                            0 };
  IO::TreeBranch<Float_t> RichCharge                            { "RichCharge",                           999.0 };
  IO::TreeBranch<Float_t> RichNExpectedPhotoElectrons           { "RichNExpectedPhotoElectrons",          999.0 };
  IO::TreeBranch<Float_t> RichNCollectedPhotoElectrons          { "RichNCollectedPhotoElectrons",         999.0 };
  IO::TreeBranch<Float_t> RichNPhotoElectrons                   { "RichNPhotoElectrons",                99999.0 };

  IO::TreeBranch<Bool_t> RichIsGood                         { "RichIsGood",                     false };
  IO::TreeBranch<Bool_t> RichIsClean                        { "RichIsClean",                    false };
  IO::TreeBranch<Short_t> RichNumberOfPmts                  { "RichNumberOfPmts",                  -5 };
  IO::TreeBranch<Bool_t> RichIsNaF                          { "RichIsNaF",                      false };
  IO::TreeBranch<Bool_t> RichIsClearNaF                     { "RichIsClearNaF",                 false };
  IO::TreeBranch<Float_t> RichBetaConsistency               { "RichBetaConsistency",          99999.0 };
  IO::TreeBranch<Float_t> RichRingX                         { "RichRingX",                    99999.0 };
  IO::TreeBranch<Float_t> RichRingY                         { "RichRingY",                    99999.0 };
  IO::TreeBranch<Float_t> RichRingProbability               { "RichRingProbability",          99999.0 };
  IO::TreeBranch<UChar_t> RichTileIndex                     { "RichTileIndex",                    255 };
  IO::TreeBranch<Float_t> RichDistanceToTileBorder          { "RichDistanceToTileBorder",     99999.0 };
  IO::TreeBranch<Float_t> RichAtNaFTopPositionX             { "RichAtNaFTopPositionX",        99999.0 };
  IO::TreeBranch<Float_t> RichAtNaFTopPositionY             { "RichAtNaFTopPositionY",        99999.0 };
  IO::TreeBranch<Float_t> RichAtPmtsTopPositionX            { "RichAtPmtsTopPositionX",       99999.0 };
  IO::TreeBranch<Float_t> RichAtPmtsTopPositionY            { "RichAtPmtsTopPositionY",       99999.0 };
 
  
  // Custom accessors
  bool IsMC() const { return McGeneratedMomentum() > 0.0f; }
  void SetFillAmsVariables(bool value) {fFillAmsVariables = value;}
  bool HasPhysicsTrigger() const                                                        { return (TriggerFlags() & 0x3e) != 0; }
  bool HasOnlyUnbiasedTrigger() const                                                   { return (TriggerFlags() & 0x3f) == 1; }

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new NucleiAnalysisTree; }
  virtual const IO::TreeBranch<UInt_t>* CurrentEventTime() const { return &Time; }
  virtual const IO::TreeBranch<Double_t>* CurrentWeight() const { return &McEventWeight; }
  virtual const IO::TreeBranch<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
  virtual float CalculateGeneratedMomentumAtGivenZ(const Analysis::Event& event, const float zPosition /*Tof or Rich position from AMSGeometry*/, double zCoordinateTolerance = 1e-1 /* 0.1cm */);

  bool fFillAmsVariables = false;
};
#endif 
