#ifndef MultiTrackChargeConfusionStudyTree_hh
#define MultiTrackChargeConfusionStudyTree_hh

#include <cassert>
#include <cmath>

#include "TreeInterface.hh"
#include "AnalysisSettings.hh"

namespace Analysis {
class Event;
}

class MultiTrackChargeConfusionStudyTree : public IO::TreeInterface {
public:
  MultiTrackChargeConfusionStudyTree();

  // Event information
  IO::TreeBranch<UInt_t>   Run                { "Run",                  0 };
  IO::TreeBranch<UInt_t>   Event              { "Event",                0 };

  // MC information
  IO::TreeBranch<Float_t>  McGeneratedMomentum       { "McGeneratedMomentum",          0.0 };

  // TRD information
  IO::TreeBranch<Float_t> TrdEcalDeltaPosY   { "TrdEcalDeltaPosY",   -999.0 };
  IO::TreeBranch<Float_t> TrdEcalDeltaNegY   { "TrdEcalDeltaNegY",   -999.0 };

  // Tracker information
  IO::TreeBranch<UShort_t> TrackerNumberOfTracks   { "TrackerNumberOfTracks",   0 };

  IO::TreeBranch<UShort_t> Trk1TrdHitsOnTrack      { "Trk1TrdHitsOnTrack",      0 };
  IO::TreeBranch<UShort_t> Trk1YHitCount           { "Trk1YHitCount",           0 };
  IO::TreeBranch<Float_t>  Trk1Rig                 { "Trk1Rig",               0.0 };
  IO::TreeBranch<Float_t>  Trk1ChiSquareX          { "Trk1ChiSquareX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk1ChiSquareY          { "Trk1ChiSquareY",     -999.0 };
  IO::TreeBranch<Float_t>  Trk1TofBeta             { "Trk1TofBeta",        -999.0 };
  IO::TreeBranch<Float_t>  Trk1EoverR              { "Trk1EoverR",         -999.0 };
  IO::TreeBranch<Float_t>  Trk1EcalDeltaX          { "Trk1EcalDeltaX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk1EcalDeltaY          { "Trk1EcalDeltaY",     -999.0 };

  IO::TreeBranch<UShort_t> Trk2TrdHitsOnTrack      { "Trk2TrdHitsOnTrack",      0 };
  IO::TreeBranch<UShort_t> Trk2YHitCount           { "Trk2YHitCount",           0 };
  IO::TreeBranch<Float_t>  Trk2Rig                 { "Trk2Rig",               0.0 };
  IO::TreeBranch<Float_t>  Trk2ChiSquareX          { "Trk2ChiSquareX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk2ChiSquareY          { "Trk2ChiSquareY",     -999.0 };
  IO::TreeBranch<Float_t>  Trk2TofBeta             { "Trk2TofBeta",        -999.0 };
  IO::TreeBranch<Float_t>  Trk2EoverR              { "Trk2EoverR",         -999.0 };
  IO::TreeBranch<Float_t>  Trk2EcalDeltaX          { "Trk2EcalDeltaX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk2EcalDeltaY          { "Trk2EcalDeltaY",     -999.0 };

  IO::TreeBranch<UShort_t> Trk3TrdHitsOnTrack      { "Trk3TrdHitsOnTrack",      0 };
  IO::TreeBranch<UShort_t> Trk3YHitCount           { "Trk3YHitCount",           0 };
  IO::TreeBranch<Float_t>  Trk3Rig                 { "Trk3Rig",               0.0 };
  IO::TreeBranch<Float_t>  Trk3ChiSquareX          { "Trk3ChiSquareX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk3ChiSquareY          { "Trk3ChiSquareY",     -999.0 };
  IO::TreeBranch<Float_t>  Trk3TofBeta             { "Trk3TofBeta",        -999.0 };
  IO::TreeBranch<Float_t>  Trk3EoverR              { "Trk3EoverR",         -999.0 };
  IO::TreeBranch<Float_t>  Trk3EcalDeltaX          { "Trk3EcalDeltaX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk3EcalDeltaY          { "Trk3EcalDeltaY",     -999.0 };

  IO::TreeBranch<UShort_t> Trk4TrdHitsOnTrack      { "Trk4TrdHitsOnTrack",      0 };
  IO::TreeBranch<UShort_t> Trk4YHitCount           { "Trk4YHitCount",           0 };
  IO::TreeBranch<Float_t>  Trk4Rig                 { "Trk4Rig",               0.0 };
  IO::TreeBranch<Float_t>  Trk4ChiSquareX          { "Trk4ChiSquareX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk4ChiSquareY          { "Trk4ChiSquareY",     -999.0 };
  IO::TreeBranch<Float_t>  Trk4TofBeta             { "Trk4TofBeta",        -999.0 };
  IO::TreeBranch<Float_t>  Trk4EoverR              { "Trk4EoverR",         -999.0 };
  IO::TreeBranch<Float_t>  Trk4EcalDeltaX          { "Trk4EcalDeltaX",     -999.0 };
  IO::TreeBranch<Float_t>  Trk4EcalDeltaY          { "Trk4EcalDeltaY",     -999.0 };

  // Ecal information
  IO::TreeBranch<UShort_t> EcalNumberOfShowers            { "EcalNumberOfShowers",               0 };
  IO::TreeBranch<Float_t> EcalEnergyDeposited             { "EcalEnergyDeposited",             0.0 };
  IO::TreeBranch<Float_t> EcalEnergyElectron              { "EcalEnergyElectron",              0.0 };

  IO::TreeBranch<Float_t> EcalCentreOfGravityX            { "EcalCentreOfGravityX",         -999.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityY            { "EcalCentreOfGravityY",         -999.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityZ            { "EcalCentreOfGravityZ",         -999.0 };

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new MultiTrackChargeConfusionStudyTree; }
  virtual const IO::TreeBranchBase<UInt_t>* CurrentEventTime() const { return nullptr; }
  virtual const IO::TreeBranchBase<Double_t>* CurrentWeight() const { return nullptr; }
  virtual const IO::TreeBranchBase<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
};

#endif
