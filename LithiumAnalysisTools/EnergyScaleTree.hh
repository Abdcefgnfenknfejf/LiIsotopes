#ifndef EnergyScaleTree_hh
#define EnergyScaleTree_hh

#include "TreeInterface.hh"

class EnergyScaleTree : public IO::TreeInterface {
public:
  EnergyScaleTree();

  // General information
  IO::TreeBranch<UInt_t> Time  { "Time",   0  };
  IO::TreeBranch<Int_t>  Run   { "Run",   -1  };
  IO::TreeBranch<Int_t>  Event { "Event", -1  };

  // Beam-test information
  IO::TreeBranch<UChar_t>  BtParticleId      { "BtParticleId",      IO::ValueLimitMode::HighestValue };
  IO::TreeBranch<Float_t>  BtNominalMomentum { "BtNominalMomentum",                              0.0 };

  // Monte-Carlo generator information
  IO::TreeBranch<UChar_t>  McParticleId        { "McParticleId",        IO::ValueLimitMode::HighestValue };
  IO::TreeBranch<Float_t>  McGeneratedMomentum { "McGeneratedMomentum",                              0.0 };

  // Tracker information
  IO::TreeBranch<Float_t> TrackerRigidity           { "TrackerRigidity",              0.0 };
  IO::TreeBranch<Char_t>  TrackerPattern            { "TrackerPattern",                -1 };
  IO::TreeBranch<Float_t> TrackerTrackAtEcalTopX    { "TrackerTrackAtEcalTopX",    -999.0 };
  IO::TreeBranch<Float_t> TrackerTrackAtEcalTopY    { "TrackerTrackAtEcalTopY",    -999.0 };
  IO::TreeBranch<Float_t> TrackerTrackAtEcalBottomX { "TrackerTrackAtEcalBottomX", -999.0 };
  IO::TreeBranch<Float_t> TrackerTrackAtEcalBottomY { "TrackerTrackAtEcalBottomY", -999.0 };

  // Ecal information
  IO::TreeBranch<std::vector<Float_t>> EcalEnergyDepositedPerLayer { "EcalEnergyDepositedPerLayer", IO::TreeVectorSize(18) };
  IO::TreeBranch<Float_t> EnergyFractionInLastTwoLayers  { "EnergyFractionInLastTwoLayers",   0.0 };
  IO::TreeBranch<Float_t> EcalEnergyDepositedXRaw        { "EcalEnergyDepositedXRaw",         0.0 };
  IO::TreeBranch<Float_t> EcalEnergyDepositedYRaw        { "EcalEnergyDepositedYRaw",         0.0 };
  IO::TreeBranch<Float_t> EcalEnergyDepositedX           { "EcalEnergyDepositedX",            0.0 };
  IO::TreeBranch<Float_t> EcalEnergyDepositedY           { "EcalEnergyDepositedY",            0.0 };
  IO::TreeBranch<Float_t> EcalEnergyDeposited            { "EcalEnergyDeposited",             0.0 };
  IO::TreeBranch<Float_t> EcalEnergyReconstructed        { "EcalEnergyReconstructed",         0.0 };
  IO::TreeBranch<Float_t> EcalEnergyElectron             { "EcalEnergyElectron",              0.0 };
  IO::TreeBranch<Float_t> EcalEnergyElectronNew          { "EcalEnergyElectronNew",           0.0 };
  IO::TreeBranch<Float_t> EcalShowerPhi                  { "EcalShowerPhi",                -999.0 };
  IO::TreeBranch<Float_t> EcalShowerTheta                { "EcalShowerTheta",              -999.0 };
  IO::TreeBranch<Float_t> EcalBDT                        { "EcalBDT",                        -2.0 };
  IO::TreeBranch<Float_t> EcalShowerMaximum              { "EcalShowerMaximum",              -1.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityX           { "EcalCentreOfGravityX",         -999.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityY           { "EcalCentreOfGravityY",         -999.0 };

  // TRD information
  IO::TreeBranch<UChar_t> TrdPActiveLayersTracker { "TrdPActiveLayersTracker",    0 };
  IO::TreeBranch<Float_t> TrdPLRElecProt          { "TrdPLRElecProt",          -1.0 };
  IO::TreeBranch<Float_t> TrdPLRHeliElec          { "TrdPLRHeliElec",          -1.0 };

  // In-memory tree branches
  IO::InMemoryTreeBranch<Float_t> TrackerTrackThetaAtEcalUpper { "TrackerTrackThetaAtEcalUpper", IO::ValueLimitMode::HighestValue };

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new EnergyScaleTree; }
  virtual const IO::TreeBranch<UInt_t>* CurrentEventTime() const { return &Time; }
  virtual const IO::TreeBranch<Double_t>* CurrentWeight() const { return nullptr; }
  virtual const IO::TreeBranch<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
};

#endif
