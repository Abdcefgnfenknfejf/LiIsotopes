#ifndef MCGeneratorTree_hh
#define MCGeneratorTree_hh

#include <cassert>
#include <cmath>

#include "TreeInterface.hh"
#include "AnalysisSettings.hh"

namespace AC {
class MC;
class MCEventGenerator;
class MCParticle;
}

namespace Analysis {
class Event;
}

class MCGeneratorTree : public IO::TreeInterface {
public:
  MCGeneratorTree();

  // Event information
  IO::TreeBranch<UInt_t>   Run                { "Run",                  0 };
  IO::TreeBranch<UInt_t>   Event              { "Event",                0 };

  // MC information
  IO::TreeBranch<Float_t>  McPrimaryMomentum       { "McPrimaryMomentum",          0.0 };
  IO::TreeBranch<Float_t>  McPrimaryTheta          { "McPrimaryTheta",          -999.0 };
  IO::TreeBranch<Float_t>  McPrimaryPhi            { "McPrimaryPhi",            -999.0 };
  IO::TreeBranch<Float_t>  McPrimaryX              { "McPrimaryX",              -999.0 };
  IO::TreeBranch<Float_t>  McPrimaryY              { "McPrimaryY",              -999.0 };
  IO::TreeBranch<Float_t>  McPrimaryZ              { "McPrimaryZ",              -999.0 };

  IO::TreeBranch<Float_t>  McAllBSPhotonMomenta    { "McAllBSPhotonMomenta",       0.0 };

  IO::TreeBranch<Float_t>  McBSPhotonMomentum      { "McBSPhotonMomentum",         0.0 };
  IO::TreeBranch<Float_t>  McBSPhotonTheta         { "McBSPhotonTheta",         -999.0 };
  IO::TreeBranch<Float_t>  McBSPhotonPhi           { "McBSPhotonPhi",           -999.0 };
  IO::TreeBranch<Float_t>  McBSPhotonX             { "McBSPhotonX",             -999.0 };
  IO::TreeBranch<Float_t>  McBSPhotonY             { "McBSPhotonY",             -999.0 };
  IO::TreeBranch<Float_t>  McBSPhotonZ             { "McBSPhotonZ",             -999.0 };

  IO::TreeBranch<Float_t>  McSecondaryElecMomentum { "McSecondaryElecMomentum",    0.0 };
  IO::TreeBranch<Float_t>  McSecondaryElecTheta    { "McSecondaryElecTheta",    -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryElecPhi      { "McSecondaryElecPhi",      -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryElecX        { "McSecondaryElecX",        -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryElecY        { "McSecondaryElecY",        -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryElecZ        { "McSecondaryElecZ",        -999.0 };

  IO::TreeBranch<Float_t>  McSecondaryPosiMomentum { "McSecondaryPosiMomentum",    0.0 };
  IO::TreeBranch<Float_t>  McSecondaryPosiTheta    { "McSecondaryPosiTheta",    -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryPosiPhi      { "McSecondaryPosiPhi",      -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryPosiX        { "McSecondaryPosiX",        -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryPosiY        { "McSecondaryPosiY",        -999.0 };
  IO::TreeBranch<Float_t>  McSecondaryPosiZ        { "McSecondaryPosiZ",        -999.0 };

  // Tracker information
  IO::TreeBranch<UChar_t>              TrackerNumberOfTracks                                 { "TrackerNumberOfTracks",                                   0 };
  IO::TreeBranch<Float_t>              TrackerCharge                                         { "TrackerCharge",                                        0.0f };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoMaxSpanRigidity                    { "TrackerTrackChoutkoMaxSpanRigidity",                   0.0f };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoMaxSpanChiSquareX                  { "TrackerTrackChoutkoMaxSpanChiSquareX",                 0.0f };
  IO::TreeBranch<Float_t>              TrackerTrackChoutkoMaxSpanChiSquareY                  { "TrackerTrackChoutkoMaxSpanChiSquareY",                 0.0f };
  IO::TreeBranch<Char_t>               TrackerPattern                                        { "TrackerPattern",                                         -1 };

  // Ecal information
  IO::TreeBranch<UChar_t> EcalNumberOfShowers             { "EcalNumberOfShowers",               0 };
  IO::TreeBranch<Float_t> EcalEnergyDeposited             { "EcalEnergyDeposited",             0.0 };
  IO::TreeBranch<Float_t> EcalEnergyElectronNew           { "EcalEnergyElectronNew",           0.0 };

  IO::TreeBranch<Float_t> EcalCentreOfGravityX            { "EcalCentreOfGravityX",         -999.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityY            { "EcalCentreOfGravityY",         -999.0 };
  IO::TreeBranch<Float_t> EcalCentreOfGravityZ            { "EcalCentreOfGravityZ",         -999.0 };

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new MCGeneratorTree; }
  virtual const IO::TreeBranchBase<UInt_t>* CurrentEventTime() const { return nullptr; }
  virtual const IO::TreeBranchBase<Double_t>* CurrentWeight() const { return nullptr; }
  virtual const IO::TreeBranchBase<UChar_t>* CurrentTriggerFlags() const { return nullptr; }

  bool DoesBremsstrahlungPhotonConvert(const AC::MC& mc, const AC::MCParticle* bsPhoton) const;
  const AC::MCParticle* FindBremsstrahlungPhotonOriginatingFromPrimaryParticleThatConvertLater(const AC::MC& mc) const;
  const std::vector<const AC::MCEventGenerator*> FindAllBremsstrahlungPhotonsOriginatingFromPrimaryParticle(const AC::MC& mc) const;
  std::pair<const AC::MCEventGenerator*, const AC::MCEventGenerator*> FindElecPosiPairOriginatingFromBremsstrahlungPhoton(const AC::MC& mc, const AC::MCParticle* bsPhoton) const;
};

#endif
