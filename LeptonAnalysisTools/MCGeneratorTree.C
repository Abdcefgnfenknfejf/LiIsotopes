#include "MCGeneratorTree.hh"

#include "AMSGeometry.h"
#include "AnalysisEvent.hh"
#include "Clamping.hh"
#include "Event.h"
#include "ParticleId.hh"

#define INFO_OUT_TAG "MCGeneratorTree"
#include "debugging.hh"

MCGeneratorTree::MCGeneratorTree()
  : IO::TreeInterface("MCGeneratorTree", "MC generator tree") {

  RegisterBranches();
}

bool MCGeneratorTree::DoesBremsstrahlungPhotonConvert(const AC::MC& mc, const AC::MCParticle* bsPhoton) const {

  for (const AC::MCParticle& particle : mc.MCParticles()) {
    // Only consider particles whose parent is a bremsstrahlung photon.
    if (particle.MotherParticle() != bsPhoton->TrackID())
      continue;

    if (particle.Process() != AC::PairProduction)
      continue;

    if (int(particle.ParticleID()) == ParticleId::Id(ParticleId::Electron) || int(particle.ParticleID()) == ParticleId::Id(ParticleId::Positron))
      return true;
  }

  return false;
}

const std::vector<const AC::MCEventGenerator*> MCGeneratorTree::FindAllBremsstrahlungPhotonsOriginatingFromPrimaryParticle(const AC::MC& mc) const {

  std::vector<const AC::MCEventGenerator*> list;
  for (const AC::MCParticle& particle : mc.MCParticles()) {
    // Only consider direct children of the primary particle.
    if (particle.MotherParticle() != 1)
      continue;

    if (particle.Process() == AC::Bremsstrahlung)
      list.push_back(&particle.EventGenerators().front());
  }

  return list;
}

const AC::MCParticle* MCGeneratorTree::FindBremsstrahlungPhotonOriginatingFromPrimaryParticleThatConvertLater(const AC::MC& mc) const {

  for (const AC::MCParticle& particle : mc.MCParticles()) {
    // Only consider direct children of the primary particle.
    if (particle.MotherParticle() != 1)
      continue;

    if (particle.Process() == AC::Bremsstrahlung) {
      if (DoesBremsstrahlungPhotonConvert(mc, &particle))
        return &particle;
    }
  }

  return nullptr;
}

std::pair<const AC::MCEventGenerator*, const AC::MCEventGenerator*> MCGeneratorTree::FindElecPosiPairOriginatingFromBremsstrahlungPhoton(const AC::MC& mc, const AC::MCParticle* bsPhoton) const {

  const AC::MCEventGenerator* secondaryElectron = nullptr;
  const AC::MCEventGenerator* secondaryPositron = nullptr;

  for (const AC::MCParticle& particle : mc.MCParticles()) {
    // Only consider particles whose parent is a bremsstrahlung photon.
    if (particle.MotherParticle() != bsPhoton->TrackID())
      continue;

    if (particle.Process() != AC::PairProduction)
      continue;

    if (!secondaryElectron && int(particle.ParticleID()) == ParticleId::Id(ParticleId::Electron))
      secondaryElectron = &particle.EventGenerators().front();

    if (!secondaryPositron && int(particle.ParticleID()) == ParticleId::Id(ParticleId::Positron))
      secondaryPositron = &particle.EventGenerators().front();

    if (secondaryElectron && secondaryPositron)
      break;
  }

  return std::make_pair(secondaryElectron, secondaryPositron);
}

void MCGeneratorTree::Fill(const Analysis::Event& event) {

  const Analysis::Particle* particle = event.PrimaryParticle();

  // General
  Run = clampTo<UInt_t>(event.Run());
  Event = clampTo<UInt_t>(event.EventNumber());

  // MC
  if (event.IsMC()) {
    const auto& mc = event.RawEvent()->MC();
    const AC::MCEventGenerator* primaryGenerator = mc.PrimaryEventGenerator();
    assert(primaryGenerator);

    McPrimaryMomentum = primaryGenerator->Momentum();
    McPrimaryTheta = primaryGenerator->Theta();
    McPrimaryPhi = primaryGenerator->Phi();
    McPrimaryX = primaryGenerator->X0();
    McPrimaryY = primaryGenerator->Y0();
    McPrimaryZ = primaryGenerator->Z0();

    const auto& listOfBSPhotons = FindAllBremsstrahlungPhotonsOriginatingFromPrimaryParticle(mc);
    double momentumSum = 0.0;
    for (const auto* generator : listOfBSPhotons)
      momentumSum += generator->Momentum();
    McAllBSPhotonMomenta = momentumSum;

    const AC::MCParticle* bsPhoton = FindBremsstrahlungPhotonOriginatingFromPrimaryParticleThatConvertLater(mc);
    if (bsPhoton) {
      const AC::MCEventGenerator& bsPhotonGenerator = bsPhoton->EventGenerators().front();
      McBSPhotonMomentum = bsPhotonGenerator.Momentum();
      McBSPhotonTheta = bsPhotonGenerator.Theta();
      McBSPhotonPhi = bsPhotonGenerator.Phi();
      McBSPhotonX = bsPhotonGenerator.X0();
      McBSPhotonY = bsPhotonGenerator.Y0();
      McBSPhotonZ = bsPhotonGenerator.Z0();

      auto pair = FindElecPosiPairOriginatingFromBremsstrahlungPhoton(mc, bsPhoton);
      if (pair.first && pair.second) {
        McSecondaryElecMomentum = pair.first->Momentum();
        McSecondaryElecTheta = pair.first->Theta();
        McSecondaryElecPhi = pair.first->Phi();
        McSecondaryElecX = pair.first->X0();
        McSecondaryElecY = pair.first->Y0();
        McSecondaryElecZ = pair.first->Z0();

        McSecondaryPosiMomentum = pair.second->Momentum();
        McSecondaryPosiTheta = pair.second->Theta();
        McSecondaryPosiPhi = pair.second->Phi();
        McSecondaryPosiX = pair.second->X0();
        McSecondaryPosiY = pair.second->Y0();
        McSecondaryPosiZ = pair.second->Z0();
      }
    }
  }

  // Tracker
  TrackerNumberOfTracks = clampTo<UChar_t>(event.NumberOfTrackerTracks());

  const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
  const AC::ParticleHypothesis particleHypothesis = AC::DefaultMass;

  if (particle && particle->HasTrackerTrack()) {
    const AC::TrackerTrack* associatedTrackerTrack = particle->TrackerTrack();
    assert(associatedTrackerTrack);

    TrackerCharge = associatedTrackerTrack->GetChargeAndError(3, 1 /* electron mass */).Charge();
    TrackerPattern = clampTo<Char_t>(particle->TrackerLayerPatternClassification());

    // Choutko max-span
    int choutkoMaxSpanIndex = associatedTrackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
    if (choutkoMaxSpanIndex >= 0) {
      auto& trackFit = associatedTrackerTrack->TrackFits().at(choutkoMaxSpanIndex);
      TrackerTrackChoutkoMaxSpanRigidity = trackFit.Rigidity();
      TrackerTrackChoutkoMaxSpanChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoMaxSpanChiSquareY = trackFit.ChiSquareNormalizedY();
    }
  }

  // ECAL
  EcalNumberOfShowers = clampTo<UChar_t>(event.NumberOfEcalShower());

  if (particle && particle->HasEcalShower()) {
    const AC::ECALShower* shower = particle->EcalShower();
    assert(shower);
    Vector3 sTrackerTrackPoint;
    Vector3 sTrackerTrackDirection;

    if (const Analysis::SplineTrack* spline = particle->GetSplineTrack())
      spline->CalculateLocalPositionAndDirection(AC::AMSGeometry::ZECALUpper, sTrackerTrackPoint, sTrackerTrackDirection);
    else {
      sTrackerTrackPoint.SetXYZ(0, 0, 0);
      sTrackerTrackDirection.SetXYZ(0, 0, 0);
    }

    EcalEnergyDeposited = shower->DepositedEnergyCorrectedForAnodeEfficiencyInMeV();
    EcalEnergyElectronNew = shower->ReconstructedEnergyElectron2017();
    
    EcalCentreOfGravityX = shower->X();
    EcalCentreOfGravityY = shower->Y();
    EcalCentreOfGravityZ = shower->Z();
  }
}

void MCGeneratorTree::UpdateInMemoryBranches() {

}
