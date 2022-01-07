#include "EnergyScaleTree.hh"

#include "AMSGeometry.h"
#include "AnalysisEvent.hh"
#include "BeamtestSettings.hh"
#include "EcalLongitudinalShowerFit.hh"
#include "Environment.hh"
#include "Event.h"
#include "SlowControlLookup.hh"
#include "TrackerTrack.h"
#include "TrdLikelihoodCalculation.hh"
#include "Utilities.hh"

#include <cmath>

#include <TFile.h>
#include <TGraphErrors.h>

#define INFO_OUT_TAG "EnergyScaleTree"
#include "debugging.hh"

EnergyScaleTree::EnergyScaleTree()
  : IO::TreeInterface("EnergyScaleTree", "Energy scale study tree") {

  RegisterBranches();
}

void EnergyScaleTree::Fill(const Analysis::Event& event) {

  // Following assertions are guarenteed by the selection.cfg selector.
  const Analysis::Particle* particle = event.PrimaryParticle();
  assert(particle);

  const AC::ECALShower* associatedEcalShower = particle->EcalShower();
  assert(associatedEcalShower);

  const AC::TrackerTrack* associatedTrackerTrack = particle->TrackerTrack();
  assert(associatedTrackerTrack);

  const Analysis::SplineTrack* splineTrack = particle->GetSplineTrack();
  assert(splineTrack);

  // General information
  Time = event.TimeStamp().GetSec();
  Run = event.Run();
  Event = event.EventNumber();

  // Beam-test information
  if (event.IsBeamTest()) {
    const std::vector<Utilities::BeamDescription>& beamDescriptions = Utilities::BeamtestSettings::Self().BeamDescriptionsForEvent(event);
    if (beamDescriptions.empty())
      FATAL_OUT_WITH_EVENT << "No beam test description available for this run. Aborting!" << std::endl;

    const Utilities::BeamDescription& beamDescription = beamDescriptions.front();
    BtParticleId = ParticleId::Id(Utilities::BeamtestSettings::Self().CherenkovParticleId(event));
    BtNominalMomentum = beamDescription.momentum;
  }

  // Monte-Carlo generator information
  McParticleId = event.McParticleId();
  McGeneratedMomentum = event.McMomentum();

  // Tracker information
  const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
  int choutkoMaxSpanIndex = associatedTrackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, AC::DefaultMass);
  if (choutkoMaxSpanIndex >= 0) {
    auto& trackFit = associatedTrackerTrack->TrackFits().at(choutkoMaxSpanIndex);
    TrackerRigidity = trackFit.Rigidity();
  }

  TrackerPattern = particle->TrackerLayerPatternClassification();

  const Vector3& pointAtEcalTop = splineTrack->PositionAtZ(AC::AMSGeometry::ZECALUpper);
  TrackerTrackAtEcalTopX = pointAtEcalTop.X();
  TrackerTrackAtEcalTopY = pointAtEcalTop.Y();

  const Vector3& pointAtEcalBottom = splineTrack->PositionAtZ(AC::AMSGeometry::ZECALLower);
  TrackerTrackAtEcalBottomX = pointAtEcalBottom.X();
  TrackerTrackAtEcalBottomY = pointAtEcalBottom.Y();

  // Ecal information
  Vector3 trackerTrackPoint;
  Vector3 trackerTrackDirection;
  splineTrack->CalculateLocalPositionAndDirection(AC::AMSGeometry::ZECALUpper, trackerTrackPoint, trackerTrackDirection);

  EcalShowerPhi = particle->EcalShower()->PhiStraightLineFitMethod();
  EcalShowerTheta = particle->EcalShower()->ThetaStraightLineFitMethod();
  // not so sure about this one, was originally EnergyReconstruction.EnergyFractionInLastTwoLayers
  EnergyFractionInLastTwoLayers = particle->EcalShower()->RelativeRearLeak();
  for (unsigned int layer = 0; layer < 18; ++layer)
    EcalEnergyDepositedPerLayer().emplace_back(particle->EcalShower()->DepositedEnergyInLayerInMeV(layer + 1));
  EcalEnergyDepositedXRaw = particle->EcalShower()->DepositedEnergyForProjectionInMeV(AC::MeasurementMode::XZMeasurement);
  EcalEnergyDepositedYRaw = particle->EcalShower()->DepositedEnergyForProjectionInMeV(AC::MeasurementMode::YZMeasurement);
  EcalEnergyDepositedX = particle->EcalShower()->DepositedEnergyForProjectionCorrectedForAnodeEfficiencyInMeV(AC::MeasurementMode::XZMeasurement);
  EcalEnergyDepositedY = particle->EcalShower()->DepositedEnergyForProjectionCorrectedForAnodeEfficiencyInMeV(AC::MeasurementMode::YZMeasurement);
  EcalEnergyDeposited = particle->EcalShower()->DepositedEnergyCorrectedForAnodeEfficiencyInMeV();
  EcalEnergyReconstructed = particle->EcalShower()->ReconstructedEnergy();
  EcalEnergyElectron = particle->EcalShower()->ReconstructedEnergyElectron();
  EcalEnergyElectronNew = particle->EcalShower()->ReconstructedEnergyElectron2017();
  EcalBDT = particle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD);

  EcalShowerMaximum = particle->EcalShower()->ShowerMaximumLongitudinalShowerFit();

  EcalCentreOfGravityX = associatedEcalShower->X();
  EcalCentreOfGravityY = associatedEcalShower->Y();

  // TRD information
  bool pXeOk = false;
  double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(event.TimeStamp(), pXeOk);
  if (!pXeOk && event.IsISS()) {
    WARN_OUT << "Can not query xenon partial ressure for time stamp=" << event.TimeStamp().AsDouble() << std::endl;
  }

  auto& trdHitsFromTrackerTrack = particle->TrdHitsFromTrackerTrack();
  Analysis::TrdLikelihoodCalculation likelihoodCalculatorTracker(trdHitsFromTrackerTrack, pXe);
  TrdPActiveLayersTracker = likelihoodCalculatorTracker.NumberOfActiveLayers();

  float trdPLikelihoodTrackerHitsElectron = likelihoodCalculatorTracker.ComputeTrdLikelihood(ParticleId::Electron, TrackerRigidity());
  float trdPLikelihoodTrackerHitsProton = likelihoodCalculatorTracker.ComputeTrdLikelihood(ParticleId::Proton, TrackerRigidity());
  float trdPLikelihoodTrackerHitsHelium = likelihoodCalculatorTracker.ComputeTrdLikelihood(ParticleId::Helium, TrackerRigidity());

  TrdPLRElecProt = trdPLikelihoodTrackerHitsElectron     > 0 ? -std::log(trdPLikelihoodTrackerHitsElectron     / (trdPLikelihoodTrackerHitsElectron     + trdPLikelihoodTrackerHitsProton))     : -1.0f;
  TrdPLRHeliElec = trdPLikelihoodTrackerHitsHelium       > 0 ? -std::log(trdPLikelihoodTrackerHitsHelium       / (trdPLikelihoodTrackerHitsHelium       + trdPLikelihoodTrackerHitsElectron))   : -1.0f;
}

void EnergyScaleTree::UpdateInMemoryBranches() {

  // TrackerTrackThetaAtEcalUpper
  {
    double dx = TrackerTrackAtEcalTopX() - TrackerTrackAtEcalBottomX();
    double dy = TrackerTrackAtEcalTopY() - TrackerTrackAtEcalBottomY();
    double dz = AC::AMSGeometry::ZECALUpper - AC::AMSGeometry::ZECALLower;
    TrackerTrackThetaAtEcalUpper = std::atan2(std::sqrt(dx * dx + dy * dy), dz) * RadToDeg();
  }
}
