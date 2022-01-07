#include "MultiTrackChargeConfusionStudyTree.hh"

#include "AMSGeometry.h"
#include "AnalysisEvent.hh"
#include "Clamping.hh"
#include "Event.h"
#include "TrackFactory.hh"
#include "TrdTracking.hh"

#define INFO_OUT_TAG "MultiTrackChargeConfusionStudyTree"
#include "debugging.hh"

MultiTrackChargeConfusionStudyTree::MultiTrackChargeConfusionStudyTree()
  : IO::TreeInterface("MultiTrackChargeConfusionStudyTree", "Multi-track CC study tree") {

  RegisterBranches();
}

const AC::ECALShower* ChooseMaximumEcalShower(const AC::Event* event) {

  const AC::ECALShower* useEcalShower = nullptr;

  double highestEnergy = 0.0;
  for (const AC::ECALShower& shower : event->ECAL().Showers()) {
    if (shower.DepositedEnergyInMeV() > highestEnergy) {
      highestEnergy = shower.DepositedEnergyInMeV();
      useEcalShower = &shower;
    }
  }

  return useEcalShower;
}

void PropagateTrdTrackThroughMagneticField(const Analysis::TrdTrack& trdTrack, double z, double chargeTimesMomentum, Vector3& extrapolatedPoint, Vector3& extrapolatedDirection) {

  static const double sDeltaZ = 10.0;                  // Arbitary delta in Z direction (10 cm).
  static const double sFieldStrength = 0.14;           // B-Field in Tesla
  static const double sMagnetLengthHalf = 100.0 / 2.0; // Magnet length in cm

  double charge = chargeTimesMomentum > 0 ? 1.0 : -1.0;

  // Non-bending plane (X-Z)
  double x = trdTrack.X(z);
  double x2 = trdTrack.X(z + sDeltaZ);

  // Bending plane (Y-Z)
  // outside the magnet use straight line:
  // y = y01 + m1 * z above the magnet
  // y = y03 + m3 * z below the magnet
  // y = y0 + Q * std::sqrt(std::pow(R, 2) - std::pow(z - z0, 2)) inside the magnet
  double m1 = trdTrack.DerivativeYZ(z);
  double y01 = trdTrack.Y(0);

  double momentumAlongZ = std::abs(chargeTimesMomentum) * trdTrack.DirectionZ();
  double curvatureSquared = std::pow(momentumAlongZ / (0.3 * sFieldStrength) * 100.0 /* m -> cm */, 2);

  double z0 = charge * m1 * std::sqrt(curvatureSquared / (std::pow(m1, 2) + 1.0)) + sMagnetLengthHalf;
  double y0 = y01 + m1 * sMagnetLengthHalf - charge * std::sqrt(curvatureSquared - std::pow(sMagnetLengthHalf - z0, 2));
  double m3 = charge * (z0 + sMagnetLengthHalf) / std::sqrt(curvatureSquared - std::pow(-sMagnetLengthHalf - z0, 2));
  double y03 = y0 + charge * std::sqrt(curvatureSquared - std::pow((-sMagnetLengthHalf) - z0, 2)) - m3 * (-sMagnetLengthHalf);

  double y = 0;
  double dyOverDz = 0;
  if (z > sMagnetLengthHalf) {
    y = y01 + m1 * z;
    dyOverDz = m1;
  }
  else if (z <= sMagnetLengthHalf && z >= -sMagnetLengthHalf) {
    y = y0 + charge * std::sqrt(curvatureSquared - std::pow(z - z0, 2));
    dyOverDz = -charge * (z - z0) / std::sqrt(curvatureSquared - std::pow(z - z0, 2));
  }
  else {
    y = y03 + m3 * z;
    dyOverDz = m3;
  }

  double y2 = y + dyOverDz * sDeltaZ;
  extrapolatedPoint.SetXYZ(x, y, z);

  extrapolatedDirection.SetXYZ(x2 - x, y2 - y, sDeltaZ);
  extrapolatedDirection.Normalize();
}

void MultiTrackChargeConfusionStudyTree::Fill(const Analysis::Event& event) {

  const AC::ECALShower* ecalShower = ChooseMaximumEcalShower(event.RawEvent());

  // General
  Run = clampTo<UInt_t>(event.Run());
  Event = clampTo<UInt_t>(event.EventNumber());

  if (event.IsMC()) {
    const AC::MCEventGenerator* primaryGenerator = event.RawEvent()->MC().PrimaryEventGenerator();
    assert(primaryGenerator);

    McGeneratedMomentum = primaryGenerator->Momentum();
  }

  // TRD tracking
  if (ecalShower) {
    static Analysis::TrdTracking sTrdTracking(&Utilities::ConfigHandler::GetGlobalInstance());
    std::vector<Analysis::TrdTrack> trdTracks;
    sTrdTracking.CreateTrdTracksFromSegments_TrdStandalone(trdTracks, event);

    if (!trdTracks.empty()) {
      float energy = std::abs(ecalShower->ReconstructedEnergyElectron());

      Vector3 extrapolatedPointAssumingPosRig, extrapolatedPointAssumingNegRig, extrapolatedDirection;
      PropagateTrdTrackThroughMagneticField(trdTracks.front(), ecalShower->Z(), energy, extrapolatedPointAssumingPosRig, extrapolatedDirection);
      PropagateTrdTrackThroughMagneticField(trdTracks.front(), ecalShower->Z(), -energy, extrapolatedPointAssumingNegRig, extrapolatedDirection);

      TrdEcalDeltaPosY = ecalShower->Y() - extrapolatedPointAssumingPosRig.Y();
      TrdEcalDeltaNegY = ecalShower->Y() - extrapolatedPointAssumingNegRig.Y();
    }
  }

  // Tracker
  static Analysis::TrackFactory sTrackFactory;
  static Analysis::SplineTrack sSplineTrack;

  const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
  const AC::ParticleHypothesis particleHypothesis = AC::DefaultMass;

  TrackerNumberOfTracks = event.NumberOfTrackerTracks();
  if (event.NumberOfTrackerTracks() >= 1) {
    const AC::TrackerTrack& trackOne = event.RawEvent()->Tracker().Tracks()[0];
    for (const AC::TOFBeta& tofBeta : event.RawEvent()->TOF().Betas()) {
      if (tofBeta.TrackerTrackIndex() == 0)
        Trk1TofBeta = tofBeta.Beta();
    }

    Trk1TrdHitsOnTrack = trackOne.TrdKNumberOfHitsForLikelihoods(); 
    Trk1YHitCount = trackOne.NumberOfHitsY();
    int choutkoMaxSpanIndex = trackOne.GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
    if (choutkoMaxSpanIndex >= 0) {
      auto& trackFit = trackOne.TrackFits().at(choutkoMaxSpanIndex);
      Trk1Rig = trackFit.Rigidity();
      Trk1ChiSquareX = trackFit.ChiSquareNormalizedX();
      Trk1ChiSquareY = trackFit.ChiSquareNormalizedY();
      if (ecalShower) {
        sSplineTrack.Clear();
        sTrackFactory.CreateSplineTrackFrom(trackOne.TrackFitCoordinates(), sSplineTrack);
        assert(!sSplineTrack.IsEmpty());

        Trk1EcalDeltaX = sSplineTrack.X(ecalShower->Z()) - ecalShower->X();
        Trk1EcalDeltaY = sSplineTrack.Y(ecalShower->Z()) - ecalShower->Y();
        Trk1EoverR = std::abs(ecalShower->DepositedEnergyInMeV() / (1000.0 * trackFit.Rigidity()));
      }
    }
  }

  if (event.NumberOfTrackerTracks() >= 2) {
    const AC::TrackerTrack& trackTwo = event.RawEvent()->Tracker().Tracks()[1];
    for (const AC::TOFBeta& tofBeta : event.RawEvent()->TOF().Betas()) {
      if (tofBeta.TrackerTrackIndex() == 1)
        Trk2TofBeta = tofBeta.Beta();
    }

    Trk2TrdHitsOnTrack = trackTwo.TrdKNumberOfHitsForLikelihoods(); 
    Trk2YHitCount = trackTwo.NumberOfHitsY();
    int choutkoMaxSpanIndex = trackTwo.GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
    if (choutkoMaxSpanIndex >= 0) {
      auto& trackFit = trackTwo.TrackFits().at(choutkoMaxSpanIndex);
      Trk2Rig = trackFit.Rigidity();
      Trk2ChiSquareX = trackFit.ChiSquareNormalizedX();
      Trk2ChiSquareY = trackFit.ChiSquareNormalizedY();
      if (ecalShower) {
        sSplineTrack.Clear();
        sTrackFactory.CreateSplineTrackFrom(trackTwo.TrackFitCoordinates(), sSplineTrack);
        assert(!sSplineTrack.IsEmpty());

        Trk2EcalDeltaX = sSplineTrack.X(ecalShower->Z()) - ecalShower->X();
        Trk2EcalDeltaY = sSplineTrack.Y(ecalShower->Z()) - ecalShower->Y();
        Trk2EoverR = std::abs(ecalShower->DepositedEnergyInMeV() / (1000.0 * trackFit.Rigidity()));
      }
    }
  }

  if (event.NumberOfTrackerTracks() >= 3) {
    const AC::TrackerTrack& trackThree = event.RawEvent()->Tracker().Tracks()[2];
    for (const AC::TOFBeta& tofBeta : event.RawEvent()->TOF().Betas()) {
      if (tofBeta.TrackerTrackIndex() == 2)
        Trk3TofBeta = tofBeta.Beta();
    }

    Trk3TrdHitsOnTrack = trackThree.TrdKNumberOfHitsForLikelihoods(); 
    Trk3YHitCount = trackThree.NumberOfHitsY();
    int choutkoMaxSpanIndex = trackThree.GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
    if (choutkoMaxSpanIndex >= 0) {
      auto& trackFit = trackThree.TrackFits().at(choutkoMaxSpanIndex);
      Trk3Rig = trackFit.Rigidity();
      Trk3ChiSquareX = trackFit.ChiSquareNormalizedX();
      Trk3ChiSquareY = trackFit.ChiSquareNormalizedY();
      if (ecalShower) {
        sSplineTrack.Clear();
        sTrackFactory.CreateSplineTrackFrom(trackThree.TrackFitCoordinates(), sSplineTrack);
        assert(!sSplineTrack.IsEmpty());

        Trk3EcalDeltaX = sSplineTrack.X(ecalShower->Z()) - ecalShower->X();
        Trk3EcalDeltaY = sSplineTrack.Y(ecalShower->Z()) - ecalShower->Y();
        Trk3EoverR = std::abs(ecalShower->DepositedEnergyInMeV() / (1000.0 * trackFit.Rigidity()));
      }
    }
  }

  if (event.NumberOfTrackerTracks() >= 4) {
    const AC::TrackerTrack& trackFour = event.RawEvent()->Tracker().Tracks()[3];
    for (const AC::TOFBeta& tofBeta : event.RawEvent()->TOF().Betas()) {
      if (tofBeta.TrackerTrackIndex() == 3)
        Trk4TofBeta = tofBeta.Beta();
    }

    Trk4TrdHitsOnTrack = trackFour.TrdKNumberOfHitsForLikelihoods(); 
    Trk4YHitCount = trackFour.NumberOfHitsY();
    int choutkoMaxSpanIndex = trackFour.GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
    if (choutkoMaxSpanIndex >= 0) {
      auto& trackFit = trackFour.TrackFits().at(choutkoMaxSpanIndex);
      Trk4Rig = trackFit.Rigidity();
      Trk4ChiSquareX = trackFit.ChiSquareNormalizedX();
      Trk4ChiSquareY = trackFit.ChiSquareNormalizedY();
      if (ecalShower) {
        sSplineTrack.Clear();
        sTrackFactory.CreateSplineTrackFrom(trackFour.TrackFitCoordinates(), sSplineTrack);
        assert(!sSplineTrack.IsEmpty());

        Trk4EcalDeltaX = sSplineTrack.X(ecalShower->Z()) - ecalShower->X();
        Trk4EcalDeltaY = sSplineTrack.Y(ecalShower->Z()) - ecalShower->Y();
        Trk4EoverR = std::abs(ecalShower->DepositedEnergyInMeV() / (1000.0 * trackFit.Rigidity()));
      }
    }
  }

  EcalNumberOfShowers = event.NumberOfEcalShower();
  if (ecalShower) {
    EcalEnergyDeposited = ecalShower->DepositedEnergyInMeV() / 1000.0;
    EcalEnergyElectron = ecalShower->ReconstructedEnergyElectron();
    EcalCentreOfGravityX = ecalShower->X();
    EcalCentreOfGravityY = ecalShower->Y();
    EcalCentreOfGravityZ = ecalShower->Z();
  }
}

void MultiTrackChargeConfusionStudyTree::UpdateInMemoryBranches() {

}
