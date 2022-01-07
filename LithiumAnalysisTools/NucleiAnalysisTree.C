#include "NucleiAnalysisTree.hh"
#include "AMSGeometry.h"
#include "AnalysisEvent.hh"
#include "AnalysisParticle.hh"
#include "Event.h"
#include "MCProcess.h"
#include "ParticleId.hh"
#include "StraightLineTrack.hh"
#include "TofChargeCalculator.hh"
#include "TOF.h"
#include "TrackerTrack.h"
#include "TrdTrackTofMatching.hh"
#include "KinematicVariable.hh"
#include "RTIReader.hh"

#include <TFile.h>
#include <TH1.h>

#include <cmath>

#ifdef HAVE_AMS_SUPPORT
#ifndef __ROOTSHAREDLIBRARY__
#define __ROOTSHAREDLIBRARY__
#endif
#include "DisableWarnings.h"
#include "root.h"
#include "Tofrec02_ihep.h"
#include "TrCharge.h"
#include "TrExtAlignDB.h"
#include "EnableWarnings.h"
#endif
#define INFO_OUT_TAG "NucleiAnalysisTree"
#include "debugging.hh"

NucleiAnalysisTree::NucleiAnalysisTree()
  : IO::TreeInterface("NucleiAnalysisTree", "Beryllium analysis tree") {

  RegisterBranches();
}


void NucleiAnalysisTree::Fill(const Analysis::Event& event) {

  const Analysis::Particle* primaryParticle = event.PrimaryParticle();

  // Event information
  Run = event.Run();
  Event = event.EventNumber();
  Time = event.TimeStamp().GetSec();
  TriggerFlags = event.TriggerFlags();

  // RTI
  if (primaryParticle && !event.IsMC() && !event.IsBeamTest()) {
    IGRFCutOff25PN = event.IGRFMaxCutOff(RTI::CutOffMode::CutOff25PN);  //cut off in 25 deg field-of-view for positively charged particles.
    IGRFCutOff30PN = event.IGRFMaxCutOff(RTI::CutOffMode::CutOff30PN);  
    IGRFCutOff35PN = event.IGRFMaxCutOff(RTI::CutOffMode::CutOff35PN);
    IGRFCutOff40PN = event.IGRFMaxCutOff(RTI::CutOffMode::CutOff40PN);
  }

  // RTI
  if (primaryParticle && !event.IsMC() && !event.IsBeamTest()) {
    GeomagneticCutOff25PN = event.GeomagneticMaxCutOff(RTI::CutOffMode::CutOff25PN);
    GeomagneticCutOff30PN = event.GeomagneticMaxCutOff(RTI::CutOffMode::CutOff30PN);
    GeomagneticCutOff35PN = event.GeomagneticMaxCutOff(RTI::CutOffMode::CutOff35PN);
    GeomagneticCutOff40PN = event.GeomagneticMaxCutOff(RTI::CutOffMode::CutOff40PN);
  }

  McEventWeight = event.McWeight(); // If it is not MC it will be filled with 1.
  // Monte-Carlo generator information

  //  RTI::Reader* rti = RTI::Reader::Self();
  if (event.IsMC()) {

    //    const RTI::Record* record = rti->RTIRecordForTimeStamp(event.TimeStamp().GetSec());
    //  McSimulateRigidityCutOff = record->IGRFCutOff(RTI::CutOffMode::CutOff40PN);

    McParticleId = event.McParticleId();
    McGeneratedMomentum = event.McMomentum();
    McGeneratedBeta = Utilities::KinematicVariable::ConvertFromMomentum(McGeneratedMomentum(), Utilities::KinematicVariable::Beta, ParticleId::ToSpecies(McParticleId()));
    static const float sUpperTOFAcceptancePlane = 62.875; // cm //CHECK, see const double AMSG4SteppingAction::facc_pl in gbatch/CC/geant4.C
    McGeneratedMomentumAtUpperTof = CalculateGeneratedMomentumAtGivenZ(event, sUpperTOFAcceptancePlane);
    static const float sLowerTOFAcceptancePlane = -62.875; // cm
    McGeneratedMomentumAtLowerTof = CalculateGeneratedMomentumAtGivenZ(event, sLowerTOFAcceptancePlane);
    static const float sRichAcceptancePlane = -69.975; // cm
    McGeneratedMomentumAtRich = CalculateGeneratedMomentumAtGivenZ(event, sRichAcceptancePlane);
    McGeneratedBetaAtRich = Utilities::KinematicVariable::ConvertFromMomentum(McGeneratedMomentumAtRich(), Utilities::KinematicVariable::Beta, ParticleId::ToSpecies(McParticleId()));
    const AC::MCEventGenerator* primaryGenerator = event.RawEvent()->MC().PrimaryEventGenerator();
    assert(primaryGenerator);
    McGeneratedStartPosX = primaryGenerator->X0();
    McGeneratedStartPosY = primaryGenerator->Y0();
    McGeneratedStartPosZ = primaryGenerator->Z0();
    McGeneratedDirectionTheta = primaryGenerator->Theta();
    McGeneratedDirectionPhi = primaryGenerator->Phi();

    if (event.GetMcPrimarySplineTrack()) {
      // Here GetXmin() really refers to Z-direction as the spline track is defined in different coordinate system.
      double zTopFromZY = event.GetMcPrimarySplineTrack()->SplineZY().GetXmin();
      double zTopFromZX = event.GetMcPrimarySplineTrack()->SplineZX().GetXmin();
      assert(fuzzyCompare(zTopFromZY, zTopFromZX));
      const double mcSplineEndPositionZ = zTopFromZY;
      const Vector3& mcSplineEndposition = event.GetMcPrimarySplineTrack()->PositionAtZ(zTopFromZY);

      McPrimaryParticleEndPosX = mcSplineEndposition.X();
      McPrimaryParticleEndPosY = mcSplineEndposition.Y();
      McPrimaryParticleEndPosZ = mcSplineEndPositionZ;

      const AC::MCParticle* mcPrimaryParticle = event.RawEvent()->MC().PrimaryParticle();
      assert(mcPrimaryParticle);
      assert(mcPrimaryParticle->TrackID() == 1);

      const AC::MC::MCParticlesVector& mcParticles = event.RawEvent()->MC().MCParticles();
      bool foundPrimaryProcess = false;
      // Set Decay process type to unknown
      McPrimaryParticleDecayProcess = static_cast<UChar_t>(AC::MCProcess::Unknown);
      DEBUG_OUT << "Event No.: " << event.EventNumber() << std::endl;
      for (const auto& mcParticle : mcParticles) {
	// Only continue if particle originates from primary particle
	if (mcParticle.MotherParticle() != 1)
	  continue;
	const AC::MC::EventGeneratorsVector& eventGenerators = mcParticle.EventGenerators();
	for (const auto& eventGenerator : eventGenerators) {

	  if(!fuzzyCompare(static_cast<float>(eventGenerator.X0()), static_cast<float>(mcSplineEndposition.X())))
	    continue;
	  if(!fuzzyCompare(static_cast<float>(eventGenerator.Y0()), static_cast<float>(mcSplineEndposition.Y())))
	    continue;
	  if(!fuzzyCompare(static_cast<float>(eventGenerator.Z0()), static_cast<float>(mcSplineEndPositionZ)))
	    continue;

	  UChar_t furtherProcessType = mcParticle.Process();
	  if(foundPrimaryProcess && furtherProcessType != McPrimaryParticleDecayProcess()){
	    WARN_OUT_ON_MASTER_ONCE << "Second found primary process of type "
				    << AC::ProcessAsString(static_cast<AC::MCProcess>(furtherProcessType))
				    << " does not match to first found one of type "
				    << AC::ProcessAsString(static_cast<AC::MCProcess>(McPrimaryParticleDecayProcess()))
				    << std::endl;
	  }
	  McPrimaryParticleDecayProcess = furtherProcessType;
	  foundPrimaryProcess = true;
	  break;
	}
      }
      DEBUG_OUT << "Decay Process: " << AC::ProcessAsString(static_cast<AC::MCProcess>(McPrimaryParticleDecayProcess()))
		<< std::endl;
    }
    else {
      WARN_OUT << "MC Spline Primary spline track is empty!" << std::endl;
    }
  }

  // TOF information
  TofNumberOfLayers = event.TofNumberOfLayers();
  if (primaryParticle && primaryParticle->HasTofBeta()) {
    const AC::TOFBeta* tofBeta = primaryParticle->TofBeta();
    assert(tofBeta);
    HasTofBeta = primaryParticle->HasTofBeta();
    Analysis::TofChargeCalculator chargeCalculator(event.RawEvent()->TOF(), *tofBeta);
    TofBeta = primaryParticle->BetaTof();
    TofInverseBeta = tofBeta->BetaConverted();
    HasTofCharge = primaryParticle->HasTofCharge();
    TofCharge = chargeCalculator.GetChargeAndError(Analysis::TofChargeCalculator::GoodLayers).charge;
    TofLowerCharge = chargeCalculator.GetChargeAndError(Analysis::TofChargeCalculator::GoodLowerTOFLayers).charge;
    TofUpperCharge = chargeCalculator.GetChargeAndError(Analysis::TofChargeCalculator::GoodUpperTOFLayers).charge;
    TofBetaNumberOfLayers = primaryParticle->TofBetaNumberOfLayers(); //FIXME: don't use particle accessor.
    TofBetaNormalizedChi2Coo = tofBeta->NormalizedChi2Coo();
    TofBetaNormalizedChi2Time = tofBeta->NormalizedChi2Time();
    TofInverseBetaError = tofBeta->InverseBetaError();

    //  if (AC::CurrentACQtVersion() >= AC::ACQtVersion(7, 7)) {
    //  TofRawCharge = chargeCalculator.GetRawChargeAndError(Analysis::TofChargeCalculator::GoodLayers).charge;
    //  TofRawLowerCharge = chargeCalculator.GetRawChargeAndError(Analysis::TofChargeCalculator::GoodLowerTOFLayers).charge;
    //  TofRawUpperCharge = chargeCalculator.GetRawChargeAndError(Analysis::TofChargeCalculator::GoodUpperTOFLayers).charge;
    //  }
    
    if (primaryParticle->GetTrdTrack()) {
      Analysis::TrdTrackTofMatching& match = Analysis::TrdTrackTofMatching::SharedInstance();
      match.ApplyMatching(event, primaryParticle->GetTrdTrack());
      TofTimeDifference = match.TimeDifference();
      TofTrdMatchingNorm = match.MatchingNorm();
    }
  }

  // ACC information
  AccNumberOfHits = event.RawEvent()->Trigger().ACCFlagsBitset().count();
  AccNumberOfStoredClusters = event.NumberOfAccHits();

  // Tracker information
  TrackerNumberOfTracks = event.NumberOfTrackerTracks();
  if (primaryParticle && primaryParticle->HasTrackerTrack()) {
    const AC::TrackerTrack* trackerTrack = primaryParticle->TrackerTrack();
    assert(trackerTrack);
    HasTrackerTrack = primaryParticle->HasTrackerTrack();
    HasTrackerCharge = primaryParticle->HasTrackerCharge();
    Rigidity = primaryParticle->Rigidity(); //FixMe: is this overlap with other varibales
    
    TrackerNumberOfHitsX = trackerTrack->NumberOfHitsX(); // Hits in L1 - L9.
    TrackerNumberOfHitsY = trackerTrack->NumberOfHitsY();
    TrackerNumberOfHitsYInner = trackerTrack->NHyInner(); // Only for L2 - L8.
    TrackerNumberOfHitsXYInner = trackerTrack->NHxyInner();

    InnerTrackerCharge = trackerTrack->GetChargeAndError(3, 0).Charge(); // 3 means inner tracker.
    assert(TrackerCharges().empty());
    assert(TrackerChargesX().empty());
    assert(TrackerChargesY().empty());
    for (int layerIndex = 1; layerIndex <= 9; ++layerIndex) {
      TrackerCharges().push_back(trackerTrack->GetLayerCharge(layerIndex, 3, 0));  // take Y if X is not good
      TrackerChargesX().push_back(trackerTrack->GetLayerCharge(layerIndex, 1, 0)); // take X cluster
      TrackerChargesY().push_back(trackerTrack->GetLayerCharge(layerIndex, 2, 0)); // take Y cluster
    }

    if (AC::CurrentACQtVersion() >= AC::ACQtVersion(7, 7)) {
      assert(trackerTrack->RawCharges().size() > 0);
      const AC::TrackerCharge& innerTrackerRawCharge = trackerTrack->RawCharges().at(0);
      TrackerRawChargeFromDeDx = innerTrackerRawCharge.Charge(); // by choice the inner charge, RawCharges().at(1) would return the allTrackerRawCharge.

      TrackerChargeHL = trackerTrack->GetChargeAndErrorHuLiu(3).Charge();
    }

    const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
    const AC::ParticleHypothesis defaultParticleHypothesis = AC::DefaultMass;
    int trackFitIndex = trackerTrack->GetFitFuzzy(AC::Choutko, AC::Inner, refitPattern, defaultParticleHypothesis);
    if (trackFitIndex >= 0) {
      const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
      HasTrackerTrackChoutkoInnerFit = true;
      TrackerTrackChoutkoInnerRigidity = trackFit.Rigidity();
      TrackerTrackChoutkoInnerRigidityRelError = std::abs(trackFit.InverseRigidityError() * trackFit.Rigidity());
      TrackerTrackChoutkoInnerChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoInnerChiSquareY = trackFit.ChiSquareNormalizedY();
    }

    trackFitIndex = trackerTrack->GetFitFuzzy(AC::Choutko, AC::InnerPlusL1, refitPattern, defaultParticleHypothesis);
    if (trackFitIndex >= 0) {
      const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
      HasTrackerTrackChoutkoInnerPlusL1Fit = true;
      TrackerTrackChoutkoInnerPlusL1Rigidity = trackFit.Rigidity();
      TrackerTrackChoutkoInnerPlusL1RigidityRelError = std::abs(trackFit.InverseRigidityError() * trackFit.Rigidity());
      TrackerTrackChoutkoInnerPlusL1ChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoInnerPlusL1ChiSquareY = trackFit.ChiSquareNormalizedY();
    }

    trackFitIndex = trackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, defaultParticleHypothesis);
    if (trackFitIndex >= 0) {
      const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
      HasTrackerTrackChoutkoAllFit = true;
      TrackerTrackChoutkoAllRigidity = trackFit.Rigidity();
      TrackerTrackChoutkoAllRigidityRelError = std::abs(trackFit.InverseRigidityError() * trackFit.Rigidity());
      TrackerTrackChoutkoAllChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoAllChiSquareY = trackFit.ChiSquareNormalizedY();
    }

    trackFitIndex = trackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, AC::BetaFromTOF);
    if (trackFitIndex >= 0) {
      const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
      TrackerTrackChoutkoAllBetaTofRigidity = trackFit.Rigidity();
      TrackerTrackChoutkoAllBetaTofChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoAllBetaTofChiSquareY = trackFit.ChiSquareNormalizedY();
    }

    trackFitIndex = trackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, AC::BetaFromRICH);
    if (trackFitIndex >= 0) {
      const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
      TrackerTrackChoutkoAllBetaRichRigidity = trackFit.Rigidity();
      TrackerTrackChoutkoAllBetaRichChiSquareX = trackFit.ChiSquareNormalizedX();
      TrackerTrackChoutkoAllBetaRichChiSquareY = trackFit.ChiSquareNormalizedY();
    }

    if (AC::CurrentACQtVersion() >= AC::ACQtVersion(7, 7)) {
      const int refitPatternKalman = AC::PGMA + AC::RefitIfNeeded;
      trackFitIndex = trackerTrack->GetFitFuzzy(AC::Kalman, AC::All, refitPatternKalman, defaultParticleHypothesis);
      if (trackFitIndex >= 0) {
	const auto& trackFit = trackerTrack->TrackFits().at(trackFitIndex);
	HasTrackerTrackKalmanAllFit = true;
	TrackerTrackKalmanAllRigidity = trackFit.Rigidity();
	Rigidity = trackFit.Rigidity();
	TrackerTrackKalmanAllRigidityRelError = std::abs(trackFit.InverseRigidityError() * trackFit.Rigidity());
	TrackerTrackKalmanAllChiSquareX = trackFit.ChiSquareNormalizedX();
	TrackerTrackKalmanAllChiSquareY = trackFit.ChiSquareNormalizedY();
      }
    }

    const unsigned int numberOfTrackerLayers = 9;
    assert(trackerTrack->LayerYPatternBitset().size() == numberOfTrackerLayers+1);
    assert(trackerTrack->LayerXYPatternBitset().size() == numberOfTrackerLayers+1);

    for (unsigned int i = 0; i < numberOfTrackerLayers; ++i) {
      TrackerLayerYPattern().push_back(trackerTrack->LayerYPatternBitset().test(i+1));
      TrackerLayerXYPattern().push_back(trackerTrack->LayerXYPatternBitset().test(i+1));
      TrackerFitResidualsVector().push_back(-1000.0);
      TrackerFitResidualsVectorInner().push_back(-1000.0);
    }

    for (const AC::TrackerReconstructedHit& reconstructedHit : trackerTrack->ReconstructedHits()) {
      TrackerFitResidualsVector().at(int(reconstructedHit.Layer()) - 1) = reconstructedHit.ResidualYInMicrons();
      TrackerFitResidualsVectorInner().at(int(reconstructedHit.Layer()) - 1) = reconstructedHit.ResidualYWithRespectToInnerOnlyInMicrons();
    }

    TrackerPattern = primaryParticle->TrackerLayerPatternClassification();
    TrackerTrackIsNotInSolarArrayShadow = !primaryParticle->AmsParticle()->IsInSolarArrayShadow();

    // Store the position and direction of the TrackerTrack in the middle of AMS.
    const Analysis::SplineTrack* spline = primaryParticle->GetSplineTrack();
    assert(spline);

    const double positionZ = 0.0;
    Analysis::StraightLineTrack straightLineTrack = spline->TangentAtZ(positionZ);
    TrackerTrackX = straightLineTrack.X(positionZ);
    TrackerTrackY = straightLineTrack.Y(positionZ);
    TrackerTrackZ = positionZ;
    TrackerTrackTheta = straightLineTrack.Theta() * RadToDeg();
    TrackerTrackPhi = straightLineTrack.Phi() * RadToDeg();

    // Store the position of the TrackerTrack at the upper and lower ToF
    TrackerTrackUpTofX = spline->PositionAtZ(AC::AMSGeometry::ZTOFUpper).X();
    TrackerTrackUpTofY = spline->PositionAtZ(AC::AMSGeometry::ZTOFUpper).Y();
    TrackerTrackUpTofTheta = spline->PositionAtZ(AC::AMSGeometry::ZTOFUpper).Theta();
    TrackerTrackUpTofPhi = spline->PositionAtZ(AC::AMSGeometry::ZTOFUpper).Phi();

    TrackerTrackLowTofX = spline->PositionAtZ(AC::AMSGeometry::ZTOFLower).X();
    TrackerTrackLowTofY = spline->PositionAtZ(AC::AMSGeometry::ZTOFLower).Y();
    TrackerTrackLowTofTheta = spline->PositionAtZ(AC::AMSGeometry::ZTOFLower).Theta();
    TrackerTrackLowTofPhi = spline->PositionAtZ(AC::AMSGeometry::ZTOFLower).Phi();

    // Store the position of the TrackerTrack at the Agl and NaF radiators
    TrackerTrackAglTopX = spline->PositionAtZ(AC::AMSGeometry::ZRICHradiatorAglTop).X();
    TrackerTrackAglTopY = spline->PositionAtZ(AC::AMSGeometry::ZRICHradiatorAglTop).Y();

    TrackerTrackNaFTopX = spline->PositionAtZ(AC::AMSGeometry::ZRICHradiatorNaFTop).X();
    TrackerTrackNaFTopY = spline->PositionAtZ(AC::AMSGeometry::ZRICHradiatorNaFTop).Y();

    // Check that the track is in the fiducial volume. The values for the radii and Y come from QY's nuclei twiki page.
    const float fiducialRadii[8] = {62, 62, 46, 46, 46, 46, 46, 46}; // 8 and not 9 because for L1 there is a condition on X and not on R.
    const float fiducialY[9] = {47, 40, 44, 44, 36, 36, 44, 44, 29}; //
    const float fiducialX = 43;                                      // only one value instead of 9, simply because so are the cuts by QY.

    for (unsigned int i = 0; i < numberOfTrackerLayers; i++) {

      float splineX = spline->X(AC::AMSGeometry::ZTrackerLayer[i]);
      float splineY = spline->Y(AC::AMSGeometry::ZTrackerLayer[i]);
      float radius = sqrt(splineX * splineX + splineY * splineY);

      if (i < numberOfTrackerLayers - 1) {
	if (radius < fiducialRadii[i] && std::abs(splineY) < fiducialY[i])
	  TrackerTrackIsInFiducialLayer().push_back(true);
	else
	  TrackerTrackIsInFiducialLayer().push_back(false);
      }
      else {
	if (std::abs(splineX) < fiducialX && std::abs(splineY) < fiducialY[i])
	  TrackerTrackIsInFiducialLayer().push_back(true);
	else
	  TrackerTrackIsInFiducialLayer().push_back(false);
      }
    }
  }

  // TRD information
  TrdNumberOfSegmentsXZ = event.TrdSegmentsXZ().size();
  TrdNumberOfSegmentsYZ = event.TrdSegmentsYZ().size();

  // Count number of TRD vertices from both projections that match in Z position, within errors.
  UShort_t numberOfTrdVertices = 0;
  for (unsigned int i = 0; i < event.TrdVerticesXZ().size(); ++i) {
    const double xzSegmentZPosition = event.TrdVerticesXZ()[i].Z();
    const double xzSegmentZPositionUncertainty = event.TrdVerticesXZ()[i].ErrorZ();

    for (unsigned int j = 0; j < event.TrdVerticesYZ().size(); ++j) {
      const double yzSegmentZPosition = event.TrdVerticesYZ()[j].Z();
      const double yzSegmentZPositionUncertainty = event.TrdVerticesYZ()[j].ErrorZ();
      if (std::abs(xzSegmentZPosition - yzSegmentZPosition) > 2.0 * std::abs(xzSegmentZPositionUncertainty + yzSegmentZPositionUncertainty))
	continue;
      ++numberOfTrdVertices;
    }
  }
  TrdNumberOfVertices = numberOfTrdVertices;
  TrdNumberOfVerticesXZ = event.TrdVerticesXZ().size();
  TrdNumberOfVerticesYZ = event.TrdVerticesYZ().size();
  if (primaryParticle && primaryParticle->HasTrackerTrack()) {
    TrdNumberOfHitsOnTrack = primaryParticle->TrdHitsFromTrackerTrack().size();
    TrdNumberOfActiveLayers = primaryParticle->NumberOfTrdActiveLayers();
  }
  TrdNumberOfRawHits = event.NumberOfTrdRawHits();

  // RICH information
  RichNumberOfRings = event.NumberOfRichRings();
  if (primaryParticle && primaryParticle->RichRing()) {
    const auto* richRing = primaryParticle->RichRing();
    assert(richRing);
    HasRichRing = primaryParticle->HasRichRing();
    RichBeta = richRing->Beta();
    TofRichBetaMatch = TMath::Abs(richRing->Beta() -  primaryParticle->BetaTof()) / primaryParticle->BetaTof();
    RichBetaError = richRing->BetaError();
    RichCharge = richRing->ChargeEstimate();
    RichNExpectedPhotoElectrons = richRing->NExpectedPhotoElectrons();
    RichNCollectedPhotoElectrons = richRing->NCollectedPhotoElectrons();
    RichNPhotoElectrons = richRing->NPhotoElectrons();
    RichIsGood = richRing->IsGood();
    RichIsClean = (richRing->Status()&1) == 0;
    RichNumberOfPmts = richRing->PMTs();
    RichNumberOfHits = richRing->NumberOfHits();
    RichIsNaF = richRing->IsNaF();
    RichBetaConsistency = richRing->BetaConsistency();
    RichRingX = richRing->X();
    RichRingY= richRing->Y();
    RichRingProbability = richRing->Probability();
    RichTileIndex = richRing->TileIndex();
    RichDistanceToTileBorder = richRing->DistanceToTileBorder();

    // Area of NaF: 35x35 cm2.
    if (primaryParticle->HasTrackerTrack()) {

      const Analysis::SplineTrack* spline = primaryParticle->GetSplineTrack();
      assert(spline);

      // FIXME: Check that z-coordinates are precise!
      const Vector3& splinePositionAtRichNaFTop = spline->PositionAtZ(AC::AMSGeometry::ZRICHradiatorNaFTop);
      RichAtNaFTopPositionX = splinePositionAtRichNaFTop.X();
      RichAtNaFTopPositionY = splinePositionAtRichNaFTop.Y();
      const Vector3& splinePositionAtRichPmtsTop = spline->PositionAtZ(AC::AMSGeometry::ZRICHpmtsTop);
      RichAtPmtsTopPositionX = splinePositionAtRichPmtsTop.X();
      RichAtPmtsTopPositionY = splinePositionAtRichPmtsTop.Y();
    }
  }

#ifdef HAVE_AMS_SUPPORT
  if (!fFillAmsVariables)
    return;

  // Code only for AMS ROOT support
  AMSEventR* amsEvent = const_cast<AMSEventR*>(event.RawEvent()->AssociatedAMSEvent());
  if (!amsEvent) {
    WARN_OUT_WITH_EVENT << "Could not load AMS Event!" << std::endl;
    return;
  }

  if (event.IsMC()) {
    TrExtAlignDB::SmearExtAlign();
    TRCLFFKEY.UseSensorAlign = 0;
    TRCLFFKEY.ClusterCofGOpt = 1;
    TRFITFFKEY.Zshift = -1;

    // Charge dependent resolution in track fits, but only for newer MC.
    if (amsEvent->Version() >= 1106)
      TRFITFFKEY.ErcHeY = 0.0f;
  }
  else {
    // Activate eta uniformity correction for ISS and BT for all events.
    TrLinearEtaDB::SetLinearCluster();

    // Zshift correction, in particular for ions
    TRFITFFKEY.Zshift = 2;

    // Charge dependent resolution in track fits
    TRFITFFKEY.ErcHeY = 0.0f;
  }

  ParticleR* amsParticle = nullptr;
  if (amsEvent->nParticle() > 0)
    amsParticle = amsEvent->pParticle(0);

  BetaHR* amsBetaH = nullptr;
  if (amsParticle)
    amsBetaH = amsParticle->pBetaH();

  TrTrackR* amsTrackerTrack = nullptr;
  if (amsParticle)
    amsTrackerTrack = amsParticle->pTrTrack();

  if (amsParticle && amsBetaH) {
    for (int iLayer = 0; iLayer < 4; ++iLayer) {
      TofGoodQPath().push_back(amsBetaH->IsGoodQPathL(iLayer));
      TofRawChargesFromDeDx().push_back(amsBetaH->GetQL(iLayer, 2, (TofRecH::kThetaCor|TofRecH::kBirkCor|TofRecH::kReAttCor|TofRecH::kQ2Q)));
      TofBetaFromDeDx().push_back(amsBetaH->GetQBetaL(iLayer, 4.0));
    }
  }

  if (amsParticle && amsTrackerTrack) {
    // For ISS data
    if (TrTrackR::DefaultChargeVersion == 2) {
      TrackerRawChargeFromDeDx = TrCharge::GetCombinedMean(TrCharge::kInner|TrCharge::kTruncMean|TrCharge::kSqrt, amsTrackerTrack, 1,-1,
							   TrClusterR::kAsym|TrClusterR::kGain|TrClusterR::kAngle|
							   TrClusterR::kLoss|TrClusterR::kMIP, -1, 0).Mean;
    }
    // For B1106+ MC
    else if (TrTrackR::DefaultChargeVersion == 3) {
      TrackerRawChargeFromDeDx = TrCharge::GetCombinedMean(TrCharge::kInner|TrCharge::kTruncMean|TrCharge::kSqrt, amsTrackerTrack, 1,-1,
							   TrClusterR::kTotSign2017|TrClusterR::kSimAsym|TrClusterR::kSimSignal|
							   TrClusterR::kLoss|TrClusterR::kAngle, -1, 0).Mean;
    }
    else {
      WARN_OUT_WITH_EVENT << "TrTrackR::DefaultChargeVersion is unexpected: " << TrTrackR::DefaultChargeVersion << std::endl;
    }

    //    TrackerBetaFromDeDx = 0.94 * std::pow(1.0 / (TrackerRawChargeFromDeDx() / 4.0), 1.1);

    float charge = std::floor(amsTrackerTrack->GetInnerQ_all().Mean + 0.5);
    float mass = ParticleId::Mass(ParticleId::Alpha) * (charge / 2.0);
    float beta = 999.0f;

    int iTrTrackPar = amsTrackerTrack->iTrTrackPar(6, AC::All, AC::RebuildFromTDV + AC::PGMA, mass, charge, beta);
    if (iTrTrackPar > 0) {
      HasTrackerTrackKalmanAllFit = true;
      TrackerTrackKalmanAllRigidity = amsTrackerTrack->GetRigidity(iTrTrackPar);
      TrackerTrackKalmanAllChiSquareX = amsTrackerTrack->GetNormChisqX(iTrTrackPar);
      TrackerTrackKalmanAllChiSquareY = amsTrackerTrack->GetNormChisqY(iTrTrackPar);
      TrackerTrackKalmanAllRigidityRelError = std::abs(amsTrackerTrack->GetErrRinv(iTrTrackPar) * amsTrackerTrack->GetRigidity(iTrTrackPar));
    }
  }
#endif
}



float NucleiAnalysisTree::CalculateGeneratedMomentumAtGivenZ(const Analysis::Event& event, const float zPosition, double zCoordinateTolerance) {

  /*
    e.g. for Rich
    zPosition = -69.975; // cm
    see const double AMSG4SteppingAction::facc_pl in gbatch/CC/geant4.C
   */
  const AC::MCParticle* mcPrimaryParticle = event.RawEvent()->MC().PrimaryParticle();
  assert(mcPrimaryParticle);
  assert(mcPrimaryParticle->TrackID() == 1);

  float generatedMomentumAtZ = 0;

  for (const auto& generator : mcPrimaryParticle->EventGenerators()) {
    if (std::abs(generator.Z0() - zPosition) >= zCoordinateTolerance)
      continue;
    generatedMomentumAtZ = generator.Momentum();
  }
  return generatedMomentumAtZ;
}





void NucleiAnalysisTree::UpdateInMemoryBranches() {
  //  double rigidity = TrackerTrackChoutkoAllRigidity();

  
  
}


//  const float degToRad = M_PI / 180;
//  DeltaYOverDeltaZ = std::abs(sin(TrackerTrackPhi() * degToRad) * sin(TrackerTrackTheta() * degToRad) / cos(TrackerTrackTheta() * degToRad));

