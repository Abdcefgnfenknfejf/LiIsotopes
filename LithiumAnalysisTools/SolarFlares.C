#include "SolarFlares.hh"
#include "Environment.hh"

#include <cassert>
#include <cmath>

#include <TFile.h>
#include <TGraph.h>

bool IsTimeDuringSolarFlare(double timeStamp, double magnitudeTreshold, double excludeSecondsBefore, double excludeSecondsAfter) {

  static TFile* sSolarXRayFluxFile = 0;
  if (!sSolarXRayFluxFile) {
    std::string solarXRayFluxFileName = "${LEPTONANALYSIS}/LookupFiles/SolarXRayFlux.root";
    Environment::ExpandEnvironmentVariables(solarXRayFluxFileName);
    sSolarXRayFluxFile = TFile::Open(solarXRayFluxFileName.c_str(), "OPEN");
  }

  static TGraph* sSolarFlareDurationGraph = 0;
  if (!sSolarFlareDurationGraph) {
    sSolarFlareDurationGraph = dynamic_cast<TGraph*>(sSolarXRayFluxFile->Get("grSolarFlareDuration"));
    assert(sSolarFlareDurationGraph);
  }

  static TGraph* sSolarFlareAmplitudeGraph = 0;
  if (!sSolarFlareAmplitudeGraph) {
    sSolarFlareAmplitudeGraph = dynamic_cast<TGraph*>(sSolarXRayFluxFile->Get("grSolarFlareAmplitude"));
    assert(sSolarFlareAmplitudeGraph);
  }

  assert(sSolarFlareDurationGraph->GetN() == sSolarFlareAmplitudeGraph->GetN());
  for (int point = 0; point < sSolarFlareDurationGraph->GetN(); ++point) {
    double xDuration = 0;
    double yDuration = 0;
    sSolarFlareDurationGraph->GetPoint(point, xDuration, yDuration);

    double xAmplitude = 0;
    double yAmplitude = 0;
    sSolarFlareAmplitudeGraph->GetPoint(point, xAmplitude, yAmplitude);

    assert(std::abs(xDuration - xAmplitude) < 1e-5);
    if (xDuration > timeStamp)
      break;

    double magnitude = yDuration * yAmplitude;
    if (magnitude < magnitudeTreshold)
      continue;

    if (timeStamp >= (xDuration - excludeSecondsBefore) && timeStamp < (xDuration + excludeSecondsAfter))
      return true;
  }

  return false;
}
