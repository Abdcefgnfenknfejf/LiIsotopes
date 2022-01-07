#include "EcalBDTEfficiencyParameterization.hh"

#include <TF1.h>

#include <cassert>
#include <cmath>

static const double sEfficiencyUntilFirstBreak = 0.99;
static const double sFirstBreakEnergy = 200.0;

static const double sSecondBreakEnergy = 500.0;
static const double sEfficiencyAtSecondBreak = 0.94;

double EcalBDTEfficiencyParameterizationFunction(double* x, double* par) {

  double firstBreakEnergy = par[0];
  double efficiencyUntilFirstBreak = par[1];
  double secondBreakEnergy = par[2];
  double efficiencyAtSecondBreak = par[3];

  if (x[0] <= firstBreakEnergy)
    return efficiencyUntilFirstBreak;

  double slope = (efficiencyAtSecondBreak - efficiencyUntilFirstBreak) / (std::log10(secondBreakEnergy) - std::log10(firstBreakEnergy));
  double offset = 0.5 * ((efficiencyUntilFirstBreak + efficiencyAtSecondBreak) - slope * (std::log10(secondBreakEnergy) + std::log10(firstBreakEnergy)));
  return offset + slope * std::log10(x[0]);
}

TF1* EcalBDTEfficiencyParameterization() {

  static TF1* sEcalBDTCut = nullptr;
  if (!sEcalBDTCut) {
    sEcalBDTCut = new TF1("sEcalBDTCut", EcalBDTEfficiencyParameterizationFunction, 0.5, 2000.0, 4);
    sEcalBDTCut->SetParameters(sFirstBreakEnergy, sEfficiencyUntilFirstBreak, sSecondBreakEnergy, sEfficiencyAtSecondBreak);
    sEcalBDTCut->SetNpx(1e5);
  }
  return sEcalBDTCut;
}
