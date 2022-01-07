#include "FluxModel.hh"

#include "ModelFunctions.hh"
#include "ParticleId.hh"

#include <cmath>

// Parameters and functions from fit to AMS-02 data performed by S. Schael and H. Gast.
const double FluxModel::Cminus             =  4.21318e+02;
const double FluxModel::gammaMinus         =  3.26895e+00;
const double FluxModel::deltaGamma         = -5.74649e-01;
const double FluxModel::lambdaMinus        =  2.50819e-01;
const double FluxModel::inverseBreakEnergy =  2.89656e-02;
const double FluxModel::Cplus              =  9.15804e+01;
const double FluxModel::gammaPlus          =  3.75985e+00;
const double FluxModel::Csource            =  1.20153e+00;
const double FluxModel::gammaSource        =  2.41208e+00;
const double FluxModel::lambdaSource       =  1.85000e-03;
const double FluxModel::phiMinus           =  1.40036e+00;
const double FluxModel::phiPlus            =  9.90855e-01;
const double FluxModel::referenceEnergy    = 1.;
const double FluxModel::integrationFactor  = 1.001;

double FluxModel::GetPositronDiffuseFlux(double energy) {
  return Modelling::PowerLawSolarMod(energy, Cplus, gammaPlus, referenceEnergy, phiPlus, ParticleId::Mass(ParticleId::Positron));
}

double FluxModel::GetElectronDiffuseFlux(double energy) {
  return Modelling::SmoothlyBrokenPowerLawSolarModLambda(energy, Cminus, gammaMinus, deltaGamma, lambdaMinus, inverseBreakEnergy, referenceEnergy, phiMinus, ParticleId::Mass(ParticleId::Electron));
}

double FluxModel::GetPositronSourceFlux(double energy) {
  return Modelling::ExpCutoffPowerLawSolarMod(energy, Csource, gammaSource, lambdaSource, referenceEnergy, phiPlus, ParticleId::Mass(ParticleId::Positron));
}

double FluxModel::GetElectronSourceFlux(double energy) {
  return Modelling::ExpCutoffPowerLawSolarMod(energy, Csource, gammaSource, lambdaSource, referenceEnergy, phiMinus, ParticleId::Mass(ParticleId::Electron));
}

double FluxModel::GetPositronDiffuseFluxIntegral(double energyStart, double energyEnd) {
  double Estart = energyStart;
  double Eend   = energyStart * integrationFactor;
  double integral = 0.;
  while (Eend < energyEnd) {
    double phi1 = GetPositronDiffuseFlux(Estart);
    double phi2 = GetPositronDiffuseFlux(Eend);
    integral += (phi1 + phi2) / 2. * (Eend - Estart);
    Estart = Eend;
    Eend *= integrationFactor;
  }
  return integral;
}

double FluxModel::GetElectronDiffuseFluxIntegral(double energyStart, double energyEnd) {
  double Estart = energyStart;
  double Eend   = energyStart * integrationFactor;
  double integral = 0.;
  while (Eend < energyEnd) {
    double phi1 = GetElectronDiffuseFlux(Estart);
    double phi2 = GetElectronDiffuseFlux(Eend);
    integral += (phi1 + phi2) / 2. * (Eend - Estart);
    Estart = Eend;
    Eend *= integrationFactor;
  }
  return integral;
}

double FluxModel::GetPositronSourceFluxIntegral(double energyStart, double energyEnd) {
  double Estart = energyStart;
  double Eend   = energyStart * integrationFactor;
  double integral = 0.;
  while (Eend < energyEnd) {
    double phi1 = GetPositronSourceFlux(Estart);
    double phi2 = GetPositronSourceFlux(Eend);
    integral += (phi1 + phi2) / 2. * (Eend - Estart);
    Estart = Eend;
    Eend *= integrationFactor;
  }
  return integral;
}

double FluxModel::GetElectronSourceFluxIntegral(double energyStart, double energyEnd) {
  double Estart = energyStart;
  double Eend   = energyStart * integrationFactor;
  double integral = 0.;
  while (Eend < energyEnd) {
    double phi1 = GetElectronSourceFlux(Estart);
    double phi2 = GetElectronSourceFlux(Eend);
    integral += (phi1 + phi2) / 2. * (Eend - Estart);
    Estart = Eend;
    Eend *= integrationFactor;
  }
  return integral;
}
