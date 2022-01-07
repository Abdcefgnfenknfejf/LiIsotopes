#include "ElectronModelTwoPowerLawsWithSolarModulation.hh"

#include <TF1.h>

#include "FitFunction.hh"
#include "HistogramDataset.hh"
#include "ModelAnalysis.hh"
#include "ModelFunctions.hh"
#include "ModellingData.hh"
#include "Utilities.hh"

ElectronModelTwoPowerLawsWithSolarModulation::ElectronModelTwoPowerLawsWithSolarModulation(double _Epower)
  : Modelling::Model()
  , Epower(_Epower) {

  DefineParameter(phiMinus = Modelling::ModelParameter("#phi_{-}", 0.87, 0.01));

  DefineParameter(E_t = Modelling::ModelParameter("E_{t}", 3.94, 0.1));
  DefineParameter(DeltaGamma_t = Modelling::ModelParameter("#Delta#gamma_{t}", -2.14, 0.01));

  DefineParameter(C_a = Modelling::ModelParameter("C_{a}", 1.13e-2, 0.01));
  DefineParameter(E_a = Modelling::ModelParameter("E_{a}", 20.0, 0.1));
  DefineParameter(Gamma_a = Modelling::ModelParameter("#gamma_{a}", -4.31, 0.01));

  DefineParameter(C_b = Modelling::ModelParameter("C_{b}", 3.96e-6, 0.01));
  DefineParameter(E_b = Modelling::ModelParameter("E_{b}", 300.0, 0.1));
  DefineParameter(Gamma_b = Modelling::ModelParameter("#gamma_{b}", -3.14, 0.01));

  C_a.SetUnit("[\\meter\\squared~\\steradian~\\second~\\GeV]^{-1}");
  C_b.SetUnit("[\\meter\\squared~\\steradian~\\second~\\GeV]^{-1}");
  E_t.SetUnit("\\GeV");
  E_a.SetUnit("\\GeV");
  E_b.SetUnit("\\GeV");
  phiMinus.SetUnit("\\GeV");

  E_a.Fix();
  E_b.Fix();

  double Emin = 0.5;
  double Emax = 1500.;

  ElecFlux = StyleTF1(new TF1("ElecFlux", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxF"), kGreen + 3);
  ElecFluxPowerLawA = StyleTF1(new TF1("ElecFluxPowerLawA", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawAF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxPowerLawAF"), kGray);
  ElecFluxPowerLawB = StyleTF1(new TF1("ElecFluxPowerLawB", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawBF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxPowerLawBF"), kCyan - 4);
  ElecFluxEpower = StyleTF1(new TF1("ElecFluxEpower", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxEpowerF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxEpowerF"), kGreen + 3);
  ElecFluxPowerLawAEpower = StyleTF1(new TF1("ElecFluxPowerLawAEpower", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawAEpowerF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxPowerLawAEpowerF"), kGray);
  ElecFluxPowerLawBEpower = StyleTF1(new TF1("ElecFluxPowerLawBEpower", this, &ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawBEpowerF, Emin, Emax, 0, "ElectronModelTwoPowerLawsWithSolarModulation", "ElecFluxPowerLawBEpowerF"), kCyan - 4);

  ElecFlux->SetNpx(1e4);
  ElecFluxPowerLawA->SetNpx(1e4);
  ElecFluxPowerLawB->SetNpx(1e4);
  ElecFluxEpower->SetNpx(1e4);
  ElecFluxPowerLawAEpower->SetNpx(1e4);
  ElecFluxPowerLawBEpower->SetNpx(1e4);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxF(double* x, double* par) {

  return ElecFluxPowerLawAF(x, par) + ElecFluxPowerLawBF(x, par);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxLowEnergyTransition(double* x, double*) {

  double E = x[0];
  double Ehat = E + phiMinus;
  return 1.0 / (1.0 + std::pow(Ehat / E_t, DeltaGamma_t));
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawAF(double* x, double* par) {

  double E = x[0];
  double Ehat = E + phiMinus;
  return ElecFluxLowEnergyTransition(x, par) * Modelling::SolarModulationTerm(E, phiMinus, M) * C_a * std::pow(Ehat / E_a, Gamma_a);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawBF(double* x, double* par) {

  double E = x[0];
  double Ehat = E + phiMinus;
  return ElecFluxLowEnergyTransition(x, par) * Modelling::SolarModulationTerm(E, phiMinus, M) * C_b * std::pow(Ehat / E_b, Gamma_b);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * ElecFluxF(x, par);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawAEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * ElecFluxPowerLawAF(x, par);
}

double ElectronModelTwoPowerLawsWithSolarModulation::ElecFluxPowerLawBEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * ElecFluxPowerLawBF(x, par);
}
