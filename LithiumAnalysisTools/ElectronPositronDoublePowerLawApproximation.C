#include "ElectronPositronDoublePowerLawApproximation.hh"

#include <TF1.h>

#include "FitFunction.hh"
#include "HistogramDataset.hh"
#include "ModelAnalysis.hh"
#include "ModelFunctions.hh"
#include "ModellingData.hh"
#include "Utilities.hh"

ElectronPositronDoublePowerLawApproximation::ElectronPositronDoublePowerLawApproximation(double _Epower, StartValue startValue)
  : Modelling::Model()
  , Epower(_Epower) {

  if (startValue == ElectronLowEnergy) {
    DefineParameter(C = Modelling::ModelParameter("C", 2.33e-2, 0.01));
    DefineParameter(E_0 = Modelling::ModelParameter("E_{0}", 42.1, 0.1));
    DefineParameter(E_n = Modelling::ModelParameter("E_{n}", 21.13, 0.2));
    DefineParameter(Gamma = Modelling::ModelParameter("#gamma", -3.28, 0.01));
    DefineParameter(Delta_Gamma = Modelling::ModelParameter("#delta#gamma", 0.098, 0.01));
  } else if (startValue == PositronLowEnergy) {
    DefineParameter(C = Modelling::ModelParameter("C", 7.6e-5, 0.01));
    DefineParameter(E_0 = Modelling::ModelParameter("E_{0}", 25.2, 0.1));
    DefineParameter(E_n = Modelling::ModelParameter("E_{n}", 55.58, 0.2));
    DefineParameter(Gamma = Modelling::ModelParameter("#gamma", -2.97, 0.01));
    DefineParameter(Delta_Gamma = Modelling::ModelParameter("#delta#gamma", 0.13, 0.01));
  } else if (startValue == PositronHighEnergy) {
    DefineParameter(C = Modelling::ModelParameter("C", 8.8e-5, 0.01));
    DefineParameter(E_0 = Modelling::ModelParameter("E_{0}", 284.0, 0.1));
    DefineParameter(E_n = Modelling::ModelParameter("E_{n}", 63.02, 0.2));
    DefineParameter(Gamma = Modelling::ModelParameter("#gamma", -2.80, 0.01));
    DefineParameter(Delta_Gamma = Modelling::ModelParameter("#delta#gamma", -0.60, 0.01));
  }

  C.SetUnit("[\\meter\\squared~\\steradian~\\second~\\GeV]^{-1}");
  E_0.SetUnit("\\GeV");
  E_n.SetUnit("\\GeV");
  E_n.Fix();

  double Emin = 0.5;
  double Emax = 1500.;

  Flux = StyleTF1(new TF1("Flux", this, &ElectronPositronDoublePowerLawApproximation::FluxF, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxF"), kGreen + 3);
  FluxPowerLaw1 = StyleTF1(new TF1("FluxPowerLaw1", this, &ElectronPositronDoublePowerLawApproximation::FluxPowerLaw1F, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxPowerLaw1F"), kGray);
  FluxPowerLaw2 = StyleTF1(new TF1("FluxPowerLaw2", this, &ElectronPositronDoublePowerLawApproximation::FluxPowerLaw2F, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxPowerLaw2F"), kCyan - 4);
  FluxEpower = StyleTF1(new TF1("FluxEpower", this, &ElectronPositronDoublePowerLawApproximation::FluxEpowerF, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxEpowerF"), kGreen + 3);
  FluxPowerLaw1Epower = StyleTF1(new TF1("FluxPowerLaw1Epower", this, &ElectronPositronDoublePowerLawApproximation::FluxPowerLaw1EpowerF, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxPowerLaw1EpowerF"), kGray);
  FluxPowerLaw2Epower = StyleTF1(new TF1("FluxPowerLaw2Epower", this, &ElectronPositronDoublePowerLawApproximation::FluxPowerLaw2EpowerF, Emin, Emax, 0, "ElectronPositronDoublePowerLawApproximation", "FluxPowerLaw2EpowerF"), kCyan - 4);

  Flux->SetNpx(1e4);
  FluxPowerLaw1->SetNpx(1e4);
  FluxPowerLaw2->SetNpx(1e4);
  FluxEpower->SetNpx(1e4);
  FluxPowerLaw1Epower->SetNpx(1e4);
  FluxPowerLaw2Epower->SetNpx(1e4);
}

double ElectronPositronDoublePowerLawApproximation::FluxF(double* x, double* par) {

  double E = x[0];
  return E <= E_0 ? FluxPowerLaw1F(x, par) : FluxPowerLaw2F(x, par);
}

double ElectronPositronDoublePowerLawApproximation::FluxPowerLaw1F(double* x, double*) {

  double E = x[0];
  return C * std::pow(E / E_n, Gamma);
}

double ElectronPositronDoublePowerLawApproximation::FluxPowerLaw2F(double* x, double*) {

  double E = x[0];
  return C * std::pow(E / E_n, Gamma) * std::pow(E / E_0, Delta_Gamma);
}

double ElectronPositronDoublePowerLawApproximation::FluxEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * FluxF(x, par);
}

double ElectronPositronDoublePowerLawApproximation::FluxPowerLaw1EpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * FluxPowerLaw1F(x, par);
}

double ElectronPositronDoublePowerLawApproximation::FluxPowerLaw2EpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * FluxPowerLaw2F(x, par);
}
