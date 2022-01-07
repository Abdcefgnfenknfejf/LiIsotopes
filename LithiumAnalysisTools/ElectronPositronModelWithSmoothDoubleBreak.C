#include "ElectronPositronModelWithSmoothDoubleBreak.hh"

#include <TF1.h>

#include "FitFunction.hh"
#include "HistogramDataset.hh"
#include "ModelAnalysis.hh"
#include "ModelFunctions.hh"
#include "ModellingData.hh"
#include "Utilities.hh"

ElectronPositronModelWithSmoothDoubleBreak::ElectronPositronModelWithSmoothDoubleBreak(double _Epower)
  : Modelling::Model()
  , Epower(_Epower) {

  // start values
  double Nposi = 0.216 * std::pow(E0 / 5., -3.75);
  double normS = 3.28e-2 * std::pow(10., -2.4) * std::pow(E1 / 50., -2.4);

  DefineParameter(Cplus = Modelling::ModelParameter("C_{+}", Nposi, 0.1 * Nposi, 0., 120.));
  DefineParameter(gammaPlus = Modelling::ModelParameter("#gamma_{+}", 3.83, 0.05));
  DefineParameter(Csource = Modelling::ModelParameter("C_{S}", normS, 0.1 * normS, 0., 100.));
  DefineParameter(gammaSource = Modelling::ModelParameter("#gamma_{S}", 2.53, 0.05));
  DefineParameter(lambdaSource = Modelling::ModelParameter("#lambda_{S}", 1.02e-3, 1.e-4));
  DefineParameter(phiMinus = Modelling::ModelParameter("#phi_{-}", 0.8, 0.01));
  DefineParameter(phiPlus = Modelling::ModelParameter("#phi_{+}", 1.02, 0.05));

  DefineParameter(C_d = Modelling::ModelParameter("C_{-}", 3.9, 0.03));
  DefineParameter(Gamma_d = Modelling::ModelParameter("#gamma_{-}", -2.56, 0.01));
  DefineParameter(DeltaG1 = Modelling::ModelParameter("#delta#gamma_{1}", -0.7, 0.01));
  DefineParameter(S1 = Modelling::ModelParameter("S1", -0.3, 0.001));
  DefineParameter(Eb1 = Modelling::ModelParameter("Eb1", 2.0, 0.1));
  DefineParameter(DeltaG2 = Modelling::ModelParameter("#delta#gamma_{2}", 0.5, 0.01));
  DefineParameter(S2 = Modelling::ModelParameter("S2", -0.285, 0.001));
  DefineParameter(Eb2 = Modelling::ModelParameter("Eb2", 30.0, 0.1));

  S1.Fix();
  Eb1.Fix();
  Eb2.Fix();

  double Emin = 0.5;
  double Emax = 1500.;

  Color_t posiColor = kRed;
  Color_t elecColor = kBlue;

  PosiFlux = StyleTF1(new TF1("PosiFlux", this, &ElectronPositronModelWithSmoothDoubleBreak::PosiFluxF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "PosiFluxF"), posiColor);
  ElecFlux = StyleTF1(new TF1("ElecFlux", this, &ElectronPositronModelWithSmoothDoubleBreak::ElecFluxF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "ElecFluxF"), elecColor);

  Posfrac = StyleTF1(new TF1("Posfrac", this, &ElectronPositronModelWithSmoothDoubleBreak::PosfracF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "PosfracF"), kRed);
  Ratio = StyleTF1(new TF1("Ratio", this, &ElectronPositronModelWithSmoothDoubleBreak::RatioF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "RatioF"), kBlack);
  AllElec = StyleTF1(new TF1("AllElec", this, &ElectronPositronModelWithSmoothDoubleBreak::AllElecF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "AllElecF"), kBlack);
  PosiFluxEpower = StyleTF1(new TF1("PosiFluxEpower", this, &ElectronPositronModelWithSmoothDoubleBreak::PosiFluxEpowerF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "PosiFluxEpowerF"), kBlack);
  ElecFluxEpower = StyleTF1(new TF1("ElecFluxEpower", this, &ElectronPositronModelWithSmoothDoubleBreak::ElecFluxEpowerF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "ElecFluxEpowerF"), kBlack);
  AllElecEpower = StyleTF1(new TF1("AllElecEpower", this, &ElectronPositronModelWithSmoothDoubleBreak::AllElecEpowerF, Emin, Emax, 0, "ElectronPositronModelWithSmoothDoubleBreak", "AllElecEpowerF"), kBlack);
}

double ElectronPositronModelWithSmoothDoubleBreak::PosiFluxF(double* x, double*) {

  double E = x[0];

  return Modelling::PowerLawSolarMod(E, Cplus, gammaPlus, E0, phiPlus, M) +
         Modelling::ExpCutoffPowerLawSolarMod(E, Csource, gammaSource, lambdaSource, E1, phiPlus, M);
}

double ElectronPositronModelWithSmoothDoubleBreak::ElecFluxF(double* x, double*) {

  double E = x[0];
  double Ehat = E + phiMinus;

  double Diffuse    = C_d * std::pow(Ehat/E0,Gamma_d);
  double One = 1.0;
  if (DeltaG2!=0 && S2<0) One = std::pow(1.0 + std::pow(Ehat/Eb2,DeltaG2/S2),S2);
  if (S1<0) Diffuse *= std::pow(1.0 + std::pow(Ehat/Eb1*One,DeltaG1/S1),S1);

  return Modelling::SolarModulationTerm(E, phiMinus, M) * Diffuse +
         Modelling::ExpCutoffPowerLawSolarMod(x[0], Csource, gammaSource, lambdaSource, E1, phiMinus, M);
}

double ElectronPositronModelWithSmoothDoubleBreak::PosfracF(double* x, double* par) {

  return PosiFluxF(x, par) / (PosiFluxF(x, par) + ElecFluxF(x, par));
}

double ElectronPositronModelWithSmoothDoubleBreak::RatioF(double* x, double* par) {

  return PosiFluxF(x, par) / ElecFluxF(x, par);
}

double ElectronPositronModelWithSmoothDoubleBreak::AllElecF(double* x, double* par) {

  return PosiFluxF(x, par) + ElecFluxF(x, par);
}

double ElectronPositronModelWithSmoothDoubleBreak::AllElecEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * AllElecF(x, par);
}

double ElectronPositronModelWithSmoothDoubleBreak::PosiFluxEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * PosiFluxF(x, par);
}

double ElectronPositronModelWithSmoothDoubleBreak::ElecFluxEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * ElecFluxF(x, par);
}
