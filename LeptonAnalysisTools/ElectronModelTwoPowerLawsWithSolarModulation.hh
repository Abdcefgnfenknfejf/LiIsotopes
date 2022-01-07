#ifndef ElectronModelTwoPowerLawsWithSolarModulation_hh
#define ElectronModelTwoPowerLawsWithSolarModulation_hh

#include <utility>

#include "Model.hh"
#include "ModelParameter.hh"

class TF1;

class ElectronModelTwoPowerLawsWithSolarModulation : public Modelling::Model {

public:
  ElectronModelTwoPowerLawsWithSolarModulation(double _Epower = 3.0);
  virtual ~ElectronModelTwoPowerLawsWithSolarModulation() {}

  virtual TF1* GetPredictionForIdentifier(int id) const {

    if (id == 0)
      return ElecFlux;
    return 0;
  }

  double ElecFluxF(double* x, double* par);
  double ElecFluxLowEnergyTransition(double* x, double* par);
  double ElecFluxPowerLawAF(double* x, double* par);
  double ElecFluxPowerLawBF(double* x, double* par);
  double ElecFluxEpowerF(double* x, double* par);
  double ElecFluxPowerLawAEpowerF(double* x, double* par);
  double ElecFluxPowerLawBEpowerF(double* x, double* par);

  TF1* ElecFlux;
  TF1* ElecFluxPowerLawA;
  TF1* ElecFluxPowerLawB;
  TF1* ElecFluxEpower;
  TF1* ElecFluxPowerLawAEpower;
  TF1* ElecFluxPowerLawBEpower;

  double Epower = 3.0;
  double M = 0.511e-3;

  Modelling::ModelParameter phiMinus;

  Modelling::ModelParameter E_t;
  Modelling::ModelParameter DeltaGamma_t;

  Modelling::ModelParameter C_a;
  Modelling::ModelParameter E_a;
  Modelling::ModelParameter Gamma_a;

  Modelling::ModelParameter C_b;
  Modelling::ModelParameter E_b;
  Modelling::ModelParameter Gamma_b;
};

#endif
