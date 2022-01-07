#ifndef PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff_hh
#define PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff_hh

#include <utility>

#include "Model.hh"
#include "ModelParameter.hh"

class TF1;

class PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff : public Modelling::Model {

public:
  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff(double _Epower = 3.0);
  virtual ~PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff() {}

  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff* DeepCopy() const;

  virtual TF1* GetPredictionForIdentifier(int id) const {

    if (id == 0)
      return PosiFlux;
    return 0;
  }

  double PosiFluxF(double* x, double* par);
  double PosiFluxPowerLawDF(double* x, double* par);
  double PosiFluxPowerLawSF(double* x, double* par);
  double PosiFluxEpowerF(double* x, double* par);
  double PosiFluxPowerLawDEpowerF(double* x, double* par);
  double PosiFluxPowerLawSEpowerF(double* x, double* par);

  TF1* PosiFlux;
  TF1* PosiFluxPowerLawD;
  TF1* PosiFluxPowerLawS;
  TF1* PosiFluxEpower;
  TF1* PosiFluxPowerLawDEpower;
  TF1* PosiFluxPowerLawSEpower;

  double Epower = 3.0;
  double M = 0.511e-3;

  Modelling::ModelParameter phiPlus;

  Modelling::ModelParameter C_d;
  Modelling::ModelParameter E_1;
  Modelling::ModelParameter Gamma_d;

  Modelling::ModelParameter C_s;
  Modelling::ModelParameter E_2;
  Modelling::ModelParameter invE_s;
  Modelling::ModelParameter Gamma_s;
};

#endif
