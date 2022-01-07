#ifndef ElectronPositronDoublePowerLawApproximation_hh
#define ElectronPositronDoublePowerLawApproximation_hh

#include <utility>

#include "Model.hh"
#include "ModelParameter.hh"

class TF1;

class ElectronPositronDoublePowerLawApproximation : public Modelling::Model {
public:
  enum StartValue {
    ElectronLowEnergy,
    PositronLowEnergy,
    PositronHighEnergy
  };

  ElectronPositronDoublePowerLawApproximation(double _Epower = 3.0, StartValue startValue = ElectronLowEnergy);
  virtual ~ElectronPositronDoublePowerLawApproximation() {}

  virtual TF1* GetPredictionForIdentifier(int id) const {

    if (id == 0)
      return Flux;
    return 0;
  }

  double FluxF(double* x, double* par);
  double FluxPowerLaw1F(double* x, double* par);
  double FluxPowerLaw2F(double* x, double* par);
  double FluxEpowerF(double* x, double* par);
  double FluxPowerLaw1EpowerF(double* x, double* par);
  double FluxPowerLaw2EpowerF(double* x, double* par);

  TF1* Flux;
  TF1* FluxPowerLaw1;
  TF1* FluxPowerLaw2;
  TF1* FluxEpower;
  TF1* FluxPowerLaw1Epower;
  TF1* FluxPowerLaw2Epower;

  double Epower = 3.0;

  Modelling::ModelParameter C;
  Modelling::ModelParameter E_0;
  Modelling::ModelParameter E_n;
  Modelling::ModelParameter Gamma;
  Modelling::ModelParameter Delta_Gamma;
};

#endif
