#ifndef ElectronPositronModelWithSmoothDoubleBreak_hh
#define ElectronPositronModelWithSmoothDoubleBreak_hh

#include <utility>

#include "Model.hh"
#include "ModelParameter.hh"

class TF1;
class TH1D;

class ElectronPositronModelWithSmoothDoubleBreak : public Modelling::Model {

public:
  ElectronPositronModelWithSmoothDoubleBreak(double _Epower = 3.0);
  virtual ~ElectronPositronModelWithSmoothDoubleBreak() {}

  virtual TF1* GetPredictionForIdentifier(int id) const {

    if (id == 0)
      return ElecFlux;
    if (id == 1)
      return PosiFlux;
    if (id == 2)
      return AllElec;
    if (id == 3)
      return Posfrac;
    if (id == 4)
      return Ratio;
    return 0;
  }

  double PosiFluxF(double* x, double* par);
  double ElecFluxF(double* x, double* par);

  double PosfracF(double* x, double* par);
  double RatioF(double* x, double* par);
  double AllElecF(double* x, double* par);

  double PosiFluxEpowerF(double* x, double* par);
  double ElecFluxEpowerF(double* x, double* par);
  double AllElecEpowerF(double* x, double* par);

  TF1* PosiFlux;
  TF1* ElecFlux;

  TF1* Posfrac;
  TF1* Ratio;
  TF1* AllElec;

  TF1* PosiFluxEpower;
  TF1* ElecFluxEpower;
  TF1* AllElecEpower;

  double E0 = 5.0;
  double E1 = 60.0;
  double M = 0.511e-3;

  double Epower = 3.0;

  Modelling::ModelParameter Cplus;
  Modelling::ModelParameter gammaPlus;
  Modelling::ModelParameter Csource;
  Modelling::ModelParameter gammaSource;
  Modelling::ModelParameter lambdaSource;
  Modelling::ModelParameter phiMinus;
  Modelling::ModelParameter phiPlus;

  Modelling::ModelParameter C_d;
  Modelling::ModelParameter Gamma_d;
  Modelling::ModelParameter DeltaG1;
  Modelling::ModelParameter S1;
  Modelling::ModelParameter Eb1;
  Modelling::ModelParameter DeltaG2;
  Modelling::ModelParameter S2;
  Modelling::ModelParameter Eb2;
};

#endif
