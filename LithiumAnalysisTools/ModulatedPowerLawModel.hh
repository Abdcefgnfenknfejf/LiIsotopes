#ifndef ModulatedPowerLawModel_hh
#define ModulatedPowerLawModel_hh

#include <Model.hh>
#include <ModelParameter.hh>

class TF1;

class ModulatedPowerLawModel : public Modelling::Model {
public:
  ModulatedPowerLawModel(double _E0, double cStart, double gammaStart, double phiStart);
  virtual ~ModulatedPowerLawModel();

  virtual TF1* GetPredictionForIdentifier(int) const { return Flux; }

  double FluxF(double* x, double*);
  TF1* Flux;

  double FluxE3F(double* x, double*);
  TF1* FluxE3;

  double E0;
  double M;

  Modelling::ModelParameter C;
  Modelling::ModelParameter gamma;
  Modelling::ModelParameter phi;
};

#endif
