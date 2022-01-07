#ifndef SimplePowerLawModel_hh
#define SimplePowerLawModel_hh

#include <Model.hh>
#include <ModelParameter.hh>

class TF1;

class SimplePowerLawModel : public Modelling::Model {
public:
  SimplePowerLawModel(double _E0, double cStart, double gammaStart);
  virtual ~SimplePowerLawModel();

  virtual TF1* GetPredictionForIdentifier(int) const { return Flux; }

  double FluxF(double* x, double*);
  TF1* Flux;

  double E0;

  Modelling::ModelParameter C;
  Modelling::ModelParameter gamma;
};

#endif
