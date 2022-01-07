#include "SimplePowerLawModel.hh"

#include <TF1.h>

#include <cmath>

#define INFO_OUT_TAG "SimplePowerLawModel"
#include "debugging.hh"

SimplePowerLawModel::SimplePowerLawModel(double _E0, double cStart, double gammaStart)
  : Modelling::Model()
  , E0(_E0) {

  DefineParameter(C = Modelling::ModelParameter("C", cStart, 0.005 * cStart));
  DefineParameter(gamma = Modelling::ModelParameter("#gamma", gammaStart, 0.02));

  Flux = StyleTF1(new TF1("SimplePowerLawFlux", this, &SimplePowerLawModel::FluxF, 0.01, 10000.0, 0, "SimplePowerLawModel", "FluxF"), kRed);
}

SimplePowerLawModel::~SimplePowerLawModel() {

  delete Flux;
}

double SimplePowerLawModel::FluxF(double* x, double*) {

  return C * std::pow(x[0] / E0, -gamma);
}
