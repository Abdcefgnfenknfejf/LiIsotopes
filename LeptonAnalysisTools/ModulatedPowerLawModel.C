#include "ModulatedPowerLawModel.hh"

#include "ParticleId.hh"

#include <ModelFunctions.hh>

#include <TF1.h>

#include <cmath>

#define INFO_OUT_TAG "ModulatedPowerLawModel"
#include "debugging.hh"

ModulatedPowerLawModel::ModulatedPowerLawModel(double _E0, double cStart, double gammaStart, double phiStart)
  : Modelling::Model()
  , E0(_E0) {

  DefineParameter(C = Modelling::ModelParameter("C", cStart, 0.005 * cStart));
  DefineParameter(gamma = Modelling::ModelParameter("#gamma", gammaStart, 0.02));
  DefineParameter(phi = Modelling::ModelParameter("#phi", phiStart, 0.01 * phiStart));

  Flux = StyleTF1(new TF1("ModulatedPowerLawModel", this, &ModulatedPowerLawModel::FluxF, 0.5, 1000.0, 0, "ModulatedPowerLawModel", "FluxF"), kRed);
  FluxE3 = StyleTF1(new TF1("ModulatedPowerLawModelE3", this, &ModulatedPowerLawModel::FluxE3F, 0.5, 1000.0, 0, "ModulatedPowerLawModel", "FluxE3F"), kRed);
}

ModulatedPowerLawModel::~ModulatedPowerLawModel() {

  delete Flux;
  delete FluxE3;
}

double ModulatedPowerLawModel::FluxF(double* x, double*) {

  return Modelling::PowerLawSolarMod(x[0], C, gamma, E0, phi, ParticleId::Mass(ParticleId::Electron));
}

double ModulatedPowerLawModel::FluxE3F(double* x, double* par) {

  return std::pow(x[0], 3) * FluxF(x, par);
}
