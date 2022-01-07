#ifndef FluxModel_hh
#define FluxModel_hh

class FluxModel {
private:
  static const double Cminus;
  static const double gammaMinus;
  static const double deltaGamma;
  static const double lambdaMinus;
  static const double inverseBreakEnergy;
  static const double Cplus;
  static const double gammaPlus;
  static const double Csource;
  static const double gammaSource;
  static const double lambdaSource;
  static const double phiMinus;
  static const double phiPlus;
  static const double referenceEnergy;
  static const double integrationFactor;

public:
  static double GetPositronDiffuseFlux(double energy);
  static double GetElectronDiffuseFlux(double energy);
  static double GetPositronSourceFlux(double energy);
  static double GetElectronSourceFlux(double energy);
  
  static double GetPositronDiffuseFluxIntegral(double energyStart, double energyEnd);
  static double GetElectronDiffuseFluxIntegral(double energyStart, double energyEnd);
  static double GetPositronSourceFluxIntegral(double energyStart, double energyEnd);
  static double GetElectronSourceFluxIntegral(double energyStart, double energyEnd);
};

#endif
