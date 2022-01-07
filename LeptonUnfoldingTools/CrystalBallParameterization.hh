#ifndef CrystalBallParameterization_hh
#define CrystalBallParameterization_hh

#include <string>

#include "MigrationParameterization.hh"

class TCanvas;
class TF1;
class TGraph;

class CrystalBallParameterization : public MigrationParameterization {

public:
  CrystalBallParameterization(const std::string& paramFilename);
  virtual ~CrystalBallParameterization();

  virtual TCanvas* DrawParameterCanvas(const std::string& nameSuffix = "");

  virtual TF1* MakeFunction(double E) const;
  virtual bool IsPdf() const { return false; }

  double GaussianEnergyResolution(double E) const;

  double Alpha(double E) const;
  double N(double E) const;
  double Amp(double E) const;

private:

  double fGaussianStatisticTerm = 0.15;
  double fGaussianConstantTerm = 0.02;

  TGraph* fCrystalBallAlpha = nullptr;
  TGraph* fCrystalBallN = nullptr;
  TGraph* fCrystalBallAmp = nullptr;

};

#endif
