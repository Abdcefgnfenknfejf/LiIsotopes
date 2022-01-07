#ifndef CompositeParameterization_hh
#define CompositeParameterization_hh

#include <string>

#include "MigrationParameterization.hh"

class TCanvas;
class TF1;
class TGraph;

class CompositeParameterization : public MigrationParameterization {

public:
  CompositeParameterization(const std::string& paramFilename);
  CompositeParameterization(const CompositeParameterization& p);
  virtual ~CompositeParameterization();

  virtual TCanvas* DrawParameterCanvas(const std::string& nameSuffix = "");

  virtual void Dump() const;

  virtual TF1* MakeFunction(double E) const;
  virtual bool IsPdf() const { return true; }

  virtual double UnderflowProbability(double E) const;
  virtual double OverflowProbability(double E) const;

  virtual void Randomize(unsigned int seed = 0, bool randomize = true, double nSigma = 0.0);

  TGraph* GraphLandauMPV() const;
  TGraph* GraphLandauSigma() const;
  TGraph* GraphCbAlpha() const;
  TGraph* GraphCbN() const;
  TGraph* GraphCbMu() const;
  TGraph* GraphCbRelSigma() const;
  TGraph* GraphLandauFraction() const;
  TGraph* GraphProbOverflow() const;
  TGraph* GraphProbUnderflow() const;

  static TGraph* B1091_EnergyResolution();

private:

  TGraph* grLandauMPV = nullptr;
  TGraph* grLandauSigma = nullptr;
  TGraph* grCbAlpha = nullptr;
  TGraph* grCbN = nullptr;
  TGraph* grCbMu = nullptr;
  TGraph* grCbRelSigma = nullptr;
  TGraph* grLandauFraction = nullptr;
  TGraph* grProbUnderflow = nullptr;
  TGraph* grProbOverflow = nullptr;

  double fEmin = 0.5;
  double fEmax = 1.e3;

};

#endif
