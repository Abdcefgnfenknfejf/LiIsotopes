#ifndef TemplateFit2DFitFunction_hh
#define TemplateFit2DFitFunction_hh

#include <TMinuitMinimizer.h>

class TH2;

class TemplateFit2DFitFunction : public ROOT::Math::IMultiGenFunction {
public:
  enum class ContributingCategories {
    Electrons   = 1 << 1,
    Positrons   = 1 << 2,
    CCElectrons = 1 << 3,
    CCPositrons = 1 << 4,
    Protons     = 1 << 5,
    CCProtons   = 1 << 6
  };

  TemplateFit2DFitFunction(const TH2* dataDistribution,
                           const TH2* electronTemplate, const TH2* positronTemplate,
                           const TH2* ccElectronTemplate, const TH2* ccPositronTemplate,
                           const TH2* protonTemplate, const TH2* ccProtonTemplate);

  /** Number of parameters */
  virtual unsigned int NDim() const { return 5; }

  /** Defines the actual function call. */
  virtual double DoEval(const double* parameters) const;

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const { return new TemplateFit2DFitFunction(*this); }

  double PdfValue(int xbin, int ybin) const;
  TH2* CalculateHistogram(int contributingCategories) const;
  TH2* CalculateTotalHistogram() const;
  TH2* DrawZeroProbabilityEvents() const;

  void UpdateParameters(const double*);

private:
  const TH2* fDataDistribution { nullptr };

  const TH2* fElectronTemplate { nullptr };
  const TH2* fPositronTemplate { nullptr };
  const TH2* fCCElectronTemplate { nullptr };
  const TH2* fCCPositronTemplate { nullptr };
  const TH2* fProtonTemplate { nullptr };
  const TH2* fCCProtonTemplate { nullptr };

  double fNumberOfEvents { 0.0 };
  double fProtonFraction { 0.0 };
  double fCCProtonFraction { 0.0 };
  double fElectronFraction { 0.0 };
  double fPositronFraction { 0.0 };
  double fChargeConfusionProbability { 0.0 };
};

#endif
