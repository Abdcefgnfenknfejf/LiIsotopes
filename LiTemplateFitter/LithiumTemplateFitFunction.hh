#ifndef LithiumTemplateFitFunction_hh
#define LithiumTemplateFitFunction_hh

#include <TMinuitMinimizer.h>

#include "Quantity.hh"

class TH1;

enum FitType {FitRelativeFractions, FitRatios};

/* FitRelativeFraction: NTotalEvents, n7, n9, n10 */
/* FitRatios: NTotalEvents, n7/n9, n9, n10/n9 */
/* FitSumAndRatio: NTotalEvents, n7, n9 + n10, n10/n9 */

class LithiumTemplateFitFunction : public ROOT::Math::IMultiGenFunction {
public:
  LithiumTemplateFitFunction(const TH1* dataDistribution,
			       const TH1* li6Template, const TH1* li7Template,
			       double rangeMinimum, double rangeMaximum, FitType fitType);

  /** Number of parameters */
  virtual unsigned int NDim() const { return 2; }

  /** Defines the actual function call. */
  virtual double DoEval(const double* parameters) const;

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const { return new LithiumTemplateFitFunction(*this); }

  double PdfValue(int xbin) const;

  TH1* DataHistogram() const;
  TH1* CalculateHistogram(std::string species) const;
  TH1* CalculateTotalHistogram() const;
  TH1* DrawZeroProbabilityEvents() const;

  void UpdateParameters(const double*);
  void UpdateErrors(const double*);

  int NumberOfEvents() const;

  /** FitMode returns the value of fitRatios. TO BE CHANGED! */
  std::string FitMode() const;

  Utilities::Quantity GetLi6Fraction() const;
  Utilities::Quantity GetLi7Fraction() const;
  Utilities::Quantity GetLi6Li7Ratio() const;

private:
  double GetLi6Events() const;
  double GetLi7Events() const;

  const TH1* fDataDistribution { nullptr };

  const TH1* fLi6Template { nullptr };
  const TH1* fLi7Template { nullptr };

  double fNumberOfEvents { 0.0 };

  double fLi6Fraction { 0.0 };
  double fLi7Fraction { 0.0 };
  double fLi6FractionError { 0.0 };
  double fLi7FractionError { 0.0 };

  double fLi6Li7Ratio { 0.0 };
  double fLi6Li7RatioError { 0.0 };

  double fRangeMinimum {0.0};
  double fRangeMaximum {0.0};
  int fBinMinimum {-1};
  int fBinMaximum {-1};
  bool fFitRatios { true };
};

#endif
