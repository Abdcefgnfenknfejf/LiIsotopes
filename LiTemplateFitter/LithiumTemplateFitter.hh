#ifndef LithiumTemplateFitter_hh
#define LithiumTemplateFitter_hh

#include <string>
#include <vector>

#include "Quantity.hh"

#include "LithiumTemplateFitFunction.hh"

#include "TCanvas.h"
#include <Math/Minimizer.h>

class TH1;

class LithiumTemplateFitter {
public:
  LithiumTemplateFitter(FitType fitType);

  struct FitResults {
    std::string fitType;                 // if true the "main" result vector is relativeFractions, otherwise ratios.
    int numberEvents = -1.0;
    std::vector<Utilities::Quantity> relativeFractions;        // Li6fraction/TotalNumber, Be9fraction/TotalNumber, Be10fraction/TotalNumber.
    std::vector<Utilities::Quantity> absoluteNumberEvents;     // Li6fraction*TotalNumber, Be9fraction*TotalNumber, Be10fraction*TotalNumber.
    std::vector<Utilities::Quantity> ratios;                   // Li6/Li7, Li7/total number events

    TH1* dataHisto = nullptr;
    std::vector<TH1*> templateHistos;
    TH1* resultHisto = nullptr;
    int ndf = -1;
    double chiSquare = -3.0;
    int detectorIndex = -1;
  };

  struct FitNormParameters {
    double fractionLi6 = -1.0;
    double fractionLi7 = -1.0;
  };

  struct FitRatioParameters {
    double ratioLi6Li7 = -2.0;
    double fractionLi7 = -2.0;
  };


  void SetTemplateHistogram(TH1* templateHistogram);
  void SetDataHistogram(TH1* dataHistogram);
  void SetFitRange(double fitRangeMinimum, double fitRangeMaximum);
  bool PerformTemplateFit(const FitNormParameters& startParameters, FitResults& fitResult);
  bool PerformTemplateFit(const FitRatioParameters& startParameters, FitResults& fitResult);
  bool PerformTemplateFitterFit(FitResults& fitResult); // use ACSoft template fitter
  void FillResults(ROOT::Math::Minimizer* minimizer, LithiumTemplateFitFunction& fitFunction, FitResults& fitResult);
  TCanvas* CreateResultDrawing(FitResults& fitResult, std::string canvasName);
  int TemplateColor(int i) const;


private:
  //  bool fFitRatios;
  FitType fFitType;
  double fFitRangeMinimum = 0.0;
  double fFitRangeMaximum = 0.0;
  int fFitBinMinimum = -1.0;
  int fFitBinMaximum = -1.0;
  TH1* fDataHisto;
  std::vector<TH1*> fTemplateHistos;
  const std::vector<std::string> fTemplateNames = {"Li6", "Li7"};
  int fNDF = -1; // effective number of degrees of freedom for chi2 calculation (# bins - number parameters)

  bool CheckHistogramConsistency();
  void ResetFitResults(FitResults& fitResults);
  ROOT::Math::Minimizer* SetupMinimizer(LithiumTemplateFitFunction& fitFunction, FitNormParameters& fitParameters, double toleranceFactor, double stepSizeFactor);
  ROOT::Math::Minimizer* SetupMinimizer(LithiumTemplateFitFunction& fitFunction, FitRatioParameters& fitParameters, double toleranceFactor, double stepSizeFactor);
  double CalculateChi2(TH1* dataHisto, TH1* modelHisto, int freeParameters);
  bool IsPhysicalResult(FitResults& fitResults);

};

#endif
