#ifndef TemplateFitTools_hh
#define TemplateFitTools_hh

#include "Quantity.hh"
#include "TemplateFitter.hh"

class TCanvas;
class TH1D;

class TemplateFitTools {
public:
  TemplateFitTools();

  bool PerformTemplateFit(TH1D* dataDistribution, TH1D* signalTemplate, TH1D* backgroundTemplate);
  TCanvas* DrawResults(const std::string& canvasName, const std::string& canvasTitle);

  double ChiSquarePerNDF() const { return fChiSquare / double(fNDF); }
  const Utilities::Quantity& SignalResults() const { return fSignalResults; }
  const Utilities::Quantity& BackgroundResults() const { return fBackgroundResults; }

private:
  TH1D* fDataDistribution;
  TH1D* fSignalTemplate;
  TH1D* fBackgroundTemplate;

  Utilities::TemplateFitter fTemplateFitter;
  Utilities::Quantity fSignalResults;
  Utilities::Quantity fBackgroundResults;
  double fChiSquare;
  int fNDF;
};

#endif
