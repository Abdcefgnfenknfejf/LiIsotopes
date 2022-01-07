#include "GraphTools.hh"

#include "AnalysisSettings.hh"
#include "FluxTools.hh"

#include <TGraphAsymmErrors.h>

#include <cassert>
#include <cmath>
#include <iostream>

TGraphAsymmErrors* CorrectBinCentersLaffertyWyatt(TGraphAsymmErrors* spectrum, const Binning::Definition& binning, double gamma) {

  assert(spectrum->GetN() == int(binning.NumberOfBins()));
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(binning.NumberOfBins());
  graph->SetTitle(spectrum->GetTitle());

  int nPoint = 0;
  for (int i = 0; i < spectrum->GetN(); ++i) {
    bool binningIdenticalCaseOne = std::abs(spectrum->GetX()[i] - binning.LowEdge(i + 1));
    bool binningIdenticalCaseTwo = std::abs(spectrum->GetX()[i] - binning.Value(i + 1));
    assert(binningIdenticalCaseOne || binningIdenticalCaseTwo);
    double xTrue = Modelling::LaffertyWyatt(binning.LowEdge(i + 1), binning.LowEdge(i + 2), gamma);

    double x, y;
    spectrum->GetPoint(i, x, y);
    double exLow = spectrum->GetEXlow()[i];
    double exHigh = spectrum->GetEXhigh()[i];
    double eyLow = spectrum->GetEYlow()[i];
    double eyHigh = spectrum->GetEYhigh()[i];

    graph->SetPoint(nPoint, xTrue, y);
    graph->SetPointError(nPoint, exLow, exHigh, eyLow, eyHigh);
    ++nPoint;
  }

  return graph;
}

TGraphAsymmErrors* MultiplyByRigidityOrEnergy(TGraphAsymmErrors* spectrum, double power) {

  TGraphAsymmErrors* multiplied = new TGraphAsymmErrors(spectrum->GetN());

  for (int i = 0; i < spectrum->GetN(); ++i) {
    double x = 0;
    double y = 0;
    spectrum->GetPoint(i, x, y);

    double exLow = spectrum->GetErrorXlow(i);
    double exHigh = spectrum->GetErrorXhigh(i);
    double eyLow = spectrum->GetErrorYlow(i);
    double eyHigh = spectrum->GetErrorYhigh(i);

    multiplied->SetPoint(i, x, y * std::pow(x, power));
    multiplied->SetPointError(i, exLow, exHigh, eyLow * std::pow(x, power), eyHigh * std::pow(x, power));
  }

  return multiplied;
}

void RemoveHorizontalErrorsFromGraph(TGraphAsymmErrors* graph) {

  for (int i = 0; i < graph->GetN(); ++i)
    graph->SetPointError(i, 0, 0, graph->GetErrorYlow(i), graph->GetErrorYhigh(i));
}

void DivideFluxUncertaintyByBinWidth(TGraph* fluxUncertainty) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();
  for (int i = 0; i < fluxUncertainty->GetN(); ++i) {
    double xFlux, yFlux;
    fluxUncertainty->GetPoint(i, xFlux, yFlux);

    double binWidth = binning.LowEdge(i + 2) - binning.LowEdge(i + 1);
    fluxUncertainty->SetPoint(i, xFlux, yFlux / binWidth);
  }
}

void DivideFluxByBinWidth(TGraphAsymmErrors* flux) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();
  for (int i = 0; i < flux->GetN(); ++i) {
    double xFlux, yFlux;
    flux->GetPoint(i, xFlux, yFlux);
    double eyFluxLow = flux->GetErrorYlow(i);
    double eyFluxHigh = flux->GetErrorYhigh(i);

    if (yFlux > 0) {
      double binWidth = binning.LowEdge(i + 2) - binning.LowEdge(i + 1);
      flux->SetPoint(i, xFlux, yFlux / binWidth);
      flux->SetPointError(i, 0, 0, eyFluxLow / binWidth, eyFluxHigh / binWidth);
    }
  }
}
