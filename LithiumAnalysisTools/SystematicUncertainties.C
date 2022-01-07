#include "SystematicUncertainties.hh"

#include <cassert>
#include <cmath>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TStyle.h>

#include "AnalysisSettings.hh"
#include "FluxTools.hh"
#include "GraphTools.hh"
#include "PositronRatioTools.hh"
#include "ResolutionModels.hh"
#include "SmoothingTools.hh"
#include "Utilities.hh"

#define INFO_OUT_TAG "SystematicUncertainties"
#include "debugging.hh"

void AddAdditionalRelativeSystematicUncertaintyToGraph(TGraphAsymmErrors* graph, TGraph* additionalSystematicUncertainty) {

  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);

    double xSystematic, ySystematic;
    additionalSystematicUncertainty->GetPoint(i, xSystematic, ySystematic);

    assert(std::abs(x - xSystematic) < 1e-5);
    assert(std::abs(graph->GetErrorYlow(i) - graph->GetErrorYhigh(i)) < 1e-5);
    assert(!std::isinf(ySystematic));
    assert(!std::isnan(ySystematic));

    double eySystematic = ySystematic * y;
    double ey = std::sqrt(std::pow(graph->GetErrorY(i), 2) + std::pow(eySystematic, 2));
    assert(!std::isinf(ey));
    assert(!std::isnan(ey));

    graph->SetPointError(i, 0, 0, ey, ey);
  }
}

void DumpTableEntry(unsigned int bin, double x, double y, double ey) {

  INFO_OUT << "bin "  <<  std::setw(2) << std::right << bin          << std::fixed << std::setprecision(2)
           << " ("    <<  std::setw(7) << std::right << x            << " GeV) "
           << "y="    << std::setw(16) << std::right << y            << " "
           << "+/-"   << std::setw(14) << std::right << ey           << " "
           << "---> " <<  std::setw(8) << std::right << ey / y * 100 << " %" << std::endl;
}

void DumpTableEntryInvalid(unsigned int bin, double x, double y, double ey, const std::string& problem) {

  WARN_OUT << "bin "  <<  std::setw(2) << std::right << bin          << std::fixed << std::setprecision(2)
           << " ("    <<  std::setw(7) << std::right << x            << " GeV) "
           << "y="    << std::setw(16) << std::right << y            << " "
           << "+/-"   << std::setw(14) << std::right << ey           << " "
           << "---> " << problem       << std::endl;
}

std::string SpeciesTitleForType(SpeciesType type) {

  switch (type) {
  case ElectronFlux:
    return "Electron flux";
  case PositronFlux:
    return "Positron flux";
  case AllElectronFlux:
    return "All-electron flux";
  case PositronFraction:
    return "Positron fraction";
   case PositronElectronRatio:
    return "Positron/electron ratio";
  }

  assert(false);
  return "";
}

std::string SpeciesNameForType(SpeciesType type) {

  switch (type) {
  case ElectronFlux:
    return "ElectronFlux";
  case PositronFlux:
    return "PositronFlux";
  case AllElectronFlux:
    return "AllElectronFlux";
  case PositronFraction:
    return "PositronFraction";
   case PositronElectronRatio:
    return "PositronElectronRatio";
  }

  assert(false);
  return "";
}

template<typename T>
TH1D* HistogramFromGraph(T* graph, const Binning::Definition& binning) {

  TH1D* histogram = Make<TH1D>(Form("%sHistogram", graph->GetName()), graph->GetTitle(), binning);
  histogram->SetLineColor(graph->GetMarkerColor());
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    histogram->SetBinContent(i + 1, y ? (graph->GetY()[i]) : -1.0);
  }
  return histogram;
}

TH1D* DrawUncertainty(SpeciesType speciesType, const Binning::Definition& energyBinning, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper, double xLow, double xHigh) {

  TGraph* uncertaintyGraphClone = dynamic_cast<TGraph*>(uncertaintyGraph->Clone());

  if (speciesType == ElectronFlux)
    RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
  else if (speciesType == PositronFlux || speciesType == PositronFraction || speciesType == PositronElectronRatio)
    RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
  else {
    assert(speciesType == AllElectronFlux);
    RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
  }

  if (!same) {
    auto* frame = DrawEnergyFrame(xLow, xHigh, yLower, yUpper, "Energy / GeV", "Relative uncertainty");
    if (xLow == gAMSTimeDependentStartShowEnergy && xHigh == gAMSTimeDependentStopShowEnergy)
      frame->GetXaxis()->SetLabelSize(frame->GetXaxis()->GetLabelSize() * 0.9);
  }

  TH1D* histogram = HistogramFromGraph(uncertaintyGraphClone, energyBinning);
  histogram->GetXaxis()->SetRangeUser(xLow, xHigh);
  histogram->SetLineWidth(4);
  histogram->Draw("][.SAME");
  return histogram;
}

TH1D* DrawUncertainty(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper) {

  const Binning::Definition& energyBinning = Binning::Predefined::AbsoluteEnergyBinning();
  return DrawUncertainty(speciesType, energyBinning, uncertaintyGraph, same, yLower, yUpper, gAMSStartShowEnergy, gAMSStopShowEnergy);
}

TH1D* DrawUncertaintyInBartelsRotation(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper) {

  Binning::Definition energyBinning = Binning::Tools::SubRange(Binning::Predefined::AbsoluteEnergyBinning(),
                                                               gAMSTimeDependentStartShowEnergy, gAMSTimeDependentStopShowEnergy, true);

  return DrawUncertainty(speciesType, energyBinning, uncertaintyGraph, same, yLower, yUpper, gAMSTimeDependentStartShowEnergy, gAMSTimeDependentStopShowEnergy);
}

void DrawSingleUncertainty(TGraph* uncertainty, SpeciesType speciesType, const std::string& systName, const std::string& systTitle, std::string fileSuffix, double yLower, double yUpper) {

  std::string speciesTitle = SpeciesTitleForType(speciesType);
  std::string speciesName = SpeciesNameForType(speciesType);

  TCanvas* canvas = new TCanvas(Form("canvasUncertainty_%s_%s", speciesName.c_str(), systName.c_str()),
                                Form("%s) Single uncertainty '%s'", speciesTitle.c_str(), systTitle.c_str()));
  canvas->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D* uncertaintyHistogram = nullptr;
  if (!fileSuffix.empty())
    uncertaintyHistogram = DrawUncertaintyInBartelsRotation(speciesType, uncertainty, false, yLower, yUpper);
  else
    uncertaintyHistogram = DrawUncertainty(speciesType, uncertainty, false, yLower, yUpper);

  TLegend* legend = new TLegend(gPad->GetLeftMargin() + 0.03, gPad->GetBottomMargin() + 0.05, 1.0 - gPad->GetRightMargin() - 0.03, gPad->GetBottomMargin() + 0.15, nullptr, "brNDC");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(uncertaintyHistogram, systTitle.c_str(), "l");

  legend->Draw();
  if (!fileSuffix.empty()) {
    canvas->SaveAs(Form("%s_%s.png", canvas->GetName(), fileSuffix.c_str()));
    canvas->SaveAs(Form("%s_%s.pdf", canvas->GetName(), fileSuffix.c_str()));
  } else {
    canvas->SaveAs(Form("%s.png", canvas->GetName()));
    canvas->SaveAs(Form("%s.pdf", canvas->GetName()));
  }
}

void DrawTwoUncertainties(TGraph* uncertainty1, TGraph* uncertainty2, SpeciesType speciesType1, SpeciesType speciesType2, const std::string& systName, const std::string& systTitle, std::string fileSuffix, double yLower, double yUpper) {

  std::string speciesTitle1 = SpeciesTitleForType(speciesType1);
  std::string speciesName1 = SpeciesNameForType(speciesType1);

  std::string speciesTitle2 = SpeciesTitleForType(speciesType2);
  std::string speciesName2 = SpeciesNameForType(speciesType2);

  TCanvas* canvas = new TCanvas(Form("canvasUncertainties_%s_and_%s_%s", speciesName1.c_str(), speciesName2.c_str(), systName.c_str()),
                                Form("%s and %s) Single uncertainty '%s'", speciesTitle1.c_str(), speciesTitle2.c_str(), systTitle.c_str()));
  canvas->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D* uncertaintyHistogram1 = nullptr;
  TH1D* uncertaintyHistogram2 = nullptr;
  if (!fileSuffix.empty()) {
    uncertaintyHistogram1 = DrawUncertaintyInBartelsRotation(speciesType1, uncertainty1, false, yLower, yUpper);
    uncertaintyHistogram2 = DrawUncertaintyInBartelsRotation(speciesType2, uncertainty2, true, yLower, yUpper);
  } else {
    uncertaintyHistogram1 = DrawUncertainty(speciesType1, uncertainty1, false, yLower, yUpper);
    uncertaintyHistogram2 = DrawUncertainty(speciesType2, uncertainty2, true, yLower, yUpper);
  }
  uncertaintyHistogram2->SetLineStyle(2);
  uncertaintyHistogram2->SetLineWidth(2);

  TLegend* legend = new TLegend(gPad->GetLeftMargin() + 0.03, gPad->GetBottomMargin() + 0.05, 1.0 - gPad->GetRightMargin() - 0.03, gPad->GetBottomMargin() + 0.30, nullptr, "brNDC");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(uncertaintyHistogram1, Form("#splitline{%s}{for '%s'}", systTitle.c_str(), speciesTitle1.c_str()), "l");
  legend->AddEntry(uncertaintyHistogram2, Form("#splitline{%s}{for '%s'}", systTitle.c_str(), speciesTitle2.c_str()), "l");

  legend->Draw();
  if (!fileSuffix.empty()) {
    canvas->SaveAs(Form("%s_%s.png", canvas->GetName(), fileSuffix.c_str()));
    canvas->SaveAs(Form("%s_%s.pdf", canvas->GetName(), fileSuffix.c_str()));
  } else {
    canvas->SaveAs(Form("%s.png", canvas->GetName()));
    canvas->SaveAs(Form("%s.pdf", canvas->GetName()));
  }
}

TGraph* ExtractTotalRelativeUncertaintyFromGraph(TGraphAsymmErrors* graph) {

  TGraph* uncertaintyGraph = new TGraph(graph->GetN());
  for (int point = 0; point < graph->GetN(); ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    assert(std::abs(graph->GetErrorYlow(point) - graph->GetErrorYhigh(point)) < 1e-5);
    double ey = graph->GetErrorY(point);
    uncertaintyGraph->SetPoint(point, x, y ? ey / y : 0.0);
  }

  return uncertaintyGraph;
}

TGraph* CombineRelativeUncertainties(TGraphAsymmErrors* graph, const std::vector<TGraph*>& uncertainties) {

  assert(uncertainties.size() > 0);
  TGraph* uncertaintyGraph = new TGraph(uncertainties[0]->GetN());
  for (int point = 0; point < uncertainties[0]->GetN(); ++point) {
    double x = 0.0;
    double ySquareSum = 0.0;
    for (unsigned int i = 0; i < uncertainties.size(); ++i) {
      double xNew, y;
      uncertainties[i]->GetPoint(point, xNew, y);
      ySquareSum += std::pow(y, 2);

      if (i > 0)
        assert(std::abs(x - xNew) < 1e-5);
      else
        x = xNew;
    }

    // If the graph value is 0.0 (invalid graph point), set all uncertainties also to 0.0, by convention.
    double ey = std::sqrt(ySquareSum);
    if (graph->GetY()[point] == 0.0)
      ey = 0.0;

    uncertaintyGraph->SetPoint(point, x, ey);
  }

  return uncertaintyGraph;
}

void VerifyRelativeUncertaintyGraphIsNull(TGraph* graph, const std::string& description) {

  for (int i = 0; i < graph->GetN(); ++i) {
    double x1, y1;
    graph->GetPoint(i, x1, y1);

    // Warn about any inconsistencies.
    if (y1 != 0.0)
      FATAL_OUT << "bin=" << i << " x=" << x1 << " (GeV) -> y should be null, but is " << y1 << " in " << description << "." << std::endl;
  }
}

void VerifyRelativeUncertaintyGraphsAreEqual(TGraph* graph1, TGraph* graph2, const std::string& description, double tolerance) {

  assert(graph1->GetN() == graph2->GetN());
  for (int i = 0; i < graph1->GetN(); ++i) {
    double x1, y1;
    graph1->GetPoint(i, x1, y1);

    double x2, y2;
    graph2->GetPoint(i, x2, y2);
    assert(std::abs(x1 - x2) < 1e-5);

    // Warn about any inconsistencies.
    if (std::abs(y1 - y2) > tolerance)
      FATAL_OUT << "bin=" << i << " x=" << x1 << " (GeV) -> Deviation of " << std::abs(y1 - y2) * 1000.0 << " per-mille in " << description << "." << std::endl;
  }
}

double EnergyScaleUncertainty(double energy) {

  // Calibrated in test beam from 10GeV - 290GeV, uncertainty 2%
  static const double sSystError = 0.02;

  // Below 10 GeV, systematic uncertainties from threshold effects, energy lost in absorber material.
  if (energy < 10.0) {
    double x1 = std::log(0.5);
    double y1 = 0.05;             // systematic error at 0.5 GeV: 5%
    double x2 = std::log(10.0);
    double y2 = sSystError;
    double m = (y2 - y1) / (x2 - x1);
    double b = y2 - m * x2;
    return m * std::log(energy) + b;
  }

  // Above 290 GeV, systematic uncertainties from leakage und non-linearities in electronic.
  if (energy > 290.0) {
    double x1 = std::log(290.0);
    double y1 = sSystError;
    double x2 = std::log(700.0);
    double y2 = 0.04;             // systematic error at 700 GeV: 4%
    double m = (y2 - y1) / (x2 - x1);
    double b = y2 - m * x2;
    return m * std::log(energy) + b;
  }

  return sSystError;
}

double FunEnergyScaleUncertainty(double* x, double*) {

  return EnergyScaleUncertainty(x[0]);
}

TF1* EnergyScaleRelativeUncertaintyFunction() {

  const auto& energyBinning = Binning::Predefined::AbsoluteEnergyBinning();

  TF1* energyScaleErrorFunction = new TF1("energyScaleErrorFunction", FunEnergyScaleUncertainty, energyBinning.Min(), energyBinning.Max());
  energyScaleErrorFunction->SetNpx(10000);
  return energyScaleErrorFunction;
}

TGraph* EnergyScaleRelativeUncertaintyGraph() {

  TF1* energyScaleErrorFunction = EnergyScaleRelativeUncertaintyFunction();

  const auto& energyBinning = Binning::Predefined::AbsoluteEnergyBinning();
  TGraph* energyScaleRelativeUncertainty = new TGraph;
  for (unsigned int i = 0; i < energyBinning.NumberOfBins(); ++i) {
    double xTrue = Modelling::LaffertyWyatt(energyBinning.LowEdge(i + 1), energyBinning.UpEdge(i + 1), gAMSLaffertyWyattSpectralIndex);
    energyScaleRelativeUncertainty->SetPoint(i, xTrue, energyScaleErrorFunction->Eval(xTrue));
  }

  delete energyScaleErrorFunction;
  return energyScaleRelativeUncertainty;
}

TGraph* EnergyScaleRelativeUncertaintyGraphForFlux(TGraph* spectralIndexGraph, TGraphAsymmErrors* fluxGraph) {

  assert(spectralIndexGraph);
  assert(fluxGraph);
  assert(spectralIndexGraph->GetN() == fluxGraph->GetN());

  auto* energyScaleRelativeUncertaintyGraph = EnergyScaleRelativeUncertaintyGraph();

  TGraph* uncertaintyGraph = new TGraph(spectralIndexGraph->GetN());
  for (int i = 0; i < spectralIndexGraph->GetN(); ++i) {
    double x = spectralIndexGraph->GetX()[i];
    double xFlux = fluxGraph->GetX()[i];
    assert(std::abs(x - xFlux) < 1e-5);

    double sigmaEnergyOverEnergy = energyScaleRelativeUncertaintyGraph->Eval(x);
    assert(sigmaEnergyOverEnergy > 0);
    assert(!std::isnan(sigmaEnergyOverEnergy));

    double gamma = -spectralIndexGraph->GetY()[i]; // The graph contains -gamma.
    assert(!std::isnan(gamma));

    double ey = std::abs(std::pow(1.0 + sigmaEnergyOverEnergy, -gamma + 1.0) - 1.0);
    assert(ey > 0);
    assert(!std::isnan(ey));

    uncertaintyGraph->SetPoint(i, x, ey);
  }

  delete energyScaleRelativeUncertaintyGraph;
  return uncertaintyGraph;
}

TGraph* EnergyScaleRelativeUncertaintyGraphForFlux(double sigmaEnergyOverEnergy, TGraph* spectralIndexGraph, TGraphAsymmErrors* fluxGraph) {

  assert(sigmaEnergyOverEnergy > 0);
  assert(!std::isnan(sigmaEnergyOverEnergy));

  assert(spectralIndexGraph);
  assert(fluxGraph);
  assert(spectralIndexGraph->GetN() == fluxGraph->GetN());

  TGraph* uncertaintyGraph = new TGraph(spectralIndexGraph->GetN());
  for (int i = 0; i < spectralIndexGraph->GetN(); ++i) {
    double x = spectralIndexGraph->GetX()[i];
    double xFlux = fluxGraph->GetX()[i];
    assert(std::abs(x - xFlux) < 1e-5);

    double xSpectralIndex = x;
    if (xSpectralIndex > gAMSTimeDependentFluxesMaxEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMaxEnergy;
    else if (xSpectralIndex < gAMSTimeDependentFluxesMinEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMinEnergy;

    double gamma = -spectralIndexGraph->Eval(xSpectralIndex); // The graph contains -gamma.
    assert(!std::isnan(gamma));

    double ey = std::abs(std::pow(1.0 + sigmaEnergyOverEnergy, -gamma + 1.0) - 1.0);
    assert(ey > 0);
    assert(!std::isnan(ey));

    uncertaintyGraph->SetPoint(i, x, ey);
  }

  return uncertaintyGraph;
}

TGraph* EnergyScaleRelativeUncertaintyGraphForPositronFraction(TGraph* electronSpectralIndexGraph, TGraph* positronSpectralIndexGraph,
                                                               TGraphAsymmErrors* positronElectronRatioGraph, TGraphAsymmErrors* positronFractionGraph) {

  assert(electronSpectralIndexGraph);
  assert(positronSpectralIndexGraph);
  assert(positronElectronRatioGraph);
  assert(positronFractionGraph);

  assert(electronSpectralIndexGraph->GetN() == positronSpectralIndexGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronElectronRatioGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronFractionGraph->GetN());

  auto* energyScaleRelativeUncertaintyGraph = EnergyScaleRelativeUncertaintyGraph();

  TGraph* uncertaintyGraph = new TGraph(electronSpectralIndexGraph->GetN());
  for (int i = 0; i < electronSpectralIndexGraph->GetN(); ++i) {
    double xElectron = electronSpectralIndexGraph->GetX()[i];
    double xPositron = positronSpectralIndexGraph->GetX()[i];
    assert(std::abs(xElectron - xPositron) < 1e-5);

    double xPositronFraction = positronFractionGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronFraction) < 1e-5);

    double xPositronElectronRatio = positronElectronRatioGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronElectronRatio) < 1e-5);

    double sigmaEnergyOverEnergy = energyScaleRelativeUncertaintyGraph->Eval(xElectron);
    assert(sigmaEnergyOverEnergy > 0);
    assert(!std::isnan(sigmaEnergyOverEnergy));

    double gammaElectron = -electronSpectralIndexGraph->GetY()[i]; // The graph contains -gamma.
    assert(!std::isnan(gammaElectron));

    double gammaPositron = -positronSpectralIndexGraph->GetY()[i]; // The graph contains -gamma.
    assert(!std::isnan(gammaPositron));

    double y = positronFractionGraph->GetY()[i];

    double ey = 0.0;
    if (y > 0) {
      double ratio = positronElectronRatioGraph->GetY()[i];
      assert(ratio > 0);
      ey = std::abs(1.0 / (y * (1.0 + 1.0 / ratio * std::pow(1.0 + sigmaEnergyOverEnergy, gammaPositron - gammaElectron))) - 1.0);

      assert(ey > 0);
      assert(!std::isnan(ey));
    }

    uncertaintyGraph->SetPoint(i, xElectron, ey);
  }

  delete energyScaleRelativeUncertaintyGraph;
  return uncertaintyGraph;
}

TGraph* EnergyScaleRelativeUncertaintyGraphForPositronFraction(double sigmaEnergyOverEnergy, TGraph* electronSpectralIndexGraph, TGraph* positronSpectralIndexGraph,
                                                               TGraphAsymmErrors* positronElectronRatioGraph, TGraphAsymmErrors* positronFractionGraph) {

  assert(sigmaEnergyOverEnergy > 0);
  assert(!std::isnan(sigmaEnergyOverEnergy));

  assert(electronSpectralIndexGraph);
  assert(positronSpectralIndexGraph);
  assert(positronElectronRatioGraph);
  assert(positronFractionGraph);

  assert(electronSpectralIndexGraph->GetN() == positronSpectralIndexGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronElectronRatioGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronFractionGraph->GetN());

  TGraph* uncertaintyGraph = new TGraph(electronSpectralIndexGraph->GetN());
  for (int i = 0; i < electronSpectralIndexGraph->GetN(); ++i) {
    double xElectron = electronSpectralIndexGraph->GetX()[i];
    double xPositron = positronSpectralIndexGraph->GetX()[i];
    assert(std::abs(xElectron - xPositron) < 1e-5);

    double xPositronFraction = positronFractionGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronFraction) < 1e-5);

    double xPositronElectronRatio = positronElectronRatioGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronElectronRatio) < 1e-5);

    double xSpectralIndex = xElectron;
    if (xSpectralIndex > gAMSTimeDependentFluxesMaxEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMaxEnergy;
    else if (xSpectralIndex < gAMSTimeDependentFluxesMinEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMinEnergy;

    double gammaElectron = -electronSpectralIndexGraph->Eval(xSpectralIndex); // The graph contains -gamma.
    assert(!std::isnan(gammaElectron));

    double gammaPositron = -positronSpectralIndexGraph->Eval(xSpectralIndex); // The graph contains -gamma.
    assert(!std::isnan(gammaPositron));

    double y = positronFractionGraph->GetY()[i];

    double ey = 0.0;
    if (y > 0) {
      double ratio = positronElectronRatioGraph->GetY()[i];
      assert(ratio > 0);
      ey = std::abs(1.0 / (y * (1.0 + 1.0 / ratio * std::pow(1.0 + sigmaEnergyOverEnergy, gammaPositron - gammaElectron))) - 1.0);

      assert(ey > 0);
      assert(!std::isnan(ey));
    }

    uncertaintyGraph->SetPoint(i, xElectron, ey);
  }

  return uncertaintyGraph;
}


TGraph* EnergyScaleRelativeUncertaintyGraphForPositronElectronRatio(TGraph* electronSpectralIndexGraph, TGraph* positronSpectralIndexGraph, TGraphAsymmErrors* positronElectronRatioGraph) {

  assert(electronSpectralIndexGraph);
  assert(positronSpectralIndexGraph);
  assert(positronElectronRatioGraph);
  assert(electronSpectralIndexGraph->GetN() == positronSpectralIndexGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronElectronRatioGraph->GetN());

  auto* energyScaleRelativeUncertaintyGraph = EnergyScaleRelativeUncertaintyGraph();

  TGraph* uncertaintyGraph = new TGraph(electronSpectralIndexGraph->GetN());
  for (int i = 0; i < electronSpectralIndexGraph->GetN(); ++i) {
    double xElectron = electronSpectralIndexGraph->GetX()[i];
    double xPositron = positronSpectralIndexGraph->GetX()[i];
    assert(std::abs(xElectron - xPositron) < 1e-5);

    double xPositronElectronRatio = positronElectronRatioGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronElectronRatio) < 1e-5);

    double sigmaEnergyOverEnergy = energyScaleRelativeUncertaintyGraph->Eval(xElectron);
    assert(sigmaEnergyOverEnergy > 0);
    assert(!std::isnan(sigmaEnergyOverEnergy));

    double gammaElectron = -electronSpectralIndexGraph->GetY()[i]; // The graph contains -gamma.
    assert(!std::isnan(gammaElectron));

    double gammaPositron = -positronSpectralIndexGraph->GetY()[i]; // The graph contains -gamma.
    assert(!std::isnan(gammaPositron));

    double ey = std::abs(std::pow(1.0 + sigmaEnergyOverEnergy, gammaElectron - gammaPositron) - 1.0);
    assert(ey > 0);
    assert(!std::isnan(ey));

    uncertaintyGraph->SetPoint(i, xElectron, ey);
  }

  delete energyScaleRelativeUncertaintyGraph;
  return uncertaintyGraph;
}

TGraph* EnergyScaleRelativeUncertaintyGraphForPositronElectronRatio(double sigmaEnergyOverEnergy, TGraph* electronSpectralIndexGraph, TGraph* positronSpectralIndexGraph, TGraphAsymmErrors* positronElectronRatioGraph) {

  assert(sigmaEnergyOverEnergy > 0);
  assert(!std::isnan(sigmaEnergyOverEnergy));

  assert(electronSpectralIndexGraph);
  assert(positronSpectralIndexGraph);
  assert(positronElectronRatioGraph);
  assert(electronSpectralIndexGraph->GetN() == positronSpectralIndexGraph->GetN());
  assert(electronSpectralIndexGraph->GetN() == positronElectronRatioGraph->GetN());

  TGraph* uncertaintyGraph = new TGraph(electronSpectralIndexGraph->GetN());
  for (int i = 0; i < electronSpectralIndexGraph->GetN(); ++i) {
    double xElectron = electronSpectralIndexGraph->GetX()[i];
    double xPositron = positronSpectralIndexGraph->GetX()[i];
    assert(std::abs(xElectron - xPositron) < 1e-5);

    double xPositronElectronRatio = positronElectronRatioGraph->GetX()[i];
    assert(std::abs(xElectron - xPositronElectronRatio) < 1e-5);

    double xSpectralIndex = xElectron;
    if (xSpectralIndex > gAMSTimeDependentFluxesMaxEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMaxEnergy;
    else if (xSpectralIndex < gAMSTimeDependentFluxesMinEnergy)
      xSpectralIndex = gAMSTimeDependentFluxesMinEnergy;

    double gammaElectron = -electronSpectralIndexGraph->Eval(xSpectralIndex); // The graph contains -gamma.
    assert(!std::isnan(gammaElectron));

    double gammaPositron = -positronSpectralIndexGraph->Eval(xSpectralIndex); // The graph contains -gamma.
    assert(!std::isnan(gammaPositron));

    double ey = std::abs(std::pow(1.0 + sigmaEnergyOverEnergy, gammaElectron - gammaPositron) - 1.0);
    assert(ey > 0);
    assert(!std::isnan(ey));

    uncertaintyGraph->SetPoint(i, xElectron, ey);
  }

  return uncertaintyGraph;
}
