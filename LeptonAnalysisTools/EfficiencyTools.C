#include "EfficiencyTools.hh"

#include "AnalysisSettings.hh"
#include "Cut.hh"
#include "DynamicSelector.hh"
#include "EfficiencyTemplateFits.hh"
#include "EfficiencyHistograms.hh"
#include "GraphTools.hh"
#include "Utilities.hh"

#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <cassert>
#include <cmath>
#include <iomanip>

#define INFO_OUT_TAG "EfficiencyTools"
#include "debugging.hh"

static const double sTagAndProbeSignificantCorrectionTresholdInSigma = 2.0;
static const double sTagAndProbeSignificantCorrectionTresholdAbsolute = 0.005;

TGraphErrors* EfficiencyComparisionData::CombineUncorrelatedEfficiencies(const std::vector<TGraphErrors*>& efficiencies) {

  if (efficiencies.empty())
    return nullptr;

  int points = efficiencies.front()->GetN();
  for (unsigned int k = 1; k < efficiencies.size(); ++k)
    assert(efficiencies[k]->GetN() == points);

  TGraphErrors* combinedEfficiency = new TGraphErrors(points);
  for (int i = 0; i < points; ++i) {
    double xCombined = 0.0;
    double yCombined = 0.0;
    efficiencies.front()->GetPoint(i, xCombined, yCombined);

    assert(yCombined >= 0);
    double eySquared = yCombined > 0 ? std::pow(efficiencies.front()->GetErrorY(i) / yCombined, 2) : 0.0;

    for (unsigned int k = 1; k < efficiencies.size(); ++k) {
      double x = 0.0;
      double y = 0.0;
      efficiencies[k]->GetPoint(i, x, y);
      assert(std::abs(x - xCombined) < 1e-5);

      yCombined *= y;
      assert(y >= 0);
      if (y > 0)
        eySquared += std::pow(efficiencies[k]->GetErrorY(i) / y, 2);
    }

    combinedEfficiency->SetPoint(i, xCombined, yCombined);
    combinedEfficiency->SetPointError(i, efficiencies.front()->GetErrorX(i), std::sqrt(eySquared) * yCombined);
  }

  return combinedEfficiency;
}

TGraphAsymmErrors* EfficiencyComparisionData::CombineUncorrelatedEfficiencies(const std::vector<TGraphAsymmErrors*>& efficiencies) {

  if (efficiencies.empty())
    return nullptr;

  int points = efficiencies.front()->GetN();
  for (unsigned int k = 1; k < efficiencies.size(); ++k)
    assert(efficiencies[k]->GetN() == points);

  TGraphAsymmErrors* combinedEfficiency = new TGraphAsymmErrors(points);
  for (int i = 0; i < points; ++i) {
    double xCombined = 0.0;
    double yCombined = 0.0;
    efficiencies.front()->GetPoint(i, xCombined, yCombined);

    assert(yCombined >= 0);
    double eyLowSquared = yCombined > 0 ? std::pow(efficiencies.front()->GetErrorYlow(i) / yCombined, 2) : 0.0;
    double eyHighSquared = yCombined > 0 ? std::pow(efficiencies.front()->GetErrorYhigh(i) / yCombined, 2) : 0.0;

    for (unsigned int k = 1; k < efficiencies.size(); ++k) {
      double x = 0.0;
      double y = 0.0;
      efficiencies[k]->GetPoint(i, x, y);
      assert(std::abs(x - xCombined) < 1e-5);

      yCombined *= y;
      assert(y >= 0);
      if (y > 0) {
        eyLowSquared += std::pow(efficiencies[k]->GetErrorYlow(i) / y, 2);
        eyHighSquared += std::pow(efficiencies[k]->GetErrorYhigh(i) / y, 2);
      }
    }

    combinedEfficiency->SetPoint(i, xCombined, yCombined);
    combinedEfficiency->SetPointEXlow(i, efficiencies.front()->GetErrorXlow(i));
    combinedEfficiency->SetPointEXhigh(i, efficiencies.front()->GetErrorXhigh(i));
    combinedEfficiency->SetPointEYlow(i, std::sqrt(eyLowSquared) * yCombined);
    combinedEfficiency->SetPointEYhigh(i, std::sqrt(eyHighSquared) * yCombined);
  }

  return combinedEfficiency;
}

static void DrawControlPlot(TH2* histogram, bool fillsPad = true) {

  assert(histogram);
  TString cutName = histogram->GetTitle();
  cutName.Replace(0, cutName.First(":") + 2, "");

  histogram->GetYaxis()->SetTitle(cutName.Data());
  histogram->GetYaxis()->SetTitleOffset(fillsPad ? 0.95 : 1.14);
  histogram->GetYaxis()->SetTitleSize(0.04);
  histogram->GetYaxis()->SetLabelSize(0.04);
  histogram->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
  histogram->GetXaxis()->SetTitleOffset(fillsPad ? 1.22 : 1.14);
  histogram->GetXaxis()->SetTitleSize(0.04);
  histogram->GetXaxis()->SetLabelSize(0.04);
  histogram->GetXaxis()->SetLabelOffset(fillsPad ? 0.004 : -0.005);
  Utilities::NormalizeHistogramXSlices(histogram);

  histogram->Draw("col.z");
  MoveTitleBoxToCenter();
}

static void StyleEfficiencyGraph(TGraphAsymmErrors* graph) {

  assert(graph);
  graph->GetYaxis()->SetTitle("Efficiency");
  graph->GetYaxis()->SetNoExponent(kTRUE);
  graph->GetYaxis()->SetRangeUser(-0.05, 1.09);
  graph->GetYaxis()->SetTitleOffset(0.6);
  graph->GetYaxis()->SetTitleSize(0.1);
  graph->GetYaxis()->SetLabelSize(0.1);
  graph->GetYaxis()->SetLabelOffset(0.01);

  graph->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
  graph->GetXaxis()->SetTitleOffset(1.22);
  graph->GetXaxis()->SetTitleSize(0.04);
  graph->GetXaxis()->SetLabelSize(0.04);
  graph->GetXaxis()->SetLabelOffset(0.004);
}

void SingleEfficiencyData::Initialize(TFile* effectiveAcceptanceFileElectron, TFile* effectiveAcceptanceFileProton, TFile* issTagAndProbeEfficiencyFile, TFile* mcTagAndProbeEfficiencyFile, const std::string& selectorName, unsigned int cutIndex) {

  const Binning::Definition* binning = &Binning::Predefined::AbsoluteEnergyBinning();
  if (selectorName == "EnergyScale") // Uses rigidity scale.
    binning = &Binning::Predefined::AbsoluteRigidityBinning();

  if (issTagAndProbeEfficiencyFile) {
    issControlPlotElectron = dynamic_cast<TH2F*>(issTagAndProbeEfficiencyFile->Get(Form("Electron/tagAndProbeEfficiencyControlPlot_for_cut_%i", cutIndex)));
    issControlPlotProton = dynamic_cast<TH2F*>(issTagAndProbeEfficiencyFile->Get(Form("Proton/tagAndProbeEfficiencyControlPlot_for_cut_%i", cutIndex)));
    issTagAndProbeElectron = dynamic_cast<TGraphAsymmErrors*>(issTagAndProbeEfficiencyFile->Get(Form("Electron/tagAndProbeEfficiency_for_cut_%i", cutIndex)));
    issTagAndProbeProton = dynamic_cast<TGraphAsymmErrors*>(issTagAndProbeEfficiencyFile->Get(Form("Proton/tagAndProbeEfficiency_for_cut_%i", cutIndex)));

    if (issTagAndProbeElectron) {
      issTagAndProbeElectron = CorrectBinCentersLaffertyWyatt(issTagAndProbeElectron, *binning, gAMSLaffertyWyattSpectralIndex);
      RemoveHorizontalErrorsFromGraph(issTagAndProbeElectron);
    }

    if (issTagAndProbeProton) {
      issTagAndProbeProton = CorrectBinCentersLaffertyWyatt(issTagAndProbeProton, *binning, gAMSLaffertyWyattSpectralIndex);
      RemoveHorizontalErrorsFromGraph(issTagAndProbeProton);
    }

    Cuts::Selector* issSelector = dynamic_cast<Cuts::Selector*>(issTagAndProbeEfficiencyFile->Get(Form("Selectors/%s", selectorName.c_str())));
    assert(issSelector);

    Cuts::Cut* issCut = issSelector->GetCut(cutIndex);
    assert(issCut);

    if (auto* efficiencyTemplateFitAttachment = const_cast<Cuts::EfficiencyTemplateFits*>(issCut->FindAttachment<Cuts::EfficiencyTemplateFits>())) {
      issTagAndProbeElectronPassedHistogram = efficiencyTemplateFitAttachment->PassedHistogramForCategory("Electron");
      issTagAndProbeElectronTotalHistogram = efficiencyTemplateFitAttachment->TotalHistogramForCategory("Electron");
    }
  }

  if (mcTagAndProbeEfficiencyFile) {
    mcControlPlotElectron = dynamic_cast<TH2F*>(mcTagAndProbeEfficiencyFile->Get(Form("Electron/tagAndProbeEfficiencyControlPlot_for_cut_%i", cutIndex)));
    mcTagAndProbeElectron = dynamic_cast<TGraphAsymmErrors*>(mcTagAndProbeEfficiencyFile->Get(Form("Electron/tagAndProbeEfficiency_for_cut_%i", cutIndex)));
    mcTagAndProbeProton = dynamic_cast<TGraphAsymmErrors*>(mcTagAndProbeEfficiencyFile->Get(Form("Proton/tagAndProbeEfficiency_for_cut_%i", cutIndex)));

    if (mcTagAndProbeElectron) {
      mcTagAndProbeElectron = CorrectBinCentersLaffertyWyatt(mcTagAndProbeElectron, *binning, gAMSLaffertyWyattSpectralIndex);
      RemoveHorizontalErrorsFromGraph(mcTagAndProbeElectron);
    }

    if (mcTagAndProbeProton) {
      mcTagAndProbeProton = CorrectBinCentersLaffertyWyatt(mcTagAndProbeProton, *binning, gAMSLaffertyWyattSpectralIndex);
      RemoveHorizontalErrorsFromGraph(mcTagAndProbeProton);
    }

    Cuts::Selector* mcSelector = dynamic_cast<Cuts::Selector*>(mcTagAndProbeEfficiencyFile->Get(Form("Selectors/%s", selectorName.c_str())));
    assert(mcSelector);

    Cuts::Cut* mcCut = mcSelector->GetCut(cutIndex);
    assert(mcCut);

    if (auto* efficiencyTemplateFitAttachment = const_cast<Cuts::EfficiencyTemplateFits*>(mcCut->FindAttachment<Cuts::EfficiencyTemplateFits>())) {
      mcTagAndProbeElectronPassedHistogram = efficiencyTemplateFitAttachment->PassedHistogramForCategory("Electron");
      mcTagAndProbeElectronTotalHistogram = efficiencyTemplateFitAttachment->TotalHistogramForCategory("Electron");
    }
  }

  if (effectiveAcceptanceFileElectron) {
    fLowerCutFunction = dynamic_cast<TF1*>(effectiveAcceptanceFileElectron->Get(Form("%s/lowerCutValueFunction_for_cut_%i", selectorName.c_str(), cutIndex)));
    fUpperCutFunction = dynamic_cast<TF1*>(effectiveAcceptanceFileElectron->Get(Form("%s/upperCutValueFunction_for_cut_%i", selectorName.c_str(), cutIndex)));
    mcCountingElectronGenerated = dynamic_cast<TGraphAsymmErrors*>(effectiveAcceptanceFileElectron->Get(Form("%s/mcCountingEfficiency_after_cut_%i", selectorName.c_str(), cutIndex)));
    RemoveHorizontalErrorsFromGraph(mcCountingElectronGenerated);
  }

  if (effectiveAcceptanceFileProton) {
    mcCountingProtonGenerated = dynamic_cast<TGraphAsymmErrors*>(effectiveAcceptanceFileProton->Get(Form("%s/mcCountingEfficiency_after_cut_%i", selectorName.c_str(), cutIndex)));
    RemoveHorizontalErrorsFromGraph(mcCountingProtonGenerated);
  }

  fSelectorName = selectorName;
  fCutIndex = cutIndex;
}

void SingleEfficiencyData::SetControlPlotYRange(float yMin, float yMax) {

  fRangeYMin = yMin;
  fRangeYMax = yMax;
}

std::string SingleEfficiencyData::CutName() const {

  assert(mcCountingElectronGenerated || issControlPlotElectron);
  TString cutName = mcCountingElectronGenerated ? mcCountingElectronGenerated->GetTitle() : issControlPlotElectron->GetTitle();
  cutName.Replace(0, cutName.First(":") + 2, "");
  return cutName.Data();
}

// Electrons
void SingleEfficiencyData::DrawControlPlotElectron() {

  assert(issControlPlotElectron);

  gPad->SetRightMargin(0.1);
  DrawControlPlot(issControlPlotElectron);

  if (fRangeYMin != fRangeYMax)
    issControlPlotElectron->GetYaxis()->SetRangeUser(fRangeYMin, fRangeYMax);

  if (fLowerCutFunction) {
    fLowerCutFunction->SetLineStyle(2);
    fLowerCutFunction->SetLineWidth(4);
    fLowerCutFunction->DrawCopy("l.same");
  }

  if (fUpperCutFunction) {
    fUpperCutFunction->SetLineStyle(2);
    fUpperCutFunction->SetLineWidth(4);
    fUpperCutFunction->DrawCopy("l.same");
  }
}

void SingleEfficiencyData::DrawMonteCarloCountingElectron() {

  assert(mcCountingElectronGenerated);
  mcCountingElectronGenerated->SetMarkerColor(kViolet);
  mcCountingElectronGenerated->SetMarkerStyle(24);
  mcCountingElectronGenerated->Draw("AP");
  StyleEfficiencyGraph(mcCountingElectronGenerated);
}

void SingleEfficiencyData::DrawTagAndProbeElectron(bool drawSame) {

  if (mcTagAndProbeElectron && mcTagAndProbeElectron->GetN() > 0) {
    mcTagAndProbeElectron->SetMarkerColor(kBlue);
    mcTagAndProbeElectron->SetMarkerStyle(kOpenCircle);
    if (drawSame)
      mcTagAndProbeElectron->Draw("P.SAME");
    else {
      drawSame = true;
      mcTagAndProbeElectron->Draw("AP");
      StyleEfficiencyGraph(mcTagAndProbeElectron);
    }
  }

  if (mcTagAndProbeElectronTemplateFit && mcTagAndProbeElectronTemplateFit->GetN() > 0) {
    mcTagAndProbeElectronTemplateFit->SetMarkerColor(kAzure);
    mcTagAndProbeElectronTemplateFit->SetMarkerStyle(kOpenSquare);
    if (drawSame)
      mcTagAndProbeElectronTemplateFit->Draw("P.SAME");
    else {
      drawSame = true;
      mcTagAndProbeElectronTemplateFit->Draw("AP");
      StyleEfficiencyGraph(mcTagAndProbeElectronTemplateFit);
    }
  }

  if (issTagAndProbeElectron && issTagAndProbeElectron->GetN() > 0) {
    issTagAndProbeElectron->SetMarkerColor(kRed);
    issTagAndProbeElectron->SetMarkerStyle(kFullCircle);
    if (drawSame)
      issTagAndProbeElectron->Draw("P.SAME");
    else {
      issTagAndProbeElectron->Draw("AP");
      StyleEfficiencyGraph(issTagAndProbeElectron);
    }
  }

  if (issTagAndProbeElectronTemplateFit && issTagAndProbeElectronTemplateFit->GetN() > 0) {
    issTagAndProbeElectronTemplateFit->SetMarkerSize(0.9);
    issTagAndProbeElectronTemplateFit->SetMarkerColor(kGreen + 2);
    issTagAndProbeElectronTemplateFit->SetMarkerStyle(kFullSquare);
    if (drawSame)
      issTagAndProbeElectronTemplateFit->Draw("P.SAME");
    else {
      issTagAndProbeElectronTemplateFit->Draw("AP");
      StyleEfficiencyGraph(issTagAndProbeElectronTemplateFit);
    }
  }
}

void SingleEfficiencyData::DrawLegendElectron() {

  TLegend* legendElectron = new TLegend(0.2088237, 0.1744418, 0.7992088, 0.5763274, nullptr, "brNDC");
  legendElectron->SetFillColor(kWhite);
  legendElectron->SetTextSize(0.042);
  if (mcCountingElectronGenerated && mcCountingElectronGenerated->GetN() > 0)
    legendElectron->AddEntry(mcCountingElectronGenerated, "MC counting (E_{gen})", "pl");
  if (mcTagAndProbeElectron && mcTagAndProbeElectron->GetN() > 0)
    legendElectron->AddEntry(mcTagAndProbeElectron, "MC tag & probe", "pl");
  if (mcTagAndProbeElectronTemplateFit && mcTagAndProbeElectronTemplateFit->GetN() > 0)
    legendElectron->AddEntry(mcTagAndProbeElectronTemplateFit, "MC tag & probe (temp. fit)", "pl");
  if (issTagAndProbeElectron && issTagAndProbeElectron->GetN() > 0)
    legendElectron->AddEntry(issTagAndProbeElectron, "ISS tag & probe", "pl");
  if (issTagAndProbeElectronTemplateFit && issTagAndProbeElectronTemplateFit->GetN() > 0)
    legendElectron->AddEntry(issTagAndProbeElectronTemplateFit, "ISS tag & probe (temp. fit)", "pl");
  legendElectron->Draw();
}

// Protons
void SingleEfficiencyData::DrawControlPlotProton() {

  assert(issControlPlotProton);

  gPad->SetRightMargin(0.1);
  DrawControlPlot(issControlPlotProton);

  if (fRangeYMin != fRangeYMax)
    issControlPlotProton->GetYaxis()->SetRangeUser(fRangeYMin, fRangeYMax);

  if (fLowerCutFunction) {
    fLowerCutFunction->SetLineStyle(2);
    fLowerCutFunction->SetLineWidth(4);
    fLowerCutFunction->DrawCopy("l.same");
  }

  if (fUpperCutFunction) {
    fUpperCutFunction->SetLineStyle(2);
    fUpperCutFunction->SetLineWidth(4);
    fUpperCutFunction->DrawCopy("l.same");
  }
}

void SingleEfficiencyData::DrawMonteCarloCountingProton() {

  assert(mcCountingProtonGenerated);
  mcCountingProtonGenerated->SetMarkerColor(kViolet);
  mcCountingProtonGenerated->SetMarkerStyle(24);
  mcCountingProtonGenerated->Draw("AP");
  StyleEfficiencyGraph(mcCountingProtonGenerated);
}

void SingleEfficiencyData::DrawTagAndProbeProton(bool drawSame) {

  if (mcTagAndProbeProton && mcTagAndProbeProton->GetN() > 0) {
    mcTagAndProbeProton->SetMarkerColor(kBlue);
    mcTagAndProbeProton->SetMarkerStyle(22);
    if (drawSame)
      mcTagAndProbeProton->Draw("P.SAME");
    else {
      drawSame = true;
      mcTagAndProbeProton->Draw("AP");
      StyleEfficiencyGraph(mcTagAndProbeProton);
    }
  }

  if (issTagAndProbeProton && issTagAndProbeProton->GetN() > 0) {
    issTagAndProbeProton->SetMarkerColor(kRed);
    issTagAndProbeProton->SetMarkerStyle(20);
    if (drawSame)
      issTagAndProbeProton->Draw("P.SAME");
    else {
      issTagAndProbeProton->Draw("AP");
      StyleEfficiencyGraph(issTagAndProbeProton);
    }
  }
}

void SingleEfficiencyData::DrawLegendProton() {

  TLegend* legendProton = new TLegend(0.2088237, 0.1744418, 0.7992088, 0.5763274, nullptr, "brNDC");
  legendProton->SetFillColor(kWhite);
  legendProton->SetTextSize(0.042);
  if (mcCountingProtonGenerated && mcCountingProtonGenerated->GetN() > 0)
    legendProton->AddEntry(mcCountingProtonGenerated, "MC counting (E_{gen})", "pl");
  if (mcTagAndProbeProton && mcTagAndProbeProton->GetN() > 0)
    legendProton->AddEntry(mcTagAndProbeProton, "MC tag & probe", "pl");
  if (issTagAndProbeProton && issTagAndProbeProton->GetN() > 0)
    legendProton->AddEntry(issTagAndProbeProton, "ISS tag & probe", "pl");
  legendProton->Draw();
}

void SingleEfficiencyData::DrawControlPlotPlaceholder(const std::string& species) {

  TPaveText* pt = new TPaveText(0.03982597,0.3251637,0.9413487,0.6271537,"br");
  pt->SetFillColor(kWhite);
  pt->SetTextColor(kRed);
  pt->SetTextSize(0.07);
  pt->SetTextAlign(22);
  pt->AddText("No ISS control plot");
  pt->AddText(Form("for %s", species.c_str()));
  pt->Draw();
}

void SingleEfficiencyData::DrawEfficiencyPlaceholder(const std::string& species) {

  TPaveText* pt = new TPaveText(0.03982597,0.3251637,0.9413487,0.6271537,"br");
  pt->SetFillColor(kWhite);
  pt->SetTextColor(kRed);
  pt->SetTextSize(0.07);
  pt->SetTextAlign(22);
  pt->AddText("Neither MC counting, nor");
  pt->AddText("tag & probe data available");
  pt->AddText(Form("for %s", species.c_str()));
  pt->Draw();
}

void SingleEfficiencyData::DrawOverview() {

  assert(mcCountingElectronGenerated);
  std::string cutName = CutName();

  bool hasProtonControlPlot = issControlPlotProton;
  bool hasProtonMc = mcCountingProtonGenerated ? mcCountingProtonGenerated->GetN() > 0 : false;
  bool hasProtonTagAndProbe = (mcTagAndProbeProton ? mcTagAndProbeProton->GetN() > 0 : false) || (issTagAndProbeProton ? issTagAndProbeProton->GetN() > 0 : false);

  bool hasElectronControlPlot = issControlPlotElectron;
  bool hasElectronMc = mcCountingElectronGenerated ? mcCountingElectronGenerated->GetN() > 0 : false;
  bool hasElectronTagAndProbe = (mcTagAndProbeElectron ? mcTagAndProbeElectron->GetN() > 0 : false) || (issTagAndProbeElectron ? issTagAndProbeElectron->GetN() > 0 : false);

  TCanvas* canvas = new TCanvas(Form("canvas_%s_%i", fSelectorName.c_str(), fCutIndex), Form("%s: Efficiency %i - %s", fSelectorName.c_str(), fCutIndex + 1, cutName.c_str()));
  canvas->Divide(2, 2, 1e-5, 1e-5);

  canvas->cd(1);
  if (hasElectronControlPlot) {
    gPad->SetLogx();
    gPad->SetLogz();
    DrawControlPlotElectron();
  } else
    DrawControlPlotPlaceholder("electrons");

  canvas->cd(2);
  if (hasProtonControlPlot) {
    gPad->SetLogx();
    gPad->SetLogz();
    DrawControlPlotProton();
  } else
    DrawControlPlotPlaceholder("protons");

  canvas->cd(3);
  if (hasElectronMc || hasElectronTagAndProbe) {
    gPad->SetGrid();
    gPad->SetLogx();

    if (hasElectronMc)
      DrawMonteCarloCountingElectron();
    if (hasElectronTagAndProbe)
      DrawTagAndProbeElectron(hasElectronMc);

    MoveTitleBoxToCenter();
    DrawLegendElectron();
  } else
    DrawEfficiencyPlaceholder("electrons");

  canvas->cd(4);
  if (hasProtonMc || hasProtonTagAndProbe) {
    gPad->SetGrid();
    gPad->SetLogx();

    if (hasProtonMc)
      DrawMonteCarloCountingProton();
    if (hasProtonTagAndProbe)
      DrawTagAndProbeProton(hasProtonMc);

    MoveTitleBoxToCenter();
    DrawLegendProton();
  } else
    DrawEfficiencyPlaceholder("protons");

  canvas->Modified();
  canvas->Update();

  canvas->SaveAs(Form("%s_%s.png", canvas->GetName(), cutName.c_str()));
  canvas->SaveAs(Form("%s_%s.pdf", canvas->GetName(), cutName.c_str()));
}

static double StraightLineInLogScaleParameterization(double* x, double* par) {

  double a0 = par[0];
  double a1 = par[1];
  return a0 + a1 * std::log10(x[0]);
}

double StraightLineWithBreakInLogScaleParameterization(double* x, double* par) {

  double a0 = par[0];
  double a1 = par[1];
  double b1 = par[2];
  double xBreak = par[3];
  double b0 = a0 + a1 * std::log10(xBreak) - b1 * std::log10(xBreak);
  if (x[0] < xBreak)
    return a0 + a1 * std::log10(x[0]);
  return b0 + b1 * std::log10(x[0]);
}

double StraightLineWithBreakInLogScaleEndingInAConstantParameterization(double* x, double* par) {

  double a0 = par[0];
  double a1 = par[1];
  double b1 = par[2];
  double xBreak1 = par[3];
  double xBreak2 = par[4];
  double b0 = a0 + a1 * std::log10(xBreak1) - b1 * std::log10(xBreak1);
  if (x[0] < xBreak1)
    return a0 + a1 * std::log10(x[0]);
  if (x[0] > xBreak2)
    return b0 + b1 * std::log10(xBreak2);

  return b0 + b1 * std::log10(x[0]);
}

static TF1* ConstantFunction() {

  static TF1* sFunction = nullptr;
  if (!sFunction) {
    sFunction = new TF1("ConstantParameterization", "pol0", 0.0, 2000.0);
    sFunction->SetParNames("p_{0}");
    sFunction->SetNpx(1e6);
  }
  return sFunction;
}

static TF1* StraightLineInLogScaleFunction() {

  static TF1* sFunction = nullptr;
  if (!sFunction) {
    sFunction = new TF1("StraightLineInLogScaleParameterization", StraightLineInLogScaleParameterization, 0.0, 2000.0, 2);
    sFunction->SetNpx(1e6);
    sFunction->SetParNames("p_{0}", "p_{1}");
    sFunction->SetParLimits(1, -0.5, 0.5);
  }
  return sFunction;
}

static TF1* StraightLineWithBreakInLogScaleFunction() {

  static TF1* sFunction = nullptr;
  if (!sFunction) {
    sFunction = new TF1("StraightLineWithBreakInLogScaleParameterization", StraightLineWithBreakInLogScaleParameterization, 0.0, 2000.0, 4);
    sFunction->SetNpx(1e6);
    sFunction->SetParameter(3, 30);
    sFunction->SetParLimits(1, -0.5, 0.5);
    sFunction->SetParLimits(2, -0.5, 0.5);
    sFunction->SetParNames("p_{0}", "p_{1}", "p_{2}", "p_{3}");
  }
  return sFunction;
}

static TF1* StraightLineWithBreakInLogScaleEndingInAConstantFunction() {

  static TF1* sFunction = nullptr;
  if (!sFunction) {
    sFunction = new TF1("StraightLineWithBreakInLogScaleEndingInAConstantParameterization", StraightLineWithBreakInLogScaleEndingInAConstantParameterization, 0.0, 2000.0, 5);
    sFunction->SetNpx(1e6);
    sFunction->SetParameter(3, 30);
    sFunction->SetParameter(4, 150);
    sFunction->SetParLimits(1, -0.5, 0.5);
    sFunction->SetParLimits(2, -0.5, 0.5);
    sFunction->SetParNames("p_{0}", "p_{1}", "p_{2}", "p_{3}", "p_{4}");
  }
  return sFunction;
}

double RelativeTagAndProbeRatioErrorAtEnergy(TGraphErrors* tagAndProbeRatioWithAllErrors, double energy) {

  // Only an approximation, we need this only to dump to the console.
  int tagAndProbePoint = 0;
  for (; tagAndProbePoint < tagAndProbeRatioWithAllErrors->GetN(); ++tagAndProbePoint) {
    if (tagAndProbeRatioWithAllErrors->GetX()[tagAndProbePoint] > energy) {
      if (tagAndProbePoint > 0)
        --tagAndProbePoint;
      break;
    }
  }

  return tagAndProbeRatioWithAllErrors->GetEY()[tagAndProbePoint] / tagAndProbeRatioWithAllErrors->GetY()[tagAndProbePoint];
}

TGraphErrors* SingleEfficiencyData::FitTagAndProbeRatio(TagAndProbeRatioFitMode mode, float xMin, float xMax, float yMin, float yMax, float break1Minimum, float break1Maximum, float fixedBreak1, float fixedBreak2, float keepRelativeErrorConstantAboveEnergy, TagAndProbeRatioFitTrustSource trust) {

  assert(mcCountingElectronGenerated);
  std::string cutName = CutName();

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  std::cout << std::endl;
  INFO_OUT << "Executing tag & probe ratio fit for cut='" << cutName << std::fixed << std::setprecision(2) << "' from xMin=" << xMin << " xMax=" << xMax << ":" << std::endl;
  TCanvas* canvas = new TCanvas(Form("canvas_tag_probe_%s_%i", fSelectorName.c_str(), fCutIndex), Form("%s: Tag & Probe ratio %i - %s", fSelectorName.c_str(), fCutIndex + 1, cutName.c_str()));
  canvas->Divide(1, 2, 1e-5, 1e-5);
  canvas->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(1e-5);
  gPad->SetPad(0, 0.33, 1, 1);
  gPad->SetLogx();

  TGraphAsymmErrors* mcTagAndProbeElectronClone = dynamic_cast<TGraphAsymmErrors*>(mcTagAndProbeElectron->Clone());
  mcTagAndProbeElectronClone->SetMarkerColor(kBlue);
  mcTagAndProbeElectronClone->SetMarkerStyle(kOpenCircle);
  mcTagAndProbeElectronClone->SetMarkerSize(2.2);
  StyleEfficiencyGraph(mcTagAndProbeElectronClone);

  TGraphAsymmErrors* mcTagAndProbeElectronTemplateFitClone = nullptr;
  if (mcTagAndProbeElectronTemplateFit) {
    mcTagAndProbeElectronTemplateFitClone = dynamic_cast<TGraphAsymmErrors*>(mcTagAndProbeElectronTemplateFit->Clone());
    mcTagAndProbeElectronTemplateFitClone->SetMarkerColor(kBlue);
    mcTagAndProbeElectronTemplateFitClone->SetMarkerStyle(kOpenSquare);
    mcTagAndProbeElectronTemplateFitClone->SetMarkerSize(2.2);
    StyleEfficiencyGraph(mcTagAndProbeElectronTemplateFitClone);
    mcTagAndProbeElectronTemplateFitClone->Draw("AP");
  } else
    mcTagAndProbeElectronClone->Draw("AP");

  TGraphAsymmErrors* issTagAndProbeElectronClone = dynamic_cast<TGraphAsymmErrors*>(issTagAndProbeElectron->Clone());
  issTagAndProbeElectronClone->SetMarkerColor(kRed);
  issTagAndProbeElectronClone->SetMarkerStyle(kFullCircle);
  issTagAndProbeElectronClone->SetMarkerSize(2.2);
  TruncateGraphEnergyRange(issTagAndProbeElectronClone, 0.0, xMax);

  TGraphAsymmErrors* issTagAndProbeElectronTemplateFitClone = nullptr;
  if (issTagAndProbeElectronTemplateFit) {
    issTagAndProbeElectronTemplateFitClone = dynamic_cast<TGraphAsymmErrors*>(issTagAndProbeElectronTemplateFit->Clone());
    issTagAndProbeElectronTemplateFitClone->SetMarkerColor(kRed);
    issTagAndProbeElectronTemplateFitClone->SetMarkerStyle(kFullSquare);
    issTagAndProbeElectronTemplateFitClone->SetMarkerSize(2.2);
    TruncateGraphEnergyRange(issTagAndProbeElectronTemplateFitClone, 0.0, xMax);
    issTagAndProbeElectronTemplateFitClone->Draw("P.SAME");
  } else
    issTagAndProbeElectronClone->Draw("P.SAME");

  TGraphAsymmErrors* mcTagAndProbeElectronCorrected = dynamic_cast<TGraphAsymmErrors*>(mcTagAndProbeElectron->Clone("mcTagAndProbeElectronCorrected"));
  mcTagAndProbeElectronCorrected->SetMarkerColor(kOrange);
  mcTagAndProbeElectronCorrected->SetMarkerStyle(kFullTriangleUp);
  mcTagAndProbeElectronCorrected->SetMarkerSize(2.2);
  MoveTitleBoxToCenter();

  canvas->cd(2);
  gPad->SetPad(0, 0, 1, 0.33);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.0);
  gPad->SetBottomMargin(0.4);
  gPad->SetLogx();

  if (mcTagAndProbeElectron->GetN() != issTagAndProbeElectron->GetN())
    FATAL_OUT << "mcTagAndProbeElectron('" << mcTagAndProbeElectron->GetName() << "')->GetN()=" << mcTagAndProbeElectron->GetN() << " "
              << "issTagAndProbeElectron('" << issTagAndProbeElectron->GetName() << "')->GetN()=" << issTagAndProbeElectron->GetN() << ". Aborting!" << std::endl;

  if (mcTagAndProbeElectronTemplateFit) {
    if (mcTagAndProbeElectron->GetN() != mcTagAndProbeElectronTemplateFit->GetN())
      FATAL_OUT << "mcTagAndProbeElectron('" << mcTagAndProbeElectron->GetName() << "')->GetN()=" << mcTagAndProbeElectron->GetN() << " "
                << "mcTagAndProbeElectronTemplateFit('" << mcTagAndProbeElectronTemplateFit->GetName() << "')->GetN()=" << mcTagAndProbeElectronTemplateFit->GetN() << ". Aborting!" << std::endl;
  }

  if (issTagAndProbeElectronTemplateFit) {
    if (issTagAndProbeElectron->GetN() != issTagAndProbeElectronTemplateFit->GetN())
      FATAL_OUT << "issTagAndProbeElectron('" << issTagAndProbeElectron->GetName() << "')->GetN()=" << issTagAndProbeElectron->GetN() << " "
                << "issTagAndProbeElectronTemplateFit('" << issTagAndProbeElectronTemplateFit->GetName() << "')->GetN()=" << issTagAndProbeElectronTemplateFit->GetN() << ". Aborting!" << std::endl;
  }

  TF1* fitFunction = nullptr;
  switch (mode) {
  case FitNothing:
    break;
  case FitConstant:
    fitFunction = ConstantFunction();
    break;
  case FitStraightLineInLogScale:
    fitFunction = StraightLineInLogScaleFunction();
    break;
  case FitStraightLineWithBreakInLogScale:
    fitFunction = StraightLineWithBreakInLogScaleFunction();
    if (fixedBreak1 != 0.0)
      fitFunction->FixParameter(3, fixedBreak1);
    if (break1Minimum != 0.0 && break1Maximum != 0.0)
      fitFunction->SetParLimits(3, break1Minimum, break1Maximum);
    break;
  case FitStraightLineWithBreakInLogScaleEndingInAConstant:
    fitFunction = StraightLineWithBreakInLogScaleEndingInAConstantFunction();
    if (fixedBreak1 != 0.0)
      fitFunction->FixParameter(3, fixedBreak1);
    if (fixedBreak2 != 0.0)
      fitFunction->FixParameter(4, fixedBreak2);
    if (break1Minimum != 0.0 && break1Maximum != 0.0)
      fitFunction->SetParLimits(3, break1Minimum, break1Maximum);
    break;
  }

  if (fitFunction) {
    fitFunction->SetLineColor(kMagenta + 2);
    fitFunction->SetLineWidth(0);
  }

  TGraphAsymmErrors* tagAndProbeRatio = new TGraphAsymmErrors(mcTagAndProbeElectron->GetN());

  unsigned int iteration = 0;
  float relativeSystematicUncertainty = 0.0001;
  do {
    for (int point = 0; point < mcTagAndProbeElectron->GetN(); ++point) {
      float x = issTagAndProbeElectron->GetX()[point];
      float y1 = issTagAndProbeElectronTemplateFit ? issTagAndProbeElectronTemplateFit->GetY()[point] : issTagAndProbeElectron->GetY()[point];
      float y2 = mcTagAndProbeElectronTemplateFit ? mcTagAndProbeElectronTemplateFit->GetY()[point] : mcTagAndProbeElectron->GetY()[point];

      float y1ErrorLow = issTagAndProbeElectronTemplateFit ? issTagAndProbeElectronTemplateFit->GetErrorYlow(point) : issTagAndProbeElectron->GetErrorYlow(point);
      float y2ErrorLow = mcTagAndProbeElectronTemplateFit ? mcTagAndProbeElectronTemplateFit->GetErrorYlow(point) : mcTagAndProbeElectron->GetErrorYlow(point);
      float y1ErrorHigh = issTagAndProbeElectronTemplateFit ? issTagAndProbeElectronTemplateFit->GetErrorYhigh(point) : issTagAndProbeElectron->GetErrorYhigh(point);
      float y2ErrorHigh = mcTagAndProbeElectronTemplateFit ? mcTagAndProbeElectronTemplateFit->GetErrorYhigh(point) : mcTagAndProbeElectron->GetErrorYhigh(point);

      float ratio = y2 ? y1 / y2 : 0.0;
      float ratioUncertaintyLow  = y2 ? std::sqrt(std::pow(y1ErrorLow / y2, 2)  + std::pow(-(y1 * y2ErrorLow)  / (y2 * y2), 2)) : 1.0;
      float ratioUncertaintyHigh = y2 ? std::sqrt(std::pow(y1ErrorHigh / y2, 2) + std::pow(-(y1 * y2ErrorHigh) / (y2 * y2), 2)) : 1.0;

      ratioUncertaintyLow = std::sqrt(std::pow(ratioUncertaintyLow, 2) + std::pow(relativeSystematicUncertainty * ratio, 2));
      ratioUncertaintyHigh = std::sqrt(std::pow(ratioUncertaintyHigh, 2) + std::pow(relativeSystematicUncertainty * ratio, 2));

      tagAndProbeRatio->SetPoint(point, x, ratio);
      tagAndProbeRatio->SetPointError(point, 0, 0, ratioUncertaintyLow, ratioUncertaintyHigh);
    }

    if (!fitFunction)
      break;

    ++iteration;
    tagAndProbeRatio->Fit(fitFunction, "EX0QRN", "", xMin, xMax);
    double chiSquarePerNdf = fitFunction->GetChisquare() / fitFunction->GetNDF();
    if (chiSquarePerNdf <= 1.02)
      break;

    if (iteration > 100)
      break;

    relativeSystematicUncertainty += 0.0001;
  } while (true);

  tagAndProbeRatio->SetMarkerColor(kGreen + 2);
  tagAndProbeRatio->SetMarkerStyle(kOpenCircle);
  tagAndProbeRatio->SetMarkerSize(1.3);
  tagAndProbeRatio->GetXaxis()->SetLabelOffset(0.008);
  tagAndProbeRatio->GetXaxis()->SetLabelSize(0.19);
  tagAndProbeRatio->GetXaxis()->SetTitleOffset(1.00);
  tagAndProbeRatio->GetXaxis()->SetTitleSize(0.19);
  tagAndProbeRatio->GetYaxis()->SetLabelOffset(0.01);
  tagAndProbeRatio->GetYaxis()->SetLabelSize(0.15);
  tagAndProbeRatio->GetYaxis()->SetTitleOffset(0.4);
  tagAndProbeRatio->GetYaxis()->SetTitleSize(0.15);
  tagAndProbeRatio->GetYaxis()->SetTitle("1 + #delta^{c}(E)");
  tagAndProbeRatio->GetYaxis()->SetNdivisions(205);
  tagAndProbeRatio->GetYaxis()->CenterTitle();
  tagAndProbeRatio->GetXaxis()->SetMoreLogLabels();
  tagAndProbeRatio->GetXaxis()->SetNoExponent(kTRUE);
  tagAndProbeRatio->GetYaxis()->SetRangeUser(yMin, yMax);
  tagAndProbeRatio->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
  if (fSelectorName == "EnergyScale") // Uses rigidity scale.
    tagAndProbeRatio->GetXaxis()->SetTitle("Rigidity / GV");
  else
    tagAndProbeRatio->GetXaxis()->SetTitle("ECAL Energy / GeV");

  tagAndProbeRatio->Draw("AP");

  std::cout << std::endl;
  std::cout << " Minuit output:" << std::endl;
  int fitStatus = 4000;
  if (fitFunction)
    fitStatus = tagAndProbeRatio->Fit(fitFunction, "EX0RM", "", xMin, xMax);
  std::cout << std::endl;

  TF1* fitFunctionClone = dynamic_cast<TF1*>(fitFunction->Clone());
  fitFunctionClone->SetLineWidth(4);

  double chiSquarePerNdf = fitFunction ? (fitFunction->GetChisquare() / fitFunction->GetNDF()) : 1.0;
  std::cout << " Done! Added " << relativeSystematicUncertainty * 100.0 << " % additional systematic uncertainty." << std::endl;
  std::cout << " Needed iterations=" << std::setw(3) << iteration << " Chi^2/ndf=" << std::setw(7) << chiSquarePerNdf << " status=" << fitStatus << std::endl;
  std::cout << std::endl;
  tagAndProbeRatio->SetTitle(Form("Additional MC systematic uncertainty: %.2f %%", relativeSystematicUncertainty * 100.0));

  if (iteration > 1000 || fitStatus != 4000)
    FATAL_OUT << " -> Tag & probe ratio fit did not converge for cut '" << cutName << "'. Aborting!" << std::endl;

  MoveTitleBoxToCenter();
  delete dynamic_cast<TPaveText*>(gPad->FindObject("title"));

  // Figure out whether the correction is significant.
  bool significantCorrection = true;
  switch (mode) {
  case FitNothing:
    significantCorrection = false;
    break;
  case FitConstant: {
    double offset = fitFunction->GetParameter(0);
    double offsetUncertainty = fitFunction->GetParError(0);
    if (std::abs(offset - 1.0) < sTagAndProbeSignificantCorrectionTresholdInSigma * offsetUncertainty
     || std::abs(offset - 1.0) < sTagAndProbeSignificantCorrectionTresholdAbsolute)
      significantCorrection = false;
    break;
  }
  case FitStraightLineInLogScale: {
    double offset = fitFunction->GetParameter(0);
    double offsetUncertainty = fitFunction->GetParError(0);
    double slope = fitFunction->GetParameter(1);
    double slopeUncertainty = fitFunction->GetParError(1);
    if (std::abs(offset - 1.0) < sTagAndProbeSignificantCorrectionTresholdInSigma * offsetUncertainty &&
        std::abs(slope) < sTagAndProbeSignificantCorrectionTresholdInSigma * slopeUncertainty)
      significantCorrection = false;

    // Additionaly check the minimum/maximum of the fit function, and it's deviation from unity.
    if (significantCorrection) {
      double fitFunctionMinimum = fitFunction->GetMinimum(xMin, xMax);
      double fitFunctionMaximum = fitFunction->GetMaximum(xMin, xMax);
      if (std::abs(fitFunctionMaximum - 1.0) < sTagAndProbeSignificantCorrectionTresholdAbsolute &&
          std::abs(fitFunctionMinimum - 1.0) < sTagAndProbeSignificantCorrectionTresholdAbsolute)
        significantCorrection = false;
    }

    break;
  }
  default:
    break;
  }

  if (fitFunction) {
    if (trust == TrustMeanDataMC) {
      std::cout << " ----> Significant correction? " << (significantCorrection ? "yes" : "no") << " (trust mode=MeanDataMC)" << std::endl
                << "      @   1 GeV: " << std::fixed << std::setprecision(2) << (1.0 - (fitFunction->Eval(1)   + 0.5 * (1.0 - fitFunction->Eval(1))))   * 100.0 << " %" << std::endl
                << "      @  10 GeV: " << std::fixed << std::setprecision(2) << (1.0 - (fitFunction->Eval(10)  + 0.5 * (1.0 - fitFunction->Eval(10))))  * 100.0 << " %" << std::endl
                << "      @ 100 GeV: " << std::fixed << std::setprecision(2) << (1.0 - (fitFunction->Eval(100) + 0.5 * (1.0 - fitFunction->Eval(100)))) * 100.0 << " %" << std::endl
                << "      @ 500 GeV: " << std::fixed << std::setprecision(2) << (1.0 - (fitFunction->Eval(500) + 0.5 * (1.0 - fitFunction->Eval(500)))) * 100.0 << " %" << std::endl;
    } else {
      std::cout << " ----> Significant correction? " << (significantCorrection ? "yes" : "no") << std::endl
                << "      @   1 GeV: " << std::fixed << std::setprecision(2) << (1.0 - fitFunction->Eval(1))   * 100.0 << " %" << std::endl
                << "      @  10 GeV: " << std::fixed << std::setprecision(2) << (1.0 - fitFunction->Eval(10))  * 100.0 << " %" << std::endl
                << "      @ 100 GeV: " << std::fixed << std::setprecision(2) << (1.0 - fitFunction->Eval(100)) * 100.0 << " %" << std::endl
                << "      @ 500 GeV: " << std::fixed << std::setprecision(2) << (1.0 - fitFunction->Eval(500)) * 100.0 << " %" << std::endl;
    }
  }

  gPad->Update();

  // Move statistics box to pad one.
  TPaveStats* statisticsBox = dynamic_cast<TPaveStats*>(tagAndProbeRatio->FindObject("stats"));
  TPaveStats* statisticsBoxClone = nullptr;
  if (statisticsBox) {
    statisticsBoxClone = dynamic_cast<TPaveStats*>(statisticsBox->Clone());
    statisticsBoxClone->SetX1NDC(0.50);
    statisticsBoxClone->SetY1NDC(0.036);
    statisticsBoxClone->SetX2NDC(0.936);
    statisticsBoxClone->SetY2NDC(0.58);
    statisticsBoxClone->SetTextFont(62);
    statisticsBoxClone->SetTextAlign(11);
    statisticsBoxClone->SetTextSize(0.07);
    statisticsBoxClone->SetBorderSize(0);

    TString chiSquareText = statisticsBoxClone->GetLine(0)->GetTitle();
    chiSquareText.ReplaceAll("ndf","dof");
    statisticsBoxClone->GetLine(0)->SetTitle(chiSquareText);
  }

  // Construct 1-sigma confidence intervals from fit results.
  TGraphErrors* tagAndProbeRatioWithAllErrors = new TGraphErrors(tagAndProbeRatio->GetN());
  for (int k = 0; k < tagAndProbeRatio->GetN(); ++k)
    tagAndProbeRatioWithAllErrors->SetPoint(k, tagAndProbeRatio->GetX()[k], 1.0);

  if (fitFunction)
    TVirtualFitter::GetFitter()->GetConfidenceIntervals(tagAndProbeRatioWithAllErrors, 0.683);

  // Respect the deviation of the ratio from unity as systematic error.
  for (int k = 0; k < tagAndProbeRatio->GetN(); ++k) {
    double x = tagAndProbeRatio->GetX()[k];

    if (!fitFunction) {
      tagAndProbeRatioWithAllErrors->SetPointError(k, 0, 0);
      continue;
    }

    double correctionValue = fitFunction->Eval(x);
    assert(std::abs(tagAndProbeRatioWithAllErrors->GetY()[k] - correctionValue) < 1e-5);
    if (trust == TrustMeanDataMC) {
      // Preserve relative uncertainty from fit when changing the correction from correctionValue -> 0.5 * (1 - correctionValue)
      double relativeUncertaintyFromFit = tagAndProbeRatioWithAllErrors->GetErrorY(k) / tagAndProbeRatioWithAllErrors->GetY()[k];
      correctionValue += 0.5 * (1.0 - correctionValue);
      tagAndProbeRatioWithAllErrors->SetPoint(k, x, correctionValue);
      tagAndProbeRatioWithAllErrors->SetPointError(k, 0, correctionValue * relativeUncertaintyFromFit);
    } else {
      // Preserve relative uncertainty from fit when changing the correction from correctionValue -> 0.5 * (1 - correctionValue)
      double relativeUncertaintyFromFit = tagAndProbeRatioWithAllErrors->GetErrorY(k) / tagAndProbeRatioWithAllErrors->GetY()[k];
      correctionValue += 0.5 * (1.0 - correctionValue);
      tagAndProbeRatioWithAllErrors->SetPointError(k, 0, (correctionValue + 0.5 * (1.0 - correctionValue)) * relativeUncertaintyFromFit);
    }

    double uncertaintyFromFit = tagAndProbeRatioWithAllErrors->GetErrorY(k);
    if (std::isnan(uncertaintyFromFit) || std::isinf(uncertaintyFromFit)) {
      if (std::isnan(uncertaintyFromFit))
        FATAL_OUT << "x=" << x << " error from tag & probe fit is NaN!" << std::endl;
      else
        FATAL_OUT << "x=" << x << " error from tag & probe fit is Inf!" << std::endl;
    }

    double absoluteUncertaintyFromModelFit = uncertaintyFromFit;
    double absoluteUncertaintyFromAdditionalSystematicsGoodnessOfFit = relativeSystematicUncertainty * tagAndProbeRatioWithAllErrors->GetY()[k];
    double absoluteUncertaintyFromDeviationToUnity = std::abs(correctionValue - 1.0);

    double absoluteUncertainty = std::sqrt(std::pow(absoluteUncertaintyFromModelFit, 2) + std::pow(absoluteUncertaintyFromAdditionalSystematicsGoodnessOfFit, 2) + std::pow(absoluteUncertaintyFromDeviationToUnity, 2));
    assert(!std::isnan(absoluteUncertainty));
    assert(!std::isinf(absoluteUncertainty));
    tagAndProbeRatioWithAllErrors->SetPointError(k, 0, absoluteUncertainty);
  }

  // For some cuts we might want to keep the relative error at a higher level then needed.
  if (keepRelativeErrorConstantAboveEnergy > 0.0) {
    double relativeUncertaintyReference = 0.0;
    for (int k = 0; k < tagAndProbeRatioWithAllErrors->GetN(); ++k) {
      double x, y;
      tagAndProbeRatioWithAllErrors->GetPoint(k, x, y);
      if (x < keepRelativeErrorConstantAboveEnergy)
        continue;

      if (relativeUncertaintyReference == 0.0) {
        relativeUncertaintyReference = tagAndProbeRatioWithAllErrors->GetErrorY(k) / y;
        continue;
      }

      tagAndProbeRatioWithAllErrors->SetPointError(k, 0.0, relativeUncertaintyReference * y);
    }
  }

  std::cout << std::endl;
  std::cout << " ----> Final error including deviation from unity as additional systematic uncertainty:" << std::endl
            << "      @   1 GeV: " << std::fixed << std::setprecision(2) << RelativeTagAndProbeRatioErrorAtEnergy(tagAndProbeRatioWithAllErrors,   1) * 100.0 << " %" << std::endl
            << "      @  10 GeV: " << std::fixed << std::setprecision(2) << RelativeTagAndProbeRatioErrorAtEnergy(tagAndProbeRatioWithAllErrors,  10) * 100.0 << " %" << std::endl
            << "      @ 100 GeV: " << std::fixed << std::setprecision(2) << RelativeTagAndProbeRatioErrorAtEnergy(tagAndProbeRatioWithAllErrors, 100) * 100.0 << " %" << std::endl
            << "      @ 500 GeV: " << std::fixed << std::setprecision(2) << RelativeTagAndProbeRatioErrorAtEnergy(tagAndProbeRatioWithAllErrors, 500) * 100.0 << " %" << std::endl;

  // Only truncate energy range, after the fit was performed!
  TruncateGraphEnergyRange(tagAndProbeRatio, 0.0, xMax);
  tagAndProbeRatioWithAllErrors->SetTitle(cutName.c_str());
  tagAndProbeRatioWithAllErrors->SetFillColor(kOrange);
  tagAndProbeRatioWithAllErrors->SetFillStyle(1001);
  tagAndProbeRatioWithAllErrors->SetLineColor(kBlack);
  tagAndProbeRatioWithAllErrors->SetLineWidth(4);
  tagAndProbeRatioWithAllErrors->SetLineStyle(7);

  tagAndProbeRatioWithAllErrors->Draw("E3.same");
  fitFunctionClone->Draw("L.same");
  tagAndProbeRatioWithAllErrors->Draw("LX.same");
  tagAndProbeRatio->Draw("P.same");

  // Create legend & statistics box on upper pad, remove statistics box from lower pad.
  canvas->cd(1);
  if (statisticsBoxClone)
    statisticsBoxClone->Draw();

  TLegend* legend = new TLegend(0.15, 0.036, 0.50, 0.58, nullptr, "brNDC");
  legend->SetTextSize(0.07);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  if (mcTagAndProbeElectronTemplateFitClone)
    legend->AddEntry(mcTagAndProbeElectronTemplateFitClone, "MC tag & probe", "p");
  else
    legend->AddEntry(mcTagAndProbeElectronClone, "MC tag & probe", "p");
  if (issTagAndProbeElectronTemplateFitClone)
    legend->AddEntry(issTagAndProbeElectronTemplateFitClone, "ISS tag & probe", "p");
  else
    legend->AddEntry(issTagAndProbeElectronClone, "ISS tag & probe", "p");
  legend->AddEntry(tagAndProbeRatio, "ISS / MC ratio", "p");
  legend->AddEntry(fitFunctionClone, "Ratio fit", "fl");
  legend->AddEntry(tagAndProbeRatioWithAllErrors, "Correction", "l");
  legend->Draw();

  gPad->Update();
  canvas->cd(2);
  gStyle->SetOptFit(0);
  if (statisticsBox)
    tagAndProbeRatio->GetListOfFunctions()->Remove(statisticsBox);
  gPad->Update();

  if (significantCorrection)
    tagAndProbeRatioWithAllErrors->SetName(Form("cut_%i_tag_and_probe_correction_significant", fCutIndex));
  else
    tagAndProbeRatioWithAllErrors->SetName(Form("cut_%i_tag_and_probe_correction", fCutIndex));

#if 0
  // Eventually draw MC tag&probe corrected result to upper pad, for debugging purpose.
   canvas->cd(1);

  // Replace errors with the ones from the tag&probe fit result.
  for (int i = 0; i < mcTagAndProbeElectronCorrected->GetN(); ++i)
    mcTagAndProbeElectronCorrected->SetPointError(i, 0, 0, 0, 0);
  EfficiencyComparisionData::CorrectEfficiencyForTagAndProbeRatio(mcTagAndProbeElectronCorrected, tagAndProbeRatioWithAllErrors);
  mcTagAndProbeElectronCorrected->Draw("P.SAME");
#endif

  canvas->SaveAs(Form("%s_%s.png", canvas->GetName(), cutName.c_str()));
  canvas->SaveAs(Form("%s_%s.pdf", canvas->GetName(), cutName.c_str()));
  return tagAndProbeRatioWithAllErrors;
}

TGraphErrors* SingleEfficiencyData::FitBartelsRotationToAverageDataRatio(const SingleEfficiencyData& averageData, const SingleEfficiencyData& bartelsRotationData, float yMin, float yMax, const std::string& bartelsRotationOutputFileSuffix) {

  assert(averageData.issTagAndProbeElectron);
  assert(averageData.issControlPlotElectron);

  TGraphAsymmErrors* issTagAndProbeElectronAverage = averageData.issTagAndProbeElectron;
  TGraphAsymmErrors* issTagAndProbeElectronBartelsRotation = bartelsRotationData.issTagAndProbeElectron;

  if (!issTagAndProbeElectronBartelsRotation) {
    WARN_OUT << "No time dependent tag & probe information available. Skipping!" << std::endl;
    issTagAndProbeElectronBartelsRotation = dynamic_cast<TGraphAsymmErrors*>(issTagAndProbeElectronAverage->Clone());
  }

  std::string cutName = averageData.CutName();

  // If the BR is known to have no/small events, assume no correction for this BR (tagAndProbeRatio 100% uncertainty).
  if (IsBlacklistedBartelsRotation(bartelsRotationOutputFileSuffix)) {
    WARN_OUT << "Bartels rotation '" << bartelsRotationOutputFileSuffix << "' is black-listed. Assuming tag & probe ratio for cut '" << cutName << "' is at 1.0 +/- 1.0." << std::endl;

    TGraphErrors* tagAndProbeRatio = new TGraphErrors(issTagAndProbeElectronAverage->GetN());
    for (int point = 0; point < issTagAndProbeElectronAverage->GetN(); ++point) {
      float x = issTagAndProbeElectronAverage->GetX()[point];
      tagAndProbeRatio->SetPoint(point, x, 1.0);
      tagAndProbeRatio->SetPointError(point, 0.0, 1.0);
    }

    tagAndProbeRatio->SetTitle(cutName.c_str());
    tagAndProbeRatio->SetName(Form("cut_%i_bartels_rotation_average_tag_and_probe_correction", averageData.fCutIndex));
    return tagAndProbeRatio;
  }

  // For some cuts the ratio is exactly 1.0 everywhere in the energy range used for the time-dep fluxes.
  // IN that case skip the fit, and continue assuming 1.0 with 0% uncertainty.
  bool ratioIsOneEverywhere = true;
  for (int point = 0; point < issTagAndProbeElectronAverage->GetN(); ++point) {
    float x = issTagAndProbeElectronAverage->GetX()[point];
    if (x < gAMSTimeDependentFluxesMinEnergy || x > gAMSTimeDependentFluxesMaxEnergy)
      continue;

    float y1 = issTagAndProbeElectronBartelsRotation->GetY()[point];
    float y2 = issTagAndProbeElectronAverage->GetY()[point];
    float ratio = y2 ? y1 / y2 : 0.0;
    if (std::abs(ratio - 1.0) > 1e-10) {
      ratioIsOneEverywhere = false;
      break;
    }
  }

  if (ratioIsOneEverywhere) {
    INFO_OUT << "Ratio is 1.0 everywhere, skipping fit for cut '" << cutName << "'." << std::endl;
    TGraphErrors* tagAndProbeRatio = new TGraphErrors(issTagAndProbeElectronAverage->GetN());
    for (int point = 0; point < issTagAndProbeElectronAverage->GetN(); ++point) {
      float x = issTagAndProbeElectronAverage->GetX()[point];
      tagAndProbeRatio->SetPoint(point, x, 1.0);
      tagAndProbeRatio->SetPointError(point, 0.0, 0.0);
    }

    tagAndProbeRatio->SetTitle(cutName.c_str());
    tagAndProbeRatio->SetName(Form("cut_%i_bartels_rotation_average_tag_and_probe_correction", averageData.fCutIndex));
    return tagAndProbeRatio;
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas* canvas = new TCanvas(Form("canvas_tag_probe_%s_%i_bartels_rotation_to_average", averageData.fSelectorName.c_str(), averageData.fCutIndex), Form("%s: Tag & Probe ratio %i - %s Bartels rotation/Average", averageData.fSelectorName.c_str(), averageData.fCutIndex + 1, cutName.c_str()));
  canvas->Divide(1, 2, 1e-5, 1e-5);
  canvas->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(1e-5);
  gPad->SetPad(0, 0.33, 1, 1);
  gPad->SetLogx();

  issTagAndProbeElectronAverage->SetMarkerColor(kBlue);
  issTagAndProbeElectronAverage->SetMarkerStyle(kOpenCircle);
  issTagAndProbeElectronAverage->SetMarkerSize(2.2);
  issTagAndProbeElectronAverage->Draw("AP");
  StyleEfficiencyGraph(issTagAndProbeElectronAverage);
  issTagAndProbeElectronAverage->GetXaxis()->SetLimits(gAMSTimeDependentStartShowEnergy, gAMSTimeDependentStopShowEnergy);

  issTagAndProbeElectronBartelsRotation->SetMarkerColor(kRed);
  issTagAndProbeElectronBartelsRotation->SetMarkerSize(2.2);
  issTagAndProbeElectronBartelsRotation->SetMarkerStyle(kFullCircle);
  issTagAndProbeElectronBartelsRotation->Draw("P.SAME");

  MoveTitleBoxToCenter();

  canvas->cd(2);
  gPad->SetPad(0, 0, 1, 0.33);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.0);
  gPad->SetBottomMargin(0.4);
  gPad->SetLogx();

  if (issTagAndProbeElectronAverage->GetN() != issTagAndProbeElectronBartelsRotation->GetN())
    FATAL_OUT << "issTagAndProbeElectronAverage('" << issTagAndProbeElectronAverage->GetName() << "')->GetN()=" << issTagAndProbeElectronAverage->GetN() << " "
              << "issTagAndProbeElectronBartelsRotation('" << issTagAndProbeElectronBartelsRotation->GetName() << "')->GetN()=" << issTagAndProbeElectronBartelsRotation->GetN() << ". Aborting!" << std::endl;

  TGraphAsymmErrors* tagAndProbeRatio = new TGraphAsymmErrors(issTagAndProbeElectronAverage->GetN());

  for (int point = 0; point < issTagAndProbeElectronAverage->GetN(); ++point) {
    float x = issTagAndProbeElectronAverage->GetX()[point];

    float y1 = issTagAndProbeElectronBartelsRotation->GetY()[point];
    float y2 = issTagAndProbeElectronAverage->GetY()[point];

    float y1ErrorLow = issTagAndProbeElectronBartelsRotation->GetErrorYlow(point);
    float y2ErrorLow = issTagAndProbeElectronAverage->GetErrorYlow(point);
    float y1ErrorHigh = issTagAndProbeElectronBartelsRotation->GetErrorYhigh(point);
    float y2ErrorHigh = issTagAndProbeElectronAverage->GetErrorYhigh(point);

    float ratio = y2 ? y1 / y2 : 0.0;
    float ratioUncertaintyLow  = y2 ? std::sqrt(std::pow(y1ErrorLow / y2, 2)  + std::pow(-(y1 * y2ErrorLow)  / (y2 * y2), 2)) : 1.0;
    float ratioUncertaintyHigh = y2 ? std::sqrt(std::pow(y1ErrorHigh / y2, 2) + std::pow(-(y1 * y2ErrorHigh) / (y2 * y2), 2)) : 1.0;

    // Add 1e-4 relative error on top to stability fits numerically (for cases where eyh=0, eyl=<very-small>).
    ratioUncertaintyLow = std::sqrt(std::pow(ratioUncertaintyLow / ratio, 2) + std::pow(1.0e-4, 2)) * ratio;
    ratioUncertaintyHigh = std::sqrt(std::pow(ratioUncertaintyHigh / ratio, 2) + std::pow(1.0e-4, 2)) * ratio;

    tagAndProbeRatio->SetPoint(point, x, ratio);
    tagAndProbeRatio->SetPointError(point, 0, 0, ratioUncertaintyLow, ratioUncertaintyHigh);
  }

  tagAndProbeRatio->SetMarkerColor(kGreen + 2);
  tagAndProbeRatio->SetMarkerStyle(kOpenCircle);
  tagAndProbeRatio->SetMarkerSize(1.3);
  tagAndProbeRatio->GetXaxis()->SetLabelOffset(0.008);
  tagAndProbeRatio->GetXaxis()->SetLabelSize(0.19);
  tagAndProbeRatio->GetXaxis()->SetTitleOffset(1.00);
  tagAndProbeRatio->GetXaxis()->SetTitleSize(0.19);
  tagAndProbeRatio->GetYaxis()->SetLabelOffset(0.01);
  tagAndProbeRatio->GetYaxis()->SetLabelSize(0.15);
  tagAndProbeRatio->GetYaxis()->SetTitleOffset(0.4);
  tagAndProbeRatio->GetYaxis()->SetTitleSize(0.15);
  tagAndProbeRatio->GetYaxis()->SetTitle("Ratio");
  tagAndProbeRatio->GetYaxis()->SetNdivisions(205);
  tagAndProbeRatio->GetYaxis()->CenterTitle();
  tagAndProbeRatio->GetXaxis()->SetMoreLogLabels();
  tagAndProbeRatio->GetXaxis()->SetNoExponent(kTRUE);
  tagAndProbeRatio->GetYaxis()->SetRangeUser(yMin, yMax);
  tagAndProbeRatio->GetXaxis()->SetLimits(gAMSTimeDependentStartShowEnergy, gAMSTimeDependentStopShowEnergy);
  tagAndProbeRatio->GetXaxis()->SetTitle("ECAL Energy / GeV");
 
  tagAndProbeRatio->Draw("AP");
  MoveTitleBoxToCenter();
  delete dynamic_cast<TPaveText*>(gPad->FindObject("title"));

  std::cout << std::endl;
  std::cout << " Minuit output:" << std::endl;
  TF1* fitFunction = ConstantFunction();
  fitFunction->SetLineWidth(0);

  int fitStatus = tagAndProbeRatio->Fit(fitFunction, "EX0RM", "", gAMSTimeDependentFluxesMinEnergy, gAMSTimeDependentFluxesMaxEnergy);
  std::cout << std::endl;

  if (fitStatus != 4000)
    FATAL_OUT << " -> Fit did not converge for cut '" << cutName << "'. Aborting!" << std::endl;

  TF1* fitFunctionClone = dynamic_cast<TF1*>(fitFunction->Clone());
  fitFunctionClone->SetLineWidth(4);

  // Figure out whether the correction is significant.
  bool significantCorrection = true;
  if (std::abs(1.0 - fitFunction->Eval(10)) * 100.0 < 0.56)
    significantCorrection = false;

  std::cout << "Cut='" << cutName << "' significant correction? " << (significantCorrection ? "yes" : "no") << " magnitude: "
            << std::fixed << std::setprecision(2) << (1.0 - fitFunction->Eval(10)) * 100 << " %" << std::endl;
  gPad->Update();

  // Move statistics box to pad one.
  TPaveStats* statisticsBox = dynamic_cast<TPaveStats*>(tagAndProbeRatio->FindObject("stats"));
  TPaveStats* statisticsBoxClone = dynamic_cast<TPaveStats*>(statisticsBox->Clone());
  statisticsBoxClone->SetX1NDC(0.60);
  statisticsBoxClone->SetY1NDC(0.036);
  statisticsBoxClone->SetX2NDC(0.936);
  statisticsBoxClone->SetY2NDC(0.58);
  statisticsBoxClone->SetTextFont(62);
  statisticsBoxClone->SetTextAlign(11);
  statisticsBoxClone->SetTextSize(0.07);
  statisticsBoxClone->SetBorderSize(0);

  TString chiSquareText = statisticsBoxClone->GetLine(0)->GetTitle();
  chiSquareText.ReplaceAll("ndf","dof");
  statisticsBoxClone->GetLine(0)->SetTitle(chiSquareText);

  // Construct 1-sigma confidence intervals from fit results.
  TGraphErrors* tagAndProbeRatioWithAllErrors = new TGraphErrors(tagAndProbeRatio->GetN());
  for (int k = 0; k < tagAndProbeRatio->GetN(); ++k)
    tagAndProbeRatioWithAllErrors->SetPoint(k, tagAndProbeRatio->GetX()[k], 1.0);
  TVirtualFitter::GetFitter()->GetConfidenceIntervals(tagAndProbeRatioWithAllErrors, 0.683);

  // Only truncate energy range, after the fit was performed!
  TruncateGraphEnergyRange(tagAndProbeRatio, 0.0, gAMSTimeDependentFluxesMaxEnergy);

  if (significantCorrection)
    tagAndProbeRatioWithAllErrors->SetName(Form("cut_%i_bartels_rotation_average_tag_and_probe_correction_significant", averageData.fCutIndex));
  else
    tagAndProbeRatioWithAllErrors->SetName(Form("cut_%i_bartels_rotation_average_tag_and_probe_correction", averageData.fCutIndex));

  tagAndProbeRatioWithAllErrors->SetTitle(cutName.c_str());
  tagAndProbeRatioWithAllErrors->SetFillColor(kOrange);
  tagAndProbeRatioWithAllErrors->SetFillStyle(1001);
  tagAndProbeRatioWithAllErrors->SetLineColor(kBlack);
  tagAndProbeRatioWithAllErrors->SetLineWidth(4);
  tagAndProbeRatioWithAllErrors->SetLineStyle(7);

  tagAndProbeRatioWithAllErrors->Draw("E3.same");
  fitFunctionClone->Draw("L.same");
  tagAndProbeRatioWithAllErrors->Draw("LX.same");
  tagAndProbeRatio->Draw("P.same");

  // Create legend & statistics box on upper pad, remove statistics box from lower pad.
  canvas->cd(1);
  statisticsBoxClone->Draw();

  TLegend* legend = new TLegend(0.15, 0.036, 0.50, 0.58, nullptr, "brNDC");
  legend->SetTextSize(0.07);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(issTagAndProbeElectronBartelsRotation, "Bartels rotation", "p");
  legend->AddEntry(issTagAndProbeElectronAverage, "Average", "p");
  legend->AddEntry(tagAndProbeRatio, "Bartels rotation / Average", "p");
  legend->AddEntry(tagAndProbeRatioWithAllErrors, "Ratio fit", "fl");
  legend->Draw();

  gPad->Update();
  canvas->cd(2);
  gStyle->SetOptFit(0);
  tagAndProbeRatio->GetListOfFunctions()->Remove(statisticsBox);
  gPad->Update();

  canvas->SaveAs(Form("%s_%s_%s.png", canvas->GetName(), cutName.c_str(), bartelsRotationOutputFileSuffix.c_str()));
  canvas->SaveAs(Form("%s_%s_%s.pdf", canvas->GetName(), cutName.c_str(), bartelsRotationOutputFileSuffix.c_str()));
  return tagAndProbeRatioWithAllErrors;
}

void EfficiencyComparisionData::CorrectEfficiencyForTagAndProbeRatio(TGraphAsymmErrors* efficiency, TGraphErrors* tagAndProbeRatio) {

  assert(efficiency->GetN() == tagAndProbeRatio->GetN());
  for (int i = 0; i < efficiency->GetN(); ++i) {
    double x, y;
    efficiency->GetPoint(i, x, y);

    double exl = efficiency->GetErrorXlow(i);
    double exh = efficiency->GetErrorXhigh(i);
    double eyl = efficiency->GetErrorYlow(i);
    double eyh = efficiency->GetErrorYhigh(i);

    double correctionPoint, correctionValue;
    tagAndProbeRatio->GetPoint(i, correctionPoint, correctionValue);
    double correctionValueError = tagAndProbeRatio->GetErrorY(i);
    assert(std::abs(x - correctionPoint) < 1e-3);

    eyl = std::sqrt(std::pow(eyl / y, 2) + std::pow(correctionValueError / correctionValue, 2)) * (y * correctionValue);
    eyh = std::sqrt(std::pow(eyh / y, 2) + std::pow(correctionValueError / correctionValue, 2)) * (y * correctionValue);

    efficiency->SetPoint(i, x, y * correctionValue);
    efficiency->SetPointError(i, exl, exh, eyl, eyh);
  }
}

void EfficiencyComparisionData::CorrectEfficiencyForTagAndProbeRatio(TGraphAsymmErrors* efficiency, double tagAndProbeRatio, double tagAndProbeRatioUncertainty) {

  for (int i = 0; i < efficiency->GetN(); ++i) {
    double x, y;
    efficiency->GetPoint(i, x, y);

    double exl = efficiency->GetErrorXlow(i);
    double exh = efficiency->GetErrorXhigh(i);
    double eyl = efficiency->GetErrorYlow(i);
    double eyh = efficiency->GetErrorYhigh(i);

    eyl = std::sqrt(std::pow(eyl / y, 2) + std::pow(tagAndProbeRatioUncertainty / tagAndProbeRatio, 2)) * (y * tagAndProbeRatio);
    eyh = std::sqrt(std::pow(eyh / y, 2) + std::pow(tagAndProbeRatioUncertainty / tagAndProbeRatio, 2)) * (y * tagAndProbeRatio);

    efficiency->SetPoint(i, x, y * tagAndProbeRatio);
    efficiency->SetPointError(i, exl, exh, eyl, eyh);
  }
}

void EfficiencyComparisionData::Initialize(TFile* effectiveAcceptanceFileElectron, TFile* effectiveAcceptanceFileProton, TFile* issTagAndProbeEfficiencyFile, TFile* mcTagAndProbeEfficiencyFile, const std::string& selectorName, unsigned int numberOfCuts) {

  fSelectorName = selectorName;
  fNumberOfCuts = numberOfCuts;

  for (unsigned int cutIndex = 0; cutIndex < fNumberOfCuts; ++cutIndex) {
    SingleEfficiencyData cutData;
    cutData.Initialize(effectiveAcceptanceFileElectron, effectiveAcceptanceFileProton, issTagAndProbeEfficiencyFile, mcTagAndProbeEfficiencyFile, fSelectorName, cutIndex);
    fEfficiencies.emplace_back(cutData);
  }
}

void EfficiencyComparisionData::CalculateOverallEfficiency() {

  std::vector<TGraphAsymmErrors*> mcCountingElectronGeneratedEfficiencies;
  std::vector<TGraphAsymmErrors*> mcTagAndProbeElectronEfficiencies;
  std::vector<TGraphAsymmErrors*> issTagAndProbeElectronEfficiencies;
  for (unsigned int cut = 0; cut < fNumberOfCuts; ++cut) {
    if (fEfficiencies.at(cut).mcCountingElectronGenerated)
      mcCountingElectronGeneratedEfficiencies.push_back(fEfficiencies.at(cut).mcCountingElectronGenerated);

    if (fEfficiencies.at(cut).mcTagAndProbeElectron)
      mcTagAndProbeElectronEfficiencies.push_back(fEfficiencies.at(cut).mcTagAndProbeElectron);

    if (fEfficiencies.at(cut).issTagAndProbeElectron)
      issTagAndProbeElectronEfficiencies.push_back(fEfficiencies.at(cut).issTagAndProbeElectron);
  }

  fMCCountingElectronGeneratedOverallEfficiency = CombineUncorrelatedEfficiencies(mcCountingElectronGeneratedEfficiencies);
  assert(fMCCountingElectronGeneratedOverallEfficiency);

  fMCTagAndProbeElectronOverallEfficiency = CombineUncorrelatedEfficiencies(mcTagAndProbeElectronEfficiencies);
  fISSTagAndProbeElectronOverallEfficiency = CombineUncorrelatedEfficiencies(issTagAndProbeElectronEfficiencies);
}

void EfficiencyComparisionData::SetControlPlotYRange(unsigned int cut, float yMin, float yMax) {

  assert(cut < fNumberOfCuts);
  fEfficiencies.at(cut).SetControlPlotYRange(yMin, yMax);
}

void EfficiencyComparisionData::DrawSingleCutComparisions() {

  for (unsigned int cut = 0; cut < fNumberOfCuts; ++cut)
    fEfficiencies.at(cut).DrawOverview();
}

TGraphErrors* EfficiencyComparisionData::FitTagAndProbeRatio(unsigned int cut, const TagAndProbeRatioFitSettings& settings) {

  assert(cut < fNumberOfCuts);
  return fEfficiencies.at(cut).FitTagAndProbeRatio(settings.fitMode, settings.fitXMin, settings.fitXMax, settings.fitYMin, settings.fitYMax, settings.fitBreak1Minimum, settings.fitBreak1Maximum, settings.fitFixedBreak1, settings.fitFixedBreak2, settings.fitKeepRelativeErrorConstantAboveEnergy, settings.fitTrustSource);
}

void EfficiencyComparisionData::DrawOverallEfficiencyComparision() {

  static unsigned int gCounter = 0;
  ++gCounter;

  TCanvas* canvas = new TCanvas(Form("canvas_overall_efficiency_%i", gCounter), Form("%s: Overall efficiency", fSelectorName.c_str()));
  canvas->cd();

  assert(fMCCountingElectronGeneratedOverallEfficiency);

  gPad->SetGrid();
  gPad->SetLogx();

  fMCCountingElectronGeneratedOverallEfficiency->SetMarkerColor(kViolet);
  fMCCountingElectronGeneratedOverallEfficiency->SetMarkerStyle(24);
  fMCCountingElectronGeneratedOverallEfficiency->Draw("APE0");
  StyleEfficiencyGraph(fMCCountingElectronGeneratedOverallEfficiency);

  if (fMCTagAndProbeElectronOverallEfficiency) {
    fMCTagAndProbeElectronOverallEfficiency->SetMarkerColor(kBlue);
    fMCTagAndProbeElectronOverallEfficiency->SetMarkerStyle(22);
    fMCTagAndProbeElectronOverallEfficiency->Draw("P.SAME");
  }

  if (fISSTagAndProbeElectronOverallEfficiency) {
    fISSTagAndProbeElectronOverallEfficiency->SetMarkerColor(kRed);
    fISSTagAndProbeElectronOverallEfficiency->SetMarkerStyle(20);
    fISSTagAndProbeElectronOverallEfficiency->Draw("P.SAME");
  }

  TLegend* legend = new TLegend(0.2630522, 0.182243, 0.8082329, 0.4750779, nullptr, "brNDC");
  legend->SetTextSize(0.035);
  legend->SetFillColor(kWhite);
  legend->AddEntry(fMCCountingElectronGeneratedOverallEfficiency, "MC counting (E_{gen})", "pl");
  if (fMCTagAndProbeElectronOverallEfficiency && fMCTagAndProbeElectronOverallEfficiency->GetN() > 0)
    legend->AddEntry(fMCTagAndProbeElectronOverallEfficiency, "MC tag & probe", "pl");
  if (fISSTagAndProbeElectronOverallEfficiency && fISSTagAndProbeElectronOverallEfficiency->GetN() > 0)
    legend->AddEntry(fISSTagAndProbeElectronOverallEfficiency, "ISS tag & probe", "pl");
  legend->Draw();

  canvas->SaveAs(Form("%s_overview.png", canvas->GetName()));
  canvas->SaveAs(Form("%s_overview.pdf", canvas->GetName()));
}
