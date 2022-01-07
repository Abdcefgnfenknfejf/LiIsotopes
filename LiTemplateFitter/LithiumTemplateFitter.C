#include "LithiumTemplateFitter.hh"

#include <cassert>
#include <string>
#include <vector>

#include "BinningTools.hh"
#include "Quantity.hh"
#include "TemplateFitter.hh"
#include "Utilities.hh"

#include "LithiumTemplateFitFunction.hh"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include <TColor.h>
#include <QString>
#include <QVector>

#define INFO_OUT_TAG "LithiumTemplateFitter"
#include "debugging.hh"

LithiumTemplateFitter::LithiumTemplateFitter(FitType fitType)
: fDataHisto(nullptr) {

  fFitType = fitType;

}

void LithiumTemplateFitter::SetDataHistogram(TH1* dataHisto) {

  assert(dataHisto);
  fDataHisto = dataHisto;
}

void LithiumTemplateFitter::SetTemplateHistogram(TH1* templateHisto) {

  assert(templateHisto);
  fTemplateHistos.push_back(templateHisto);
}

void LithiumTemplateFitter::SetFitRange(double fitRangeMinimum, double fitRangeMaximum) {

  fFitRangeMinimum = fitRangeMinimum;
  fFitRangeMaximum = fitRangeMaximum;
  fFitBinMinimum = std::max(fDataHisto->FindFixBin(fFitRangeMinimum), 1);
  fFitBinMaximum = std::min(fDataHisto->FindFixBin(fFitRangeMaximum), fDataHisto->GetNbinsX());
}

ROOT::Math::Minimizer* LithiumTemplateFitter::SetupMinimizer(LithiumTemplateFitFunction& fitFunction, FitNormParameters& fitParameters,
							       double toleranceFactor, double stepSizeFactor) {

  static ROOT::Math::Minimizer* sMinimizer = nullptr;
    if (!sMinimizer) {
      sMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
      sMinimizer->SetMaxFunctionCalls(100000);
      sMinimizer->SetPrintLevel(-1); // for debugging purposes, set to -1 when done.
      sMinimizer->SetStrategy(2);
      sMinimizer->SetTolerance(toleranceFactor);
    }

    sMinimizer->Clear();
    sMinimizer->SetFunction(fitFunction);
    sMinimizer->SetErrorDef(0.5);  // Likelihood minimization
    sMinimizer->SetLimitedVariable(0, "li6Fraction", fitParameters.fractionLi6, stepSizeFactor * fitParameters.fractionLi6, 0.0, 1.0);
    sMinimizer->SetLimitedVariable(1, "li7Fraction", fitParameters.fractionLi7, stepSizeFactor * fitParameters.fractionLi7, 0.0, 1.0);


    return sMinimizer;
}

ROOT::Math::Minimizer* LithiumTemplateFitter::SetupMinimizer(LithiumTemplateFitFunction& fitFunction, FitRatioParameters& fitParameters,
							       double toleranceFactor, double stepSizeFactor) {

  static ROOT::Math::Minimizer* sMinimizer = nullptr;
    if (!sMinimizer) {
      sMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
      sMinimizer->SetMaxFunctionCalls(1000000);
      sMinimizer->SetPrintLevel(-1); // for debugging purposes, set to -1 when done.
      sMinimizer->SetStrategy(2);
      sMinimizer->SetTolerance(toleranceFactor);
    }

    sMinimizer->Clear();
    sMinimizer->SetFunction(fitFunction);
    sMinimizer->SetErrorDef(0.5);  // Likelihood minimization
    sMinimizer->SetLimitedVariable(0, "li6li7Ratio", fitParameters.ratioLi6Li7, stepSizeFactor * fitParameters.ratioLi6Li7, 0.0, 5.0);
    sMinimizer->SetLimitedVariable(1, "li7Fraction", fitParameters.fractionLi7, stepSizeFactor * fitParameters.fractionLi7, 0.0, 1.2);

    return sMinimizer;
}


int LithiumTemplateFitter::TemplateColor(int i) const {
  QVector<int> firstColors;
  firstColors << kAzure - 3 << kPink << kGreen << kYellow << kOrange;
  if (i < firstColors.size())
    return firstColors.at(i);
  return Utilities::RootColor(i);
}


double LithiumTemplateFitter::CalculateChi2(TH1* dataHisto, TH1* modelHisto, int freeParameters) {

  double chi2 = 0.0;
  int effectiveNumberBins = 0.0;

  for (int iBin = fFitBinMinimum; iBin <= fFitBinMaximum; ++iBin) {
    double dataEntries = dataHisto->GetBinContent(iBin);
    if (TMath::Nint(dataEntries) == 0)
      continue;

    double modelEntries = modelHisto->GetBinContent(iBin);
    double dataError = dataHisto->GetBinError(iBin);
    ++effectiveNumberBins;

    chi2 += pow((dataEntries - modelEntries) / dataError, 2);
  }

  fNDF = effectiveNumberBins - freeParameters;
  return chi2 / (effectiveNumberBins - freeParameters);
}

bool LithiumTemplateFitter::PerformTemplateFit(const FitNormParameters& startParameters, FitResults& fitResult) {

  CheckHistogramConsistency();
  ResetFitResults(fitResult);

  LithiumTemplateFitFunction fitFunction(fDataHisto, fTemplateHistos.at(0), fTemplateHistos.at(1), fFitRangeMinimum, fFitRangeMaximum, fFitType);

  FitNormParameters initialValues = startParameters;
  const double toleranceFactor = 1.0;
  const double stepSizeFactor = 1e-2;
  ROOT::Math::Minimizer* minimizer = SetupMinimizer(fitFunction, initialValues, toleranceFactor, stepSizeFactor);
  if (minimizer->Minimize()) {
    FillResults(minimizer, fitFunction, fitResult);
    return true;
  }

  return false;
}

bool LithiumTemplateFitter::PerformTemplateFitterFit(FitResults& fitResult) {

  if (fFitType != FitRelativeFractions)
    FATAL_OUT << "In this mode you can only fit normalisation, fitRatios must be false";

  ResetFitResults(fitResult);

  Utilities::TemplateFitter templateFitter(-1);
  templateFitter.UseUnboundedFit(false);
  templateFitter.SetDataHistogram(fDataHisto);
  for (int iTemplate = 0; iTemplate < static_cast<int>(fTemplateHistos.size()); iTemplate++)
    templateFitter.AddTemplateHistogram(fTemplateHistos.at(iTemplate));
  //  templateFitter.AutoStartValues();
  fitResult.chiSquare = templateFitter.Fit(1, fFitRangeMinimum, fFitRangeMaximum);

  if (!templateFitter.HasBeenFit())
    return false;

  if (fFitType == FitRatios)
    fitResult.fitType = "Ratios"; // temporary solution, the fit types are 3 now.
  fitResult.numberEvents = fDataHisto->Integral(templateFitter.FitRangeLow(), templateFitter.FitRangeHigh());

  for (int iTemplate = 0; iTemplate < static_cast<int>(fTemplateHistos.size()); iTemplate++) {
    fitResult.relativeFractions.push_back({templateFitter.GetRelativeResult().at(iTemplate), templateFitter.GetRelativeResultError().at(iTemplate)});
    fitResult.absoluteNumberEvents.push_back({templateFitter.GetAbsoluteResult().at(iTemplate), templateFitter.GetAbsoluteResultError().at(iTemplate)});
  }

  // convert relative fractions in ratios and absolute number of events for comparison
  //To be modified
  fitResult.ratios.push_back({templateFitter.GetRelativeResult().at(0) / templateFitter.GetRelativeResult().at(1),
	templateFitter.GetRelativeResult().at(0) / templateFitter.GetRelativeResult().at(1) * sqrt(pow(templateFitter.GetRelativeResultError().at(0)/templateFitter.GetRelativeResult().at(0), 2) + pow(templateFitter.GetRelativeResultError().at(1)/templateFitter.GetRelativeResult().at(1), 2))});
  fitResult.ratios.push_back({templateFitter.GetRelativeResult().at(1), templateFitter.GetRelativeResultError().at(1)});
  //  fitResult.ratios.push_back({templateFitter.GetRelativeResult().at(1) / templateFitter.GetRelativeResult().at(0),
	//	templateFitter.GetRelativeResult().at(2) / templateFitter.GetRelativeResult().at(1) * sqrt(pow(templateFitter.GetRelativeResultError().at(2)/templateFitter.GetRelativeResult().at(2), 2) + pow(templateFitter.GetRelativeResultError().at(1)/templateFitter.GetRelativeResult().at(1), 2))});


  fitResult.ndf = fTemplateHistos.size();
  fitResult.dataHisto = (TH1*)templateFitter.DataHistogram();
  fitResult.resultHisto = (TH1*)templateFitter.FitResult();
  fitResult.templateHistos = templateFitter.ResultHistograms();

  return true;
}

bool LithiumTemplateFitter::PerformTemplateFit(const FitRatioParameters& startParameters, FitResults& fitResult) {

  CheckHistogramConsistency();
  ResetFitResults(fitResult);

  LithiumTemplateFitFunction fitFunction(fDataHisto, fTemplateHistos.at(0), fTemplateHistos.at(1), fFitRangeMinimum, fFitRangeMaximum, fFitType);
  FitRatioParameters initialValues = startParameters;
  const double toleranceFactor = 1.0;
  const double stepSizeFactor = 1e-2;

  ROOT::Math::Minimizer* minimizer = SetupMinimizer(fitFunction, initialValues, toleranceFactor, stepSizeFactor);
  if (minimizer->Minimize()) {
    FillResults(minimizer, fitFunction, fitResult);
    if (IsPhysicalResult(fitResult))
      return true;
  }

  return false;
}


void LithiumTemplateFitter::ResetFitResults(FitResults& fitResults) {

  fitResults.fitType = "RelativeFractions";
  fitResults.numberEvents = -1.0;
  fitResults.relativeFractions.clear();
  fitResults.absoluteNumberEvents.clear();
  fitResults.ratios.clear();
  fitResults.dataHisto = nullptr;
  fitResults.templateHistos.clear();
  fitResults.chiSquare = -3.0;
}

bool LithiumTemplateFitter::CheckHistogramConsistency() {

  for (unsigned int iTemplate = 0; iTemplate < fTemplateHistos.size(); iTemplate++) {
    if (fTemplateHistos.at(iTemplate)->GetNbinsX() != fDataHisto->GetNbinsX())
      return false;

  if (Binning::FromTAxis(fTemplateHistos.at(iTemplate)->GetXaxis()) != Binning::FromTAxis(fDataHisto->GetXaxis()))
      return false;
  }

  return true;
}

void LithiumTemplateFitter::FillResults(ROOT::Math::Minimizer* minimizer, LithiumTemplateFitFunction& fitFunction, FitResults& fitResult) {

  const double* par = minimizer->X();
  const double* parErrors = minimizer->Errors();

  fitFunction.UpdateParameters(par);
  fitFunction.UpdateErrors(parErrors);

  fitResult.fitType = fitFunction.FitMode();
  fitResult.numberEvents = fitFunction.NumberOfEvents();

  if (fitResult.fitType == "Ratios") {
    for (int iPar = 0; iPar < static_cast<int>(fitFunction.NDim()); iPar++){
      fitResult.ratios.emplace_back(par[iPar], parErrors[iPar]);
    }
      // convert ratios in relative fractions and absolute number of events for comparison
      fitResult.relativeFractions.emplace_back(fitFunction.GetLi7Fraction().value * fitFunction.GetLi6Li7Ratio().value,
					       fitFunction.GetLi7Fraction().value * fitFunction.GetLi6Li7Ratio().value *
					       sqrt(pow(fitFunction.GetLi6Li7Ratio().uncertainty/fitFunction.GetLi6Li7Ratio().value, 2) +
						    pow(fitFunction.GetLi7Fraction().uncertainty/fitFunction.GetLi7Fraction().value, 2)));
      fitResult.relativeFractions.emplace_back(fitFunction.GetLi7Fraction().value, fitFunction.GetLi7Fraction().uncertainty);
      
      fitResult.absoluteNumberEvents.emplace_back(fitFunction.NumberOfEvents() * fitFunction.GetLi7Fraction().value * fitFunction.GetLi6Li7Ratio().value,
						  fitFunction.NumberOfEvents() * fitFunction.GetLi7Fraction().value * fitFunction.GetLi6Li7Ratio().value *
						  sqrt(pow(fitFunction.GetLi6Li7Ratio().uncertainty/fitFunction.GetLi6Li7Ratio().value, 2) +
						       pow(fitFunction.GetLi7Fraction().uncertainty/fitFunction.GetLi7Fraction().value, 2)));
      fitResult.absoluteNumberEvents.emplace_back(fitFunction.NumberOfEvents() * fitFunction.GetLi7Fraction().value,
					      fitFunction.NumberOfEvents() * fitFunction.GetLi7Fraction().uncertainty);
      
  }
  else {
    for (int iPar = 0; iPar < static_cast<int>(fitFunction.NDim()); iPar++) {
      fitResult.relativeFractions.emplace_back(par[iPar], parErrors[iPar]);
      fitResult.absoluteNumberEvents.emplace_back(fitFunction.NumberOfEvents() * par[iPar], fitFunction.NumberOfEvents() * parErrors[iPar]);
    }
    
    // convert relative fractions in ratios and absolute number of events for comparison
    fitResult.ratios.emplace_back(fitFunction.GetLi6Fraction().value / fitFunction.GetLi7Fraction().value,
				  fitFunction.GetLi6Fraction().value / fitFunction.GetLi7Fraction().value *
				  sqrt(pow(fitFunction.GetLi6Fraction().uncertainty/fitFunction.GetLi6Fraction().value, 2) + pow(fitFunction.GetLi7Fraction().uncertainty/fitFunction.GetLi7Fraction().value, 2)));
    fitResult.ratios.emplace_back(fitFunction.GetLi7Fraction().value, fitFunction.GetLi7Fraction().uncertainty);
  }
  
  fitResult.dataHisto = fitFunction.DataHistogram();
  fitResult.resultHisto = fitFunction.CalculateTotalHistogram();
  for (int iName = 0; iName < static_cast<int>(fTemplateNames.size()); iName++) {
    TH1* resultTemplate = fitFunction.CalculateHistogram(fTemplateNames.at(iName));
    fitResult.templateHistos.push_back(resultTemplate);
  }
  fitResult.chiSquare = CalculateChi2(fitFunction.DataHistogram(), fitFunction.CalculateTotalHistogram(), fitFunction.NDim());
  fitResult.ndf = fNDF;
}
  
bool LithiumTemplateFitter::IsPhysicalResult(FitResults& fitResult) {
    
  for (int iResult = 0; iResult < static_cast<int> (fitResult.relativeFractions.size()); iResult++) {
    if (std::isnan(fitResult.relativeFractions.at(iResult).value) || std::isnan(fitResult.relativeFractions.at(iResult).uncertainty) ||
	std::isnan(fitResult.ratios.at(iResult).value) || std::isnan(fitResult.ratios.at(iResult).uncertainty) ||
	fitResult.ratios.at(iResult).value < 0 || fitResult.relativeFractions.at(iResult).value < 0 || fitResult.relativeFractions.at(iResult).value > 1)
      return false;
  }

  return true;
}

TCanvas* LithiumTemplateFitter::CreateResultDrawing(FitResults& fitResult, std::string canvasName) {

  TCanvas* graphicsOutput = new TCanvas(canvasName.c_str(), "", 10, 10, 800, 600);
  graphicsOutput->cd();

  TPad* padFits = new TPad("padFits", "padFits", 0.0, 0.30, 1.0, 1.0);
  padFits->SetBottomMargin(0.01);
  padFits->Draw();
  padFits->cd();
  padFits->SetGridx();
  padFits->SetLogy();
  TLegend* legend = new TLegend(0.6, 0.4, 1 - gStyle->GetPadRightMargin(), 1 - gStyle->GetPadTopMargin(), NULL, "brNDC");
  legend->SetFillColor(kWhite);
  int colors[5] = {kBlack, kRed, kBlue, kOrange + 1, kGreen + 1};

  TH1* dataHistogram = fitResult.dataHisto->DrawCopy();
  dataHistogram->SetStats(0);
  double maximumY = dataHistogram->GetMaximum() * 1.1;
  dataHistogram->GetYaxis()->SetRangeUser(1, maximumY);
  dataHistogram->SetLineColor(colors[0]);
  dataHistogram->SetMarkerColor(colors[0]);
  dataHistogram->SetMarkerStyle(kFullCircle);
  dataHistogram->SetMarkerSize(0.8);

  legend->AddEntry(dataHistogram, Form("ISS Data (N = %i)", (int)dataHistogram->GetEntries()));

  TH1* combinedTemplateHistogram = fitResult.resultHisto->DrawCopy("hist same");
  combinedTemplateHistogram->SetStats(0);
  combinedTemplateHistogram->SetLineWidth(2);
  combinedTemplateHistogram->SetLineColor(colors[1]);
  combinedTemplateHistogram->SetLineStyle(kDashed);;
  int ndf = fitResult.ndf;
  double chiSquared = fitResult.chiSquare;
  legend->AddEntry(combinedTemplateHistogram, Form("#chi^{2}/ndf=%3.2f/%d=%3.2f", chiSquared * ndf, ndf, chiSquared), "l");

  // Individual template histograms.
  for (unsigned int histo = 0; histo < fitResult.templateHistos.size(); histo++) {
    TH1* templateHistogram = fitResult.templateHistos.at(histo)->DrawCopy("hist same");
    const int colorTemplate = TemplateColor(histo);
    templateHistogram->SetStats(0);
    templateHistogram->SetLineColor(colorTemplate);
    templateHistogram->SetFillColorAlpha(colorTemplate, 0.40);
    templateHistogram->SetFillStyle(3002);
    templateHistogram->SetMarkerColor(colorTemplate);

    TString fittedValue = Form("(%.3f #pm %.3f)%%", 100 * fitResult.relativeFractions.at(histo).value, 100 * fitResult.relativeFractions.at(histo).uncertainty);

    legend->AddEntry(templateHistogram, Form("%s fraction = %s", fTemplateNames.at(histo).c_str(), fittedValue.Data()), "f");
    if (histo != 1) {
      TString fittedRatios = Form("%.3f #pm %.3f", fitResult.ratios.at(histo).value, fitResult.ratios.at(histo).uncertainty);
      legend->AddEntry(templateHistogram, Form("Li6/Li7 = %s",fittedRatios.Data()), "l");
    }
  }

  // Redraw combined template fit and data histogram to show them on top.
  combinedTemplateHistogram->Draw("hist same");
  dataHistogram->Draw("same");

  TLine *lineFitRange = new TLine();
  lineFitRange->SetLineStyle(2);
  lineFitRange->SetLineColor(kRed);

  //  if (fFitRangeMinimum != 0.0 && fFitRangeMaximum != 0.0) { // fit range set
  //  lineFitRange->DrawLine(fFitRangeMinimum, 0, fFitRangeMinimum, maximumY);
  //   lineFitRange->DrawLine(fFitRangeMaximum, 0, fFitRangeMaximum, maximumY);
  //  }
  //  legend->AddEntry(lineFitRange, "Fit range", "l");
  legend->Draw();
  graphicsOutput->cd();

  TLine* zeroLine = new TLine(2, 0, 13, 0);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(kDashed);
  zeroLine->SetLineWidth(2);
  TPad* padPull = new TPad("padPull", "padPull", 0.0, 0.02, 1.0, 0.29);
  padPull->SetTopMargin(0.0);
  padPull->SetBottomMargin(0.2);
  padPull->SetGridx();
  padPull->Draw();
  padPull->cd();
  TH1* pullHistogram = (TH1*)dataHistogram->Clone();
  pullHistogram->Reset();

  for (int iBin = 0; iBin < dataHistogram->GetNbinsX(); iBin++) {
    pullHistogram->SetBinContent(iBin, (dataHistogram->GetBinContent(iBin) - combinedTemplateHistogram->GetBinContent(iBin)) / (dataHistogram->GetBinError(iBin) != 0 ? dataHistogram->GetBinError(iBin) : 1.0));
  }

  pullHistogram->GetXaxis()->SetLabelSize(0.09);

  pullHistogram->GetXaxis()->SetTitleSize(0.10);
  pullHistogram->GetYaxis()->SetRangeUser(-4, 4);
  pullHistogram->GetYaxis()->SetNdivisions(8);
  pullHistogram->GetYaxis()->SetLabelSize(0.09);
  pullHistogram->GetYaxis()->SetTitleSize(0.10);
  pullHistogram->GetYaxis()->SetTitle("(model - data) / #sigma");
  pullHistogram->SetMarkerStyle(kFullSquare);
  pullHistogram->SetMarkerSize(0.8);
  pullHistogram->SetLineWidth(2);
  pullHistogram->SetLineColor(kRed );
  pullHistogram->SetMarkerColor(kRed);
  pullHistogram->Draw("p e1");
  zeroLine->Draw("same");

  gPad->Modified();
  gPad->Update();

  return graphicsOutput;

}
