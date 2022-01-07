#include "TemplateFitTools.hh"

#include "AnalysisSettings.hh"

#include <sstream>
#include <iomanip>

#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>

TemplateFitTools::TemplateFitTools()
  : fDataDistribution(nullptr)
  , fSignalTemplate(nullptr)
  , fBackgroundTemplate(nullptr)
  , fTemplateFitter(-1 /* silent */)
  , fChiSquare(-1.0)
  , fNDF(-1) {

}

bool TemplateFitTools::PerformTemplateFit(TH1D* dataDistribution, TH1D* signalTemplate, TH1D* backgroundTemplate) {

  fTemplateFitter.Clear();

  fDataDistribution = dataDistribution;
  fTemplateFitter.SetDataHistogram(dataDistribution);
  fTemplateFitter.AddTemplateHistogram(signalTemplate);
  fTemplateFitter.AddTemplateHistogram(backgroundTemplate);
  fTemplateFitter.AutoStartValues();

  fTemplateFitter.Fit(1 /* likelihood fit */);
  fChiSquare = fTemplateFitter.Chi2();
  fNDF = fTemplateFitter.NDF();
  fSignalResults = Utilities::Quantity();
  fBackgroundResults = Utilities::Quantity();

  if (!fTemplateFitter.HasBeenFit()) {
    fChiSquare = -1.0;
    fNDF = -1;
    return false;
  }

  const std::vector<double>& absoluteResult = fTemplateFitter.GetAbsoluteResult();
  assert(absoluteResult.size() == 2);

  const std::vector<double>& absoluteResultError = fTemplateFitter.GetAbsoluteResultError();
  assert(absoluteResultError.size() == 2);

  fSignalResults.value = absoluteResult[0];
  fSignalResults.uncertainty = absoluteResultError[0];

  fBackgroundResults.value = absoluteResult[1];
  fBackgroundResults.uncertainty = absoluteResultError[1];

  return true;
}

TCanvas* TemplateFitTools::DrawResults(const std::string& canvasName, const std::string& canvasTitle) {

  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasTitle.c_str());
  canvas->cd();
  gPad->SetLogy();

  TLegend* legend = new TLegend(0.7204609,0.7240356,0.997996,0.9970326,NULL,"brNDC");
  legend->SetFillColor(kWhite);

  const TH1* fitResult = fTemplateFitter.FitResult();
  const std::vector<TH1*>& resultHistograms = fTemplateFitter.ResultHistograms();
  assert(resultHistograms.size() == 2);

  TH1* dataCopy = fDataDistribution->DrawCopy();
  dataCopy->SetStats(0);
  dataCopy->SetLineWidth(3);
  dataCopy->SetLineColor(kBlack);
  dataCopy->SetMarkerColor(kBlack);
  dataCopy->SetMarkerStyle(kFullCircle);
  dataCopy->SetMarkerSize(2.0);
  legend->AddEntry(dataCopy, dataCopy->GetTitle(), "LP");

  TH1* fitResultCopy = fitResult->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  fitResultCopy->SetLineWidth(5);
  fitResultCopy->SetLineColor(kGreen + 2);
  fitResultCopy->SetLineStyle(7);
  legend->AddEntry(fitResultCopy, Form("#splitline{Fit result}{#chi^{2}/ndf %5.2f / %d = %5.2f}", fChiSquare, fNDF, fChiSquare / fNDF), "L");

  TH1* signalResultHistogram = resultHistograms[0]->DrawCopy("hist.same");
  signalResultHistogram->SetStats(0);
  signalResultHistogram->SetLineWidth(2);
  signalResultHistogram->SetFillStyle(1001);
  signalResultHistogram->SetFillColorAlpha(kBlue, 0.65);
  signalResultHistogram->SetLineColor(kBlue);
  legend->AddEntry(signalResultHistogram, Form("%s #color[%i]{%.2f #pm %.2f}", signalResultHistogram->GetTitle(), signalResultHistogram->GetLineColor(), fSignalResults.value, fSignalResults.uncertainty), "F");

  TH1* backgroundResultHistogram = resultHistograms[1]->DrawCopy("hist.same");
  backgroundResultHistogram->SetStats(0);
  backgroundResultHistogram->SetLineWidth(2);
  backgroundResultHistogram->SetFillStyle(1001);
  backgroundResultHistogram->SetFillColorAlpha(kRed, 0.65);
  backgroundResultHistogram->SetLineColor(kRed);
  legend->AddEntry(backgroundResultHistogram, Form("%s #color[%i]{%.2f #pm %.2f}", backgroundResultHistogram->GetTitle(), backgroundResultHistogram->GetLineColor(), fBackgroundResults.value, fBackgroundResults.uncertainty), "F");

  // Redraw fit result and data on top of the templates
  fitResultCopy->Draw("hist.same");
  dataCopy->Draw("same");

  legend->Draw();
  return canvas;
}
