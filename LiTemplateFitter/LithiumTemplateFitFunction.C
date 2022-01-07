#include "LithiumTemplateFitFunction.hh"

#include <cassert>
#include <TH1.h>
#include <TMath.h>

#include "Quantity.hh"

#define INFO_OUT_TAG "LithiumTemplateFitFunction"
#include "debugging.hh"

LithiumTemplateFitFunction::LithiumTemplateFitFunction(const TH1* dataDistribution,
							   const TH1* li6Template, const TH1* li7Template,
							   double rangeMinimum, double rangeMaximum, FitType fitType)
: fDataDistribution(dataDistribution) {

  TH1* tempLi6Template = (TH1*)li6Template->Clone();
  tempLi6Template->Scale(1.0 / li6Template->Integral());
  tempLi6Template->Smooth();
  fLi6Template = tempLi6Template;

  TH1* tempLi7Template = (TH1*)li7Template->Clone();
  tempLi7Template->Scale(1.0 / li7Template->Integral());
  tempLi7Template->Smooth();
  fLi7Template = tempLi7Template;

  fRangeMinimum = rangeMinimum;
  fRangeMaximum = rangeMaximum;

  if (fitType == FitRelativeFractions) {
    fFitRatios = false;
  }
  else {
    fFitRatios = true;
  }

  if (rangeMinimum != 0.0 && rangeMaximum != 0.0) {
    fBinMinimum = std::max(fDataDistribution->FindFixBin(rangeMinimum), 1);
    fBinMaximum = std::min(fDataDistribution->FindFixBin(fRangeMaximum), fDataDistribution->GetNbinsX());
    fNumberOfEvents = dataDistribution->Integral(fBinMinimum, fBinMaximum);
  }
  else
    fNumberOfEvents = dataDistribution->Integral();
}

int LithiumTemplateFitFunction::NumberOfEvents() const {

  return fNumberOfEvents;
}

std::string LithiumTemplateFitFunction::FitMode() const {

  if (fFitRatios)
    return "Ratios";

  return "RelativeFractions";
}

double LithiumTemplateFitFunction::PdfValue(int xBin) const {

  double nLi6 = GetLi6Events();
  double nLi7 = GetLi7Events();

  double pdf = nLi6 * fLi6Template->GetBinContent(xBin) + nLi7 * fLi7Template->GetBinContent(xBin);

  assert(fNumberOfEvents >= 0);
  assert(!std::isinf(pdf));
  assert(!std::isnan(pdf));

  return pdf / fNumberOfEvents;
}

double LithiumTemplateFitFunction::GetLi6Events() const {

  if (fFitRatios)
    return fLi6Li7Ratio * fLi7Fraction * fNumberOfEvents;

  return fLi6Fraction * fNumberOfEvents;
}

double LithiumTemplateFitFunction::GetLi7Events() const {

  return fLi7Fraction * fNumberOfEvents;
}


Utilities::Quantity LithiumTemplateFitFunction::GetLi6Fraction() const {

  return {fLi6Fraction, fLi6FractionError};
}

Utilities::Quantity LithiumTemplateFitFunction::GetLi7Fraction() const {

  return {fLi7Fraction, fLi7FractionError};
}


Utilities::Quantity LithiumTemplateFitFunction::GetLi6Li7Ratio() const {

  return {fLi6Li7Ratio, fLi6Li7RatioError};
}


double LithiumTemplateFitFunction::DoEval(const double* parameters) const {

  const_cast<LithiumTemplateFitFunction*>(this)->UpdateParameters(parameters);

  double minusLogL = 0.0;

  int binMinimum = 1;
  int binMaximum = fDataDistribution->GetNbinsX();

  if (fRangeMinimum != 0 && fRangeMaximum != 0) {
    binMinimum = std::max(fDataDistribution->FindFixBin(fRangeMinimum), 1);
    binMaximum = std::min(fDataDistribution->FindFixBin(fRangeMaximum), fDataDistribution->GetNbinsX());
  }

  for (int iBin = binMinimum ; iBin <= binMaximum; ++iBin) {
    int binContent = TMath::Nint(fDataDistribution->GetBinContent(iBin));
    if (!binContent)
      continue;

    double pdfValue =  PdfValue(iBin);
    if (pdfValue <= 0.0)
      continue;

    minusLogL += (binContent * (-std::log(pdfValue)));
  }

  // term from extended ML, to take normalization into account
  // http://hepunx.rl.ac.uk/~adye/thesis/html/node49.html
  // or Cowan pp 84f
  minusLogL += GetLi6Events() + GetLi7Events();

  return minusLogL;
}

void LithiumTemplateFitFunction::UpdateParameters(const double* parameters) {

  if (fFitRatios) {
    fLi6Li7Ratio = parameters[0];
    fLi7Fraction = parameters[1];
  }
  else {
    fLi6Fraction = parameters[0];
    fLi7Fraction = parameters[1];
  }

}

void LithiumTemplateFitFunction::UpdateErrors(const double* errors) {

  if (fFitRatios) {
    fLi6Li7RatioError = errors[0];
    fLi7FractionError = errors[1];
  }
  else {
    fLi6FractionError = errors[0];
    fLi7FractionError = errors[1];
  }
}

TH1* LithiumTemplateFitFunction::DataHistogram() const {

  TH1* dataHisto = (TH1*)fDataDistribution->Clone();
  return dataHisto;
}

TH1* LithiumTemplateFitFunction::CalculateHistogram(std::string species) const {

  double nLi6 = GetLi6Events();
  double nLi7 = GetLi7Events();

  TH1* histogram = dynamic_cast<TH1*>(fDataDistribution->Clone(Form("hTemplateFitResult_%s", species.c_str())));
  histogram->Reset();
  histogram->SetTitle("Lithium template fit result");
  histogram->SetStats(0);

  for (int xBin = fBinMinimum; xBin <= fBinMaximum; ++xBin) {
    double prediction = 0.0;
    if (species == "Li6")
      prediction += nLi6 * fLi6Template->GetBinContent(xBin);
    else if (species == "Li7")
      prediction += nLi7 * fLi7Template->GetBinContent(xBin);

    else
      INFO_OUT << "Species name is not valid." << std::endl;

    histogram->SetBinContent(xBin, prediction);
  }

  return histogram;
}

TH1* LithiumTemplateFitFunction::CalculateTotalHistogram() const {

  double nLi6 = GetLi6Events();
  double nLi7 = GetLi7Events();

  TH1* histogram = dynamic_cast<TH1*>(fDataDistribution->Clone("hTemplateFitResult_Total"));
  histogram->Reset();
  histogram->SetTitle("Lithium template fit result");
  histogram->SetStats(0);

  for (int xBin = fBinMinimum; xBin <= fBinMaximum; ++xBin) {

    double prediction = 0.0;
    prediction += nLi6 * fLi6Template->GetBinContent(xBin);
    prediction += nLi7 * fLi7Template->GetBinContent(xBin);

    histogram->SetBinContent(xBin, prediction);
  }

  return histogram;
}

TH1* LithiumTemplateFitFunction::DrawZeroProbabilityEvents() const {

  TH1* histogram = dynamic_cast<TH1*>(fDataDistribution->Clone("hZeroProbability"));
  histogram->Reset();
  histogram->SetStats(0);
  for (int xBin = 1; xBin <= histogram->GetNbinsX(); ++xBin) {
    double pdfLi6 = fLi6Template->GetBinContent(xBin);
    double pdfLi7 = fLi7Template->GetBinContent(xBin);

    if (pdfLi6 <= 0.0 && pdfLi7 <= 0.0)
      histogram->SetBinContent(xBin, fDataDistribution->GetBinContent(xBin));
  }

  histogram->SetTitle(Form("%g events with zero probability", histogram->Integral()));
  return histogram;
}
