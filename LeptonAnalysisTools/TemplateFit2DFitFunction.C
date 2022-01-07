#include "TemplateFit2DFitFunction.hh"

#include <cassert>
#include <TH2.h>
#include <TMath.h>

TemplateFit2DFitFunction::TemplateFit2DFitFunction(const TH2* dataDistribution,
                                                   const TH2* electronTemplate, const TH2* positronTemplate,
                                                   const TH2* ccElectronTemplate, const TH2* ccPositronTemplate,
                                                   const TH2* protonTemplate, const TH2* ccProtonTemplate)
  : fDataDistribution(dataDistribution)
  , fElectronTemplate(electronTemplate)
  , fPositronTemplate(positronTemplate)
  , fCCElectronTemplate(ccElectronTemplate)
  , fCCPositronTemplate(ccPositronTemplate)
  , fProtonTemplate(protonTemplate)
  , fCCProtonTemplate(ccProtonTemplate) {

  fNumberOfEvents = dataDistribution->Integral();
}

double TemplateFit2DFitFunction::PdfValue(int xBin, int yBin) const {

  double nProtons   = fProtonFraction   * fNumberOfEvents;
  double nCCProtons = fCCProtonFraction * fNumberOfEvents;
  double nElectrons = fElectronFraction * fNumberOfEvents;
  double nPositrons = fPositronFraction * fNumberOfEvents;

  double nElectronsCorrectCharge = (1.0 - fChargeConfusionProbability) * nElectrons;
  double nElectronsWrongCharge   =        fChargeConfusionProbability  * nElectrons;

  double nPositronsCorrectCharge = (1.0 - fChargeConfusionProbability) * nPositrons;
  double nPositronsWrongCharge   =        fChargeConfusionProbability  * nPositrons;

  double pdf = nProtons                * fProtonTemplate->GetBinContent(xBin, yBin)
             + nCCProtons              * fCCProtonTemplate->GetBinContent(xBin, yBin)
             + nElectronsCorrectCharge * fElectronTemplate->GetBinContent(xBin, yBin)
             + nElectronsWrongCharge   * fCCElectronTemplate->GetBinContent(xBin, yBin)
             + nPositronsCorrectCharge * fPositronTemplate->GetBinContent(xBin, yBin)
             + nPositronsWrongCharge   * fCCPositronTemplate->GetBinContent(xBin, yBin);

  assert(fNumberOfEvents >= 0);
  assert(!std::isinf(pdf));
  assert(!std::isnan(pdf));

  return pdf / fNumberOfEvents;
}

double TemplateFit2DFitFunction::DoEval(const double* parameters) const {

  const_cast<TemplateFit2DFitFunction*>(this)->UpdateParameters(parameters);

  double minusLogL = 0.0;
  for (int i = 1 ; i <= fDataDistribution->GetNbinsX(); ++i) {
    for (int j = 1; j <= fDataDistribution->GetNbinsY(); ++j) {
      int binContent = TMath::Nint(fDataDistribution->GetBinContent(i, j));
      if (!binContent)
        continue;

      double pdfValue = PdfValue(i, j);
      if (pdfValue <= 0.0)
        continue;

      minusLogL += (binContent * (-std::log(pdfValue)));
    }
  }

  // term from extended ML, to take normalization into account
  // http://hepunx.rl.ac.uk/~adye/thesis/html/node49.html
  // or Cowan pp 84f
  minusLogL += (fProtonFraction + fCCProtonFraction + fElectronFraction + fPositronFraction) * fNumberOfEvents;

  return minusLogL;
}

TH2* TemplateFit2DFitFunction::CalculateHistogram(int contributingCategories) const {

  double nProtons   = fProtonFraction   * fNumberOfEvents;
  double nCCProtons = fCCProtonFraction * fNumberOfEvents;
  double nElectrons = fElectronFraction * fNumberOfEvents;
  double nPositrons = fPositronFraction * fNumberOfEvents;

  double nElectronsCorrectCharge = (1.0 - fChargeConfusionProbability) * nElectrons;
  double nElectronsWrongCharge   =        fChargeConfusionProbability  * nElectrons;

  double nPositronsCorrectCharge = (1.0 - fChargeConfusionProbability) * nPositrons;
  double nPositronsWrongCharge   =        fChargeConfusionProbability  * nPositrons;

  TH2* histogram = dynamic_cast<TH2*>(fDataDistribution->Clone(Form("hTemplateFitResult2D_%d", contributingCategories)));
  histogram->Reset();
  histogram->SetTitle("2D template fit result");
  histogram->SetStats(0);

  for (int xBin = 1; xBin <= histogram->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= histogram->GetNbinsY(); ++yBin) {
      double prediction = 0.0;
      if (contributingCategories & int(ContributingCategories::Protons))
        prediction += nProtons * fProtonTemplate->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::CCProtons))
        prediction += nCCProtons * fCCProtonTemplate->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::Electrons))
        prediction += nElectronsCorrectCharge * fElectronTemplate->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::CCElectrons))
        prediction += nElectronsWrongCharge * fCCElectronTemplate->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::Positrons))
        prediction += nPositronsCorrectCharge * fPositronTemplate->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::CCPositrons))
        prediction += nPositronsWrongCharge * fCCPositronTemplate->GetBinContent(xBin, yBin);

      histogram->SetBinContent(xBin, yBin, prediction);
    }
  }

  return histogram;
}

TH2* TemplateFit2DFitFunction::CalculateTotalHistogram() const {

  int total = int(ContributingCategories::Electrons)
            | int(ContributingCategories::Positrons)
            | int(ContributingCategories::CCElectrons)
            | int(ContributingCategories::CCPositrons)
            | int(ContributingCategories::Protons)
            | int(ContributingCategories::CCProtons);

  return CalculateHistogram(total);
}

TH2* TemplateFit2DFitFunction::DrawZeroProbabilityEvents() const {

  TH2* histogram = dynamic_cast<TH2*>(fDataDistribution->Clone("hZeroProbability"));
  histogram->Reset();
  histogram->SetStats(0);

  for (int xBin = 1; xBin <= histogram->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= histogram->GetNbinsY(); ++yBin) {
      double pdfProton = fProtonTemplate->GetBinContent(xBin, yBin);
      double pdfCCProton = fCCProtonTemplate->GetBinContent(xBin, yBin);
      double pdfElectron = fElectronTemplate->GetBinContent(xBin, yBin);
      double pdfPositron = fPositronTemplate->GetBinContent(xBin, yBin);
      double pdfCCElectron = fCCElectronTemplate->GetBinContent(xBin, yBin);
      double pdfCCPositron = fCCPositronTemplate->GetBinContent(xBin, yBin);

      if (pdfProton <= 0.0 && pdfCCProton <= 0.0 && pdfElectron <= 0.0 && pdfPositron <= 0.0 && pdfCCElectron <= 0.0 && pdfCCPositron <= 0.0)
        histogram->SetBinContent(xBin, yBin, fDataDistribution->GetBinContent(xBin, yBin));
    }
  }

  histogram->SetTitle(Form("%g events with zero probability", histogram->Integral()));
  return histogram;
}

void TemplateFit2DFitFunction::UpdateParameters(const double* parameters) {

  fProtonFraction = parameters[0];
  fCCProtonFraction = parameters[1];
  fElectronFraction = parameters[2];
  fPositronFraction = parameters[3];
  fChargeConfusionProbability = parameters[4];
}
