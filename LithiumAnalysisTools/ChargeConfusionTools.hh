#ifndef ChargeConfusionTools_hh
#define ChargeConfusionTools_hh

#include "FluxTools.hh"

#include <TH1D.h>
#include <TGraph.h>

#include <cassert>
#include <cmath>
#include <iostream>

void ApplyChargeConfusionCorrection(double& electrons, double& electronsError, double& electronsErrorRelativeSystematic,
                                    double& positrons, double& positronsError, double& positronsErrorRelativeSystematic,
                                    double& positronElectronRatioErrorRelativeSystematic, double& positronFractionErrorRelativeSystematic,
                                    double chargeConfusion, double chargeConfusionError) {

  if (electrons > 0.0 && positrons > 0.0) {
    double trueElectrons            = (electrons * (chargeConfusion - 1.0) + positrons * chargeConfusion) / (2.0 * chargeConfusion - 1.0);
    double dtrueElectrons_dpositive = chargeConfusion / (2.0 * chargeConfusion - 1.0);
    double dtrueElectrons_dnegative = (chargeConfusion - 1.0) / (2.0 * chargeConfusion - 1.0);
    double dtrueElectrons_dcc       = (electrons - positrons) / std::pow(2.0 * chargeConfusion - 1.0, 2);

    electronsErrorRelativeSystematic = std::abs(dtrueElectrons_dcc * chargeConfusionError) / trueElectrons;
    double trueElectronsErrorStatistical = std::sqrt(std::pow(dtrueElectrons_dpositive * positronsError, 2)  + std::pow(dtrueElectrons_dnegative * electronsError, 2));

    double truePositrons            = (electrons * chargeConfusion + positrons * (chargeConfusion - 1.0)) / (2.0 * chargeConfusion - 1.0);
    double dtruePositrons_dpositive = (chargeConfusion - 1.0) / (2.0 * chargeConfusion - 1.0);
    double dtruePositrons_dnegative = chargeConfusion / (2.0 * chargeConfusion - 1.0);
    double dtruePositrons_dcc       = -(electrons - positrons) / std::pow(2.0 * chargeConfusion - 1.0, 2);

    positronsErrorRelativeSystematic = std::abs(dtruePositrons_dcc * chargeConfusionError) / truePositrons;
    double truePositronsErrorStatistical = std::sqrt(std::pow(dtruePositrons_dpositive * positronsError, 2)  + std::pow(dtruePositrons_dnegative * electronsError, 2));

    double truePositronElectronRatio = truePositrons / trueElectrons;
    double dtruePositronElectronRatio_dcc = (std::pow(positrons, 2) - std::pow(electrons, 2)) / std::pow(electrons * (chargeConfusion - 1.0) + positrons * chargeConfusion, 2);
    positronElectronRatioErrorRelativeSystematic = std::abs(dtruePositronElectronRatio_dcc * chargeConfusionError) / truePositronElectronRatio;

    double truePositronFraction = truePositrons / (truePositrons + trueElectrons);
    double dtruePositronFraction_dcc = dtruePositrons_dcc / (electrons + positrons);
    positronFractionErrorRelativeSystematic = std::abs(dtruePositronFraction_dcc * chargeConfusionError) / truePositronFraction;

    electrons = trueElectrons;
    electronsError = trueElectronsErrorStatistical;
    positrons = truePositrons;
    positronsError = truePositronsErrorStatistical;
  }

  // If any count is negative or null, reset the result for both e+ / e- flux.
  if (electrons <= 0.0 || positrons <= 0.0) {
    electrons = 0.0;
    electronsError = 0.0;

    positrons = 0.0;
    positronsError = 0.0;

    electronsErrorRelativeSystematic = 1.0;
    positronsErrorRelativeSystematic = 1.0;
    positronElectronRatioErrorRelativeSystematic = 1.0;
    positronFractionErrorRelativeSystematic = 1.0;
  }

  // Clamp relative statistical uncertainties to 1.
  if (electronsError > electrons)
      electronsError = electrons;

  if (positronsError > positrons)
      positronsError = positrons;

  // Clamp relative systematic uncertainties to 1.
  if (electronsErrorRelativeSystematic >= 1.0)
    electronsErrorRelativeSystematic = 1.0;

  if (positronsErrorRelativeSystematic >= 1.0)
    positronsErrorRelativeSystematic = 1.0;

  if (positronElectronRatioErrorRelativeSystematic >= 1.0)
    positronElectronRatioErrorRelativeSystematic = 1.0;

  if (positronFractionErrorRelativeSystematic >= 1.0)
    positronFractionErrorRelativeSystematic = 1.0;
}

void ApplyChargeConfusionCorrectionForFluxes(TH1D* electronEventCounts, TH1D* positronEventCounts, TGraph* chargeConfusionGraph,
                                             TGraph* electronChargeConfusionRelSystUncertainty,
                                             TGraph* positronChargeConfusionRelSystUncertainty) {

  assert(electronEventCounts);
  assert(positronEventCounts);
  assert(chargeConfusionGraph);

  const auto& binning = Binning::Predefined::AbsoluteEnergyBinning();

  for (int bin = 1; bin <= electronEventCounts->GetNbinsX(); ++bin) {
    double xElectrons = electronEventCounts->GetBinLowEdge(bin);
    assert(std::abs(xElectrons - binning.LowEdge(bin)) < 1e-5);

    double electrons = electronEventCounts->GetBinContent(bin);
    double electronsError = electronEventCounts->GetBinError(bin);

    double xPositrons = positronEventCounts->GetBinLowEdge(bin);
    double positrons = positronEventCounts->GetBinContent(bin);
    double positronsError = positronEventCounts->GetBinError(bin);
    assert(std::abs(xElectrons - xPositrons) < 1e-5);

    double xTrue = Modelling::LaffertyWyatt(binning.LowEdge(bin), binning.LowEdge(bin + 1), gAMSLaffertyWyattSpectralIndex);
    double xChargeConfusion = chargeConfusionGraph->GetX()[bin - 1];
    assert(std::abs(xTrue - xChargeConfusion) < 1e-5);

    double chargeConfusion = chargeConfusionGraph->GetY()[bin - 1];
    double chargeConfusionError = chargeConfusionGraph->GetErrorY(bin - 1);

    double electronsErrorRelSystematic = 0.0;
    double positronsErrorRelSystematic = 0.0;
    double positronElectronRatioErrorRelSystematic = 0.0;
    double positronFractionErrorRelSystematic = 0.0;
    ApplyChargeConfusionCorrection(electrons, electronsError, electronsErrorRelSystematic,
                                   positrons, positronsError, positronsErrorRelSystematic,
                                   positronElectronRatioErrorRelSystematic, positronFractionErrorRelSystematic,
                                   chargeConfusion, chargeConfusionError);

    electronEventCounts->SetBinContent(bin, electrons);
    electronEventCounts->SetBinError(bin, electronsError);

    positronEventCounts->SetBinContent(bin, positrons);
    positronEventCounts->SetBinError(bin, positronsError);

    if (electronChargeConfusionRelSystUncertainty)
      electronChargeConfusionRelSystUncertainty->SetPoint(bin - 1, xTrue, electronsErrorRelSystematic);
    if (electronChargeConfusionRelSystUncertainty)
      positronChargeConfusionRelSystUncertainty->SetPoint(bin - 1, xTrue, positronsErrorRelSystematic);
  }
}

void ApplyChargeConfusionCorrectionForRatios(TH1D* electronEventCounts, TH1D* positronEventCounts, TGraph* chargeConfusionGraph,
                                             TGraph* positronElectronRatioChargeConfusionRelSystUncertainty,
                                             TGraph* positronFractionChargeConfusionRelSystUncertainty) {

  assert(electronEventCounts);
  assert(positronEventCounts);
  assert(chargeConfusionGraph);
  assert(positronElectronRatioChargeConfusionRelSystUncertainty);
  assert(positronFractionChargeConfusionRelSystUncertainty);

  const auto& binning = Binning::Predefined::AbsoluteEnergyBinning();

  for (int bin = 1; bin <= electronEventCounts->GetNbinsX(); ++bin) {
    double xElectrons = electronEventCounts->GetBinLowEdge(bin);
    assert(std::abs(xElectrons - binning.LowEdge(bin)) < 1e-5);

    double electrons = electronEventCounts->GetBinContent(bin);
    double electronsError = electronEventCounts->GetBinError(bin);

    double xPositrons = positronEventCounts->GetBinLowEdge(bin);
    double positrons = positronEventCounts->GetBinContent(bin);
    double positronsError = positronEventCounts->GetBinError(bin);
    assert(std::abs(xElectrons - xPositrons) < 1e-5);

    double xTrue = Modelling::LaffertyWyatt(binning.LowEdge(bin), binning.LowEdge(bin + 1), gAMSLaffertyWyattSpectralIndex);
    double xChargeConfusion = chargeConfusionGraph->GetX()[bin - 1];
    assert(std::abs(xTrue - xChargeConfusion) < 1e-5);

    double chargeConfusion = chargeConfusionGraph->GetY()[bin - 1];
    double chargeConfusionError = chargeConfusionGraph->GetErrorY(bin - 1);

    double electronsErrorRelSystematic = 0.0;
    double positronsErrorRelSystematic = 0.0;
    double positronElectronRatioErrorRelSystematic = 0.0;
    double positronFractionErrorRelSystematic = 0.0;
    ApplyChargeConfusionCorrection(electrons, electronsError, electronsErrorRelSystematic,
                                   positrons, positronsError, positronsErrorRelSystematic,
                                   positronElectronRatioErrorRelSystematic, positronFractionErrorRelSystematic,
                                   chargeConfusion, chargeConfusionError);

    assert(electrons >= 0.0);
    assert(electronsError >= 0.0);

    assert(positrons >= 0.0);
    assert(positronsError >= 0.0);

    if (electrons > 0.0)
      assert(electronsError / electrons <= 1.0);

    if (positrons > 0.0)
      assert(positronsError / positrons <= 1.0);

    assert(electronsErrorRelSystematic >= 0.0);
    assert(electronsErrorRelSystematic <= 1.0);

    assert(positronsErrorRelSystematic >= 0.0);
    assert(positronsErrorRelSystematic <= 1.0);

    electronEventCounts->SetBinContent(bin, electrons);
    electronEventCounts->SetBinError(bin, electronsError);

    positronEventCounts->SetBinContent(bin, positrons);
    positronEventCounts->SetBinError(bin, positronsError);

    positronElectronRatioChargeConfusionRelSystUncertainty->SetPoint(bin - 1, xTrue, positronElectronRatioErrorRelSystematic);
    positronFractionChargeConfusionRelSystUncertainty->SetPoint(bin - 1, xTrue, positronFractionErrorRelSystematic);
  }
}

#endif
