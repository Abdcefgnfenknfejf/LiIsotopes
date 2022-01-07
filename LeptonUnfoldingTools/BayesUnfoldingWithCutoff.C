#include "BayesUnfoldingWithCutoff.hh"

#include <TH1.h>
#include <TH2.h>

#include <Utilities.hh>
#include <BinningTools.hh>

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

#define INFO_OUT_TAG "BayesUnfoldingWithCutoff"
#include <debugging.hh>

BayesUnfoldingWithCutoff::BayesUnfoldingWithCutoff(const TH1* hEventCounts, const TH1* hMeasuringTime, const TH1* hAcceptance, const TH1* hTriggerEfficiency,
                                                   const TH2* hMigrationMatrix, int nBayesIter, int verbosity)
  : fMigrationMatrixAfterCutoff(static_cast<TH2*>(hMigrationMatrix->Clone("hMigrationMatrixAfterCutoff")))
  , fMigrationMatrix(hMigrationMatrix)
  , fAcceptance(hAcceptance)
  , fEventCounts(hEventCounts)
  , fMeasuringTime(hMeasuringTime)
  , fTriggerEfficiency(hTriggerEfficiency) {

  fMaximumMeasuringTime = hMeasuringTime->GetBinContent(hMeasuringTime->GetNbinsX());

  TH1* hMcTrue = fMigrationMatrix->ProjectionY();
  TH1* hMcMeas = fMigrationMatrix->ProjectionX();

  // check binning consistency
  const auto binning = Binning::Tools::FromTAxis(hEventCounts->GetXaxis());
  if (binning != Binning::Tools::FromTAxis(hMeasuringTime->GetXaxis()) ||
      binning != Binning::Tools::FromTAxis(hMigrationMatrix->GetXaxis()))
    FATAL_OUT << "Binning mismatch!" << std::endl;

  if (hAcceptance && binning != Binning::Tools::FromTAxis(hAcceptance->GetXaxis()))
    FATAL_OUT << "Binning mismatch!" << std::endl;

  if (hTriggerEfficiency && binning != Binning::Tools::FromTAxis(hTriggerEfficiency->GetXaxis()))
    FATAL_OUT << "Binning mismatch!" << std::endl;

  // dump input?
  if (verbosity < 0) {
    int w = 12;
    std::cout << std::setw(16) << "bin" << " "
              << std::setw(w) << "counts" << " +- "
              << std::setw(w) << "error" << " "
              << std::setw(w) << "meas_time" << " "
              << std::setw(w) << "acceptance" << " "
              << std::setw(w) << "trig_eff" << " "
              << std::setw(w) << "migr_x" << " "
              << std::setw(w) << "migr_y"
              << std::endl;

    for (unsigned int bin = 1; bin <= binning.NumberOfBins(); ++bin) {
      std::cout << std::setw(16) << binning.BinAsString(bin) << " "
                << std::setw(w) << hEventCounts->GetBinContent(bin) << " +- "
                << std::setw(w) << hEventCounts->GetBinError(bin) << " "
                << std::setw(w) << hMeasuringTime->GetBinContent(bin) << " "
                << std::setw(w) << (hAcceptance ? hAcceptance->GetBinContent(bin) : -1.0) << " "
                << std::setw(w) << (hTriggerEfficiency ? hTriggerEfficiency->GetBinContent(bin) : - 1.0) << " "
                << std::setw(w) << hMcMeas->GetBinContent(bin) << " "
                << std::setw(w) << hMcTrue->GetBinContent(bin)
                << std::endl;
    }
  }

  // prepare migration matrix and its projection for unfolding, taking geomagnetic cutoff as encoded in measurement time
  // properly into account
  hMcMeas->Multiply(fMeasuringTime);
  hMcMeas->Scale(1.0 / fMaximumMeasuringTime);

  for (int xBin = 0; xBin <= fMigrationMatrixAfterCutoff->GetNbinsX() + 1; ++xBin) {
    for (int yBin = 0; yBin <= fMigrationMatrixAfterCutoff->GetNbinsY() + 1; ++yBin) {
      double binContent = fMigrationMatrixAfterCutoff->GetBinContent(xBin, yBin);
      double binError = fMigrationMatrixAfterCutoff->GetBinError(xBin, yBin);
      double cutoffProb = fMeasuringTime->GetBinContent(xBin) / fMaximumMeasuringTime;

      fMigrationMatrixAfterCutoff->SetBinContent(xBin, yBin, cutoffProb * binContent);
      fMigrationMatrixAfterCutoff->SetBinError(xBin, yBin, cutoffProb * binError);
    }
  }

  fResponse = new RooUnfoldResponse(hMcMeas, hMcTrue, fMigrationMatrixAfterCutoff);
  fResponse->UseOverflow(true);

  fUnfold = new RooUnfoldBayes(fResponse, fEventCounts, nBayesIter, false);
  fUnfold->SetVerbose(std::abs(verbosity));

  delete hMcTrue;
  delete hMcMeas;
}

BayesUnfoldingWithCutoff::~BayesUnfoldingWithCutoff() {

  delete fResponse;
  delete fUnfold;
  delete fMigrationMatrixAfterCutoff;
}

TH1D* BayesUnfoldingWithCutoff::UnfoldedEventCounts() const {

  TH1D* hUnfoldedEventCounts = static_cast<TH1D*>(fUnfold->Hreco()->Clone("hUnfoldedEventCounts"));
  hUnfoldedEventCounts->Multiply(fMeasuringTime);
  hUnfoldedEventCounts->Scale(1.0 / fMaximumMeasuringTime);
  return hUnfoldedEventCounts;
}

TH1D* BayesUnfoldingWithCutoff::MakeUnfoldedFlux() const {

  TH1D* hUnfoldedFlux = static_cast<TH1D*>(fUnfold->Hreco()->Clone("hUnfoldedFlux"));
  hUnfoldedFlux->GetYaxis()->SetTitle("unfolded flux");
  hUnfoldedFlux->Divide(fAcceptance);
  hUnfoldedFlux->Divide(fTriggerEfficiency);
  hUnfoldedFlux->Scale(1.0 / fMaximumMeasuringTime);
  hUnfoldedFlux->Scale(1.0, "width");

  return hUnfoldedFlux;
}

TH2D* BayesUnfoldingWithCutoff::MakeCountsCovarianceMatrixHistogram() const {

  TH2D* hCov = Make<TH2D>("hCov", ";col;row", Binning::FromTAxis(fUnfold->Hreco()->GetXaxis()), Binning::FromTAxis(fUnfold->Hreco()->GetXaxis()));

  auto Ereco = fUnfold->Ereco();
  if (Ereco.GetNcols() != hCov->GetNbinsX() || Ereco.GetNrows() != hCov->GetNbinsY())
    FATAL_OUT << "Binning mismatch!" << std::endl;

  for (int xBin = 1; xBin <= hCov->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= hCov->GetNbinsY(); ++yBin) {
      hCov->SetBinContent(xBin, yBin, Ereco(yBin - 1, xBin - 1));
    }
  }

  return hCov;
}

TH1D* BayesUnfoldingWithCutoff::MakeFoldedUnfoldedFlux() const {

  // convolve unfolded flux with migration matrix
  // this may be compared to the raw (input) flux as a closure test

  TH1D* hMigratedFlux = static_cast<TH1D*>(fUnfold->Hreco()->Clone("hMigratedFlux"));
  hMigratedFlux->Reset();
  hMigratedFlux->SetTitle("");
  hMigratedFlux->GetYaxis()->SetTitle("Flux (m^{-2} s^{-1} sr^{-1} GeV^{-1})");

  TH2* hNormalizedMigrationMatrix = static_cast<TH2*>(fMigrationMatrix->Clone("hNormalizedMigrationMatrix"));
  Utilities::NormalizeHistogramYSlices(hNormalizedMigrationMatrix, "", true);

  for (int xBin = 1; xBin <= hMigratedFlux->GetNbinsX(); ++xBin) {
    double migratedCounts = 0.0;
    for (int i = 1; i <= hMigratedFlux->GetNbinsX(); ++i)
      migratedCounts += fUnfold->Hreco()->GetBinContent(i) * hNormalizedMigrationMatrix->GetBinContent(xBin, i);

    hMigratedFlux->SetBinContent(xBin, migratedCounts);
  }

  hMigratedFlux->Scale(1.0 / fMaximumMeasuringTime);
  hMigratedFlux->Divide(fTriggerEfficiency);
  hMigratedFlux->Divide(fAcceptance);
  hMigratedFlux->Scale(1.0, "width");

  return hMigratedFlux;
}
