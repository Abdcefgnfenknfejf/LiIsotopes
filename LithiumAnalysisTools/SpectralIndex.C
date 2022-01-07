#include "SpectralIndex.hh"

#include "AnalysisSettings.hh"
#include "ModulatedPowerLawModel.hh"
#include "SimplePowerLawModel.hh"

#include <HistogramDataset.hh>
#include <ModellingData.hh>
#include <FitFunction.hh>
#include <ModelAnalysis.hh>
#include <FluxTools.hh>

#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <cassert>
#include <cmath>

#define INFO_OUT_TAG "SpectralIndex"
#include "debugging.hh"

void AddIndividualSpectralIndexFitResult(Modelling::HistogramDataset* ds, int iBinMin, int iBinMax, TGraphErrors* grResult) {

  const TH1* dataHisto = ds->GetHistogram();

  double Emin = dataHisto->GetXaxis()->GetBinLowEdge(iBinMin);
  double Emax = dataHisto->GetXaxis()->GetBinUpEdge(iBinMax);
  double Ecorr = std::pow(10.0, 0.5*(std::log10(Emax) + std::log10(Emin)));

  DEBUG_OUT << "iMin= " << iBinMin << " iMax= " << iBinMax
            << " Emin= " << Emin << " Emax= " << Emax << " Ecorr= " << Ecorr
            << std::endl;

  double startValueAmp = dataHisto->GetBinContent(dataHisto->FindFixBin(Ecorr));
  double logF1 = std::log(dataHisto->GetBinContent(iBinMin));
  double logF2 = std::log(dataHisto->GetBinContent(iBinMax));
  double logE1 = std::log(dataHisto->GetXaxis()->GetBinCenterLog(iBinMin));
  double logE2 = std::log(dataHisto->GetXaxis()->GetBinCenterLog(iBinMax));
  double dlogF = logF2 - logF1;
  double dlogE = logE2 - logE1;
  double spectralIndexStart = -(dlogF / dlogE);
  DEBUG_OUT << "start values: amp " << startValueAmp << " index " << spectralIndexStart << std::endl;
  SimplePowerLawModel* plModel = new SimplePowerLawModel(Ecorr, startValueAmp, spectralIndexStart);
  ds->SetXRange(Emin,Emax);
  Modelling::Data dataPowerLawFit;
  dataPowerLawFit.AddDataset(0,ds);

  Modelling::FitFunction* fitfPL = new Modelling::FitFunction(&dataPowerLawFit,plModel);
  Modelling::ModelAnalysis analysisPL(fitfPL, -1);
  analysisPL.RunAnalysis();

  DEBUG_OUT << "Fit amp= " << plModel->C.Value() << " fit index: " << plModel->gamma.Value()
            << " fit chi2: " << fitfPL->Chi2() << "/" << fitfPL->Ndf() << "= " << fitfPL->Chi2()/fitfPL->Ndf()
            << std::endl;

  plModel->Flux->SetRange(Emin, Emax);

  int iP = grResult->GetN();
  grResult->SetPoint(iP, Ecorr, -plModel->gamma);
  grResult->SetPointError(iP, 0., plModel->gamma.Uncertainty());

  delete plModel;
  delete fitfPL;
}

TGraphErrors* CalculateSpectralIndexGraphUsingSlidingWindowFit(Modelling::HistogramDataset* ds, std::function<int(double)> nPointsEachSide, int graphColor) {

  const TH1* dataHisto = ds->GetHistogram();

  const int NPointsMin = 3;

  DEBUG_OUT << "Calculating spectral indices for " << dataHisto->GetTitle() << std::endl;
  TGraphErrors* grIndexFit = new TGraphErrors;
  grIndexFit->SetMarkerStyle(kFullCircle);
  grIndexFit->SetMarkerColor(graphColor);
  grIndexFit->SetLineColor(graphColor);
  grIndexFit->SetFillColor(graphColor);
  grIndexFit->SetLineWidth(2);
  grIndexFit->SetTitle(dataHisto->GetTitle());

  int iFit = 0;
  int iBin = 1;
  while (true) {

    iBin += 1;

    double Ebin = dataHisto->GetXaxis()->GetBinCenterLog(iBin);
    int NPointsEachSide = nPointsEachSide(Ebin);

    int iMin = iBin - NPointsEachSide;
    if(iMin < 1) iMin = 1;
    int iMax = iBin + NPointsEachSide;
    if(iMax>dataHisto->GetNbinsX()) iMax=dataHisto->GetNbinsX();

    int Npoints = iMax-iMin+1;
    bool atEnd = (iMax >= dataHisto->GetNbinsX());

    if (atEnd && (Npoints % 2 == 1) )
      continue;
    if (atEnd && (Npoints < NPointsMin))
      break;

    ++iFit;
    DEBUG_OUT << "iFit=" << iFit << " iBin= " << iBin << " Ebin= " << Ebin
              << " N= " << Npoints << std::endl;

    AddIndividualSpectralIndexFitResult(ds, iMin, iMax, grIndexFit);
  }

  return grIndexFit;
}

TGraphErrors* CalculateSpectralIndexGraphUsingSlidingWindowFit(Modelling::HistogramDataset* ds, const std::vector<double>& energyIntervals, int graphColor) {

  const TH1* dataHisto = ds->GetHistogram();

  DEBUG_OUT << "Calculating spectral indices for " << dataHisto->GetTitle() << std::endl;
  TGraphErrors* grIndexFit = new TGraphErrors;
  grIndexFit->SetMarkerStyle(kFullCircle);
  grIndexFit->SetMarkerColor(graphColor);
  grIndexFit->SetLineColor(graphColor);
  grIndexFit->SetFillColor(graphColor);
  grIndexFit->SetLineWidth(2);
  grIndexFit->SetTitle(dataHisto->GetTitle());

  for (unsigned int i = 0; i < energyIntervals.size() - 1; ++i) {
    auto iMin = dataHisto->GetXaxis()->FindFixBin(energyIntervals[i] + 0.01);
    auto iMax = dataHisto->GetXaxis()->FindFixBin(energyIntervals[i + 1] - 0.01);
    AddIndividualSpectralIndexFitResult(ds, iMin, iMax, grIndexFit);
  }

  return grIndexFit;
}

TGraph* CalculateSpectralIndexGraphUsingModulatedPowerLawFit(Modelling::HistogramDataset* ds, ParticleId::Species particleType, const Binning::Definition& binning, double startFitEnergy, double stopFitEnergy) {

  assert(particleType == ParticleId::Electron || particleType == ParticleId::Positron);
  const TH1* dataHisto = ds->GetHistogram();

  double yAt5 = dataHisto->GetBinContent(dataHisto->GetXaxis()->FindFixBin(5.0));

  double phiStart = particleType == ParticleId::Positron ? 1.0 : 0.5;
  auto* plModel = new ModulatedPowerLawModel(5.0, yAt5, 3.0, phiStart);
  ds->SetXRange(startFitEnergy, stopFitEnergy);

  Modelling::Data dataPowerLawFit;
  dataPowerLawFit.AddDataset(0, ds);

  Modelling::FitFunction* fitfPL = new Modelling::FitFunction(&dataPowerLawFit,plModel);
  Modelling::ModelAnalysis analysisPL(fitfPL, -1);
  if (!analysisPL.RunAnalysis()) {
    stopFitEnergy -= 1.0;

    bool converged = false;
    while (stopFitEnergy > startFitEnergy) {
      ds->SetXRange(startFitEnergy, stopFitEnergy);
      converged = analysisPL.RunAnalysis();
      if (converged)
        break;
      stopFitEnergy -= 1.0;
    }

    if (!converged)
      FATAL_OUT << "Fit did not converge. Aborting!" << std::endl;
  }

  TGraph* resultGraph = new TGraph; 
  auto* sd = new Modelling::SpectralDerivative(plModel->Flux, "index");
  for (unsigned int energyBin = 1; energyBin <= binning.NumberOfBins(); ++energyBin) {
    double xTrue = Modelling::LaffertyWyatt(binning.LowEdge(energyBin), binning.LowEdge(energyBin + 1), gAMSLaffertyWyattSpectralIndex);
    resultGraph->SetPoint(resultGraph->GetN(), xTrue, sd->SpectralIndex()->Eval(xTrue));
  }

  delete sd;
  delete plModel;
  return resultGraph;
}
