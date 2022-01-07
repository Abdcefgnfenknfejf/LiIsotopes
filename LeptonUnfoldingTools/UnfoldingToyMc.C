#include "UnfoldingToyMc.hh"

#include "MigrationParameterization.hh"
#include "MigrationRandomGenerator.hh"

#include <BinningTools.hh>
#include <BrokenPowerLawModel.hh>
#include <ProgressBar.hh>
#include <ResultTables.hh>
#include <Statistics.hh>
#include <Utilities.hh>

#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#define INFO_OUT_TAG "UnfoldingToyMc"
#include <debugging.hh>

UnfoldingToyMc::UnfoldingToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, const std::string& species_symbol, int randomSeed)
  : fEnergyBinning(energyBinning)
  , fParam(param)
  , fSpeciesSymbol(species_symbol)
{
  if (std::string(fSpeciesSymbol) == "e-")
    fSpeciesName = "Electron";
  else if (std::string(fSpeciesSymbol) == "e+")
    fSpeciesName = "Positron";
  else
    FATAL_OUT_ON_MASTER << "Invalid species symbol \"" << fSpeciesSymbol << "\"" << std::endl;

  fRandom = new TRandom3(randomSeed);
  fRandomGenerator = new MigrationRandomGenerator(fEnergyBinning, fParam);
  fRandomGenerator->SetRandomGenerator(fRandom);
}


UnfoldingToyMc::~UnfoldingToyMc() {

  DeleteHistograms();

  delete fRandom;
  delete fRandomGenerator;
}

TH1* UnfoldingToyMc::MakeParameterizedFlux() const {

  // parameterization done with parameterize_flux.C on ResultTables_LeptonFluxes_Zimmermann_05_07_2018_2D.root

  Modelling::BrokenPowerLawModel* model = nullptr;
  if (std::string(fSpeciesSymbol) == "e-") {
    model = new Modelling::BrokenPowerLawModel(ParticleId::Mass(ParticleId::Electron), 2.66,
    { {4.9, 3.54, 0.25}, {32.0, 3.22, 0.3} });
    model->C.SetValue(288.6);
    model->phi.SetValue(0.83);
  }
  else if (std::string(fSpeciesSymbol) == "e+") {
    model = new Modelling::BrokenPowerLawModel(ParticleId::Mass(ParticleId::Positron), 2.38,
    { {3.55, 3.15, 0.25}, {27.0, 2.74, 0.3}, {500.0, 3.74, 0.4} } );
    model->C.SetValue(7.7);
    model->phi.SetValue(0.33);
  }

  TH1* h = Utilities::ConvertToHistogram(model->Flux, fEnergyBinning);
  delete model;
  return h;
}


TH2* UnfoldingToyMc::GenerateToyMigrationMatrix(unsigned long long nToyResponse) {

  TH1* paramFlux = MakeParameterizedFlux();

  double xlow = fEnergyBinning.Min();
  double xup = fEnergyBinning.Max();

  TH2D* toyMigrationMatrix = Make<TH2D>("hToyMigrationMatrix", ";E_{meas};E_{truth}", fEnergyBinning, fEnergyBinning);

  Utilities::ProgressBar pb(nToyResponse);
  if (fVerbosity > 0) {
    INFO_OUT << "Fill toy migration matrix..." << std::endl;
    pb.PrintScale();
  }
  for (unsigned long long int i = 0; i < nToyResponse; ++i){
    if (fVerbosity > 0)
      pb.Update(i);

    double E = Utilities::PowerLawRandomNumber(xlow, xup, -1.0, gRandom);
    double Emeas = fRandomGenerator->GetRandomSmearedEnergy(E);

    auto ibin = fEnergyBinning.FindBin(E);
    double w = paramFlux->GetBinContent(ibin) * fEnergyBinning.BinCenterLog(ibin); // = paramFlux / E^(-1)

    toyMigrationMatrix->Fill(Emeas, E, w);
  }

  delete paramFlux;
  return toyMigrationMatrix;
}

Modelling::ResultTables* UnfoldingToyMc::SimulateMeasurement(const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance,
                                                             double fractionOfMeasuringTime, MigrationRandomGenerator* overrideRandomGenerator) {

  // clean-up from previous run (if any)
  DeleteHistograms();

  MigrationRandomGenerator* rg = overrideRandomGenerator ? overrideRandomGenerator : fRandomGenerator;

  const double w = 1.0;

  TH1* hMeasuringTimeForToymc = (TH1*)hMeasuringTime->Clone("hMeasuringTimeForToymc");
  hMeasuringTimeForToymc->Scale(fractionOfMeasuringTime);

  double meastime_max = hMeasuringTimeForToymc->GetBinContent(hMeasuringTimeForToymc->GetNbinsX());

  hToyAverageFlux = MakeParameterizedFlux();
  hToyAverageFlux->SetName("hToyAverageFlux");
  hToyAverageEventCounts = MakeParameterizedFlux();
  hToyAverageEventCounts->SetName("hToyAverageEventCounts");
  hToyAverageEventCounts->SetTitle("(before cutoff)");
  hToyAverageEventCounts->GetXaxis()->SetTitle("E_{truth} (GeV)");
  hToyAverageEventCounts->GetYaxis()->SetTitle("toy MC: true mean counts");

  hToyAverageEventCounts->Multiply(hAcceptance);
  hToyAverageEventCounts->Multiply(hTriggerEfficiency);
  hToyAverageEventCounts->Scale(meastime_max);
  Utilities::ScaleByBinWidth(*hToyAverageEventCounts);

  // set poisson errors
  for (int ibin = 1; ibin <= hToyAverageEventCounts->GetNbinsX(); ++ibin)
    hToyAverageEventCounts->SetBinError(ibin, TMath::Sqrt(hToyAverageEventCounts->GetBinContent(ibin)));

  // smear average counts to arrive at counts as a function of true energy
  // and set poisson errors
  hToyTruthEventCounts = Make<TH1D>("hToyTruthEventCounts", "(before cutoff);E_{truth} (GeV);toy MC: measured event counts", fEnergyBinning);
  for (int ibin = 1; ibin <= hToyTruthEventCounts->GetNbinsX(); ++ibin) {
    hToyTruthEventCounts->SetBinContent(ibin, fRandom->Poisson(hToyAverageEventCounts->GetBinContent(ibin)));
    hToyTruthEventCounts->SetBinError(ibin, TMath::Sqrt(hToyTruthEventCounts->GetBinContent(ibin)));
  }

  unsigned long long int nToySim = (unsigned long long int)(hToyTruthEventCounts->Integral());
  unsigned long long int iTotal = 0;

  Utilities::ProgressBar pb2(nToySim);
  if (fVerbosity > 0) {
    INFO_OUT << "Toy MC: simulate measurement, using " << nToySim << " events..." << std::endl;
    pb2.PrintScale();
  }

  TH1D* hToyMeasEventCounts = Make<TH1D>("hToyMeasEventCounts", ";E_{meas} (GeV);toy MC: measured event counts", fEnergyBinning);

  for (int ebin = 1; ebin <= hToyAverageEventCounts->GetNbinsX(); ++ebin){

    int nEventsInBin = TMath::Nint(hToyAverageEventCounts->GetBinContent(ebin));
    double E = hToyAverageEventCounts->GetXaxis()->GetBinCenterLog(ebin);

    for (int iEvent = 0; iEvent < nEventsInBin; ++iEvent) {

      if (fVerbosity > 0)
        pb2.Update(iTotal++);

      double Emeas = rg->GetRandomSmearedEnergy(E);

      double cutoffProb = hMeasuringTimeForToymc->GetBinContent(hMeasuringTimeForToymc->GetXaxis()->FindFixBin(Emeas)) / meastime_max;
      if (fRandom->Rndm() < cutoffProb)
        hToyMeasEventCounts->Fill(Emeas, w);
    }
  }

  // calculate fluxes, for easier comparison
  hToyTruthFlux = (TH1D*)hToyTruthEventCounts->Clone("hToyTruthFlux");
  hToyTruthFlux->Divide(hAcceptance);
  hToyTruthFlux->Divide(hTriggerEfficiency);
  hToyTruthFlux->Scale(1.0 / meastime_max);
  hToyTruthFlux->Scale(1.0, "width");

  TH1D* hToyMeasFlux = (TH1D*)hToyMeasEventCounts->Clone("hToyMeasFlux");
  hToyMeasFlux->Divide(hAcceptance);
  hToyMeasFlux->Divide(hTriggerEfficiency);
  hToyMeasFlux->Divide(hMeasuringTimeForToymc);
  hToyMeasFlux->Scale(1.0, "width");

  std::string fluxName = fSpeciesName + "Flux";
  std::string statErrName = fSpeciesName + "StatErr";
  std::string nevtName = fSpeciesName + "NEvents";
  StringVector tables{ fluxName, statErrName, nevtName, "MeasurementTime", "TriggerEfficiency", "FinalAcceptance" };
  Binning::Definition dummyTimeBinning = Binning::Tools::Equidistant(1, 1305800000.0, 1510527533.0);
  Modelling::ResultTables* resultTables = new Modelling::ResultTables("ToyMC", fEnergyBinning, dummyTimeBinning, tables);

  // set units
  resultTables->SetKinematicUnit("GeV");
  resultTables->SetUnit(fluxName, "1/m2/sr/s/GeV");
  resultTables->SetUnit(statErrName, "1/m2/sr/s/GeV");
  resultTables->SetUnit("MeasurementTime", "s");
  resultTables->SetUnit("FinalAcceptance", "1/m2/sr");

  resultTables->FillTableRowFromHisto(fluxName, 1, hToyMeasFlux);
  resultTables->FillTableRowFromHistoErrors(statErrName, 1, hToyMeasFlux);
  resultTables->FillTableRowFromHisto(nevtName, 1, hToyMeasEventCounts);
  resultTables->FillTableRowFromHisto("TriggerEfficiency", 1, hTriggerEfficiency);
  resultTables->FillTableRowFromHisto("MeasurementTime", 1, hMeasuringTimeForToymc); // use the scaled measurement time!
  resultTables->FillTableRowFromHisto("FinalAcceptance", 1, hAcceptance);

  delete hToyMeasFlux;
  delete hToyMeasEventCounts;
  delete hMeasuringTimeForToymc;
  return resultTables;
}

TGraphErrors* UnfoldingToyMc::MakeFluxRatioOverAverageFlux(const TH1* hTestFlux) const {

  if (!Binning::Tools::FromTAxis(hTestFlux->GetXaxis()).IsSubBinningOf(fEnergyBinning))
    FATAL_OUT << "Binning mismatch!" << std::endl;

  TGraphErrors* grRatio = new TGraphErrors;
  for (int ibinTest = 1; ibinTest <= hTestFlux->GetNbinsX(); ++ibinTest) {

    double E = hTestFlux->GetXaxis()->GetBinCenterLog(ibinTest);
    int ibinAvg = hToyAverageFlux->GetXaxis()->FindFixBin(E);
    double avgFlux = hToyAverageFlux->GetBinContent(ibinAvg);
    double testFlux = hTestFlux->GetBinContent(ibinTest);
    double testFluxErr = hTestFlux->GetBinError(ibinTest);
    double ratio = testFlux / avgFlux - 1.0;
    grRatio->SetPoint(ibinTest-1, E, ratio);
    grRatio->SetPointError(ibinTest-1, 0.0, testFluxErr / avgFlux);
  }

  return grRatio;
}


TGraphErrors* UnfoldingToyMc::MakeTruthOverAverageRatio() const {

  if (!hToyTruthFlux) {
    WARN_OUT << "Received nullptr for truth flux!" << std::endl;
    return nullptr;
  }

  TGraphErrors* grToyTruthAverageRatio = new TGraphErrors;
  grToyTruthAverageRatio->SetName("grToyTruthAverageRatio");
  for (int ibin = 1; ibin <= hToyAverageFlux->GetNbinsX(); ++ibin) {

    double E = hToyAverageFlux->GetXaxis()->GetBinCenterLog(ibin);
    double averageFlux = hToyAverageFlux->GetBinContent(ibin);
    double toyTruthAverageRatio = hToyTruthFlux->GetBinContent(ibin) / averageFlux - 1.0;
    grToyTruthAverageRatio->SetPoint(ibin-1, E, toyTruthAverageRatio);
    grToyTruthAverageRatio->SetPointError(ibin-1, 0.0, hToyTruthFlux->GetBinError(ibin) / averageFlux);
  }

  return grToyTruthAverageRatio;

}

void UnfoldingToyMc::DeleteHistograms() {

  delete hToyAverageEventCounts;
  delete hToyTruthEventCounts;
  delete hToyAverageFlux;
  delete hToyTruthFlux;
}

