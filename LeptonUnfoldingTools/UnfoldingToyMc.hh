#ifndef UnfoldingToyMc_hh
#define UnfoldingToyMc_hh

#include <BinningDefinition.hh>

namespace Modelling {
class ResultTables;
}

class MigrationRandomGenerator;
class MigrationParameterization;
class TRandom;
class TGraphErrors;
class TH1;
class TH2;

class UnfoldingToyMc {

public:
  UnfoldingToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, const std::string& species_symbol, int randomSeed);
  virtual ~UnfoldingToyMc();

  virtual void SetVerbosity(int val) { fVerbosity = val; }

  TH1* MakeParameterizedFlux() const;
  TH2* GenerateToyMigrationMatrix(unsigned long long nToyResponse);
  Modelling::ResultTables* SimulateMeasurement(const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance, double fractionOfMeasuringTime,
                                               MigrationRandomGenerator* overrideRandomGenerator = nullptr);

  TH1* ToyAverageEventCounts() { return hToyAverageEventCounts; }
  TH1* ToyTruthEventCounts() { return hToyTruthEventCounts; }

  TH1* ToyAverageFlux() { return hToyAverageFlux; }
  TH1* ToyTruthFlux() { return hToyTruthFlux; }

  Binning::Definition EnergyBinning() const { return fEnergyBinning; }

  TGraphErrors* MakeFluxRatioOverAverageFlux(const TH1* hTestFlux) const;
  TGraphErrors* MakeTruthOverAverageRatio() const;

protected:
  virtual void DeleteHistograms();

protected:

  TRandom* fRandom;
  Binning::Definition fEnergyBinning;
  const MigrationParameterization* fParam = nullptr;
  MigrationRandomGenerator* fRandomGenerator = nullptr;
  int fVerbosity = 1;
  std::string fSpeciesSymbol; // "e-" or "e+"
  std::string fSpeciesName; // "Electron" or "Positron"

  //
  // results, filled by SimulateMeasurement()
  //
  /** Average expected event counts for given measurement time, i.e. no smearing applied, before cutoff, as a function of true energy. */
  TH1* hToyAverageEventCounts = nullptr;
  /** Simulated event counts, before cutoff, as a function of true energy. */
  TH1* hToyTruthEventCounts = nullptr;

  TH1* hToyAverageFlux = nullptr;
  TH1* hToyTruthFlux = nullptr;
};

#endif
