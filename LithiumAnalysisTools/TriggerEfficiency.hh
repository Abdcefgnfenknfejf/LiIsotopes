#ifndef TriggerEfficiency_hh
#define TriggerEfficiency_hh

#include "AnalysisSettings.hh"
#include "BinningDefinition.hh"
#include "BinningTools.hh"
#include "PredefinedBinnings.hh"

auto coarseTriggerBinningFunction = []() -> Binning::Definition {

  Binning::Definition lowEnergyBinning = Binning::SubRange(Binning::Predefined::AbsoluteEnergyBinning(), 0.0, 11.80, true);
  Binning::Definition highEnergyBinning = Binning::Logarithmic(1, 11.80, gAMSMaxEnergy);
  return Binning::CombinedBinning(lowEnergyBinning, highEnergyBinning);
};

#endif
