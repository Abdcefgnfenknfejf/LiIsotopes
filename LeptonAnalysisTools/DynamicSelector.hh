#ifndef DynamicSelector_hh
#define DynamicSelector_hh

#include <map>

#include "AxisInformation.hh"
#include "BinningFunctions.hh"
#include "EfficiencyTools.hh"
#include "ConfigHandler.hh"
#include "ObjectManager.hh"
#include "Selector.hh"

namespace IO {
  class TreeReader;
}

namespace Cuts {
  class Selector;
}

class TF1;

class LeptonAnalysisTree;

class DynamicSelector {
public:
  DynamicSelector(const std::string& configFileName, const std::string& selectorName);

  void Initialize(Utilities::ObjectManager&, Utilities::ObjectManager& cutsObjectManager, Cuts::AxisInformation::Expression, LeptonAnalysisTree*, IO::TreeReader*, const std::string& directoryInFile = "");
  void Initialize(Utilities::ObjectManager&, Utilities::ObjectManager& cutsObjectManager, Cuts::AxisInformation::Expression, LeptonAnalysisTree*, IO::TreeReader*, const Binning::Functions::BinningFunction& xAxisDataBinning, const std::string& directoryInFile = "");
  void InitializeFromFile(TFile* cutFile, const std::string& directoryInFile = "Selectors");
  bool Examine(const IO::TreeReader& treeReader);

  void PrintSummary() const;

  float GetXAxisRangeStart(const std::string& cutName) const;
  float GetYAxisRangeStart(const std::string& cutName) const;
  float GetYAxisRangeStop(const std::string& cutName) const;
  const std::string& GetYAxisTitle(const std::string& cutName) const;
  const TagAndProbeRatioFitSettings& GetTagAndProbeRatioFitSettings(const std::string& cutName) const;
  Cuts::Selector* Selector() const { return fSelector; }

  std::string SelectorName() const;

private:
  void ParseAdditionalTagAndProbeInformation(Cuts::Selector*);
  void ParseAdditionalCutConstantInformation(Cuts::Selector*);
  TagAndProbeRatioFitSettings ParseTagAndProbeRatioFitSettings(const std::string& cutSection);

  Utilities::ConfigHandler fConfig;
  Cuts::Selector* fSelector;
  std::string fSelectorName;
  std::string fConfigFileName;
  std::map<std::string, TagAndProbeRatioFitSettings> fTagAndProbeRatioFitSettings;
  std::map<std::string, std::string> fCutYAxisTitles;
  std::map<std::string, float> fCutXRangeStart;
  std::map<std::string, float> fCutYRangeStart;
  std::map<std::string, float> fCutYRangeStop;
};


#endif
