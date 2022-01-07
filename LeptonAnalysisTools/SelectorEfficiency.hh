#ifndef SelectorEfficiency_hh
#define SelectorEfficiency_hh

#include <string>
#include <vector>
#include "EfficiencyTools.hh"

namespace Analysis {
  class Event;
}

namespace Cuts {
  class Selector;
}

namespace Utilities {
  class ObjectManager;
}

class CutEfficiency;
class LeptonAnalysisTree;

class SelectorEfficiency {
public:
  SelectorEfficiency(Cuts::Selector*, unsigned int flags);

  void Initialize(Utilities::ObjectManager&, const char* efficiencyCategory);

  void AddManualCut(const std::string& cutDescription, const std::vector<double>& cutValueBinning);
  void UpdateManualCutStatus(unsigned int cutIndex, float cutValue, bool cutCondition, bool tagAndProbeConditionElectron, bool tagAndProbeConditionProton, bool controlPlotConditionElectron, bool controlPlotConditionProton);

  void UpdateCutStatus(unsigned int cutIndex, float cutValue, bool tagAndProbeConditionElectron, bool tagAndProbeConditionProton, bool controlPlotConditionElectron, bool controlPlotConditionProton);
  bool PassesSelector(const Analysis::Event&, unsigned short& firstFailedCut);
  bool EvaluateSelector(const Analysis::Event&);

  void RecordPassedSelector(const LeptonAnalysisTree&);
  void RecordFailedSelector(const LeptonAnalysisTree&, unsigned int failedCutIndex);

  Utilities::ObjectManager& ObjectManager() const { return *fObjectManager; }
  Cuts::Selector* Selector() const { return fSelector; }
  bool IsMC() const { return fFlags & IsMonteCarloData; }
  bool UseEcalEnergy() const { return fFlags & UseEcalEnergyScale; }
  bool ShouldRecordTagAndProbeEfficiencies() const { return fFlags & RecordTagAndProbeEfficiencies; }

private:
  Utilities::ObjectManager* fObjectManager;
  std::string fDirectoryPrefix;
  Cuts::Selector* fSelector;
  unsigned int fFlags;
  std::vector<CutEfficiency*> fCutEfficiencies;
  std::vector<float> fLastCutValues;
  std::vector<bool> fPassesElectronTagAndProbeCuts;
  std::vector<bool> fPassesProtonTagAndProbeCuts;
  std::vector<bool> fShouldFillElectronControlPlots;
  std::vector<bool> fShouldFillProtonControlPlots;
};

#endif // SelectorEfficiency_hh
