#include "SelectorEfficiency.hh"

#include "Cut.hh"
#include "CutEfficiency.hh"
#include "LeptonAnalysisTree.hh"
#include "Selector.hh"

SelectorEfficiency::SelectorEfficiency(Cuts::Selector* selector, unsigned int flags)
  : fObjectManager(0)
  , fSelector(selector)
  , fFlags(flags) {

  assert(fSelector);
}

void SelectorEfficiency::Initialize(Utilities::ObjectManager& objectManager, const char* directoryPrefix) {

  fObjectManager = &objectManager;
  fDirectoryPrefix = directoryPrefix;

  unsigned int cuts = fSelector->NumberOfCuts();
  if (!cuts)
    return;

  // If there are already cuts registered in the selector, create CutEfficiency objects for them.
  fCutEfficiencies.reserve(cuts);
  fLastCutValues.reserve(cuts);
  fPassesElectronTagAndProbeCuts.reserve(cuts);
  fPassesProtonTagAndProbeCuts.reserve(cuts);
  fShouldFillElectronControlPlots.reserve(cuts);
  fShouldFillProtonControlPlots.reserve(cuts);

  for (unsigned int i = 0; i < cuts; ++i) {
    CutEfficiency* efficiency = new CutEfficiency(this, i, directoryPrefix);
    fCutEfficiencies.push_back(efficiency);

    fLastCutValues.push_back(Cuts::Cut::gUnsetValue);
    fPassesElectronTagAndProbeCuts.push_back(false);
    fPassesProtonTagAndProbeCuts.push_back(false);
    fShouldFillElectronControlPlots.push_back(false);
    fShouldFillProtonControlPlots.push_back(false);
  }
}

void SelectorEfficiency::AddManualCut(const std::string& cutDescription, const std::vector<double>& cutValueBinning) {

  int cutIndex = fSelector->NumberOfCuts();
  fSelector->RegisterCut(new Cuts::ManualCut(cutDescription), &cutValueBinning);
  fCutEfficiencies.push_back(new CutEfficiency(this, cutIndex, fDirectoryPrefix.c_str()));

  fLastCutValues.push_back(Cuts::Cut::gUnsetValue);
  fPassesElectronTagAndProbeCuts.push_back(false);
  fPassesProtonTagAndProbeCuts.push_back(false);
  fShouldFillElectronControlPlots.push_back(false);
  fShouldFillProtonControlPlots.push_back(false);
}

void SelectorEfficiency::UpdateCutStatus(unsigned int cutIndex, float cutValue, bool tagAndProbeConditionElectron, bool tagAndProbeConditionProton, bool controlPlotConditionElectron, bool controlPlotConditionProton) {

  Cuts::Cut* cut = fSelector->GetCut(cutIndex);
  assert(cut);

  fLastCutValues.at(cutIndex) = cutValue;

  fPassesElectronTagAndProbeCuts.at(cutIndex) = tagAndProbeConditionElectron;
  fPassesProtonTagAndProbeCuts.at(cutIndex) = tagAndProbeConditionProton;

  fShouldFillElectronControlPlots.at(cutIndex) = controlPlotConditionElectron;
  fShouldFillProtonControlPlots.at(cutIndex) = controlPlotConditionProton;
}

void SelectorEfficiency::UpdateManualCutStatus(unsigned int cutIndex, float cutValue, bool cutCondition, bool tagAndProbeConditionElectron, bool tagAndProbeConditionProton, bool controlPlotConditionElectron, bool controlPlotConditionProton) {

  Cuts::ManualCut* manualCut = fSelector->GetManualCut(cutIndex);
  assert(manualCut);

  manualCut->Examine(cutCondition);
  manualCut->SetLastCutValue(cutValue);

  fLastCutValues.at(cutIndex) = cutValue;

  fPassesElectronTagAndProbeCuts.at(cutIndex) = tagAndProbeConditionElectron;
  fPassesProtonTagAndProbeCuts.at(cutIndex) = tagAndProbeConditionProton;

  fShouldFillElectronControlPlots.at(cutIndex) = controlPlotConditionElectron;
  fShouldFillProtonControlPlots.at(cutIndex) = controlPlotConditionProton;
}

bool SelectorEfficiency::PassesSelector(const Analysis::Event& event, unsigned short& firstFailedCut) {

  firstFailedCut = 0;
  if (fSelector->Passes(event))
    return true;

  const std::vector<bool>& lastCutStatus = fSelector->LastCutStatus();
  for (unsigned int cut = 0; cut < fSelector->NumberOfCuts(); ++cut) {
    if (!lastCutStatus.at(cut)) {
      firstFailedCut = cut + 1;
      break;
    }
  }

  assert(firstFailedCut >= 1);
  return false;
}

bool SelectorEfficiency::EvaluateSelector(const Analysis::Event& event) {

  return fSelector->Evaluate(event);
}

void SelectorEfficiency::RecordPassedSelector(const LeptonAnalysisTree& tree) {

  unsigned int cuts = fSelector->NumberOfCuts();
  for (unsigned int i = 0; i < cuts; ++i) {
    fCutEfficiencies.at(i)->ProcessEvent(true, fLastCutValues.at(i), tree,
                                         fPassesElectronTagAndProbeCuts.at(i), fPassesProtonTagAndProbeCuts.at(i),
                                         fShouldFillElectronControlPlots.at(i), fShouldFillProtonControlPlots.at(i));
  }
}

void SelectorEfficiency::RecordFailedSelector(const LeptonAnalysisTree& tree, unsigned int failedCutIndex) {

  for (unsigned int i = 0; i < failedCutIndex; ++i) {
    fCutEfficiencies.at(i)->ProcessEvent(true, fLastCutValues.at(i), tree,
                                         fPassesElectronTagAndProbeCuts.at(i), fPassesProtonTagAndProbeCuts.at(i),
                                         fShouldFillElectronControlPlots.at(i), fShouldFillProtonControlPlots.at(i));
  }

  fCutEfficiencies.at(failedCutIndex)->ProcessEvent(false, fLastCutValues.at(failedCutIndex), tree,
                                                    fPassesElectronTagAndProbeCuts.at(failedCutIndex), fPassesProtonTagAndProbeCuts.at(failedCutIndex),
                                                    fShouldFillElectronControlPlots.at(failedCutIndex), fShouldFillProtonControlPlots.at(failedCutIndex));
}
