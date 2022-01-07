#include "EfficiencyTemplateFits.hh"

#include "AnalysisEvent.hh"
#include "BinningTools.hh"
#include "Cut.hh"
#include "CutAttachmentRegistry.hh"
#include "EfficiencyTools.hh"
#include "ObjectManager.hh"
#include "Selector.hh"
#include "StringTools.hh"
#include "TagSelector.hh"
#include "TemporaryChange.hh"
#include "Utilities.hh"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TLegend.h>

#define INFO_OUT_TAG "EfficiencyTemplateFits"
#include "debugging.hh"

namespace Cuts {

REGISTER_CUT_ATTACHMENT(Cuts::EfficiencyTemplateFits)

CutAttachment::ActivationMode EfficiencyTemplateFits::sActivationState = CutAttachment::ActivationMode::Selective;
bool EfficiencyTemplateFits::sUseWeights = true;

EfficiencyTemplateFits::EfficiencyTemplateFits()
  : CutAttachment(&sActivationState, CutAttachment::TypeSupport::Cuts, CutAttachment::CutValueInformation::Required)
  , fTagAndProbeHistos()
  , fAssertionCheckPerformedForPassedCategory(false) {

}

EfficiencyTemplateFits::~EfficiencyTemplateFits() {

  fTagAndProbeHistos.DeleteArray();
}

void EfficiencyTemplateFits::CreateObjects() {

  Utilities::TemporaryChange<bool> addDirectoryChange(TH1::AddDirectoryStatus, TH1::AddDirectory, false);

  if (AttachedToCut() && fCut->IsStandaloneCut()) {
    FATAL_OUT << "EfficiencyTemplateFits can't be attached to standalone cut: " << fCut->Description() << std::endl;
  }

  if (!fSelector->CommonXAxisInformation()) {
    FATAL_OUT << "SetupCommonXAxisInformation has not been called yet for selector: " << fSelector->GetName() << std::endl;
  }

  const Binning::Definition& xBinning = fSelector->CommonXAxisInformation()->AxisBinning();
  std::string xTitle = fSelector->CommonXAxisInformation()->AxisTitle();

  assert(fCut);
  if (!fCut->TagSelectors().GetEntriesFast())
    return;

  // Register histograms for each tag & probe category
  auto categories = fCut->GetTagSelectorCategories();
  assert(!categories.empty());
  assert(categories.size() * 2 == fTagAndProbeHistos.ActualSize());
  assert(fTagAndProbeHistos.ExpectedSize() == fTagAndProbeHistos.ActualSize());

  unsigned int tagAndProbeHistoCounter = 0;
  for (auto category : categories) {
    auto selectorsForCategory = fCut->GetTagSelectorsByCategory(category);
    assert(!selectorsForCategory.empty());

    std::stringstream titlePassed;
    titlePassed << "Passed events in category '" << category << "';" << xTitle << ";events";
    fTagAndProbeHistos.SetEntry(tagAndProbeHistoCounter, Make<TH2F>(Form("passedHistoTemplateFit_%s", category), titlePassed.str(), xBinning, fDistributionBinning));
    ++tagAndProbeHistoCounter;

    std::stringstream titleTotal;
    titleTotal << "Total events in category '" << category << "';" << xTitle << ";events";
    fTagAndProbeHistos.SetEntry(tagAndProbeHistoCounter, Make<TH2F>(Form("totalHistoTemplateFit_%s", category), titleTotal.str(), xBinning, fDistributionBinning));
    ++tagAndProbeHistoCounter;
  }
}

void EfficiencyTemplateFits::AddEmptyObjectsForNewTagCategory(unsigned int) {

  // Called by the cut for each new tag & probe category that should be analyzed.
  fTagAndProbeHistos.ReserveEntries(2);
}

void EfficiencyTemplateFits::AddObjectsToMerge(std::vector<TObject**>& list) {

  fTagAndProbeHistos.CollectListOfObjectsForMerging(list);
}

void EfficiencyTemplateFits::Fill(const Analysis::Event& event) {

  assert(fCut);
  assert(fSelector);

  // Update x-axis value
  double xValue = Cut::gUnsetValue;
  fSelector->CommonXAxisInformation()->GetAxisFunction()(event, xValue);

  double weight = (!sGloballyDisableWeights && sUseWeights) ? event.Weight() : 1.0;

  // Fill tag & probe histograms
  const auto& categories = fCut->GetTagSelectorCategories();
  if (categories.empty())
    return;

  // Enable Sumw2 (in case weights are used)
  if (weight != 1.0) {
    for (unsigned int i = 0; i < fTagAndProbeHistos.ActualSize(); ++i)
      Utilities::EnableWeightsForHistogram(fTagAndProbeHistos.GetEntry(i));
  }

  // Evaluate the tag & probe selectors.
  assert(fTagAndProbeHistos.ExpectedSize() == fTagAndProbeHistos.ActualSize());

  unsigned int categoryIndex = 0;
  for (const auto& category : categories) {
    auto selectorsForCategory = fCut->GetTagSelectorsByCategory(category);
    assert(!selectorsForCategory.empty());

    // Fill 'fTagAndProbeHistos'.
    bool passedCategory = true;
    for (auto tagSelector : selectorsForCategory) {
      if (std::find(fExcludedTagCuts.begin(), fExcludedTagCuts.end(), tagSelector->GetName()) == fExcludedTagCuts.end()) {
        if (!tagSelector->Passed()) {
          passedCategory = false;
          break;
        }
      }
    }

    if (passedCategory) {
      assert(fDistributionFunction);
      double distributionValue = fDistributionFunction();

      TH2F* passedTagAndProbeHisto = fTagAndProbeHistos.GetEntry(2 * categoryIndex);
      assert(passedTagAndProbeHisto);
      if (!fAssertionCheckPerformedForPassedCategory)
        assert(!strcmp(passedTagAndProbeHisto->GetName(), Form("passedHistoTemplateFit_%s", category)));

      TH2F* totalTagAndProbeHisto = fTagAndProbeHistos.GetEntry(2 * categoryIndex + 1);
      assert(totalTagAndProbeHisto);
      if (!fAssertionCheckPerformedForPassedCategory)
        assert(!strcmp(totalTagAndProbeHisto->GetName(), Form("totalHistoTemplateFit_%s", category)));

      fAssertionCheckPerformedForPassedCategory = true;

      totalTagAndProbeHisto->Fill(xValue, distributionValue, weight);
      if (fCut->Passed())
        passedTagAndProbeHisto->Fill(xValue, distributionValue, weight);
    }

    ++categoryIndex;
  }
}

void EfficiencyTemplateFits::SetDistributionFunctionAndBinning(DistributionFunction function, const Binning::Definition& binning) {

  assert(!fDistributionFunction);
  fDistributionFunction = function;
  fDistributionBinning = binning;
}

void EfficiencyTemplateFits::SetExcludedTagCuts(const std::set<std::string>& excludedTagCuts) {

  assert(fExcludedTagCuts.empty());
  fExcludedTagCuts = excludedTagCuts;
}

TH2F* EfficiencyTemplateFits::TotalHistogramForCategory(const std::string& category) const {

  unsigned int i = 0;
  std::string searchString = "passedHistoTemplateFit_" + category;
  for (; i < fTagAndProbeHistos.ActualSize(); i += 2) {
    std::string thisCategory = fTagAndProbeHistos.GetEntry(i)->GetName();
    if (thisCategory == searchString)
      break;
  }

  if (i >= fTagAndProbeHistos.ActualSize()) {
    WARN_OUT << "Could not find tag and probe histograms for category " << category << "!" << std::endl;
    return nullptr;
  }

  return fTagAndProbeHistos.GetEntry(i + 1);
}

TH2F* EfficiencyTemplateFits::PassedHistogramForCategory(const std::string& category) const {

  unsigned int i = 0;
  std::string searchString = "passedHistoTemplateFit_" + category;
  for (; i < fTagAndProbeHistos.ActualSize(); i += 2) {
    std::string thisCategory = fTagAndProbeHistos.GetEntry(i)->GetName();
    if (thisCategory == searchString)
      break;
  }

  if (i >= fTagAndProbeHistos.ActualSize()) {
    WARN_OUT << "Could not find tag and probe histograms for category " << category << "!" << std::endl;
    return nullptr;
  }

  return fTagAndProbeHistos.GetEntry(i);
}

} // namespace Cuts
