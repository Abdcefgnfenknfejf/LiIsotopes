#include "DynamicSelector.hh"

#include "AnalysisSettings.hh"
#include "Cut.hh"
#include "Environment.hh"
#include "FunctionWrapper.hh"
#include "NucleiAnalysisTree.hh"
#include "TreeReader.hh"
#include "SelectionParser.hh"
#include "StringTools.hh"
#include "PredefinedBinnings.hh"
#include "LithiumBinning.hh"

#include "EfficiencyTemplateFits.hh"
#define INFO_OUT_TAG "DynamicSelector"
#include "debugging.hh"

DynamicSelector::DynamicSelector(const std::string& configFileName, const std::string& selectorName)
  : fSelector(nullptr)
  , fSelectorName(selectorName) {

  fConfigFileName = configFileName;
  Environment::ExpandEnvironmentVariables(fConfigFileName);

  std::string tagAndProbeConfigFileName = "${MYANALYSIS}/Configuration/NucleiSelectionTag.cfg";
  Environment::ExpandEnvironmentVariables(tagAndProbeConfigFileName);
  fConfig.Read(tagAndProbeConfigFileName);

  fConfig.Read(fConfigFileName);
}


void DynamicSelector::InitializeFromFile(TFile* file, const std::string& directoryInFile) {

  assert(file);

  assert(!fSelector);
  if (directoryInFile.empty())
    fSelector = dynamic_cast<Cuts::Selector*>(file->Get(fSelectorName.c_str()));
  else
    fSelector = dynamic_cast<Cuts::Selector*>(file->Get(Form("%s/%s", directoryInFile.c_str(), fSelectorName.c_str())));

  if (!fSelector) {
    FATAL_OUT << "Can't load selector named '" << fSelectorName << "' from file. Aborting!" << std::endl;
  }

  ParseAdditionalTagAndProbeInformation(fSelector);
  ParseAdditionalCutConstantInformation(fSelector);
}

void DynamicSelector::Initialize(Utilities::ObjectManager& objectManager, Utilities::ObjectManager& cutsObjectManager, Cuts::AxisInformation::Expression , NucleiAnalysisTree* treeInterface, IO::TreeReader* treeReader, const Binning::Functions::BinningFunction& , const std::string& directoryInFile) {

  Cuts::SelectionParser selectionParser(fConfig, treeReader, &cutsObjectManager);
  if (!selectionParser.HasSelector(fSelectorName)) {
    FATAL_OUT << "Couldn't find selector '" << fSelectorName << "' in cut config file '" << fConfigFileName << "'. Aborting!" << std::endl;
  }

  fSelector = objectManager.Add(selectionParser.GetSelector(fSelectorName), directoryInFile);
  assert(fSelector->GetName() == fSelectorName);

  Binning::Predefined::SetDefaultAbsoluteRigidityBinning(LithiumRigidityBinning());
  auto xRigidityAxisValue = [treeInterface] (const Analysis::Event&, double& lastBaseValue) {
    if (treeInterface->Rigidity() != 0)
      lastBaseValue = std::abs(treeInterface->Rigidity());
  };

  auto xRigidityAxisValueMc = [treeInterface] (const Analysis::Event&, double& lastBaseValue) {
    if (treeInterface->McGeneratedMomentum() != 0)
      lastBaseValue = treeInterface->McGeneratedMomentum()/3.0;  //For Li Charge = 3.0                                                                                                                      
  };

  fSelector->SetupCommonXAxisInformation(xRigidityAxisValue, "|Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  fSelector->SetupCommonXAxisInformationMc(xRigidityAxisValueMc,  "|true Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  ParseAdditionalTagAndProbeInformation(fSelector);

}

void DynamicSelector::Initialize(Utilities::ObjectManager& objectManager, Utilities::ObjectManager& cutsObjectManager, Cuts::AxisInformation::Expression mode, NucleiAnalysisTree* treeInterface, IO::TreeReader* treeReader, const std::string& directoryInFile) {

  return Initialize(objectManager, cutsObjectManager, mode, treeInterface, treeReader, Binning::Predefined::AbsoluteRigidityBinning, directoryInFile);
}

void DynamicSelector::ParseAdditionalTagAndProbeInformation(Cuts::Selector* selector) {

  assert(selector);
  assert(!fSelector || (fSelector == selector));
  fSelector = selector;

  unsigned int cuts = fSelector->NumberOfCuts();
  for (unsigned int cutIndex = 0; cutIndex < cuts; ++cutIndex) {
    const Cuts::Cut* cut = fSelector->GetCut(cutIndex);
    const std::string& cutSection = cut->GetName();
    if (!fConfig.SectionExists(cutSection))
      continue;
    fTagAndProbeRatioFitSettings[cutSection] = ParseTagAndProbeRatioFitSettings(cutSection);
  }
}

void DynamicSelector::ParseAdditionalCutConstantInformation(Cuts::Selector* selector) {

  assert(selector);
  assert(fSelector == selector);

  unsigned int cuts = fSelector->NumberOfCuts();
  for (unsigned int cutIndex = 0; cutIndex < cuts; ++cutIndex) {
    const Cuts::Cut* cut = fSelector->GetCut(cutIndex);
    const std::string& cutSection = cut->GetName();
    if (!fConfig.SectionExists(cutSection))
      continue;

    std::string yAxisTitle = fConfig.GetStringValue(cutSection, "CutYAxisTitle");
    if (!yAxisTitle.empty()) {
      findAndReplace(yAxisTitle, "*", "#"); // # is not allowed in ConfigHandler values, but needed for ROOT + LaTeX
      fCutYAxisTitles[cutSection] = yAxisTitle;
    }

    const std::string& xAxisRangeStart = fConfig.GetStringValue(cutSection, "CutXAxisRangeStart");
    if (!xAxisRangeStart.empty())
      fCutXRangeStart[cutSection] = StringToNumber<float>(xAxisRangeStart);

    const std::string& yAxisRangeStart = fConfig.GetStringValue(cutSection, "CutYAxisRangeStart");
    if (!yAxisRangeStart.empty())
      fCutYRangeStart[cutSection] = StringToNumber<float>(yAxisRangeStart);

    const std::string& yAxisRangeStop = fConfig.GetStringValue(cutSection, "CutYAxisRangeStop");
    if (!yAxisRangeStop.empty())
      fCutYRangeStop[cutSection] = StringToNumber<float>(yAxisRangeStop);

    const std::string& lowerConstantString = fConfig.GetStringValue(cutSection, "CutLowerConstantForVisualization");
    if (!lowerConstantString.empty()) {
      auto* wrapper = cut->GetFunctions();
      std::string lowerCutFunctionName = Form("lowerCutFunction_%s", cutSection.c_str());
      if (wrapper->LowerCutValueFunction()) {
        assert(wrapper->LowerCutValueFunction()->GetName() == lowerCutFunctionName);
        assert(wrapper->LowerCutValueFunction()->GetTitle() == lowerConstantString);
      } else
        wrapper->SetLowerCutValueFunction(new TF1(lowerCutFunctionName.c_str(), lowerConstantString.c_str(), 0.0, 1.0e6));
    }

    const std::string& upperConstantString = fConfig.GetStringValue(cutSection, "CutUpperConstantForVisualization");
    if (!upperConstantString.empty()) {
      auto* wrapper = cut->GetFunctions();
      std::string upperCutFunctionName = Form("upperCutFunction_%s", cutSection.c_str());
      if (wrapper->UpperCutValueFunction()) {
        assert(wrapper->UpperCutValueFunction()->GetName() == upperCutFunctionName);
        assert(wrapper->UpperCutValueFunction()->GetTitle() == upperConstantString);
      } else
        wrapper->SetUpperCutValueFunction(new TF1(Form("upperCutFunction_%s", cutSection.c_str()), upperConstantString.c_str(), 0.0, 1.0e6));
    }
  }
}

TagAndProbeRatioFitSettings DynamicSelector::ParseTagAndProbeRatioFitSettings(const std::string& cutSection) {

  TagAndProbeRatioFitSettings settings;
  const std::string& parameterString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioFitParameters");
  if (parameterString.empty())
    return settings;

  const std::string& trustString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioFitTrustSource");
  if (!trustString.empty()) {
    static std::string trustDataString = "TrustData";
    static std::string trustMeanDataMCString = "TrustMeanDataMC";
    bool trustData = StringStartsWith(trustString, trustDataString);
    bool trustMeanDataMC = StringStartsWith(trustString, trustMeanDataMCString);
    if (!trustData && !trustMeanDataMC) {
      FATAL_OUT << "'TagAndProbeRatioFitTrustSource' is set to an invalid value. It has to begin with either 'TrustData' or 'TrustMeanDataMC'." << std::endl;
    }

    if (trustData)
      settings.fitTrustSource = TrustData;
    else
      settings.fitTrustSource = TrustMeanDataMC;
  }

  const std::string& breakMinimumString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioBreakMinimum");
  if (!breakMinimumString.empty())
    settings.fitBreak1Minimum = StringToNumber<float>(breakMinimumString);

  const std::string& breakMaximumString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioBreakMaximum");
  if (!breakMaximumString.empty())
    settings.fitBreak1Maximum = StringToNumber<float>(breakMaximumString);

  const std::string& fixedBreakString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioFitFixedBreak");
  if (!fixedBreakString.empty())
    settings.fitFixedBreak1 = StringToNumber<float>(fixedBreakString);

  const std::string& fixedConstantAboveString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioFitFixedConstantAbove");
  if (!fixedConstantAboveString.empty())
    settings.fitFixedBreak2 = StringToNumber<float>(fixedConstantAboveString);

  const std::string& keepRelativeErrorString = fConfig.GetStringValue(cutSection, "TagAndProbeRatioFitKeepRelativeErrorConstantAboveEnergy");
  if (!keepRelativeErrorString.empty())
    settings.fitKeepRelativeErrorConstantAboveEnergy = StringToNumber<float>(keepRelativeErrorString);

  static std::string straightLineWithBreakEndingInConstantString = "FitStraightLineWithBreakInLogScaleEndingInAConstant";
  static std::string straightLineWithBreakString = "FitStraightLineWithBreakInLogScale";
  static std::string straightLineString = "FitStraightLineInLogScale";
  static std::string constantString = "FitConstant";
  static std::string nothingString = "FitNothing";
  bool isStraightLineWithBreakEndingInConstant = StringStartsWith(parameterString, straightLineWithBreakEndingInConstantString);
  bool isStraightLineWithBreak = StringStartsWith(parameterString, straightLineWithBreakString);
  bool isStraightLine = StringStartsWith(parameterString, straightLineString);
  bool isConstant = StringStartsWith(parameterString, constantString);
  bool isNothing = StringStartsWith(parameterString, nothingString);
  if (!isStraightLineWithBreakEndingInConstant && !isStraightLineWithBreak && !isStraightLine && !isConstant && !isNothing) {
    FATAL_OUT << "'TagAndProbeRatioFitParameters' is set to an invalid value. "
                "It has to begin with either 'FitStraightLineWithBreakInLogScaleEndingInAConstant' or 'FitStraightLineWithBreakInLogScale' or 'FitStraightLineInLogScale' or 'FitConstant' or 'FitNothing'." << std::endl;
  }

  std::vector<std::string> parameterParts;
  StringTokenize(parameterString, '|', parameterParts);
  if (parameterParts.size() != 4 && !isNothing) {
    FATAL_OUT << "'TagAndProbeRatioFitParameters' syntax is wrong. It must look like eg. FitConstant(xMin|xMax|yMin|yMax)." << std::endl;
  }

  if (isStraightLineWithBreakEndingInConstant) {
    StringReplace(parameterParts[0], straightLineWithBreakEndingInConstantString, "");
    settings.fitMode = FitStraightLineWithBreakInLogScaleEndingInAConstant;
  } else if (isStraightLineWithBreak) {
    StringReplace(parameterParts[0], straightLineWithBreakString, "");
    settings.fitMode = FitStraightLineWithBreakInLogScale;
  } else if (isStraightLine) {
    StringReplace(parameterParts[0], straightLineString, "");
    settings.fitMode = FitStraightLineInLogScale;
  } else if (isConstant) {
    StringReplace(parameterParts[0], constantString, "");
    settings.fitMode = FitConstant;
  }

  if (!isNothing) {
    StringReplace(parameterParts[0], "(", "");
    StringReplace(parameterParts[3], ")", "");

    settings.fitXMin = StringToNumber<float>(parameterParts[0]);
    settings.fitXMax = StringToNumber<float>(parameterParts[1]);
    settings.fitYMin = StringToNumber<float>(parameterParts[2]);
    settings.fitYMax = StringToNumber<float>(parameterParts[3]);
  }

std::string excludedTagCutsString;
fConfig.GetValue(cutSection, "TagAndProbeRatioTemplateFitExcludedTags", excludedTagCutsString, "Comma-separated list of tag cuts to be excluded for template fit distribution filling.");

std::vector<std::string> listOfExcludedTagCuts = split(excludedTagCutsString, ",");
for (std::string& cutSection : listOfExcludedTagCuts)
  cutSection = strip(cutSection);

settings.excludedTagCutsForTemplateFit = std::set<std::string>(listOfExcludedTagCuts.begin(), listOfExcludedTagCuts.end());
return settings;
}

bool DynamicSelector::Examine(const IO::TreeReader& treeReader) {

  return fSelector->Passes(treeReader.DummyEvent());
}

void DynamicSelector::PrintSummary() const {

  assert(fSelector);
  fSelector->PrintSummary();
}

float DynamicSelector::GetXAxisRangeStart(const std::string& cutName) const {

  auto it = fCutXRangeStart.find(cutName);
  if (it == fCutXRangeStart.end())
    return gAMSStartShowEnergy;
  return it->second;
}

float DynamicSelector::GetYAxisRangeStart(const std::string& cutName) const {

  auto it = fCutYRangeStart.find(cutName);
  if (it == fCutYRangeStart.end())
    return 1.0e-5;
  return it->second;
}

float DynamicSelector::GetYAxisRangeStop(const std::string& cutName) const {

  auto it = fCutYRangeStop.find(cutName);
  if (it == fCutYRangeStop.end())
    return 2.0e-1;
  return it->second;
}

const std::string& DynamicSelector::GetYAxisTitle(const std::string& cutName) const {

  auto it = fCutYAxisTitles.find(cutName);
  assert(it != fCutYAxisTitles.end());
  return it->second;
}

const TagAndProbeRatioFitSettings& DynamicSelector::GetTagAndProbeRatioFitSettings(const std::string& cutName) const {

  auto it = fTagAndProbeRatioFitSettings.find(cutName);
  assert(it != fTagAndProbeRatioFitSettings.end());
  return it->second;
}
