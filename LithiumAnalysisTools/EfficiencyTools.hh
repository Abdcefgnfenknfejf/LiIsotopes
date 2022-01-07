#ifndef EfficiencyTools_hh
#define EfficiencyTools_hh

#include <string>
#include <set>

#include <TF1.h>
#include <TFile.h>
#include <TH2F.h>

namespace Cuts {
  class Cut;
};

enum TagAndProbeRatioFitMode {
  FitNothing,
  FitConstant,
  FitStraightLineInLogScale,
  FitStraightLineWithBreakInLogScale,
  FitStraightLineWithBreakInLogScaleEndingInAConstant
};

enum TagAndProbeRatioFitTrustSource {
  TrustMeanDataMC,
  TrustData
};

class DynamicSelector;
class TGraphAsymmErrors;
class TGraphErrors;

struct TagAndProbeRatioFitSettings {
  TagAndProbeRatioFitSettings() = default;

  TagAndProbeRatioFitMode fitMode { FitNothing };
  TagAndProbeRatioFitTrustSource fitTrustSource { TrustData };
  float fitXMin { 0.0f };
  float fitXMax { 0.0f };
  float fitYMin { 0.8f };
  float fitYMax { 1.2f };
  float fitBreak1Minimum { 0.0f };
  float fitBreak1Maximum { 0.0f };
  float fitFixedBreak1 { 0.0f };
  float fitFixedBreak2 { 0.0f };
  float fitKeepRelativeErrorConstantAboveEnergy { 0.0f };
  std::set<std::string> excludedTagCutsForTemplateFit;
};

struct SingleEfficiencyData {
  SingleEfficiencyData()
    : fCutIndex(0)
    , fLowerCutFunction(nullptr)
    , fUpperCutFunction(nullptr)
    , fRangeYMin(0.0)
    , fRangeYMax(0.0)
    , issControlPlotElectron(nullptr)
    , issControlPlotProton(nullptr)
    , mcControlPlotElectron(nullptr)
    , mcTagAndProbeElectronPassedHistogram(nullptr)
    , mcTagAndProbeElectronTotalHistogram(nullptr)
    , issTagAndProbeElectronPassedHistogram(nullptr)
    , issTagAndProbeElectronTotalHistogram(nullptr)
    , mcCountingElectronGenerated(nullptr)
    , mcTagAndProbeElectron(nullptr)
    , mcTagAndProbeElectronTemplateFit(nullptr)
    , issTagAndProbeElectron(nullptr)
    , issTagAndProbeElectronTemplateFit(nullptr)
    , mcCountingProtonGenerated(nullptr)
    , mcTagAndProbeProton(nullptr)
    , issTagAndProbeProton(nullptr) {

  }

  void Initialize(TFile* effectiveAcceptanceFileElectron, TFile* effectiveAcceptanceFileProton, TFile* issTagAndProbeEfficiencyFile, TFile* mcTagAndProbeEfficiencyFile, const std::string& selectorName, unsigned int cutIndex);
  void DrawOverview();
  void DrawTagAndProbeRatio();
  TGraphErrors* FitTagAndProbeRatio(TagAndProbeRatioFitMode, float xMin, float xMax, float yMin, float yMax, float break1Minimum, float break1Maximum, float fixedBreak1, float fixedBreak2, float keepRelativeErrorConstantAboveEnergy, TagAndProbeRatioFitTrustSource);
  static TGraphErrors* FitBartelsRotationToAverageDataRatio(const SingleEfficiencyData& averageData, const SingleEfficiencyData& monthlyData, float yMin, float yMax, const std::string& monthOutputFileSuffix);

  void SetControlPlotYRange(float yMin, float yMax);

  std::string CutName() const;
  TF1* LowerCutFunction() const { return fLowerCutFunction; }
  TF1* UpperCutFunction() const { return fUpperCutFunction; }

private:
  void DrawControlPlotElectron();
  void DrawMonteCarloCountingElectron();
  void DrawTagAndProbeElectron(bool drawSame);
  void DrawLegendElectron();

  void DrawControlPlotProton();
  void DrawMonteCarloCountingProton();
  void DrawTagAndProbeProton(bool drawSame);
  void DrawLegendProton();

  void DrawControlPlotPlaceholder(const std::string& species);
  void DrawEfficiencyPlaceholder(const std::string& species);

  unsigned int fCutIndex;
  std::string fSelectorName;
  TF1* fLowerCutFunction;
  TF1* fUpperCutFunction;
  float fRangeYMin;
  float fRangeYMax;

public:
  TH2F* issControlPlotElectron;
  TH2F* issControlPlotProton;
  TH2F* mcControlPlotElectron;
  TH2F* mcTagAndProbeElectronPassedHistogram;
  TH2F* mcTagAndProbeElectronTotalHistogram;
  TH2F* issTagAndProbeElectronPassedHistogram;
  TH2F* issTagAndProbeElectronTotalHistogram;

  TGraphAsymmErrors* mcCountingElectronGenerated;
  TGraphAsymmErrors* mcTagAndProbeElectron;
  TGraphAsymmErrors* mcTagAndProbeElectronTemplateFit;
  TGraphAsymmErrors* issTagAndProbeElectron;
  TGraphAsymmErrors* issTagAndProbeElectronTemplateFit;

  TGraphAsymmErrors* mcCountingProtonGenerated;
  TGraphAsymmErrors* mcTagAndProbeProton;
  TGraphAsymmErrors* issTagAndProbeProton;
};

struct EfficiencyComparisionData {
  EfficiencyComparisionData()
    : fMCCountingElectronGeneratedOverallEfficiency(0)
    , fMCTagAndProbeElectronOverallEfficiency(0)
    , fISSTagAndProbeElectronOverallEfficiency(0) {

  }

  void Initialize(TFile* effectiveAcceptanceFileElectron, TFile* effectiveAcceptanceFileProton, TFile* issTagAndProbeEfficiencyFile, TFile* mcTagAndProbeEfficiencyFile, const std::string& selectorName, unsigned int numberOfCuts);
  void SetControlPlotYRange(unsigned int cutIndex, float yMin, float yMax);
  void DrawSingleCutComparisions();

  void CalculateOverallEfficiency();
  void SetOverallLegendPosition(float x1, float y1, float x2, float y2);
  void DrawOverallEfficiencyComparision();

  TGraphErrors* FitTagAndProbeRatio(unsigned int cutIndex, const TagAndProbeRatioFitSettings& settings);

  unsigned int NumberOfCuts() const { return fNumberOfCuts; }
  const SingleEfficiencyData& SingleCutEfficiency(unsigned int cutIndex) const { return fEfficiencies.at(cutIndex); }
  static void CorrectEfficiencyForTagAndProbeRatio(TGraphAsymmErrors* efficiency, TGraphErrors* tagAndProbeRatio);
  static void CorrectEfficiencyForTagAndProbeRatio(TGraphAsymmErrors* efficiency, double tagAndProbeRatio, double tagAndProbeRatioUncertainty);

  static TGraphErrors* CombineUncorrelatedEfficiencies(const std::vector<TGraphErrors*>& efficiencies);
  static TGraphAsymmErrors* CombineUncorrelatedEfficiencies(const std::vector<TGraphAsymmErrors*>& efficiencies);

private:
  std::string fSelectorName;
  unsigned int fNumberOfCuts;

  TGraphAsymmErrors* fMCCountingElectronGeneratedOverallEfficiency;
  TGraphAsymmErrors* fMCTagAndProbeElectronOverallEfficiency;
  TGraphAsymmErrors* fISSTagAndProbeElectronOverallEfficiency;
  std::vector<SingleEfficiencyData> fEfficiencies;
};

#endif // EfficiencyTools_hh
