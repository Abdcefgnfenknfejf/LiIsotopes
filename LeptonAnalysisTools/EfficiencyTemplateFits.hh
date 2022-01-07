#ifndef EfficiencyTemplateFits_hh
#define EfficiencyTemplateFits_hh

#include "BinningDefinition.hh"
#include "CutAttachment.hh"
#include "TrackedObjectArray.hh"

#include <TH2F.h>

#include <set>

#include <functional>

namespace Analysis {
class Event;
};

namespace Binning {
class Definition;
}

namespace Cuts {
class Cut;

using DistributionFunction = std::function<double()>;

class EfficiencyTemplateFits : public CutAttachment {
public:
  EfficiencyTemplateFits();
  virtual ~EfficiencyTemplateFits();

  void SetDistributionFunctionAndBinning(DistributionFunction function, const Binning::Definition& binning);
  void SetExcludedTagCuts(const std::set<std::string>& excludedTagCuts);

  virtual const char* GetName() const { return "Cut efficiency (template fits)"; }

  TH2F* TotalHistogramForCategory(const std::string& category) const;
  TH2F* PassedHistogramForCategory(const std::string& category) const;

public:
  static CutAttachment::ActivationMode sActivationState;
  static bool sUseWeights;

private:
  virtual void CreateObjects();
  virtual void Fill(const Analysis::Event& event);
  virtual void AddObjectsToMerge(std::vector<TObject**>& list);
  virtual void AddEmptyObjectsForNewTagCategory(unsigned int numberOfTagCutsInThisCategory);

private:
  Cuts::TrackedObjectArray<TH2F> fTagAndProbeHistos; ///< Internal vector of histograms, stored as TObjArray, to be browsable in an interactive ROOT session.
  bool fAssertionCheckPerformedForPassedCategory;    //!<! Internal bool to check wether time-consuming assertions were already performed.
  DistributionFunction fDistributionFunction;        //!<! Internal
  Binning::Definition fDistributionBinning;          //!<! Internal
  std::set<std::string> fExcludedTagCuts;            //!<! Internal

  ClassDef(Cuts::EfficiencyTemplateFits, 1)
};

};

#endif
