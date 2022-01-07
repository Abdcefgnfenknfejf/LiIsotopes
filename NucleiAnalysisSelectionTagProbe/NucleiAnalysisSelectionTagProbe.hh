#ifndef NucleiAnalysisSelectionTagProbe_hh
#define NucleiAnalysisSelectionTagProbe_hh

#include "TreeReader.hh"

class DynamicSelector;
class TEfficiency;
class NucleiAnalysisTree;
class TH1;

namespace Cuts {
class Selector;
}

namespace IO {
class TreeWriter;
}

namespace Utilities {
  class ConfigHandler;
  class ObjectManager;
}

class NucleiAnalysisSelectionTagProbe : public IO::TreeReader {
public:
  NucleiAnalysisSelectionTagProbe();
  virtual ~NucleiAnalysisSelectionTagProbe();

private:
  ClassDef(NucleiAnalysisSelectionTagProbe, 1)  //What does 0 here mean?

  virtual IO::TreeInterface* Initialize(Utilities::ObjectManager&);
  virtual void ParseOptions(Utilities::ConfigHandler&);
  virtual void TreeChanged(TTree*);
  virtual bool ProcessEvent();
  virtual void Finish();

private:
  NucleiAnalysisTree* fTree;                    //! Interface to our tree
  Utilities::ObjectManager* fCutsObjectManager;  
  IO::TreeWriter* fTreeTofSelection;     //! Tree writer to write reduced tree after RichNaF selection
  IO::TreeWriter* fTreeRichNaFSelection;     //! Tree writer to write reduced tree after RichNaF selection
  IO::TreeWriter* fTreeRichAglSelection;     //! Tree writer to write reduced tree after RichAgl selection                                                                                                                            
  DynamicSelector* fSelectorTof;           //! Selection cuts
  DynamicSelector* fSelectorRichNaF;           //! Selection cuts
  DynamicSelector* fSelectorRichAgl;           //! Selection cuts

  TEfficiency* fRichNaFSelectionEfficiencyInAcceptance; //! Overall efficiency (for statistics only!
  TEfficiency* fRichAglSelectionEfficiencyInAcceptance; //! Overall efficiency (for statistics only!)

  
};

#endif
