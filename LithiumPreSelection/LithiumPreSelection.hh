#ifndef LithiumPreSelection_hh
#define LithiumPreSelection_hh

#include "AcceptanceManager.hh"
#include "TreeReader.hh"

class NucleiAnalysisTree;
class TH1;

namespace Cuts {
class Selector;
}

namespace IO {
class TreeWriter;
}

class LithiumPreSelection : public IO::TreeReader {
public:
  LithiumPreSelection();
  virtual ~LithiumPreSelection();

private:
  ClassDef(LithiumPreSelection, 0)

  virtual IO::TreeInterface* Initialize(Utilities::ObjectManager&);
  virtual void ParseOptions(Utilities::ConfigHandler&);
  virtual void TreeChanged(TTree*);
  virtual bool ProcessEvent();
  virtual void Finish();

private:
  NucleiAnalysisTree* fTree;                    //! Interface to our tree                                                                                                                  
  IO::TreeWriter* fTreeAfterSelectionWriter;     //! Tree writer to write reduced tree                                                                                    
  Cuts::Selector* fSelector;                     //! Selector for 'selection cuts' based on config file
  
};

#endif
