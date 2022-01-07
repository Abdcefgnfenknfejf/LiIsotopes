#ifndef NucleiAnalysisSelection_hh
#define NucleiAnalysisSelection_hh

#include "TreeReader.hh"

class NucleiAnalysisTree;
class TH1;

namespace Cuts {
class Selector;
}

namespace IO {
class TreeWriter;
}

class NucleiAnalysisSelection : public IO::TreeReader {
public:
  NucleiAnalysisSelection();
  virtual ~NucleiAnalysisSelection();

private:
  ClassDef(NucleiAnalysisSelection, 0)  //What does 0 here mean?

  virtual IO::TreeInterface* Initialize(Utilities::ObjectManager&);
  virtual void ParseOptions(Utilities::ConfigHandler&);
  virtual void TreeChanged(TTree*);
  virtual bool ProcessEvent();
  virtual void Finish();

private:
  NucleiAnalysisTree* fTree;                    //! Interface to our tree
  
  IO::TreeWriter* fTreeTofSelection;         //! Tree writer to write reduced tree after Tof selection
  IO::TreeWriter* fTreeRichNaFSelection;     //! Tree writer to write reduced tree after RichNaF selection
  IO::TreeWriter* fTreeRichAglSelection;     //! Tree writer to write reduced tree after RichAgl selection
  
  Cuts::Selector* fSelectorTof;                     //!Tof Selection 
  Cuts::Selector* fSelectorRichNaF;                     //! RichNaF selection 
  Cuts::Selector* fSelectorRichAgl;                     //! RichAgl selection 
  
};

#endif
