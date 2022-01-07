#include "NucleiAnalysisSelectionTagProbe.hh"
#include "NucleiAnalysisTree.hh"
#include "SelectionCuts.hh"

// ACsoft includes
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "Environment.hh"
#include "ObjectManager.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "Cut.hh"
#include "CutAttachment.hh"
#include "CutFactory.hh"
#include "ValueHistograms.hh"
#include "LithiumBinning.hh"
#include "EfficiencyHistograms.hh"


#include "AnalysisSettings.hh"
#include "DynamicSelector.hh"
#include "TriggerEfficiency.hh"

#include <TEfficiency.h>

// ROOT includes
#include <TFile.h>
#include <TH1.h>
  
#define INFO_OUT_TAG "NucleiAnalysisSelectionTagProbe"
#include "debugging.hh"



NucleiAnalysisSelectionTagProbe::NucleiAnalysisSelectionTagProbe() 
  : IO::TreeReader()
  , fTree(nullptr)
  , fCutsObjectManager(0)
  , fTreeTofSelection(nullptr)
  , fTreeRichNaFSelection(nullptr)
  , fTreeRichAglSelection(nullptr)
  , fSelectorTof(nullptr)
  , fSelectorRichNaF(nullptr)
  , fSelectorRichAgl(nullptr)
  , fRichNaFSelectionEfficiencyInAcceptance(0)
  , fRichAglSelectionEfficiencyInAcceptance(0){
}

NucleiAnalysisSelectionTagProbe::~NucleiAnalysisSelectionTagProbe() {

  delete fTree;
  delete fTreeTofSelection;
  delete fTreeRichNaFSelection;
  delete fTreeRichAglSelection;
  delete fSelectorTof;
  delete fSelectorRichNaF;
  delete fSelectorRichAgl;
  delete fRichNaFSelectionEfficiencyInAcceptance;
  delete fRichAglSelectionEfficiencyInAcceptance;
  delete fCutsObjectManager;
}

IO::TreeInterface* NucleiAnalysisSelectionTagProbe::Initialize(Utilities::ObjectManager& objectManager) {

  assert(!fTree);
  fTree = new NucleiAnalysisTree;
  
  assert(fCutsObjectManager);
  AnalysisSettings::Initialize();

  RegisterTreeCuts(fTree);  

  assert(!fSelectorRichNaF);
  fSelectorTof = new DynamicSelector("${MYANALYSIS}/Configuration/NucleiAnalysis_TagAndProbe.cfg", "TofSelection");
  fSelectorTof->Initialize(objectManager, *fCutsObjectManager, Cuts::AxisInformation::Expression::Rigidity, fTree, this);

  fSelectorRichNaF = new DynamicSelector("${MYANALYSIS}/Configuration/NucleiAnalysis_TagAndProbe.cfg", "RichNaFSelection");
  fSelectorRichNaF->Initialize(objectManager, *fCutsObjectManager, Cuts::AxisInformation::Expression::Rigidity, fTree, this);

  fSelectorRichAgl = new DynamicSelector("${MYANALYSIS}/Configuration/NucleiAnalysis_TagAndProbe.cfg", "RichAglSelection");
  fSelectorRichAgl->Initialize(objectManager, *fCutsObjectManager, Cuts::AxisInformation::Expression::Rigidity, fTree, this);

  Binning::Predefined::SetDefaultAbsoluteRigidityBinning(LithiumRigidityBinning());
  fRichNaFSelectionEfficiencyInAcceptance = objectManager.Add(Make<TEfficiency>("fRichNaFSelectionEfficiencyInAcceptance", "", LithiumRigidityBinning()));
  fRichAglSelectionEfficiencyInAcceptance = objectManager.Add(Make<TEfficiency>("fRichAglSelectionEfficiencyInAcceptance", "", LithiumRigidityBinning()));
  
  // Return the tree interface pointer to the TreeReader, so it knows which kind of tree we're analzing.
  return fTree;
}

void NucleiAnalysisSelectionTagProbe::TreeChanged(TTree* tree) {

  IO::TreeReader::TreeChanged(tree);

    // Book tree(s).
  // The idea is to clone the very same tree that we're analyzing, but only write out events that fulfil our selection criteria.
  if (!fTreeTofSelection) {                                                                // if statement?
    assert(fOutputFile);
    fTreeTofSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeTofSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeTof")));
  }


  if (!fTreeRichNaFSelection) {                                                                // if statement?
    assert(fOutputFile);
    fTreeRichNaFSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeRichNaFSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeRichNaF")));
  }

  if (!fTreeRichAglSelection) {                                                                // if statement?
    assert(fOutputFile);
    fTreeRichAglSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeRichAglSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeRichAgl")));
  }
}

void NucleiAnalysisSelectionTagProbe::ParseOptions(Utilities::ConfigHandler& config) {
  // You could register custom command-line options here.
  if (!fCutsObjectManager)
    fCutsObjectManager = new Utilities::ObjectManager("CutsObjectManager", &config, fOutputDirectory, fOutputSuffix);

  IO::TreeReader::ParseOptions(config);
}

bool NucleiAnalysisSelectionTagProbe::ProcessEvent() {
  //FixMe:in the BeFluxesTreeWriter if particle HasTrackerTrack = false then TrackerCharges.size = 0
   assert(fTree->TrackerCharges().size() > 0);  // We store exactly 9 values for each tracker layer in this branch.

  if (fSelectorTof->Examine(*this)) 
    fTreeTofSelection->Fill();
  
  if (fSelectorRichNaF->Examine(*this)) {
    fRichNaFSelectionEfficiencyInAcceptance->Fill(true, fTree->RichCharge());
    fTreeRichNaFSelection->Fill();
  } else {
    fRichNaFSelectionEfficiencyInAcceptance->Fill(false, fTree->RichCharge());
  }
  
  if (fSelectorRichAgl->Examine(*this)) {
    fRichAglSelectionEfficiencyInAcceptance->Fill(true, fTree->RichCharge());
    fTreeRichAglSelection->Fill();
  } else {
    fRichAglSelectionEfficiencyInAcceptance->Fill(false, fTree->RichCharge());
   }
  
  return true;
}

void NucleiAnalysisSelectionTagProbe::Finish() {

  fOutputFile->cd();

  assert(fSelectorTof);
  fSelectorTof->PrintSummary();
  if (fTreeTofSelection)
    fTreeTofSelection->Finish();

  assert(fSelectorRichNaF);
  fSelectorRichNaF->PrintSummary();
  if (fTreeRichNaFSelection)
    fTreeRichNaFSelection->Finish();

  assert(fSelectorRichAgl);
  fSelectorRichAgl->PrintSummary();
  if (fTreeRichAglSelection)
    fTreeRichAglSelection->Finish();

  if (fCutsObjectManager)
    fCutsObjectManager->Close(); 
}

