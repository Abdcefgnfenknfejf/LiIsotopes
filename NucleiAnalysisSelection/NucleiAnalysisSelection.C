#include "NucleiAnalysisSelection.hh"
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

// ROOT includes
#include <TFile.h>
#include <TH1.h>
  
#define INFO_OUT_TAG "NucleiAnalysisSelection"
#include "debugging.hh"



NucleiAnalysisSelection::NucleiAnalysisSelection()
  : IO::TreeReader()
  , fTree(nullptr)
  , fTreeTofSelection(nullptr)
  , fTreeRichNaFSelection(nullptr)
  , fTreeRichAglSelection(nullptr)
  , fSelectorTof(nullptr)
  , fSelectorRichNaF(nullptr)
  , fSelectorRichAgl(nullptr){
}

NucleiAnalysisSelection::~NucleiAnalysisSelection() {

  delete fTree;
  delete fTreeTofSelection;
  delete fTreeRichNaFSelection;
  delete fTreeRichAglSelection;
  delete fSelectorTof;
  delete fSelectorRichNaF;
  delete fSelectorRichAgl;
}

IO::TreeInterface* NucleiAnalysisSelection::Initialize(Utilities::ObjectManager& objectManager) {

  assert(!fTree);
  fTree = new NucleiAnalysisTree;
  
  // Load config file and parse it.
  std::string configFile = "${MYANALYSIS}/Configuration/NucleiAnalysisSelection.cfg";
  Environment::ExpandEnvironmentVariables(configFile);

  Utilities::ConfigHandler config;
  config.Read(configFile);


  RegisterTreeCuts(fTree);  

  // Book histogram(s).
  // Book selector(s).
  Cuts::SelectionParser parser(config, this);
  fSelectorTof = objectManager.Add(parser.GetSelector("TofSelection"));
  fSelectorRichNaF = objectManager.Add(parser.GetSelector("RichNaFSelection"));
  fSelectorRichAgl = objectManager.Add(parser.GetSelector("RichAglSelection"));
  // In order to fill any cut value histograms, the x axis value must be known.
  // We define this in a generic way using a lamdba function and pass it on
  // to Selector::SetupCommonXAxisInformation() for each selector.

  Binning::Predefined::SetDefaultAbsoluteRigidityBinning(LithiumRigidityBinning());
  auto xRigidityAxisValue = [this] (const Analysis::Event&, double& lastBaseValue) {
    if (fTree->Rigidity() != 0)
      lastBaseValue = std::abs(fTree->Rigidity());
  };
  
  auto xRigidityAxisValueMc = [this] (const Analysis::Event&, double& lastBaseValue) {
    if (fTree->McGeneratedMomentum() != 0)
      lastBaseValue = fTree->McGeneratedMomentum()/3.0;  //For Li Charge = 3.0                                                                                                                                                                                
  };

  fSelectorTof->SetupCommonXAxisInformation(xRigidityAxisValue, "|Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  fSelectorTof->SetupCommonXAxisInformationMc(xRigidityAxisValueMc,  "|true Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  
  fSelectorRichNaF->SetupCommonXAxisInformation(xRigidityAxisValue, "|Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  fSelectorRichNaF->SetupCommonXAxisInformationMc(xRigidityAxisValueMc,  "|true Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  
  fSelectorRichAgl->SetupCommonXAxisInformation(xRigidityAxisValue, "|Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  fSelectorRichAgl->SetupCommonXAxisInformationMc(xRigidityAxisValueMc,  "|true Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
  
  // Return the tree interface pointer to the TreeReader, so it knows which kind of tree we're analzing.
  return fTree;
}

void NucleiAnalysisSelection::TreeChanged(TTree* tree) {

  IO::TreeReader::TreeChanged(tree);

    // Book tree(s).
  // The idea is to clone the very same tree that we're analyzing, but only write out events that fulfil our selection criteria.
  if (!fTreeTofSelection) {                                                                // if statement?
    assert(fOutputFile);
    fTreeTofSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeTofSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeTof")));
    fTreeRichNaFSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeRichNaFSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeRichNaF")));
    fTreeRichAglSelection = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeRichAglSelection->Initialize(dynamic_cast<TDirectoryFile*>(fOutputFile->mkdir("TreeRichAgl")));    
  }
}

void NucleiAnalysisSelection::ParseOptions(Utilities::ConfigHandler& config) {
  // You could register custom command-line options here.
  IO::TreeReader::ParseOptions(config);
}

bool NucleiAnalysisSelection::ProcessEvent() {
  //FixMe:in the BeFluxesTreeWriter if particle HasTrackerTrack = false then TrackerCharges.size = 0
   assert(fTree->TrackerCharges().size() == 9); // We store exactly 9 values for each tracker layer in this branch.

   // Example call for the selector based on config file information.
  if (fSelectorTof->Passes(DummyEvent())) {
    fTreeTofSelection->Fill();
  }

  if (fSelectorRichNaF->Passes(DummyEvent())) {
    fTreeRichNaFSelection->Fill();
  }

  if (fSelectorRichAgl->Passes(DummyEvent())) {
    fTreeRichAglSelection->Fill();
  }
  
  return true;
}

void NucleiAnalysisSelection::Finish() {

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

}

