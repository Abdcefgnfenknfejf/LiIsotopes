#include "LithiumPreSelection.hh"
#include "NucleiAnalysisTree.hh"

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
#include "ValueHistograms.hh"
#include "CutAttachment.hh"
#include "AcceptanceManager.hh"
#include "EfficiencyHistograms.hh"
#include "McSpectrumScaler.hh"  
#include "LithiumBinning.hh"
#include "SelectionCuts.hh"
#include "FileManagerController.hh"


#include "TreeFormula.hh"
#include "TreeWriter.hh"
#include "Utilities.hh"
#include <TROOT.h>


// ROOT includes                                                                                                                                                                                                                                                                
#include <TFile.h>
#include <TH1.h>

#define INFO_OUT_TAG "LithiumPreSelection"
#include "debugging.hh"

LithiumPreSelection::LithiumPreSelection()
  : IO::TreeReader()
  , fTree(nullptr)
  , fTreeAfterSelectionWriter(nullptr)
  , fSelector(nullptr) {
  
}

LithiumPreSelection::~LithiumPreSelection() {

  delete fTree;
  delete fTreeAfterSelectionWriter;
  delete fSelector;
}

IO::TreeInterface* LithiumPreSelection::Initialize(Utilities::ObjectManager& objectManager) {

  assert(!fTree);
  fTree = new NucleiAnalysisTree;

  // Load config file and parse it.                                                                                  
     std::string configFile = "/home/lk974211/lithiumanalysis/Configuration/LithiumPreSelection.cfg";
     Environment::ExpandEnvironmentVariables(configFile);

     Utilities::ConfigHandler config;
     config.Read(configFile);

     RegisterTreeCuts(fTree);
     Binning::Predefined::SetDefaultAbsoluteRigidityBinning(LithiumRigidityBinning());  

     Cuts::SelectionParser parser(config, this);
     fSelector = objectManager.Add(parser.GetSelector("PreSelection"));

     // In order to fill any cut value histograms, the x axis value must be known.
     // We define this in a generic way using a lamdba function and pass it on
     // to Selector::SetupCommonXAxisInformation() for each selector.

     


     auto xRigidityAxisValue = [this] (const Analysis::Event&, double& lastBaseValue) {
				 if (fTree->Rigidity() != 0)
				   lastBaseValue = std::abs(fTree->Rigidity());
			       };
   
     auto xRigidityAxisValueMc = [this] (const Analysis::Event&, double& lastBaseValue) {
				   if (fTree->McGeneratedMomentum() != 0)
      lastBaseValue = fTree->McGeneratedMomentum()/3.0;  //For Li Charge = 3.0
				 };
     
     fSelector->SetupCommonXAxisInformation(xRigidityAxisValue, "|Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning); 
     fSelector->SetupCommonXAxisInformationMc(xRigidityAxisValueMc,  "|true Rigidity| / GV", Binning::Predefined::AbsoluteRigidityBinning);
     // Return the tree interface pointer to the TreeReader, so it knows which kind of tree we're analzing.
     return fTree;
}

void LithiumPreSelection::TreeChanged(TTree* tree) {

  IO::TreeReader::TreeChanged(tree);

  // Book tree(s).
  // The idea is to clone the very same tree that we're analyzing, but only write out events that fulfil our selection criteria.                                                                          
  if (!fTreeAfterSelectionWriter) {
    assert(fOutputFile);
    fTreeAfterSelectionWriter = new IO::TreeWriter(fTree, IO::TreeOptions::WriteInMemoryBranches);
    fTreeAfterSelectionWriter->Initialize(fOutputFile);
  }
}


void LithiumPreSelection::ParseOptions(Utilities::ConfigHandler& config) {
  // You could register custom command-line options here.                                                                                                                                           
  IO::TreeReader::ParseOptions(config);
}

bool LithiumPreSelection::ProcessEvent() {

  if (fSelector->Passes(DummyEvent())) {
    fTreeAfterSelectionWriter->Fill();
  }

  return true;
}

void LithiumPreSelection::Finish() {

  fOutputFile->cd();

  assert(fSelector);
  fSelector->PrintSummary();

  if (fTreeAfterSelectionWriter)
    fTreeAfterSelectionWriter->Finish();
}

