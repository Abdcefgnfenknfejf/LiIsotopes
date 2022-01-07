#include "TemplateFitTree.hh"

#include <cmath>

#define INFO_OUT_TAG "TemplateFitTree"
#include "debugging.hh"

TemplateFitTree::TemplateFitTree(const std::string& treeName)
  : IO::TreeInterface(treeName.c_str(), Form("Template fit tree (%s)", treeName.c_str()))
  , fInputTree(nullptr) {

  RegisterBranches();
}

void TemplateFitTree::SetInputTree(const LeptonAnalysisTree* tree) {

  assert(tree);
  fInputTree = tree;
}

void TemplateFitTree::Fill(const Analysis::Event&) {

  if (!fInputTree)
    FATAL_OUT << "fInputTree=0. You forgot to call SetInputTree()." << std::endl;

  Weight = fInputTree->McEventWeightWithCutOff();
  McGeneratedMomentum = fInputTree->McGeneratedMomentum();
  ElectronCCMVABDT = fInputTree->ElectronCCMVABDT();
  TrackerRigidity = fInputTree->TrackerRigidity();
  TrackerPattern = fInputTree->TrackerPattern();
  TrackerNumberOfTracks = fInputTree->TrackerNumberOfTracks();
  EcalChiSquareLateralNormalized = fInputTree->EcalChiSquareLateralNormalized();
  EcalBDTBest = fInputTree->EcalBDTBest();
  EcalBDTBestSmoothed = fInputTree->EcalBDTBestSmoothed();
  EcalEnergyBestMaximumShower = fInputTree->EcalEnergyBestMaximumShower();
  TrdLRElecProt_Energy_HybridHits_TrdP = fInputTree->TrdLRElecProt_Energy_HybridHits_TrdP();
  TrdLRElecProt_Energy_HybridHits_TrdP_times_q = std::abs(fInputTree->TrdLRElecProt_Energy_HybridHits_TrdP()) * (TrackerRigidity() < 0 ? -1.0 : 1.0);
}

void TemplateFitTree::UpdateInMemoryBranches() {

}
