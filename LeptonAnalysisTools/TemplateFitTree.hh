#ifndef TemplateFitTree_hh
#define TemplateFitTree_hh

#include "LeptonAnalysisTree.hh"

class TemplateFitTree : public IO::TreeInterface {
public:
  TemplateFitTree(const std::string& treeName);

  // A LeptonAnalysisTree is required to fill the template fit tree.
  // Be sure to set it before filling the template fit tree.
  void SetInputTree(const LeptonAnalysisTree* tree);

  // Tree content
  IO::TreeBranch<Double_t> Weight { "Weight", -999.0 };
  IO::TreeBranch<Float_t> McGeneratedMomentum { "McGeneratedMomentum", -999.0 };
  IO::TreeBranch<Float_t> ElectronCCMVABDT { "ElectronCCMVABDT", -999.0 };
  IO::TreeBranch<Float_t> TrackerRigidity { "TrackerRigidity", 0.0 };
  IO::TreeBranch<Char_t> TrackerPattern { "TrackerPattern", -1 };
  IO::TreeBranch<UChar_t> TrackerNumberOfTracks { "TrackerNumberOfTracks", 0 };
  IO::TreeBranch<Float_t> EcalChiSquareLateralNormalized { "EcalChiSquareLateralNormalized", -999.0 };
  IO::TreeBranch<Float_t> EcalBDTBest { "EcalBDTBest", -999.0 };
  IO::TreeBranch<Float_t> EcalBDTBestSmoothed { "EcalBDTBestSmoothed", -999.0 };
  IO::TreeBranch<Float_t> EcalEnergyBestMaximumShower { "EcalEnergyBestMaximumShower", -999.0 };
  IO::TreeBranch<Float_t> TrdLRElecProt_Energy_HybridHits_TrdP { "TrdLRElecProt_Energy_HybridHits_TrdP", -999.0 };
  IO::TreeBranch<Float_t> TrdLRElecProt_Energy_HybridHits_TrdP_times_q { "TrdLRElecProt_Energy_HybridHits_TrdP_times_q", -999.0 };

private:
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new TemplateFitTree(""); }
  virtual void Fill(const Analysis::Event&);

private:
  const LeptonAnalysisTree* fInputTree;
};

#endif
