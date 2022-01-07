#ifndef TrdEstimatorTools_hh
#define TrdEstimatorTools_hh

#include <utility>
#include <string>
#include <sstream>

#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "LeptonAnalysisTree.hh"

TH1D* CreateCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);

TH1D* CreateTrdLikelihoodRatioHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH1D* CreateTrdLikelihoodProductHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH1D* CreateTrdLikelihoodHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio);

TH2D* CreateTrdLikelihoodRatioVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH2D* CreateTrdLikelihoodProductVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH2D* CreateTrdLikelihoodVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio);

TH2D* CreateEcalBDTvsTrdLikelihoodRatioHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH2D* CreateEcalBDTvsTrdLikelihoodProductHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description);
TH2D* CreateEcalBDTvsTrdLikelihoodHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio);

#endif
