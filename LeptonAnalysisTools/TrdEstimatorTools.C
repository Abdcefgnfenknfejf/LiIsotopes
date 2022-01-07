#include "TrdEstimatorTools.hh"

#include "AnalysisSettings.hh"

TH1D* CreateCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "CCMVA - " << description;

  static const Binning::Definition sBinning = Binning::Equidistant(100, -1.1, 1.1);
  return Make<TH1D>(name.str().c_str(), title.str().c_str(), sBinning);
}

TH1D* CreateTrdLikelihoodRatioHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "TRD -log(L_{e/p}) - " << description;

  static const Binning::Definition sBinning = Binning::Equidistant(50, 0.0, 2.5);
  return Make<TH1D>(name.str().c_str(), title.str().c_str(), sBinning);
}

TH1D* CreateTrdLikelihoodProductHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "TRD -log10(L_{e}) - " << description;

  static const Binning::Definition sBinning = Binning::Equidistant(70, 2.0, 4.8);
  return Make<TH1D>(name.str().c_str(), title.str().c_str(), sBinning);
}

TH2D* CreateTrdLikelihoodRatioVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "TRD -log(L_{e/p}) - " << description;

  static const Binning::Definition sTrdBinning = Binning::Equidistant(100, -2.5, 2.5);
  static const Binning::Definition sCCMVABinning = Binning::Equidistant(100, -1.1, 1.1);
  return Make<TH2D>(name.str().c_str(), title.str().c_str(), sTrdBinning, sCCMVABinning);
}

TH2D* CreateTrdLikelihoodProductVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "TRD -log10(L_{e}) - " << description;

  static const Binning::Definition sTrdBinning = Binning::Equidistant(140, -4.8, 4.8);
  static const Binning::Definition sCCMVABinning = Binning::Equidistant(100, -1.1, 1.1);
  return Make<TH2D>(name.str().c_str(), title.str().c_str(), sTrdBinning, sCCMVABinning);
}

TH1D* CreateTrdLikelihoodHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio) {

  if (isRatio)
    return CreateTrdLikelihoodRatioHistogram(bin, namePrefix, description);
  return CreateTrdLikelihoodProductHistogram(bin, namePrefix, description);
}

TH2D* CreateTrdLikelihoodVsCCMVAHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio) {

  if (isRatio)
    return CreateTrdLikelihoodRatioVsCCMVAHistogram(bin, namePrefix, description);
  return CreateTrdLikelihoodProductVsCCMVAHistogram(bin, namePrefix, description);
}

TH2D* CreateEcalBDTvsTrdLikelihoodRatioHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "EcalBDT vs. TRD -log(L_{e/p}) - " << description;

  static const Binning::Definition sEcalBDTBinning = Binning::Equidistant(41, -1, 1);
  static const Binning::Definition sTrdBinning = Binning::Equidistant(50, 0.0, 2.5);

  return Make<TH2D>(name.str().c_str(), title.str().c_str(), sEcalBDTBinning, sTrdBinning);
}

TH2D* CreateEcalBDTvsTrdLikelihoodProductHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description) {

  std::stringstream name;
  name << namePrefix << "_bin_" << bin;

  std::stringstream title;
  title << "EcalBDT vs. TRD -log10(L_{e}) - " << description;

  static const Binning::Definition sEcalBDTBinning = Binning::Equidistant(41, -1, 1);
  static const Binning::Definition sTrdBinning = Binning::Equidistant(70, 2.0, 4.8);

  return Make<TH2D>(name.str().c_str(), title.str().c_str(), sEcalBDTBinning, sTrdBinning);
}

TH2D* CreateEcalBDTvsTrdLikelihoodHistogram(unsigned int bin, const std::string& namePrefix, const std::string& description, bool isRatio) {

  if (isRatio)
    return CreateEcalBDTvsTrdLikelihoodRatioHistogram(bin, namePrefix, description);
  return CreateEcalBDTvsTrdLikelihoodProductHistogram(bin, namePrefix, description);
}
