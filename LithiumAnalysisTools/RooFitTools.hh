#ifndef RooFitTools_hh
#define RooFitTools_hh

#include <string>
#include <vector>

class TH1D;
class TH2D;
class RooAbsPdf;
class RooAbsData;
class RooAbsRealLValue;
class RooArgList;
class RooDataSet;
class RooPlot;
class RooRealVar;
class RooWorkspace;

namespace Utilities {
class ObjectManager;
}

void ForceDontDrawEmptyBins(RooPlot* frame, const std::string& name);

TH1D* CreateTH1DFromPdf(RooAbsPdf* pdf, const std::string& name, const RooAbsRealLValue& xVar, int xBins = -1);
TH2D* CreateTH2DFromPdf(RooAbsPdf* pdf, const std::string& name, const RooAbsRealLValue& xVar, const RooAbsRealLValue& yVar, int xBins = -1, int yBins = -1);
TH1D* CreateTH1DFromDataSet(RooDataSet* dataSet, const std::string& name, const RooAbsRealLValue& xVar, int xBins = -1);
TH2D* CreateTH2DFromDataSet(RooDataSet* dataSet, const std::string& name, const RooAbsRealLValue& xVar, const RooAbsRealLValue& yVar, int xBins = -1, int yBins = -1);
RooAbsData* ObtainBinnedClone(RooAbsData* data, bool shouldDeleteOriginal = true);

void FixPDFParameters(RooWorkspace* workspace, const std::vector<std::string>& variableNames);
void FreePDFParameters(RooWorkspace* workspace, const std::vector<std::string>& variableNames);
void FixAllParametersInWorkspace(RooWorkspace* workspace);
void TransferParametersIntoWorkspace(const RooArgList& parameters, RooWorkspace* workspace, const std::string& variableSuffix = "");

enum ImportObjects {
  ImportElectronTemplate = 1 << 0,
  ImportCCProtonTemplate = 1 << 1,
  ImportProtonTemplate   = 1 << 2
};

RooWorkspace* CreateTrdAllTracksWorkspace(RooWorkspace* workspaceSingleTrack, RooWorkspace* workspaceMultiTracks, const std::string& trdEstimatorName, unsigned int importObjects);
void ImportIntoTrdAllTracksWorkspace(RooWorkspace* workspaceAllTracks, RooWorkspace* workspaceSingleTrack, RooWorkspace* workspaceMultiTracks, const std::string& trdEstimatorName, unsigned int importObjects, const std::string& importSuffix);

// 1D + 2D analysis
void SetupTrdBinningForBoth1DAnd2DFit(RooRealVar& trdEstimator, unsigned int energyBin);

// 1D analysis specific
void SetupTrdBinningFor1DFit(RooRealVar& trdEstimator, unsigned int energyBin, bool isTimeAveragedData);

// 2D analysis specific
void SetupTrdBinningFor2DFit(RooRealVar& trdEstimator, unsigned int energyBin, bool isTimeAveragedData);
void SetupCCMVABinningFor2DFit(RooRealVar& ccmvaEstimator, unsigned int energyBin);

#endif
