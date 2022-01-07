#include "RooFitTools.hh"

#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>

#if defined(__GNUC__) &&  __GNUC__ >= 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Waddress"
#pragma GCC diagnostic ignored "-Wnonnull-compare"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif

#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooHist.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooWorkspace.h>

#if defined(__GNUC__) &&  __GNUC__ >= 6
#pragma GCC diagnostic pop
#endif

#include <cassert>
#include <cmath>

#include "ObjectManager.hh"

#define INFO_OUT_TAG "RooFitTools"
#include "debugging.hh"

void ForceDontDrawEmptyBins(RooPlot* frame, const std::string& name) {

  assert(frame);

  auto* graph = frame->getHist(name.c_str());
  assert(graph);

  for (int point = 0; point < graph->GetN(); ++point) {
    if (graph->GetY()[point] == 0.0) {
      graph->SetPointEYlow(point, 0.0);
      graph->SetPointEYhigh(point, 0.0);
    }
  }
}

TH1D* CreateTH1DFromPdf(RooAbsPdf* pdf, const std::string& name, const RooAbsRealLValue& xVar, int xBins) {

  if (xBins == -1)
    xBins = xVar.getBins();

  TH1D* histogram = new TH1D(name.c_str(), "", xBins, xVar.getMin(), xVar.getMax());

  RooArgList plotVars(xVar);

  bool scaleForDensity = true;
  double scaleFactor = 1.0;
  if (pdf->extendMode() != RooAbsPdf::CanNotBeExtended) {
    scaleFactor = pdf->expectedEvents(plotVars);
    scaleForDensity = false;
  }

  // Check that the plot variables are all actually RooRealVars and print a warning if we do not
  // explicitly depend on one of them. Fill a set (not list!) of cloned plot variables.
  RooArgSet plotClones;
  for (int index = 0; index < plotVars.getSize(); ++index) {
    const RooRealVar* realVar = dynamic_cast<const RooRealVar*>(plotVars.at(index));
    assert(realVar);
    assert(pdf->dependsOn(*realVar));
    plotClones.addClone(*realVar, true); // do not complain about duplicates
  }

  // Reconnect all plotClones to each other, imported when plotting N-dim integrals with entangled parameterized ranges
  TIterator* plotClonesIterator = plotClones.createIterator();
  RooAbsArg* plotClone = nullptr;
  while ((plotClone = static_cast<RooAbsArg*>(plotClonesIterator->Next())))
    plotClone->recursiveRedirectServers(plotClones, false, false, true);
  delete plotClonesIterator;

  // Call checkObservables
  RooArgSet allDeps(plotClones);
  assert(!pdf->checkObservables(&allDeps));

  // Create a standalone projection object to use for calculating bin contents
  RooArgSet* cloneSet = nullptr;
  const RooAbsReal* projected = pdf->createPlotProjection(plotClones, nullptr, cloneSet, nullptr, nullptr);

  // Prepare to loop over the histogram bins
  RooRealVar* xVarPlot = dynamic_cast<RooRealVar*>(plotClones.find(xVar.GetName()));
  assert(xVarPlot);

  TAxis* xAxis = histogram->GetXaxis();
  assert(xAxis);

  if (scaleForDensity)
    scaleFactor *= (xAxis->GetXmax() - xAxis->GetXmin()) / xBins;

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

  int xBin = 0;
  int bins = xBins;

  for (int bin = 0; bin < bins; ++bin) {
    xBin = (xBin % xBins) + 1;
    xVarPlot->setVal(xAxis->GetBinCenter(xBin));

    double result = scaleFactor * projected->getVal();
    if (RooAbsReal::numEvalErrors() > 0) {
      WARN_OUT << "Function evaluation error(s) at coordinates [x]=" << xVarPlot->getVal() << std::endl;
      result = 0.0;
    }
    RooAbsReal::clearEvalErrorLog();

    histogram->SetBinContent(histogram->GetBin(xBin), result);
    histogram->SetBinError(histogram->GetBin(xBin), std::sqrt(result));
  }

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);

  delete cloneSet;
  return histogram;
}

TH2D* CreateTH2DFromPdf(RooAbsPdf* pdf, const std::string& name, const RooAbsRealLValue& xVar, const RooAbsRealLValue& yVar, int xBins, int yBins) {

  if (xBins == -1)
    xBins = xVar.getBins();

  if (yBins == -1)
    yBins = yVar.getBins();

  TH2D* histogram = new TH2D(name.c_str(), "", xBins, xVar.getMin(), xVar.getMax(),
                                               yBins, yVar.getMin(), yVar.getMax());

  RooArgList plotVars(xVar, yVar);

  bool scaleForDensity = true;
  double scaleFactor = 1.0;
  if (pdf->extendMode() != RooAbsPdf::CanNotBeExtended) {
    scaleFactor = pdf->expectedEvents(plotVars);
    scaleForDensity = false;
  }

  // Check that the plot variables are all actually RooRealVars and print a warning if we do not
  // explicitly depend on one of them. Fill a set (not list!) of cloned plot variables.
  RooArgSet plotClones;
  for (int index = 0; index < plotVars.getSize(); ++index) {
    const RooRealVar* realVar = dynamic_cast<const RooRealVar*>(plotVars.at(index));
    assert(realVar);
    assert(pdf->dependsOn(*realVar));
    plotClones.addClone(*realVar, true); // do not complain about duplicates
  }

  // Reconnect all plotClones to each other, imported when plotting N-dim integrals with entangled parameterized ranges
  TIterator* plotClonesIterator = plotClones.createIterator();
  RooAbsArg* plotClone = nullptr;
  while ((plotClone = static_cast<RooAbsArg*>(plotClonesIterator->Next())))
    plotClone->recursiveRedirectServers(plotClones, false, false, true);
  delete plotClonesIterator;

  // Call checkObservables
  RooArgSet allDeps(plotClones);
  assert(!pdf->checkObservables(&allDeps));

  // Create a standalone projection object to use for calculating bin contents
  RooArgSet* cloneSet = nullptr;
  const RooAbsReal* projected = pdf->createPlotProjection(plotClones, nullptr, cloneSet, nullptr, nullptr);

  // Prepare to loop over the histogram bins
  RooRealVar* xVarPlot = dynamic_cast<RooRealVar*>(plotClones.find(xVar.GetName()));
  assert(xVarPlot);

  RooRealVar* yVarPlot = dynamic_cast<RooRealVar*>(plotClones.find(yVar.GetName()));
  assert(yVarPlot);

  TAxis* xAxis = histogram->GetXaxis();
  assert(xAxis);

  TAxis* yAxis = histogram->GetYaxis();
  assert(yAxis);

  if (scaleForDensity) {
    scaleFactor *= (xAxis->GetXmax() - xAxis->GetXmin()) / xBins;
    scaleFactor *= (yAxis->GetXmax() - yAxis->GetXmin()) / yBins;
  }

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::CollectErrors);

  int xBin = 0;
  int yBin = 0;
  int bins = xBins * yBins;

  for (int bin = 0; bin < bins; ++bin) {
    if (bin % xBins == 0) {
      yBin = (yBin % yBins) + 1;
      yVarPlot->setVal(yAxis->GetBinCenter(yBin));
    }

    xBin = (xBin % xBins) + 1;
    xVarPlot->setVal(xAxis->GetBinCenter(xBin));

    double result = scaleFactor * projected->getVal();
    if (RooAbsReal::numEvalErrors() > 0) {
      WARN_OUT << "Function evaluation error(s) at coordinates [x]=" << xVarPlot->getVal() << " [y]=" << yVarPlot->getVal() << std::endl;
      result = 0.0;
    }
    RooAbsReal::clearEvalErrorLog();

    histogram->SetBinContent(histogram->GetBin(xBin, yBin), result);
    histogram->SetBinError(histogram->GetBin(xBin, yBin), std::sqrt(result));
  }

  RooAbsReal::setEvalErrorLoggingMode(RooAbsReal::PrintErrors);

  delete cloneSet;
  return histogram;
}

TH1D* CreateTH1DFromDataSet(RooDataSet* dataSet, const std::string& name, const RooAbsRealLValue& xVar, int xBins) {

  if (xBins == -1)
    xBins = xVar.getBins();

  RooAbsReal* xVarPlot = dynamic_cast<RooAbsReal*>(dataSet->get()->find(xVar.GetName()));
  assert(xVarPlot);

  TH1D* histogram = new TH1D(name.c_str(), "", xBins, xVar.getMin(), xVar.getMax());

  auto entries = dataSet->numEntries();
  for (decltype(entries) entry = 0; entry < entries; ++entry) {
    dataSet->get(entry);
    histogram->Fill(xVarPlot->getVal(), dataSet->weight());
  }

  return histogram;
}

TH2D* CreateTH2DFromDataSet(RooDataSet* dataSet, const std::string& name, const RooAbsRealLValue& xVar, const RooAbsRealLValue& yVar, int xBins, int yBins) {

  if (xBins == -1)
    xBins = xVar.getBins();

  if (yBins == -1)
    yBins = yVar.getBins();

  RooAbsReal* xVarPlot = dynamic_cast<RooAbsReal*>(dataSet->get()->find(xVar.GetName()));
  assert(xVarPlot);

  RooAbsReal* yVarPlot = dynamic_cast<RooAbsReal*>(dataSet->get()->find(yVar.GetName()));
  assert(yVarPlot);

  TH2D* histogram = new TH2D(name.c_str(), "", xBins, xVar.getMin(), xVar.getMax(),
                                               yBins, yVar.getMin(), yVar.getMax());

  auto entries = dataSet->numEntries();
  for (decltype(entries) entry = 0; entry < entries; ++entry) {
    dataSet->get(entry);
    histogram->Fill(xVarPlot->getVal(), yVarPlot->getVal(), dataSet->weight());
  }

  return histogram;
}

void FixPDFParameters(RooWorkspace* workspace, const std::vector<std::string>& variableNames) {

  for (const auto& variableName : variableNames) {
    if (auto* variable = workspace->var(variableName.c_str()))
       variable->setConstant(true);
  }
}

void FreePDFParameters(RooWorkspace* workspace, const std::vector<std::string>& variableNames) {

  for (const auto& variableName : variableNames) {
    if (auto* variable = workspace->var(variableName.c_str())) {
       variable->setConstant(false);
       variable->removeError();
    }
  }
}

RooAbsData* ObtainBinnedClone(RooAbsData* data, bool shouldDeleteOriginal) {

 auto* binnedData = dynamic_cast<RooDataSet*>(data)->binnedClone();
 if (shouldDeleteOriginal)
   delete data;
 return binnedData;
}

void TransferParametersIntoWorkspace(const RooArgList& parameters, RooWorkspace* workspace, const std::string& variableSuffix) {

  auto* it = parameters.createIterator();

  std::string prefixedVariableSuffix;
  if (!variableSuffix.empty())
    prefixedVariableSuffix = Form("_%s", variableSuffix.c_str());

  RooRealVar* parameter = nullptr;
  while (((parameter = dynamic_cast<RooRealVar*>(it->Next())))) {
    if (auto* targetVariable = workspace->var(Form("%s%s", parameter->GetName(), prefixedVariableSuffix.c_str()))) {
      assert(targetVariable->isConstant());
      targetVariable->setVal(parameter->getVal());
    }
  }

  delete it;
}

void FixAllParametersInWorkspace(RooWorkspace* workspace) {

 const auto& allVariables = workspace->allVars();

  TIterator* it = allVariables.createIterator();
  RooAbsArg* argument = nullptr;
  while ((argument = dynamic_cast<RooAbsArg*>(it->Next()))) {
    auto* variable = workspace->var(argument->GetName());
    assert(variable);
    variable->setConstant(true);
  }

  delete it;
}

void ImportIntoTrdAllTracksWorkspace(RooWorkspace* workspaceAllTracks, RooWorkspace* workspaceSingleTrack, RooWorkspace* workspaceMultiTracks, const std::string& trdEstimatorName, unsigned int importObjects, const std::string& importSuffix) {

  RooAbsPdf* electronWithCorrectRigidityModelSingleTrack = workspaceSingleTrack->pdf("electronSignalModel");
  RooAbsPdf* protonWithWrongRigidityModelSingleTrack = workspaceSingleTrack->pdf("ccProtonSignalModel");
  RooAbsPdf* protonWithCorrectRigidityModelSingleTrack = workspaceSingleTrack->pdf("protonSignalModel");
  if (!electronWithCorrectRigidityModelSingleTrack || !protonWithWrongRigidityModelSingleTrack || !protonWithCorrectRigidityModelSingleTrack)
    FATAL_OUT << "Cannot load single-track PDFs. Aborting!" << std::endl;

  RooAbsPdf* electronWithCorrectRigidityModelMultiTracks = workspaceMultiTracks->pdf("electronSignalModel");
  RooAbsPdf* protonWithWrongRigidityModelMultiTracks = workspaceMultiTracks->pdf("ccProtonSignalModel");
  RooAbsPdf* protonWithCorrectRigidityModelMultiTracks = workspaceMultiTracks->pdf("protonSignalModel");
  if (!electronWithCorrectRigidityModelMultiTracks || !protonWithWrongRigidityModelMultiTracks || !protonWithCorrectRigidityModelMultiTracks)
    FATAL_OUT << "Cannot load multi-tracks PDFs. Aborting!" << std::endl;

  std::string prefixedImportSuffix = importSuffix;
  if (!importSuffix.empty())
    prefixedImportSuffix = Form("_%s", importSuffix.c_str());

  RooRealVar fractionOfSingleTrackEvents { Form("fractionOfSingleTrackEvents%s", prefixedImportSuffix.c_str()), "", 0.5, 0.0, 1.0 };
  fractionOfSingleTrackEvents.setConstant(true);

  const std::string& singleTrackVariableSuffix = Form("SingleTrack%s", prefixedImportSuffix.c_str());
  const std::string& multiTracksVariableSuffix = Form("MultiTracks%s", prefixedImportSuffix.c_str());

  if (importObjects & ImportElectronTemplate) {
    workspaceAllTracks->import(*electronWithCorrectRigidityModelSingleTrack,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(singleTrackVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(singleTrackVariableSuffix.c_str()));

    workspaceAllTracks->import(*electronWithCorrectRigidityModelMultiTracks,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(multiTracksVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(multiTracksVariableSuffix.c_str()));

    electronWithCorrectRigidityModelSingleTrack = workspaceAllTracks->pdf(Form("electronSignalModel_%s", singleTrackVariableSuffix.c_str()));
    electronWithCorrectRigidityModelMultiTracks = workspaceAllTracks->pdf(Form("electronSignalModel_%s", multiTracksVariableSuffix.c_str()));

    auto* electronWithCorrectRigidityModel = new RooAddPdf(Form("electronSignalModel_AllTracks%s", prefixedImportSuffix.c_str()), "@0 + @1",
                                                           RooArgList(*electronWithCorrectRigidityModelSingleTrack, *electronWithCorrectRigidityModelMultiTracks),
                                                           fractionOfSingleTrackEvents);

    workspaceAllTracks->import(*electronWithCorrectRigidityModel, RooFit::RecycleConflictNodes());
  }

  if (importObjects & ImportCCProtonTemplate) {
    workspaceAllTracks->import(*protonWithWrongRigidityModelSingleTrack,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(singleTrackVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(singleTrackVariableSuffix.c_str()));

    workspaceAllTracks->import(*protonWithWrongRigidityModelMultiTracks,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(multiTracksVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(multiTracksVariableSuffix.c_str()));

    protonWithWrongRigidityModelSingleTrack = workspaceAllTracks->pdf(Form("ccProtonSignalModel_%s", singleTrackVariableSuffix.c_str()));
    protonWithWrongRigidityModelMultiTracks = workspaceAllTracks->pdf(Form("ccProtonSignalModel_%s", multiTracksVariableSuffix.c_str()));

    auto* protonWithWrongRigidityModel = new RooAddPdf(Form("ccProtonSignalModel_AllTracks%s", prefixedImportSuffix.c_str()), "@0 + @1",
                                                       RooArgList(*protonWithWrongRigidityModelSingleTrack, *protonWithWrongRigidityModelMultiTracks),
                                                       fractionOfSingleTrackEvents);

    workspaceAllTracks->import(*protonWithWrongRigidityModel, RooFit::RecycleConflictNodes());
  }

  if (importObjects & ImportProtonTemplate) {
    workspaceAllTracks->import(*protonWithCorrectRigidityModelSingleTrack,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(singleTrackVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(singleTrackVariableSuffix.c_str()));

    workspaceAllTracks->import(*protonWithCorrectRigidityModelMultiTracks,
                               RooFit::RecycleConflictNodes(),
                               RooFit::RenameAllVariablesExcept(multiTracksVariableSuffix.c_str(), trdEstimatorName.c_str()),
                               RooFit::RenameAllNodes(multiTracksVariableSuffix.c_str()));

    protonWithCorrectRigidityModelSingleTrack = workspaceAllTracks->pdf(Form("protonSignalModel_%s", singleTrackVariableSuffix.c_str()));
    protonWithCorrectRigidityModelMultiTracks = workspaceAllTracks->pdf(Form("protonSignalModel_%s", multiTracksVariableSuffix.c_str()));

    auto* protonWithCorrectRigidityModel = new RooAddPdf(Form("protonSignalModel_AllTracks%s", prefixedImportSuffix.c_str()), "@0 + @1",
                                                         RooArgList(*protonWithCorrectRigidityModelSingleTrack, *protonWithCorrectRigidityModelMultiTracks),
                                                         fractionOfSingleTrackEvents);

    workspaceAllTracks->import(*protonWithCorrectRigidityModel, RooFit::RecycleConflictNodes());
  }
}

RooWorkspace* CreateTrdAllTracksWorkspace(RooWorkspace* workspaceSingleTrack, RooWorkspace* workspaceMultiTracks, const std::string& trdEstimatorName, unsigned int importObjects) {

  auto* workspaceAllTracks = new RooWorkspace("fWorkspaceAllTracks");
  ImportIntoTrdAllTracksWorkspace(workspaceAllTracks, workspaceSingleTrack, workspaceMultiTracks, trdEstimatorName, importObjects, "");
  return workspaceAllTracks;
}

void SetupTrdBinningForBoth1DAnd2DFit(RooRealVar& trdEstimator, unsigned int energyBin) {

  // Used by Scripts/TemplateFit/1D/FitTrdLikelihoodTemplates.C
  if (energyBin < 3)
    trdEstimator.setBins(350);
  else if (energyBin < 4)
    trdEstimator.setBins(450);
  else if (energyBin < 6)
    trdEstimator.setBins(500);
  else if (energyBin < 7)
    trdEstimator.setBins(600);
  else if (energyBin < 9)
    trdEstimator.setBins(650);
  else if (energyBin < 10)
    trdEstimator.setBins(800);
  else if (energyBin < 11)
    trdEstimator.setBins(650);
  else if (energyBin < 14)
    trdEstimator.setBins(600);
  else if (energyBin < 37)
    trdEstimator.setBins(450);
  else if (energyBin < 46)
    trdEstimator.setBins(400);
  else if (energyBin < 63)
    trdEstimator.setBins(350);
  else if (energyBin < 71)
    trdEstimator.setBins(150);
  else if (energyBin < 75)
    trdEstimator.setBins(100);
  else
    trdEstimator.setBins(50);
}

void SetupTrdBinningFor1DFit(RooRealVar& trdEstimator, unsigned int energyBin, bool isTimeAveragedData) {

  // Used by Scripts/TemplateFit/1D/PerformTrdEstimatorTemplateFits.C
  if (energyBin <= 5)
    trdEstimator.setBins(200);
  else if (energyBin <= 15)
    trdEstimator.setBins(400);
  else if (energyBin <= 20)
    trdEstimator.setBins(300);
  else if (energyBin <= 40)
    trdEstimator.setBins(200);
  else if (energyBin == 75)
    trdEstimator.setBins(50);
  else
    trdEstimator.setBins(100);

  (void) isTimeAveragedData;
}

void SetupTrdBinningFor2DFit(RooRealVar& trdEstimator, unsigned int energyBin, bool isTimeAveragedData) {

  SetupTrdBinningFor1DFit(trdEstimator, energyBin, isTimeAveragedData);
  trdEstimator.setBins(trdEstimator.getBins() * 2);
}

void SetupCCMVABinningFor2DFit(RooRealVar& ccmvaEstimator, unsigned int energyBin) {

  // Used by Scripts/TemplateFit/1D/EstimateCCMVATemplates.C and Scripts/TemplateFit/2D/Prepare2DFitData.C
  if (energyBin <= 5)
    ccmvaEstimator.setBins(300);
  else if (energyBin <= 15)
    ccmvaEstimator.setBins(600);
  else if (energyBin <= 20)
    ccmvaEstimator.setBins(400);
  else if (energyBin <= 30)
    ccmvaEstimator.setBins(300);
  else if (energyBin >= 74)
    ccmvaEstimator.setBins(75);
  else
    ccmvaEstimator.setBins(200);
}
